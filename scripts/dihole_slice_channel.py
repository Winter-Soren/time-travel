"""Exact transverse slice-zero profiles and their full radial Dirac leakage.

One regular collar and azimuthal sector per q. This is not a full quantum
source, exterior scattering calculation, or atomic-payload propagation.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from balanced_dihole import symbolic_fields
from dihole_cylinder_transfer import zero_moments
from dihole_flux_scale import load_inputs
from dihole_local_source import local_first, mask_jets
from dihole_stretched_collar import exact_target_parameters, gauss, radius_profile
from project_paths import project_path


NAME = 'dihole_slice_channel'
FLUXES = (5, 25, 105, 1005, 10005)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local_fields(function, p, radius, theta):
    values, radial, angular = local_first(function, p, radius, theta)
    a, ar, at = [v[2]*p['physical_A_factor'] for v in (values, radial, angular)]
    return dict(W=np.exp(values[1]-values[0]), NP=np.exp(values[0]+values[1]),
        log_W_r=radial[1]-radial[0], log_NP_r=radial[0]+radial[1], a=a, ar=ar, at=at)


def profile_logs(function, p, u, nu, theta, inner_order=64):
    """Log unnormalized profile per q, and its exact u derivative.

    Known polar powers are factored before integrating a regular correction.
    The azimuthal factor exp(i*mu*phi)/sqrt(2*pi) is already normalized.
    """
    radius, rate, _ = radius_profile(u, order=128)
    radius, rate = float(radius), float(rate)
    theta = np.asarray(theta)
    points = local_fields(function, p, radius, theta)
    nodes, weights = gauss(inner_order)
    half = (theta-np.pi/2)/2
    t = np.pi/2+half[..., None]*(nodes+1)
    inner = local_fields(function, p, radius, t)
    correction = (inner['NP']*(nu-inner['a'])-nu+(1-np.cos(t))/2)/np.sin(t)
    derivative = inner['NP']*(inner['log_NP_r']*(nu-inner['a'])-inner['ar'])*rate/np.sin(t)
    log_profile = nu*np.log(np.sin(theta/2))+(1-nu)*np.log(np.cos(theta/2))
    log_profile += half*(correction@weights)
    return dict(theta=theta, log_profile_per_q=log_profile, log_profile_s_per_q=half*(derivative@weights),
        W=points['W'], log_W_s=points['log_W_r']*rate,
        angular_log_derivative_per_q=points['NP']*(nu-points['a'])/np.sin(theta),
        radius=radius, radius_rate=rate)


def prepare(function, p, u, nu, order, inner_order):
    nodes, weights = gauss(order)
    theta, weights = (nodes+1)*np.pi/2, weights*np.pi/2
    out = profile_logs(function, p, u, nu, theta, inner_order)
    out.update(weights=weights, u=u, nu=nu, order=order, inner_order=inner_order)
    return out


def channel(data, q):
    mu = q*data['nu']
    if type(q) is not int or q < 1 or not (.5 <= mu <= q-.5) or not math.isclose(mu%1, .5, abs_tol=1e-11):
        raise ValueError('Require integer flux q and a regular half-integer zero-mode sector.')
    log_profile = q*data['log_profile_per_q']
    amplitude = np.exp(log_profile-np.max(log_profile))
    norm = float(data['weights']@(data['W']*amplitude**2))
    profile = amplitude/math.sqrt(norm)
    probability = data['weights']*data['W']*profile**2
    log_s = q*data['log_profile_s_per_q']
    fs_over_f = log_s-probability@(log_s+data['log_W_s']/2)
    velocity = float(probability@(1/data['W']))
    velocity_s = float(2*probability@(fs_over_f/data['W']))
    r1 = 1/data['W']-velocity
    r0 = fs_over_f/data['W']-velocity_s/2
    mask = mask_jets(data['u']) if 0 < data['u'] < .05 else np.zeros(6)
    dot_r1, dot_r0 = 2*mask[0]*r1, 2*mask[0]*r0+mask[1]*r1
    square = lambda v: float(probability@(v*v))
    summary = dict(q=q, mu=mu, nu=data['nu'], u_over_D=data['u'], radius_over_D=data['radius'],
        order=data['order'], inner_order=data['inner_order'], mask=float(mask[0]), mask_derivative_D=float(mask[1]),
        normalization_error=abs(float(np.sum(probability))-1), projected_velocity=velocity,
        projected_velocity_derivative_D=velocity_s, r1_norm=math.sqrt(square(r1)),
        r0_norm_D=math.sqrt(square(r0)), dot_r1_norm=math.sqrt(square(dot_r1)),
        dot_r0_norm_D=math.sqrt(square(dot_r0)), r0_r1_inner_D=float(probability@(r0*r1)),
        projector_orthogonality_error=max(abs(float(probability@r0)), abs(float(probability@r1))),
        rms_angular_width=float(np.sqrt(probability@((data['theta']-probability@data['theta'])**2))))
    return summary, dict(profile=profile, probability=probability, fs_over_f=fs_over_f,
        r1=r1, r0=r0, dot_r1=dot_r1, dot_r0=dot_r0, mask=mask, velocity=velocity, velocity_s=velocity_s)


def localization(function, p, u, nu):
    radius, rate, _ = radius_profile(u, order=128)
    theta = brentq(lambda t: float(local_fields(function, p, radius, t)['a'])-nu, 1e-6, np.pi-1e-6, xtol=1e-14)
    local = local_fields(function, p, radius, theta)
    drift = float(-local['ar']*rate/local['at'])
    kappa = float(local['NP']*local['at']/np.sin(theta))
    if kappa <= 0:
        raise ArithmeticError('Nondegenerate localization maximum required.')
    return dict(theta_peak_limit=theta, peak_drift_D=drift, localization_kappa=kappa,
        W_at_peak=float(local['W']), tilt_tangent=float(radius*drift),
        leading_r0_over_sqrt_q=float(abs(drift)*math.sqrt(kappa/2)/local['W']),
        leading_transverse_gap_over_sqrt_q=float(math.sqrt(2*kappa)/(radius*local['W'])),
        leading_leakage_to_gap_ratio=float(abs(radius*drift)/2),
        scope='Local large-q Laplace coefficients only; no uniform finite-q or atomic-sector error bound.')


def controls(function, p, data, q):
    row, state = channel(data, q)
    tiny = 1e-25j
    points = np.array([.4, 1.2, 2.5])
    exact = profile_logs(function, p, data['u'], data['nu'], points, 96)['angular_log_derivative_per_q']
    numerical = np.imag(profile_logs(function, p, data['u'], data['nu'], points+tiny, 96)['log_profile_per_q'])/tiny.imag
    h = 1e-5
    fd = []
    for step in (h, h/2):
        states = [channel(prepare(function, p, data['u']+s*step, data['nu'], data['order'], data['inner_order']), q)
                  for s in (-1, 1)]
        derivative = (states[1][1]['profile']-states[0][1]['profile'])/(2*step)
        actual = state['profile']*state['fs_over_f']
        error = np.sqrt(data['weights']@(data['W']*(derivative-actual)**2))
        denom = np.sqrt(data['weights']@(data['W']*actual**2))
        speed = (states[1][0]['projected_velocity']-states[0][0]['projected_velocity'])/(2*step)
        fd.append(dict(step=step, normalized_profile_derivative_relative_error=float(error/denom),
            speed_derivative_relative_error=abs(speed/row['projected_velocity_derivative_D']-1)))
    mask, dm = state['mask'][:2]
    epsilon = 1e-25j
    we = np.exp(-2*epsilon*mask)*data['W']
    fe = np.exp(epsilon*mask)*state['profile']
    fse = np.exp(epsilon*mask)*state['profile']*(state['fs_over_f']+epsilon*dm)
    ve = data['weights']@(fe**2)
    vse = 2*data['weights']@(fe*fse)
    # Direct action coefficients in the common canonical Hamiltonian space.
    r1e = 1/we-ve
    r0e = (state['fs_over_f']+epsilon*dm)/we-vse/2
    errors = [float(np.sqrt(state['probability']@((np.imag(a)/epsilon.imag-b)**2)))
              for a, b in ((r1e, state['dot_r1']), (r0e, state['dot_r0']))]
    missing_half = abs(dm)*row['r1_norm']/2
    return dict(angular_ODE_absolute_error=float(np.max(abs(exact-numerical))),
        radial_derivative_checks=fd, common_H_metric_derivative_errors=errors,
        incorrect_half_connection_discrepancy=missing_half)


def aligned_leading_control(limits):
    """Local matrix cancellation only; no finite-q wave function or source bound."""
    sx = np.array([[0., 1.], [1., 0.]])
    sy = np.array([[0., -1j], [1j, 0.]])
    sz = np.diag([1., -1.])
    records = []
    for t in [-1., -.3, 0., .7, 2.]+[r['tilt_tangent'] for r in limits]:
        kappa, mask = .5, .7
        beta, width = np.arctan(t), kappa/np.sqrt(1+t*t)
        spin = np.array([np.cos(beta/2), np.sin(beta/2)])
        beta_dot = mask*t/(1+t*t)
        width_dot = -mask*kappa*t*t/(1+t*t)**1.5
        spin_dot = beta_dot*np.array([-spin[1], spin[0]])/2
        matrix = 1j*width*(sx-t*sz)-kappa*sy
        matrix_dot = 1j*(width_dot*(sx-t*sz)-width*mask*t*sz)
        wrong = 1j*kappa*(sx-t*sz)-kappa*sy
        records.append(dict(tilt=t, leading_residual=float(np.linalg.norm(matrix@spin)),
            paired_variation_residual=float(np.linalg.norm(matrix_dot@spin+matrix@spin_dot)),
            rotation_without_width_change_residual=float(np.linalg.norm(wrong@spin))))
    return dict(records=records, identity='[i G(sigma1-t sigma3)-kappa sigma2] eta=0; t=R theta0_s, G=kappa/sqrt(1+t^2), eta=(cos(beta/2),sin(beta/2)), beta=atan(t).',
        variation='t_dot=F t; beta_dot=F t/(1+t^2); G_dot=-F kappa t^2/(1+t^2)^(3/2).',
        scope='Only the leading local square-root-q term cancels. Derivatives, nonlinear potential, normalization and polar domains still require a full trial-state calculation.')


def run():
    _, inputs = load_inputs()
    for name in ('dihole_flux_scale.py', 'dihole_cylinder_transfer.py', 'spherical_weyl_leads.py', 'curved_throat_modes.py'):
        inputs[name] = sha(project_path(name))
    function, p = symbolic_fields()[0], exact_target_parameters()
    rows, finest = [], {}
    for order, inner in ((512, 48), (1024, 80)):
        for u in (.0125, .025, .0375):
            data = prepare(function, p, u, .9, order, inner)
            for q in FLUXES:
                row, _ = channel(data, q)
                rows.append(row)
            if order == 1024:
                finest[u] = data
    coarse, fine = rows[:15], rows[15:]
    metrics = ('projected_velocity', 'r0_norm_D', 'r1_norm', 'dot_r0_norm_D', 'dot_r1_norm')
    refinement = {key: max(abs(a[key]/b[key]-1) for a, b in zip(coarse, fine)) for key in metrics}
    checks = controls(function, p, finest[.0125], 105)
    central = prepare(function, p, 0., .5, 512, 80)
    central_row, _ = channel(central, 1)
    old = zero_moments(function, p, order=160, inner_order=80)
    baseline_error = abs(central_row['projected_velocity']*old['mean_W']-1)
    asymptotic = [dict(u_over_D=u, **localization(function, p, u, .9)) for u in finest]
    for limit in asymptotic:
        row = next(r for r in fine if r['q'] == FLUXES[-1] and r['u_over_D'] == limit['u_over_D'])
        limit['last_q_relative_difference'] = abs(row['r0_norm_D']/math.sqrt(row['q'])/limit['leading_r0_over_sqrt_q']-1)
    aligned = aligned_leading_control(asymptotic)
    tests = dict(exact_current_input_hashes=all(sha(project_path(n)) == h for n, h in inputs.items()),
        regular_half_integer_sectors=all(math.isclose(r['mu']%1, .5, abs_tol=1e-11) for r in rows),
        Hamiltonian_normalization=max(r['normalization_error'] for r in rows) < 2e-15,
        transverse_zero_ODE=checks['angular_ODE_absolute_error'] < 2e-12,
        normalized_radial_derivative=checks['radial_derivative_checks'][-1]['normalized_profile_derivative_relative_error'] < 2e-5,
        projected_speed_derivative=checks['radial_derivative_checks'][-1]['speed_derivative_relative_error'] < 2e-5,
        projection_orthogonality=max(r['projector_orthogonality_error'] for r in rows) < 2e-12,
        common_H_metric_derivative=max(checks['common_H_metric_derivative_errors']) < 2e-12,
        nonzero_mask_gradient_test=checks['incorrect_half_connection_discrepancy'] > 1e-8,
        angular_and_inner_quadrature=max(refinement.values()) < 2e-7,
        existing_q1_cylinder_recovered=baseline_error < 3e-13 and central_row['r0_norm_D'] == 0,
        high_flux_local_limit=max(r['last_q_relative_difference'] for r in asymptotic) < .01,
        aligned_leading_matrix=max(r['leading_residual'] for r in aligned['records']) < 2e-15
            and max(r['rotation_without_width_change_residual'] for r in aligned['records']) > .1,
        aligned_leading_metric_derivative=max(r['paired_variation_residual'] for r in aligned['records']) < 2e-15)
    tests = {name: bool(passed) for name, passed in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests, records=rows, quadrature_refinement=refinement,
        independent_controls=checks, existing_cylinder_velocity_relative_error=baseline_error,
        central_q1=central_row, local_large_q_limits=asymptotic, aligned_leading_control=aligned,
        operator='H=W^-1/2[-i sigma3 d_s+K/R]W^-1/2, K=-i sigma1 d_theta+sigma2 NP(mu-qa)/sin(theta).',
        projection='chi=sqrt(W)f e_plus, integral W f^2=1; v=integral f^2; h=-i(v d_s+v_s/2).',
        residual='(HU-Uh)psi=-i(r1 psi_s+r0 psi); r1=chi(1/W-v), r0=chi((f_s/f)/W-v_s/2).',
        metric_derivative='dot r1=2F r1; dot r0=2F r0+F_s r1; chi is unchanged in the common canonical Hamiltonian space.',
        scope=['Exact slice-zero ODE and radial/projection identities; displayed norms are numerical angular integrals with sensitivity checks.',
            'One upper Weyl zero mode at mu/q=.9, at three right-collar slices; no sum over sectors or full-source trace.',
            'The nonzero collar field allows local large-q localization; no global transverse gap follows.',
            'The frozen-slice packet center drifts along the collar. A field-aligned ansatz must rotate its spin and change its transverse profile; changing coordinates alone preserves the residual.',
            'No finite-q transverse gap enclosure, metric-source error bound, null-ring conversion or open-channel contribution is evaluated.',
            'No atomic flux is simulated here. The asymptotic coefficient has no quantified remainder at the nominal atomic flux.',
            'Norms use dimensionless D=1. In the fixed-coupling physical family D scales with q: r0 in physical units then scales as q^(-1/2), while its local gap ratio stays nonzero.',
            'No supporting material, stable opening, atom survival or time difference is established.'],
        sources=['https://arxiv.org/pdf/1807.04726', 'https://arxiv.org/pdf/hep-th/9906160'],
        full_source_bound=None, atomic_transit=None, source_sha256=sha(__file__), input_sha256=inputs)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Slice-channel controls failed; results retained.')
    plot(result)
    print(json.dumps({k: result[k] for k in ('passed', 'tests', 'quadrature_refinement', 'local_large_q_limits')}, indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.5), layout='constrained')
    for limit in result['local_large_q_limits']:
        rows = [r for r in result['records'] if r['order'] == 1024 and r['u_over_D'] == limit['u_over_D']]
        q = np.array([r['q'] for r in rows])
        label = f"Collar u / D = {limit['u_over_D']:.4f}"
        line, = axes[0].loglog(q, [r['r0_norm_D'] for r in rows], 'o-', label=label)
        axes[0].loglog(q, limit['leading_r0_over_sqrt_q']*np.sqrt(q), '--', color=line.get_color(), alpha=.55)
        axes[1].loglog(q, [r['dot_r0_norm_D'] for r in rows], 'o-', label=label)
    axes[0].set(xlabel='Integer magnetic flux q', ylabel='Radial leakage norm times D',
        title='Exact slice zero modes still couple along the collar')
    axes[1].set(xlabel='Integer magnetic flux q', ylabel='Geometry-variation leakage norm times D',
        title='The source variation retains that coupling')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('More flux alone does not validate frozen cross-section channels')
    fig.supxlabel('Dashed curves: local square-root-q limits. One mode per q; no global source bound or atom transit.', fontsize=8)
    for extension in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+extension), dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()
