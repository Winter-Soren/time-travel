"""Physical Euclidean Dirac form on compact bridge trial spinors.

Assemble the norm of the differential action BEFORE projecting its output.
The trial spinors vanish at both cuts and extend by zero to the exterior.
This is a conforming local building block, not the sewn-domain determinant
or vacuum source. The paired variation includes the changing flux density.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh

from balanced_dihole import symbolic_fields
from dihole_azimuthal_bridge import AzimuthalBridgeMatrices
from dihole_bridge_scattering import collar_mask
from dihole_stretched_collar import field_values, gauss, radius_profile, smooth_step
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


def mask_derivative(s):
    """Exact first derivative of the existing C-infinity collar mask."""
    t = (abs(float(s))-.05)/.05
    if not 0 < t < 1:
        return 0.
    step, first = smooth_step(t)
    hp = 1/t**2+1/(1-t)**2
    hpp = -2/t**3+2/(1-t)**3
    return float(np.sign(s)*first*((1-2*step)*hp+hpp/hp)/.1)


def trial_data(bridge, radial_modes=3, radial_order=64):
    """Pointwise trial values, radial/angular actions, and physical weights.

    Angular profiles are represented without sqrt(sin(theta)); quadrature
    is consequently dcos(theta). The two Weyls have opposite T actions.
    No projected inverse of P or W is substituted for multiplication.
    """
    if radial_modes < 1 or radial_order <= radial_modes or int(radial_modes) != radial_modes:
        raise ValueError('Positive integral radial modes and a larger quadrature order required.')
    nodes, weights = gauss(radial_order)
    ss, ww = [], []
    for lo, hi in ((-.1, -.05), (-.05, .05), (.05, .1)):
        ss.extend((hi+lo)/2+(hi-lo)*nodes/2)
        ww.extend((hi-lo)*weights/2)
    theta = bridge.theta
    upper, lower = bridge.profiles(theta)
    p, m, n = bridge.p, bridge.m, bridge.n
    profiles = np.zeros((len(theta), 2, n), complex)
    profiles[:, 0, :p], profiles[:, 1, p:] = upper, 1j*lower
    angular_reference = np.zeros_like(profiles)
    lam = np.sqrt(bridge.minus_labels*(bridge.minus_labels+1))
    angular_reference[:, 0, p:] = upper[:, int(bridge.a == 0):]*lam
    angular_reference[:, 1, int(bridge.a == 0):p] = 1j*lower*lam
    sigma_y_profiles = np.stack((-1j*profiles[:, 1], 1j*profiles[:, 0]), axis=1)
    rows = {key: [] for key in ('value', 'radial', 'angular', 'P', 'W', 'f', 'fp', 'weight')}
    modes = np.arange(1, radial_modes+1)
    for s, weight in zip(ss, ww):
        radius = float(radius_profile(min(.05, max(0., abs(s)-.05)))[0])
        fields = field_values(bridge.function, bridge.parameters, radius, bridge.cosine)
        lapse, meridian, potential = [np.asarray(fields[k], float) for k in
                                     ('lapse', 'meridian_factor', 'unit_A_phi')]
        delta = (lapse*meridian*(bridge.mu-potential)
                 -(bridge.a+.5*np.cos(theta)))/np.sin(theta)
        angular = (angular_reference+delta[:, None, None]*sigma_y_profiles)/radius
        radial = np.sqrt(10.)*np.sin(modes*np.pi*(s+.1)/.2)
        derivative = np.sqrt(10.)*modes*np.pi/.2*np.cos(modes*np.pi*(s+.1)/.2)
        value = (profiles[..., None]*radial).reshape(len(theta), 2, -1)
        action_f = (-1j*np.array([1., -1.])[None, :, None, None]
                    *profiles[..., None]*derivative).reshape(len(theta), 2, -1)
        action_a = (angular[..., None]*radial).reshape(len(theta), 2, -1)
        # ponytail: compact sine trials only; add sewn exterior trial functions
        # before using this assembly to approximate the physical determinant.
        zero = np.zeros_like(value)
        for key, data in (('value', value), ('radial', action_f), ('angular', action_a)):
            sign = 1 if key == 'value' else -1
            rows[key].append(np.concatenate((np.concatenate((data, zero), axis=2),
                                            np.concatenate((zero, sign*data), axis=2)), axis=1))
        for key, data in (('P', meridian), ('W', meridian/lapse),
                          ('f', np.full(len(theta), collar_mask(s))),
                          ('fp', np.full(len(theta), mask_derivative(s))),
                          ('weight', weight*bridge.weights)):
            rows[key].append(data)
    return {key: np.concatenate(value) for key, value in rows.items()}


def assemble(data, kappa, epsilon=0., mass=0.):
    """Common physical P_kappa form, analytic variation, and direct massive form.

    v_e=exp(e*f/2)v identifies L2(P_e) with baseline L2(P), P_e=exp(-e*f)P.
    L_e=P^-1/2[exp(-e*f)i*k*W + exp(e*f)(F-i*e*f'*Sigma/2)+A].
    q_e=L_e^*L_e; the massive form is q_e+m^2*n, on the same domain.
    """
    if not all(np.isfinite(x) for x in (kappa, epsilon, mass)) or mass < 0:
        raise ValueError('Finite frequency/variation and nonnegative proper mass required.')
    v, radial, angular = [data[k] for k in ('value', 'radial', 'angular')]
    p, w, f, fp = [data[k][:, None, None] for k in ('P', 'W', 'f', 'fp')]
    sigma = np.array([1., -1., -1., 1.])[None, :, None]
    temporal = 1j*kappa*w*v
    connection = -.5j*fp*sigma*v
    action = (np.exp(-epsilon*f)*temporal
              +np.exp(epsilon*f)*(radial+epsilon*connection)+angular)/np.sqrt(p)
    tangent = (-f*np.exp(-epsilon*f)*temporal
               +np.exp(epsilon*f)*(f*(radial+epsilon*connection)+connection))/np.sqrt(p)
    tangent_without_connection = (-f*temporal+f*radial)/np.sqrt(p)
    norm_values = np.sqrt(p)*v
    beta_action = action[:, [2, 3, 0, 1], :]
    weights = np.sqrt(data['weight'])[:, None, None]

    def gram(left, right):
        left, right = (weights*left).reshape(-1, v.shape[-1]), (weights*right).reshape(-1, v.shape[-1])
        return left.conj().T@right

    q, norm = gram(action, action), gram(norm_values, norm_values)
    dq = gram(tangent, action)+gram(action, tangent)
    wrong_dq = gram(tangent_without_connection, action)+gram(action, tangent_without_connection)
    weak_dirac = gram(norm_values, beta_action)
    direct_massive = gram(beta_action+mass*norm_values, beta_action+mass*norm_values)
    return dict(q=q, norm=norm, derivative=dq, without_connection=wrong_dq,
                weak_dirac=weak_dirac, direct_massive=direct_massive)


def relative_error(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-30))


def plot(result):
    """Show the prescribed coordinate bridge and the actual form diagnostic."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    s = np.linspace(-.1, .1, 501)
    radius = radius_profile(np.clip(abs(s)-.05, 0., .05))[0]
    axes[0].fill_between(s, -radius, radius, color='#dceef4')
    axes[0].plot(s, radius, color='#167b93', linewidth=2)
    axes[0].plot(s, -radius, color='#167b93', linewidth=2)
    for lo, hi in ((-.1, -.05), (.05, .1)):
        axes[0].axvspan(lo, hi, color='#f1b04a', alpha=.22)
    axes[0].text(0, 0, 'Proposed field region\nSupport still unknown', ha='center', va='center')
    axes[0].set(xlabel='Bridge coordinate s / D', ylabel='Coordinate radius R / D',
                title='Where the new form is checked', ylim=(-.065, .065))
    labels = [str(r['mu']) for r in result['records']]
    values = [100*r['projection_before_square_relative_deficit'] for r in result['records']]
    axes[1].bar(labels, values, color='#167b93', width=.55)
    for j, value in enumerate(values):
        axes[1].text(j, value+.15, f'{value:.2f}%', ha='center')
    axes[1].set(xlabel='Azimuthal sector mu', ylabel='Relative matrix norm difference (%)',
                title='Projecting before squaring loses terms', ylim=(0, 10))
    fig.suptitle('Physical Euclidean form: 10 checks pass; supporting stress remains unknown', fontsize=12)
    fig.supxlabel('Coordinate schematic, not an embedding or machine blueprint. Differences apply only to these finite trial spaces.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_euclidean_form_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    function = symbolic_fields()[0]
    parameters = exact_target_parameters()
    records = []
    for mu in (-.5, .5, 1.5):
        bridge = AzimuthalBridgeMatrices(function, parameters, 3, 48, mu=mu)
        data = trial_data(bridge)
        base = assemble(data, .7, mass=1.3)
        plus, minus = [assemble(data, .7, epsilon=e) for e in (1e-5, -1e-5)]
        refined = assemble(trial_data(bridge, radial_order=96), .7)
        angular_refined = assemble(trial_data(AzimuthalBridgeMatrices(
            function, parameters, 3, 72, mu=mu)), .7)
        projected_square = base['weak_dirac'].conj().T@np.linalg.solve(base['norm'], base['weak_dirac'])
        missed = base['q']-projected_square
        a = np.linalg.cholesky(base['norm'])
        inverse = np.linalg.inv(a)
        physical_p = inverse@base['q']@inverse.conj().T
        tangent = inverse@base['derivative']@inverse.conj().T
        mass = 1.7
        identity = np.eye(len(a))
        resolvents = [np.linalg.solve(physical_p+j*mass**2*identity, identity) for j in range(4)]
        summed = sum(c*r for c, r in zip((1, -3, 3, -1), resolvents))
        product = 6*mass**6*resolvents[0]@resolvents[1]@resolvents[2]@resolvents[3]
        pv = -.5*float(np.trace(tangent@product).real)
        variants = [inverse@row['q']@inverse.conj().T for row in (plus, minus)]
        actions = [-.5*sum(c*np.linalg.slogdet(matrix+j*mass**2*identity)[1]
                          for j, c in enumerate((1, -3, 3, -1))) for matrix in variants]
        varied_mass = assemble(data, .7, epsilon=.001, mass=1.3)
        records.append(dict(mu=mu, full_trial_dimension=len(a),
            min_generalized_eigenvalue=float(eigh(base['q'], base['norm'], eigvals_only=True)[0]),
            hermiticity_error=relative_error(base['q'], base['q'].conj().T),
            dirac_skew_error=relative_error(base['weak_dirac'], -base['weak_dirac'].conj().T),
            mass_shift_error=relative_error(base['direct_massive'], base['q']+1.3**2*base['norm']),
            varied_mass_shift_error=relative_error(varied_mass['direct_massive'],
                                                  varied_mass['q']+1.3**2*varied_mass['norm']),
            variation_fd_error=relative_error(base['derivative'], (plus['q']-minus['q'])/2e-5),
            norm_variation_error=relative_error(plus['norm'], minus['norm']),
            radial_quadrature_change=relative_error(base['q'], refined['q']),
            derivative_quadrature_change=relative_error(base['derivative'], refined['derivative']),
            angular_quadrature_change=relative_error(base['q'], angular_refined['q']),
            omitted_connection_relative_effect=relative_error(base['derivative'], base['without_connection']),
            projection_before_square_relative_deficit=float(np.linalg.norm(missed)/np.linalg.norm(base['q'])),
            projection_deficit_min_eigenvalue=float(eigh(missed, base['norm'], eigvals_only=True)[0]),
            pv_rational_error=relative_error(summed, product),
            pv_response=pv, pv_response_fd=(actions[0]-actions[1])/2e-5,
            pv_response_fd_error=float(abs(pv-(actions[0]-actions[1])/2e-5))))
    tests = dict(positive_common_form=all(r['min_generalized_eigenvalue'] > 0 for r in records),
        physical_skew_dirac=all(r['dirac_skew_error'] < 1e-9 for r in records),
        fixed_proper_mass_shift=all(max(r['mass_shift_error'], r['varied_mass_shift_error']) < 1e-11 for r in records),
        metric_variation=all(r['variation_fd_error'] < 2e-8 for r in records),
        fixed_euclidean_norm=all(r['norm_variation_error'] == 0 for r in records),
        quadrature_control=all(max(r['radial_quadrature_change'], r['derivative_quadrature_change'],
                                  r['angular_quadrature_change']) < 1e-8 for r in records),
        retained_output_matters=all(r['projection_deficit_min_eigenvalue'] > -1e-8
                                   and r['projection_before_square_relative_deficit'] > .001 for r in records),
        connection_matters=all(r['omitted_connection_relative_effect'] > .001 for r in records),
        pv_common_operator=all(r['pv_rational_error'] < 1e-9 for r in records),
        pv_differentiation=all(r['pv_response_fd_error'] < 2e-8 for r in records))
    names = ('dihole_azimuthal_bridge.py', 'dihole_stretched_collar.py', 'dihole_bridge_scattering.py',
             'balanced_dihole.py', 'complete_angular_modes.py', 'azimuthal_modes.py', 'spherical_weyl_leads.py')
    result = dict(complete=True, passed=all(tests.values()), tests=tests, records=records,
        kappa_D0=.7, proper_mass_D0=1.7, radial_modes=3, angular_paired_modes=3,
        radial_order=64, angular_order=48, seconds=perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
        scope=['Actual q1 prescribed bridge, compact H1 trial functions vanishing at both cuts.',
               'Common physical Euclidean quadratic form and metric variation; no exterior determinant.',
               'Projection deficit is a basis diagnostic, not an error bound for the physical source.',
               'Displayed PV derivative belongs only to this finite restricted form; no source sign inference.',
               'No physical angular/frequency/IR limits, local subtraction, or regulator removal performed.'])
    project_path('dihole_euclidean_form_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise RuntimeError('Physical form checks failed; evidence retained.')
    plot(result)
    return result


if __name__ == '__main__':
    run()
