"""Local large-q Dirac normal form at the balanced dihole's magnetic null.

The background is Emparan, hep-th/9906160v3, equations 8--16:
https://arxiv.org/pdf/hep-th/9906160 . The scaling and coefficient checks below
are our application, not a published scattering solution or a novelty claim.
Fixed Q/D geometry, unit-charge spinor, hbar=c=1; changing q also changes the
Maxwell/gravity coupling calibration. No quantum stress or backreaction is solved.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from balanced_dihole import symbolic_fields
from dihole_return_geometry import return_fraction, weyl_fields
from project_paths import project_path, source_matches


def coefficient_checks(parameters, geometry):
    """Independent 70-digit differentiation of the full published Weyl fields."""
    with mp.workdps(70):
        Q, k = mp.mpf(str(parameters['Q'])), mp.mpf(str(parameters['D']))/2
        M = (mp.sqrt(k*k+4*Q*k)-k)/2
        ag = mp.sqrt(k*k-M*M)
        B = M/((k+ag)*(M+k))
        h = M*(1-B*ag/2)
        rc = mp.findroot(lambda r: r**3-3*h*r*r+(M+h)*ag*ag-2*M*ag/B,
                         mp.mpf(str(geometry['critical_r'])))
        rhoc = mp.sqrt((rc-M)**2-k*k)

        def fields(rho, z):
            upper, lower = mp.sqrt(rho*rho+(z-k)**2), mp.sqrt(rho*rho+(z+k)**2)
            r, c = M+(upper+lower)/2, (lower-upper)/(2*k)
            delta, sigma = r*r-2*M*r-ag*ag, r*r-ag*ag*c*c
            H = delta+k*k*(1-c*c)
            W = (r*r-ag*ag)**2+delta*ag*ag*(1-c*c)
            N = (delta+ag*ag*(1-c*c)+2*B*M*r*ag*(1-c*c)+B*B*W*(1-c*c)/4)/sigma
            P = N*sigma*sigma/(H*H)
            a = (2*M*r*ag+B*W/2)*(1-c*c)/(2*Q*N*sigma)
            return N, P, a

        af = lambda rho, z: fields(rho, z)[2]
        vf = lambda rho, z: fields(rho, z)[0]/fields(rho, z)[1]
        bf = lambda rho, z: fields(rho, z)[0]**2/rho
        N, P, astar = fields(rhoc, mp.mpf(0))
        arr, azz = mp.diff(af, (rhoc, 0), (2, 0)), mp.diff(af, (rhoc, 0), (0, 2))
        ar, az = mp.diff(af, (rhoc, 0), (1, 0)), mp.diff(af, (rhoc, 0), (0, 1))
        denom = rc-2*M+2*B*M*ag+B*B*(rc**3-ag*ag*rc-2*M*ag*ag)/4
        a2_formula = 3*B*rc*(rc-2*h)*rhoc**2/(4*Q*denom**2*(rc-M)**2)
        maxwell = []
        # This tests the weighted Maxwell equation independently away from the null.
        for rho, z in ((rhoc, mp.mpf(0)), (mp.mpf('.55'), mp.mpf('.08')),
                       (mp.mpf('.8'), mp.mpf('-.11'))):
            br, bz = mp.diff(bf, (rho, z), (1, 0)), mp.diff(bf, (rho, z), (0, 1))
            div = (bf(rho, z)*(mp.diff(af, (rho, z), (2, 0))+mp.diff(af, (rho, z), (0, 2)))
                   +br*mp.diff(af, (rho, z), (1, 0))+bz*mp.diff(af, (rho, z), (0, 1)))
            maxwell.append(float(abs(div)))
        v, b, a2 = N/P, N*N/rhoc, arr/2
        L0 = (v/(b*a2))**(mp.mpf(1)/3)
        E0, label0 = v/L0, v/(L0*b)
        coeff = dict(rho_star=float(rhoc), critical_r=float(rc), a_star=float(astar),
                     N_star=float(N), P_star=float(P), v_star=float(v), b_star=float(b),
                     a_rhorho=float(arr), a_zz=float(azz), a2=float(a2),
                     v_rho=float(mp.diff(vf, (rhoc, 0), (1, 0))),
                     v_z=float(mp.diff(vf, (rhoc, 0), (0, 1))),
                     b_rho=float(mp.diff(bf, (rhoc, 0), (1, 0))),
                     a_rhorhorho=float(mp.diff(af, (rhoc, 0), (3, 0))),
                     a_rhozz=float(mp.diff(af, (rhoc, 0), (1, 2))),
                     width_prefactor=float(L0), energy_prefactor=float(E0),
                     label_prefactor=float(label0))
        controls = dict(digits=70, gradient_at_null=float(max(abs(ar), abs(az))),
                        Hessian_trace=float(abs(arr+azz)),
                        a2_formula_error=float(abs(a2-a2_formula)),
                        weighted_Maxwell_residuals=maxwell,
                        saved_critical_point_error=abs(float(rhoc)-geometry['null_ring_rho']))
    return coeff, controls


def physical_coefficients(function, p, rho, z):
    """Canonical flat-d(rho)dz measure: H=-i sigma_i(v d_i+v_i/2)+sigma_y b(mu-qa)."""
    N, P, a, _, _ = weyl_fields(function, p, rho, z)
    upper, lower = np.hypot(rho, z-p['k']), np.hypot(rho, z+p['k'])
    total = upper+lower
    r, c = p['M']+total/2, 2*z/total
    sr, sz = rho/upper+rho/lower, (z-p['k'])/upper+(z+p['k'])/lower
    H, S = (r-p['M'])**2-p['k']**2*c*c, r*r-p['a']**2*c*c
    v = H*H/(S*S)
    vr = v*(4*(r-p['M'])/H-4*r/S)
    vc = v*(-4*p['k']**2*c/H+4*p['a']**2*c/S)
    v_rho = vr*sr/2-vc*2*z*sr/total**2
    v_z = vr*sz/2+vc*(2/total-2*z*sz/total**2)
    return dict(v=v, v_rho=v_rho, v_z=v_z, b=N*N/rho, a=a,
                independent_v_error=float(np.max(abs(v-N/P))))


def local_comparison(function, p, c, q, grid_order=81):
    """Sampled local coefficient errors; no continuum supremum assertion."""
    L = c['width_prefactor']*q**(-1/3)
    E = c['energy_prefactor']*q**(1/3)
    label_width = E/c['b_star']
    mu = math.floor(q*c['a_star'])+.5
    eta = (mu-q*c['a_star'])/label_width
    axis = np.linspace(-1., 1., grid_order)
    s, t = np.meshgrid(axis, axis)
    fields = physical_coefficients(function, p, c['rho_star']+L*s, L*t)
    kinetic = fields['v']/c['v_star']
    connection_rho = L*fields['v_rho']/(2*c['v_star'])
    connection_z = L*fields['v_z']/(2*c['v_star'])
    mass = fields['b']*(mu-q*fields['a'])/E
    mass0 = eta-s*s+t*t
    dv = c['v_rho']/c['v_star']
    db = c['b_rho']/c['b_star']
    mass1 = mass0+L*(db*s*mass0-c['a_rhorhorho']*s**3/(6*c['a2'])
                     -c['a_rhozz']*s*t*t/(2*c['a2']))
    maxabs = lambda x: float(np.max(abs(x)))
    leading = dict(kinetic=maxabs(kinetic-1), connection_rho=maxabs(connection_rho),
                   connection_z=maxabs(connection_z), angular_mass=maxabs(mass-mass0))
    first = dict(kinetic=maxabs(kinetic-(1+L*dv*s)),
                 connection_rho=maxabs(connection_rho-L*dv/2),
                 connection_z=maxabs(connection_z), angular_mass=maxabs(mass-mass1))
    return dict(q=q, mu=mu, eta=eta, radial_width=L, energy_scale=E,
                label_width=label_width, normalized_label_width=label_width/q,
                grid_order=grid_order, box=[-1., 1.],
                leading_coefficient_errors=leading, first_order_remainder_errors=first,
                independent_v_error=fields['independent_v_error'],
                central_slice=dict(s=axis.tolist(), actual_mass=mass[grid_order//2].tolist(),
                                   normal_mass=mass0[grid_order//2].tolist()))


def checks():
    saved_path = project_path('dihole_return_geometry_results.json')
    saved = json.loads(saved_path.read_text(encoding='utf-8'))
    assert saved['passed'] and source_matches('dihole_return_geometry.py', saved['source_sha256'])
    p, geometry = saved['parameters'], saved['geometry']
    assert p['q'] == 1 and p['D'] == 1 and p['Q'] == .2
    recomputed = return_fraction(p)
    assert recomputed == geometry
    function = symbolic_fields()[0]
    coeff, independent = coefficient_checks(p, geometry)
    rows = [local_comparison(function, p, coeff, q) for q in (1000, 8000, 64000)]
    fine = [local_comparison(function, p, coeff, row['q'], 161) for row in rows]
    mesh_error = max(abs(x[k]-y[k]) for r, f in zip(rows, fine)
                     for key in ('leading_coefficient_errors', 'first_order_remainder_errors')
                     for x, y in [(r[key], f[key])] for k in x)
    steps = []
    for old, new in zip(rows, rows[1:]):
        steps.append(dict(q_pair=[old['q'], new['q']],
                          leading_mass_ratio=old['leading_coefficient_errors']['angular_mass']/new['leading_coefficient_errors']['angular_mass'],
                          first_order_mass_ratio=old['first_order_remainder_errors']['angular_mass']/new['first_order_remainder_errors']['angular_mass']))
    gradient_error = 0.
    for rho, z in ((.55, .08), (.8, -.11), (coeff['rho_star'], 0.)):
        exact = physical_coefficients(function, p, rho, z)
        eps = 1e-4
        for direction in (0, 1):
            dr, dz = (eps, 0.) if direction == 0 else (0., eps)
            vals = [float(physical_coefficients(function, p, rho+j*dr, z+j*dz)['v'])
                    for j in (-2, -1, 1, 2)]
            fd = (vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*eps)
            gradient_error = max(gradient_error, abs(fd-float(exact['v_rho' if direction == 0 else 'v_z'])))
    tests = dict(saved_critical_point=independent['saved_critical_point_error'] < 3e-14,
                 critical_gradient=independent['gradient_at_null'] < 1e-60,
                 full_Maxwell_identity=max(independent['weighted_Maxwell_residuals']) < 1e-60,
                 Hessian_saddle=independent['Hessian_trace'] < 1e-60 and coeff['a2'] > 0,
                 analytic_a2=independent['a2_formula_error'] < 1e-60,
                 independent_metric_speed=max(r['independent_v_error'] for r in rows) < 2e-14,
                 kinetic_connection_derivative=gradient_error < 2e-10,
                 physical_half_integer_labels=all((2*r['mu']) % 2 == 1 and .5 <= r['mu'] <= r['q']-.5 for r in rows),
                 sampled_mesh_refinement=mesh_error < 2e-5,
                 leading_scaling=all(1.6 < x['leading_mass_ratio'] < 2.5 for x in steps),
                 first_order_scaling=all(3. < x['first_order_mass_ratio'] < 5.5 for x in steps))
    result = dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                  parameters=p, coefficients=coeff, independent_checks=independent,
                  local_comparisons=rows, refined_comparisons=fine,
                  observed_scaling=steps, sampled_mesh_error=mesh_error,
                  kinetic_gradient_error=gradient_error,
                  equations={
                      'metric': 'ds2=-N2 dt2+P2(d_rho2+dz2)+rho2/N2 dphi2; v=N/P; b=N2/rho; a=Aphysical/q',
                      'canonical_spinor': 'u=P sqrt(rho/N) Psi; fixed-mu polar spin frame, flat d_rho dz norm',
                      'physical_H': 'H=-i sigma_x(v partial_rho+v_rho/2)-i sigma_z(v partial_z+v_z/2)+sigma_y b(mu-q a)',
                      'Maxwell': 'partial_rho(b a_rho)+partial_z(b a_z)=0; at grad(a)=0: a_rhorho+a_zz=0',
                      'a2': 'a2=a_rhorho/2=3 B r*(r*-2h) rho*2/[4 Q D0(r*)2 (r*-M)2], h=M(1-B ag/2)',
                      'scales': 'Lq=[v*/(q b* a2)]^(1/3); Eq=v*/Lq; Mq=Eq/b*; eta=(mu-q a*)/Mq',
                      'normal_form': 'H/Eq -> -i sigma_x partial_s-i sigma_z partial_t+sigma_y(eta-s2+t2)',
                      'first_order': 'v/v*=1+Lq (v_rho*/v*)s+O(Lq2); mass=eta-s2+t2+Lq[(b_rho*/b*)s(eta-s2+t2)-a_rhorhorho*s3/(6a2)-a_rhozz*s*t2/(2a2)]+O(Lq2)',
                      'widths': 'Delta_mu=O(q^(1/3)); Delta_(mu/q)=O(q^(-2/3)); Delta_rho,Delta_z=O(q^(-1/3)); Eq=O(q^(1/3)/D)'},
                  scope=[
                      'Fixed Q/D=.2 classical balanced dihole; changing q changes the coupling calibration 4piG/g2=(2Q/q)2.',
                      'Coefficient maxima are sampled on a bounded local box, not certified continuum suprema or global operator bounds.',
                      'Smooth Taylor expansion gives the local asymptotic orders for bounded eta and compactly supported test spinors; the finite numerical comparisons do not prove a scattering limit.',
                      'The axial t2 term is leading order; a one-dimensional radial barrier is not a two-dimensional return or escape probability.',
                      'No actual return matrix, quantum channel count, Casimir stress, stable medium, throat construction or backward travel is established.',
                      'At q=1 the local scale is comparable to separation, so this large-q normal form is not a controlled replacement there.'],
                  sources=[dict(url='https://arxiv.org/pdf/hep-th/9906160', attribution='Emparan, Black diholes, v3, equations 8--16: background only; normal-form calculation is ours.')],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                                for name in ('balanced_dihole.py', 'dihole_return_geometry.py', 'dihole_return_geometry_results.json', 'project_paths.py')})
    if not result['passed']:
        raise AssertionError(json.dumps(dict(tests=result['tests'], steps=steps, mesh=mesh_error), indent=2))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    c, rows = result['coefficients'], result['local_comparisons']
    colors = ['#d18b24', '#2e8e87', '#536ab3']
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), layout='constrained')
    ax = axes[0]
    for row, color in zip(rows, colors):
        L = row['radial_width']
        ax.add_patch(Rectangle((-L, -L), 2*L, 2*L, facecolor=color, alpha=.12,
                               edgecolor=color, linewidth=1.8))
        ax.plot([], [], color=color, label=f"q = {row['q']:,}")
    ax.axline((0, 0), slope=1, color='#858585', lw=1, ls='--')
    ax.axline((0, 0), slope=-1, color='#858585', lw=1, ls='--')
    ax.plot(0, 0, 'o', color='#343a45', ms=5)
    ax.set(xlim=(-.108, .108), ylim=(-.108, .108), aspect='equal',
           xlabel=r'$(\rho-\rho_*)/D$', ylabel=r'$z/D$', title='Local mixing region')
    ax.legend(loc='upper left', fontsize=8)
    ax.text(.02, -.28, r'$a-a_*\simeq a_2[(\rho-\rho_*)^2-z^2]$'+'\nDashed: tangent flux separatrix; squares: test boxes',
            transform=ax.transAxes, fontsize=9)
    ax = axes[1]
    for index, (row, color) in enumerate(zip(rows, colors)):
        width = row['normalized_label_width']
        ax.plot([c['a_star']-width, c['a_star']+width], [index, index], color=color, lw=8, alpha=.5)
        ax.plot(row['mu']/row['q'], index, 'o', color=color, ms=6)
    ax.axvline(c['a_star'], color='#555', ls='--', lw=1)
    ax.set(yticks=range(3), yticklabels=[f"q = {r['q']:,}" for r in rows],
           xlabel=r'Flux label $\mu/q$', title='Label band narrows faster')
    ax.invert_yaxis()
    ax.text(.02, -.28, r'Bars: $|\eta|\leq1$, width $\propto q^{-2/3}$'+'\nDots: actual nearest half-integer channels',
            transform=ax.transAxes, fontsize=9)
    ax = axes[2]
    q = [r['q'] for r in rows]
    leading = [r['leading_coefficient_errors']['angular_mass'] for r in rows]
    first = [r['first_order_remainder_errors']['angular_mass'] for r in rows]
    ax.loglog(q, leading, 'o-', color='#2e8e87', label='Leading normal form')
    ax.loglog(q, first, 's-', color='#536ab3', label='With first-order terms')
    ax.set(xlabel='Magnetic charge q', ylabel='Sampled angular-mass coefficient error',
           title='Actual two-dimensional coefficients')
    ax.legend(fontsize=8)
    ax.text(.02, -.28, 'Speed and connection errors checked separately.\nLocal coefficient check; no return probability.',
            transform=ax.transAxes, fontsize=9)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Where magnetic guiding breaks down | Local Dirac scaling at the null ring', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_null_scaling_progress.'+ext), dpi=175)
    plt.close(fig)


if __name__ == '__main__':
    result = checks()
    project_path('dihole_null_scaling_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result)
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'coefficients', 'observed_scaling', 'sampled_mesh_error')}, indent=2))
