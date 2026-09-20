"""Classical return topology, a frozen-slice gap, and a conditional energy gate.

The exact exterior is Emparan's black dihole. Field lines are not quantum
trajectories. The ideal channel energy is not the actual Melvin vacuum energy.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from balanced_dihole import parameters, symbolic_fields, cut_profile
from project_paths import project_path, source_matches


def equatorial(r, p):
    m, a, b, q = (p[key] for key in ('M', 'a', 'B', 'Q'))
    length = r**3-a*a*r-2*m*a*a
    denominator = r-2*m+2*b*m*a+b*b*length/4
    return (2*m*a+b*length/2)/(2*q*denominator)


def return_fraction(p):
    m, a, b = (p[key] for key in ('M', 'a', 'B'))
    h = m*(1-b*a/2)
    polynomial = lambda r: r**3-3*h*r*r+(m+h)*a*a-2*m*a/b
    hi = 2*p['rplus']
    while polynomial(hi) <= 0:
        hi *= 2
    r = brentq(polynomial, p['rplus'], hi, xtol=5e-15)
    potential = equatorial(r, p)
    return dict(Q_over_D=p['Q']/p['D'], critical_r=r,
                null_ring_rho=math.sqrt((r-m)**2-p['k']**2),
                minimum_A_over_q=potential, returning_flux_fraction=1-potential)


def weyl_fields(function, p, rho, z):
    """Actual N, P, physical A and its Weyl-coordinate gradient."""
    rho, z = np.broadcast_arrays(np.asarray(rho, float), np.asarray(z, float))
    plus, minus = np.hypot(rho, z+p['k']), np.hypot(rho, z-p['k'])
    total = plus+minus
    r, c = p['M']+total/2, 2*z/total
    raw = function(r, c, p['M'], p['a'], p['B'])
    n, meridian, ahat, _, _, ar, _, _, ac = [np.broadcast_to(v, rho.shape) for v in raw]
    sr, sz = rho/plus+rho/minus, (z+p['k'])/plus+(z-p['k'])/minus
    factor = p['physical_A_factor']
    gradient_rho = factor*(ar*sr/2-ac*2*z*sr/total**2)
    gradient_z = factor*(ar*sz/2+ac*(2/total-2*z*sz/total**2))
    return n, meridian, ahat*factor, gradient_rho, gradient_z


def trace_line(function, p, level, tolerance=2e-10):
    """Unit Euclidean arclength; optical time integrates P/N on this curve."""
    epsilon = .05*p['D']
    cosine = brentq(lambda u: float(cut_profile(function, p, epsilon, u)['Aphysical'])-level,
                    -.999999, .999999, xtol=5e-15)
    start = [epsilon*math.sqrt(1-cosine*cosine), p['k']+epsilon*cosine, 0.]

    def equation(_, state):
        n, meridian, _, ar, az = weyl_fields(function, p, state[0], state[1])
        norm = float(np.hypot(ar, az))
        if norm < 1e-13:
            raise ValueError('Trace reached magnetic null: guiding description is singular.')
        return [-float(az)/norm, float(ar)/norm, float(meridian/n)]

    def lower(_, state):
        return math.hypot(state[0], state[1]+p['k'])-epsilon

    def infinity_proxy(_, state):
        return state[1]-4*p['D']

    lower.terminal = infinity_proxy.terminal = True
    lower.direction, infinity_proxy.direction = -1, 1
    solution = solve_ivp(equation, [0, 60*p['D']], start, method='DOP853',
                         rtol=tolerance, atol=tolerance/10, dense_output=True,
                         max_step=.03*p['D'], events=(lower, infinity_proxy))
    if not solution.success or solution.status != 1:
        raise RuntimeError(solution.message)
    samples = solution.sol(np.linspace(0, solution.t[-1], 201))
    potentials = weyl_fields(function, p, samples[0], samples[1])[2]
    return dict(A_over_q=level, endpoint='lower_cut' if len(solution.t_events[0]) else 'upper_plane_z4',
                auxiliary_arclength=float(solution.t[-1]), optical_time=float(solution.y[2, -1]),
                maximum_flux_drift=float(np.max(abs(potentials-level))),
                endpoint_coordinates=solution.y[:2, -1].tolist(), points=samples[:2].T.tolist())


def midplane_gap_certificate():
    """70-digit interval/Descartes proof: inf N^2(A-1/2)/rho > .2138/D."""
    iv, old = mp.iv, mp.iv.dps
    iv.dps = 70
    try:
        q, k = iv.mpf(1)/5, iv.mpf(1)/2
        m = (iv.sqrt(k*k+4*q*k)-k)/2
        a = iv.sqrt(k*k-m*m)
        rp = m+k
        b = m/((k+a)*rp)
        beta = b*b/4
        length = [4*m*m*rp, 3*rp*rp-a*a, 3*rp, iv.mpf(1)]
        denominator = [rp-2*m+2*b*m*a+beta*length[0], 1+beta*length[1], beta*length[2], beta]
        numerator = [m*a/q+b*length[0]/(4*q)-denominator[0]/2]
        numerator += [b*length[j]/(4*q)-denominator[j]/2 for j in range(1, 4)]

        def multiply(x, y):
            out = [iv.mpf(0) for _ in range(len(x)+len(y)-1)]
            for i, xx in enumerate(x):
                for j, yy in enumerate(y):
                    out[i+j] += xx*yy
            return out

        def evaluate(poly, x):
            out = iv.mpf(0)
            for coefficient in reversed(poly):
                out = out*x+coefficient
            return out

        product = multiply(denominator, numerator)
        assert all(coefficient.a > 0 for coefficient in product)
        derivative = [(j+1)*product[j+1] for j in range(6)]
        left = multiply([0, 4*k*rp, 2*rp+4*k, 2], derivative)
        right = multiply([2*k*rp, 2*rp+10*k, 6], product)
        stationary = [u-v for u, v in zip(left, right)]
        assert len(stationary) == 9
        assert all(c.b < 0 for c in stationary[:3]) and all(c.a > 0 for c in stationary[3:])
        lo, hi = '0.4502557132414467', '0.4502557132414469'
        assert evaluate(stationary, iv.mpf(lo)).b < 0 < evaluate(stationary, iv.mpf(hi)).a
        x = iv.mpf([lo, hi])
        u = evaluate(product, x)/((rp+x)**2*iv.sqrt(x*(x+2*k)))
        assert u.a > iv.mpf('.2138')
        return dict(passed=True, digits=70, positive_stationary_roots=1,
                    shifted_r_bracket=[lo, hi], minimum_U_interval=str(u),
                    certified_gap_lower='.2138/D', stationary_coefficients=[str(c) for c in stationary])
    finally:
        iv.dps = old


def ideal_energy_check():
    """Finite integer-channel control of the analytically derived quarter rule."""
    classical, q = 1., 100
    rows = []
    for count in (22, 25, 30):
        d = np.linspace(0., 2., count)  # L_i/pi; illustrative, not dihole data.
        scaled_derivative = lambda ell: -2*classical/ell-q/24+np.sum((ell/(ell+d))**2)/6
        c = (4*count-q)/24
        row = dict(q=q, closed_channels=count, asymptotic_scaled_derivative=c)
        if c > 0:
            lo, hi = 2*classical/c, (2*classical+np.sum(d)/3)/c
            root = brentq(scaled_derivative, lo, hi)
            energy = classical/root**2+q/(24*root)-np.sum(1/(root+d))/6
            identity = -classical/root**2-np.sum(d/(root+d)**2)/6
            assert abs(energy-identity) < 2e-14 and energy < 0
            row.update(minimum_ell=root, minimum_energy=energy, bracket=[lo, hi])
        else:
            assert np.max([scaled_derivative(x) for x in np.geomspace(.01, 1e6, 70)]) < 0
        rows.append(row)
    return rows


def draw(result, function, p):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.8), layout='constrained')
    ax = axes[0]
    rho, z = np.meshgrid(np.linspace(.004, 1.4, 330), np.linspace(-1.25, 1.25, 470))
    potential = weyl_fields(function, p, rho, z)[2]
    mask = (np.hypot(rho, z-p['k']) < .05) | (np.hypot(rho, z+p['k']) < .05)
    potential = np.ma.masked_where(mask, potential)
    ax.contour(rho, z, potential, levels=[.2, .5, .75, .8, .9, .98, 1.2], colors='#c0c7cd', linewidths=.9)
    critical = result['geometry']['minimum_A_over_q']
    ax.contour(rho, z, potential, levels=[critical], colors='#c47f00', linewidths=1.5, linestyles='--')
    for line in result['field_lines']:
        points = np.array(line['points'])
        color = '#008776' if line['endpoint'] == 'lower_cut' else '#497db1'
        ax.plot(points[:, 0], points[:, 1], color=color, lw=1.5)
    for height in (-p['k'], p['k']):
        ax.add_patch(Circle((0, height), .05, facecolor='#353b45', zorder=5))
    ax.plot(result['geometry']['null_ring_rho'], 0, 'o', color='#c47f00', ms=6)
    ax.text(.04, .59, 'Upper cut', fontsize=9)
    ax.text(.04, -.65, 'Lower cut', fontsize=9)
    ax.text(.76, .07, 'Magnetic\nnull ring', fontsize=9, color='#986200')
    ax.set(xlim=(0, 1.4), ylim=(-1.25, 1.25), xlabel=r'$\rho/D$', ylabel=r'$z/D$',
           title='Actual classical field lines')
    ax.set_aspect('equal')
    ax.text(.02, -.23, 'Green: joins cuts   Blue: outward\nDashed: separatrix; no passage is sewn',
            transform=ax.transAxes, fontsize=9)

    ax = axes[1]
    radial = np.geomspace(.02, 6., 400)
    n, _, gauge, _, _ = weyl_fields(function, p, radial, 0.)
    ax.plot(radial, n*n*(gauge-.5)/radial, color='#9159a3', lw=2)
    ax.axhline(.2138, color='#9159a3', ls='--', lw=1)
    ax.text(.55, .27, r'Certified $|E_\perp|>0.2138\,\hbar c/D$', fontsize=9)
    ax.set(xlim=(0, 3), ylim=(0, 1.), xlabel=r'$\rho/D$', ylabel=r'$U D$',
           title=r'Central slice: $q=1,\ \mu=1/2$')
    ax.text(.02, -.23, 'A transverse barrier, not a bound\non full two-dimensional transmission',
            transform=ax.transAxes, fontsize=9)

    ax = axes[2]
    family = result['geometry_family']
    ax.plot([r['Q_over_D'] for r in family], [100*r['returning_flux_fraction'] for r in family], color='#008776', lw=2)
    ax.axhline(25, color='#c47f00', ls='--', label='Conditional channel threshold')
    ax.plot(.2, 100*(1-critical), 'o', color='#253346')
    ax.annotate('Current geometry: 22.1%', (.2, 100*(1-critical)), (.08, 31),
                arrowprops=dict(arrowstyle='->', color='#253346'), fontsize=9)
    ax.set(xlim=(0, 1), ylim=(16, 39), xlabel=r'Charge radius / Weyl separation $Q/D$',
           ylabel='Classical returning flux (%)', title='Changing the geometry')
    ax.legend(loc='lower right', fontsize=8)
    ax.text(.02, -.23, 'Flux is not a quantum channel count.\nCrossing 25% does not prove support.',
            transform=ax.transAxes, fontsize=9)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Return paths and the support test | No stable medium demonstrated', fontsize=15)
    for extension in ('png', 'pdf'):
        fig.savefig(project_path('dihole_return_progress.'+extension), dpi=170)
    plt.close(fig)


def run():
    saved_path = project_path('balanced_dihole_results.json')
    saved = json.loads(saved_path.read_text(encoding='utf-8'))
    assert saved['passed'] and source_matches('balanced_dihole.py', saved['source_sha256'])
    p, function = saved['parameters'], symbolic_fields()[0]
    geometry = return_fraction(p)
    lines = [trace_line(function, p, level) for level in (.2, .5, .75, .8, .9, .98)]
    refined = [trace_line(function, p, level, 2e-12) for level in (.75, .8)]
    reference = [lines[2], lines[3]]
    refinement = max(abs(a['optical_time']-b['optical_time']) for a, b in zip(reference, refined))
    samples = np.array([[.15, .7], [.5, .1], [1.2, -.8]])
    step = 2e-5
    derivative_error = 0.
    for rho, z in samples:
        values = weyl_fields(function, p, rho, z)
        for direction in (0, 1):
            dr, dz = (step, 0.) if direction == 0 else (0., step)
            near = [weyl_fields(function, p, rho+j*dr, z+j*dz)[2] for j in (-2, -1, 1, 2)]
            finite = (near[0]-8*near[1]+8*near[2]-near[3])/(12*step)
            derivative_error = max(derivative_error, float(abs(finite-values[3+direction])))
    rr = np.geomspace(.02, 6., 40)
    fields = weyl_fields(function, p, rr, 0.)
    eq_error = float(np.max(abs(fields[2]-equatorial(p['M']+np.hypot(rr, p['k']), p))))
    null = weyl_fields(function, p, geometry['null_ring_rho'], 0.)
    family = [return_fraction(parameters(x)) for x in np.linspace(.01, .99, 99)]
    crossing = brentq(lambda x: return_fraction(parameters(x))['returning_flux_fraction']-.25, .3, .5)
    gap, energy = midplane_gap_certificate(), ideal_energy_check()
    tests = dict(equatorial_reduction=eq_error < 2e-12,
                 Weyl_gradient=derivative_error < 1e-8,
                 magnetic_null=float(np.hypot(null[3], null[4])) < 2e-12,
                 conserved_flux=max(x['maximum_flux_drift'] for x in lines+refined) < 2e-8,
                 endpoint_topology=all((x['endpoint']=='lower_cut') == (x['A_over_q']>geometry['minimum_A_over_q']) for x in lines),
                 optical_time_refinement=refinement < 2e-7,
                 interval_midplane_gap=gap['passed'], ideal_energy_identity=True,
                 sampled_geometry_crossing=abs(return_fraction(parameters(crossing))['returning_flux_fraction']-.25) < 2e-13)
    result = dict(passed=all(tests.values()), tests=tests, parameters=p, geometry=geometry,
                  field_lines=lines, refined_lines=refined, optical_time_refinement=refinement,
                  derivative_error=derivative_error, equatorial_error=eq_error,
                  midplane_gap=gap, ideal_energy_controls=energy, geometry_family=family,
                  quarter_crossing_Q_over_D=crossing,
                  ideal_energy_gate=dict(criterion='f>1/4 iff one finite global minimum, fixed finite lengths and channels',
                    baseline_positive_coefficient=(1-4*geometry['returning_flux_fraction'])/24,
                    weak_interaction_return_fraction=(9-3*4**(2/3))/8),
                  scope=['Classical flux fraction is not a quantum return probability or finite-q channel fraction.',
                    'Frozen central-slice gap is not a full two-dimensional tunneling/transmission bound.',
                    'Field-line optical times use finite cuts; they are not measured quantum delays.',
                    'Quarter rule assumes fixed closed NS channels, transparent open vacua, all-channel anomaly, and fixed geometry/clock matching.',
                    'Actual scattering, Melvin vacuum stress, backreaction and stability are uncomputed.',
                    'Magnetic guiding fails at the null ring; large-q control is still needed.',
                    'The crossing near Q/D=.389 is outside a parametrically well-separated limit; it does not validate the energy ansatz.',
                    'No physical source, obtainable matter, sewn throat or backward travel is demonstrated.'],
                  sources=['https://arxiv.org/abs/hep-th/9906160', 'https://arxiv.org/html/1807.04726v3'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                                for name in ('balanced_dihole.py', 'balanced_dihole_results.json', 'project_paths.py')})
    if not result['passed']:
        raise AssertionError(tests)
    project_path('dihole_return_geometry_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result, function, p)
    print(json.dumps({k: result[k] for k in ('passed', 'tests', 'geometry', 'optical_time_refinement', 'derivative_error', 'quarter_crossing_Q_over_D')}, indent=2))
    return result


if __name__ == '__main__':
    run()
