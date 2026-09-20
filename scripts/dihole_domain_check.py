"""Illustrate the optical ends and check the domain note against actual fields.

The global proof is in docs/DIHOLE_DIRAC_DOMAIN.md. These sampled checks are
sanity checks, not replacements for completeness or essential self-adjointness.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from balanced_dihole import parameters, symbolic_fields
from dihole_return_geometry import weyl_fields
from project_paths import project_path, source_matches


def run():
    p, function = parameters(), symbolic_fields()[0]
    epsilon = np.geomspace(1e-5, .1, 180)
    angles = np.array([0., np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    # The endpoint formula uses distances, avoiding loss of the small Delta.
    rho = epsilon[:, None]*np.sin(angles)
    z = p['k']+epsilon[:, None]*np.cos(angles)
    dplus, dminus = np.hypot(rho, z-p['k']), np.hypot(rho, z+p['k'])
    r, c = p['M']+(dplus+dminus)/2, 2*z/(dplus+dminus)
    sigma = r*r-p['a']**2*c*c
    alpha = sigma**2/(dplus*dminus)**2
    ratio = alpha*epsilon[:, None]**2/p['Q']**2

    # Independently compare the published metric implementation away from ends.
    rr = np.array([.02, .15, .5, 1., 3.])[:, None]
    zz = np.array([-.8, -.2, 0., .3, .9])[None, :]
    N, P, *_ = weyl_fields(function, p, rr, zz)
    dp, dm = np.hypot(rr, zz-p['k']), np.hypot(rr, zz+p['k'])
    r0, c0 = p['M']+(dp+dm)/2, 2*zz/(dp+dm)
    direct = (r0*r0-p['a']**2*c0*c0)**2/(dp*dm)**2
    relative_error = float(np.max(abs((P/N)/direct-1)))
    axis_z = np.array([-1.2, -.8, -.3, 0., .3, .8, 1.2])
    axis_N, axis_P, axis_a, *_ = weyl_fields(function, p, 0., axis_z)
    expected_a = (abs(axis_z)<p['k']).astype(float)
    checks = dict(
        published_optical_metric=relative_error < 3e-13,
        regular_balanced_axis=float(np.max(abs(axis_N*axis_P-1))) < 3e-13,
        regular_gauge_patches=float(np.max(abs(axis_a-expected_a))) < 3e-13,
        horizon_limit_all_sampled_directions=float(np.max(abs(ratio[0]-1))) < .0002,
        sampled_optical_lower_bound=bool(np.all(direct >= 1)),
    )
    assert all(checks.values()), checks
    calibration = json.loads(project_path('dihole_null_scaling_results.json').read_text())
    assert calibration['passed'] and source_matches('dihole_null_scaling.py', calibration['source_sha256'])
    coefficients = calibration['coefficients']
    q, outer = 512000, 9.
    length = coefficients['width_prefactor']*q**(-1/3)
    result = dict(passed=True, checks=checks, parameters=p,
        exact_target='D=1,Q=1/5,k=D/2; M,a,B derived exactly as in the domain note. Floats approximate this target.',
        optical_metric_relative_error=relative_error,
        axis_NP_error=float(np.max(abs(axis_N*axis_P-1))),
        axis_gauge_error=float(np.max(abs(axis_a-expected_a))),
        epsilon=epsilon.tolist(), angles=angles.tolist(),
        scaled_optical_factor=ratio.tolist(),
        probe=dict(q=q, outer=outer, length=length, rho_star=coefficients['rho_star']),
        scope='Sampled consistency checks accompanying an analytic domain proof; no quadrature certification or traversable throat.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
            for name in ('balanced_dihole.py', 'dihole_return_geometry.py', 'dihole_null_scaling_results.json')})
    project_path('dihole_domain_check_results.json').write_text(json.dumps(result, indent=2)+'\n')
    draw(result)
    print(json.dumps({key: result[key] for key in ('passed', 'checks', 'optical_metric_relative_error', 'axis_NP_error')}))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={'width_ratios': [1.05, 1]})
    k, Q = result['parameters']['k'], result['parameters']['Q']
    left.axvline(0., color='#247b62', lw=3, label='Regular axis included')
    left.scatter([0., 0.], [-k, k], marker='o', s=90, facecolors='white', edgecolors='#a94a31', lw=2, zorder=4)
    left.annotate('Horizon end\n(infinite optical distance)', (0, k), (.13, .82),
        arrowprops={'arrowstyle': '->', 'color': '#a94a31'}, color='#853922')
    left.text(.12, -.58, 'Second horizon end', color='#853922')
    p = result['probe']
    radius = np.sqrt(2)*p['length']*p['outer']
    center = p['rho_star']
    vertices = np.array([[center-radius, 0], [center, radius], [center+radius, 0], [center, -radius]])
    left.add_patch(Polygon(vertices, facecolor='#cfe3ff', edgecolor='#235b9f', lw=2))
    left.annotate('Compact quantum probe\naway from axis and horizons', (center, 0), (.37, -.38),
        arrowprops={'arrowstyle': '->', 'color': '#235b9f'}, color='#235b9f')
    left.text(.1, .12, 'Smooth exterior', color='#454b55')
    left.set(xlim=(-.07, 1.07), ylim=(-.9, 1.04), xlabel=r'$\rho/D$', ylabel=r'$z/D$',
        title='The domain is the uncut exterior')
    left.set_aspect('equal')
    left.legend(loc='lower right', fontsize=9)
    eps = np.array(result['epsilon'])
    ratios = np.array(result['scaled_optical_factor'])
    right.loglog(eps, ratios[:, 2]*Q**2/eps**2, color='#235b9f', lw=2, label=r'Actual $P/N$, equatorial approach')
    right.loglog(eps, Q**2/eps**2, '--', color='#a94a31', label=r'Horizon limit $Q^2/\epsilon^2$')
    right.set(xlabel=r'Distance to horizon end $\epsilon/D$', ylabel=r'Optical length factor $P/N$',
        title='A finite coordinate gap has infinite optical length')
    right.text(.04, .15, r'$\int_0^{\epsilon_0}(P/N)\,d\epsilon=\infty$', transform=right.transAxes, fontsize=15)
    right.grid(alpha=.15)
    right.legend(fontsize=9, loc='upper right')
    fig.suptitle('One well-defined quantum operator on the full exterior', fontsize=15, y=.99)
    fig.text(.5, .02, 'Domain proved for this known black-hole geometry. No traversable medium or backward passage established.', ha='center', fontsize=10)
    fig.tight_layout(rect=(0, .06, 1, .95))
    for extension in ('png', 'pdf'):
        fig.savefig(project_path('dihole_domain_progress.'+extension), dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()
