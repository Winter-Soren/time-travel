"""A smooth, finite-energy spherical packet in the coupled linear theory.

The spectral window is smooth and compact, not a reflecting spatial wall.
All quantities use Lambda=1; stress is the TOTAL linear Einstein source.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.special import expit, spherical_jn

from expanding_field_completion import coefficients
from project_paths import project_path


def spectrum(order=96):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    q = np.concatenate([4*(nodes+1), nodes+9])
    dq = np.concatenate([4*weights, weights])
    window = np.ones_like(q)
    u = (q[order:]-8)/2
    window[order:] = expit(1/u-1/(1-u))
    weights = math.sqrt(2/math.pi)*dq*q*q*np.exp(-q*q/2)*window
    return q, weights/weights.sum()


def packet(h, *, order=96, epsilon=.001, rtol=2e-11, atol=2e-13, step=.07):
    if not (math.isfinite(epsilon) and 0 < epsilon < .01):
        raise ValueError('This linear test requires 0<epsilon<0.01.')
    base = coefficients(.01, h)
    P = 1e4
    K, G, W = [base[key] for key in ('kinetic_over_cutoff2', 'gradient_over_cutoff2', 'fourth_derivative')]
    c, V = base['X'], base['X']**2  # alpha=Lambda=1 in this experiment.
    Z = c*V+h*h*(P+3*c)*(P+12*c)
    J = P*P*h*(P+9*c)/Z
    radius = 1/base['jeans_wavenumber_over_cutoff']
    q, weights = spectrum(order)
    k = q/radius
    ratio = base['growth_over_expansion']
    tau = np.linspace(0., 4., 81)
    def rhs(t, y):
        f, v = y.reshape(2, -1)
        x2 = q*q*math.exp(-2*t)
        return np.concatenate((v, -3*v-4*ratio*ratio*x2*(x2-1)*f))
    solution = solve_ivp(rhs, (0., 4.), np.concatenate((np.ones_like(q), np.zeros_like(q))),
        method='DOP853', t_eval=tau, rtol=rtol, atol=atol, max_step=step)
    if not solution.success:
        raise ArithmeticError(solution.message)
    f, ftau = solution.y.reshape(2, len(q), len(tau))
    f, ftau = f.T, ftau.T
    scale = np.exp(tau)[:, None]
    p2 = (k[None, :]/scale)**2
    z = -epsilon*weights[None, :]*f
    v = -epsilon*weights[None, :]*h*ftau
    acc = -3*h*v-(W/K)*p2*(p2-G/W)*z
    d, dt = -p2*z, -p2*v+2*h*p2*z
    n = (J*v+W*d)/P
    s = (K*v-J*d)/P
    ndot, sdot = (J*acc+W*dt)/P, (K*acc-J*dt)/P
    clock, clockdot = v-h*n, acc-h*ndot
    x = np.linspace(0., 10., 401)
    r = radius*x
    kr = k[:, None]*r[None, :]
    j0, j1 = spherical_jn(0, kr), spherical_jn(1, kr)
    transverse = np.divide(j1, kr, out=np.full_like(kr, 1/3), where=kr != 0)
    shear = j0-2*transverse
    B_full = (-2*clockdot@j0+2*sdot@transverse+h*s@(-2*j0+6*transverse)
              +2*(p2*z)@(j0-transverse)-2*(p2*n)@transverse)
    B_reduced = -2*(clockdot+h*s+d)@j0
    flux = (2*P*(k[None, :]/scale)*clock)@j1
    plus, minus = P*B_full+2*flux, P*B_full-2*flux
    rho = P*(-2*d+6*h*clock-2*h*s)@j0
    rho_independent = (2*h*(3*G-K)*v-2*(G+3*h*h*W)*d)@j0
    radial_gradient = z@(-k[:, None]*j1)
    shift_over_r = s@transverse
    compactness = 2*h*(scale*r[None, :])**2*(clock@j0-shift_over_r)-2*r[None, :]*radial_gradient
    conformal, lapse, velocity = z@j0, n@j0, v@j0
    # Hamiltonian mass check: excess mass includes energy FLUX across de Sitter spheres.
    effective_density = rho-h*scale*r[None, :]*flux
    mass_integral = np.array([simpson(row*r*r, x=r) for row in effective_density])
    endpoint_mass = P*r[-1]*compactness[:, -1]/scale[:, 0]**2
    absolute_mass = np.array([simpson(abs(row)*r*r, x=r) for row in effective_density])
    initial_gaussian = -epsilon*np.exp(-x*x/2)
    analytic_d0 = epsilon/radius**2*(3-x*x)*np.exp(-x*x/2)
    analytic_compactness = (-2*G/P*epsilon*x*x*np.exp(-x*x/2)
        -2*h*h*W/P*r*r*analytic_d0)
    initial_density_integral = float(simpson(rho[0]*r*r, x=r))
    initial_absolute_density = float(simpson(abs(rho[0])*r*r, x=r))
    full_scale = max(float(np.max(abs(P*B_full))), 1e-30)
    center_B = (-2*clockdot+2*sdot/3-4*d/3-2*p2*n/3).sum(axis=1)
    angular_equation = sdot+3*h*s+d-p2*n
    return dict(base=base, epsilon=epsilon, radius=radius, q=q, weights=weights,
        tau=tau, x=x, zeta=conformal, lapse=lapse, velocity=velocity,
        radial_gradient_times_r=r[None, :]*radial_gradient, shift_over_r=shift_over_r,
        null_plus=plus, null_minus=minus, flux=flux, density=rho,
        compactness_excess=compactness, mode_amplitude=f, mode_velocity=ftau,
        controls=dict(initial_gaussian_error=float(np.max(abs(conformal[0]-initial_gaussian))),
            initial_compactness_error=float(np.max(abs(compactness[0]-analytic_compactness))),
            full_reduced_null_relative_error=float(np.max(abs(P*(B_full-B_reduced)))/full_scale),
            center_projection_relative_error=float(np.max(abs(P*(center_B-B_full[:, 0])))/full_scale),
            density_identity_relative_error=float(np.max(abs(rho-rho_independent))/max(float(np.max(abs(rho))), 1e-30)),
            angular_Einstein_equation_max=float(np.max(abs(angular_equation))),
            center_flux_max=float(np.max(abs(flux[:, 0]))),
            integrated_mass_relative_error=float(np.max(abs(mass_integral-endpoint_mass)/absolute_mass)),
            initial_density_compensation_fraction=abs(initial_density_integral)/initial_absolute_density))


def plot(result, state):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.4), layout='constrained')
    x, tau = state['x'], state['tau']
    for key, label, style in [('null_plus', 'Outward null stress', '-'),
                              ('null_minus', 'Inward null stress', '--'), ('density', 'Energy density', ':')]:
        axes[0].plot(x, state[key][0]*1e5, style, label=label)
    axes[0].axhline(0, color='black', lw=.7)
    axes[0].set(xlim=(0., 5.), xlabel='Initial radius / packet width', ylabel='Linear stress / cutoff^4 (x 100,000)',
                title='Negative core and positive surroundings')
    axes[0].legend(fontsize=7)
    axes[1].plot(tau, state['null_plus'][:, 0]*1e5, color='#205e80')
    axes[1].set(xlabel='Expansion time H t', ylabel='Central null stress / cutoff^4 (x 100,000)',
                title='The source dilutes as the packet expands')
    for index in (0, 20, 60):
        axes[2].plot(x, state['compactness_excess'][index]*1e6, label=f'H t = {tau[index]:g}')
    axes[2].axhline(0, color='black', lw=.7)
    axes[2].set(xlim=(0., 6.), xlabel='Comoving radius / initial width', ylabel='Excess compactness (x 1,000,000)',
                title='Linear compactness only')
    axes[2].legend(fontsize=8)
    fig.suptitle('A localized spherical packet in the coupled linear theory', fontsize=13)
    fig.supxlabel('First-order packet: compactness needs nonlinear validation at this amplitude. No spatial throat, constructed material or backward passage.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('localized_expanding_packet_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    baseline_path = project_path('expanding_field_completion_results.json')
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if not baseline['passed'] or sha(project_path('expanding_field_completion.py')) != baseline['source_sha256']:
        raise ArithmeticError('The coupled-action baseline provenance failed.')
    h = baseline['optima'][1]['H_over_cutoff']
    state = packet(h)
    fine = packet(h, order=128, rtol=7e-13, atol=8e-15, step=.035)
    names = ('zeta', 'lapse', 'null_plus', 'null_minus', 'density', 'flux', 'compactness_excess')
    changes = {name: float(np.max(abs(state[name]-fine[name]))/max(np.max(abs(fine[name])), 1e-30)) for name in names}
    controls = state['controls']
    # For q<=10 and r_growth<=.025, energy decay before band entry bounds
    # |f| by14.2 at all later times; |x*j0'(x)|<=2 then excludes a neck.
    analytic_mode_bound = math.exp(.025**2/3)*(1+5*(math.log(10)+1/3))
    radial_derivative_bound = 2*state['epsilon']*analytic_mode_bound
    checks = dict(initial_transform=controls['initial_gaussian_error'] < 2e-13,
        initial_metric_check=controls['initial_compactness_error'] < 2e-12,
        independent_null_projection=controls['full_reduced_null_relative_error'] < 2e-10,
        regular_center=controls['center_projection_relative_error'] < 2e-10 and controls['center_flux_max'] == 0,
        density_constraint=controls['density_identity_relative_error'] < 2e-10,
        angular_Einstein_equation=controls['angular_Einstein_equation_max'] < 2e-14,
        flux_corrected_mass=controls['integrated_mass_relative_error'] < 2e-6,
        positive_compensation=bool(controls['initial_density_compensation_fraction'] < 2e-6
                            and np.min(state['density'][0]) < 0 < np.max(state['density'][0])),
        resolution=max(changes.values()) < 3e-8,
        negative_center=bool(np.all(state['null_plus'][:, 0] < 0) and np.all(state['null_minus'][:, 0] < 0)),
        small_metric=bool(np.max(abs(state['zeta'])) < .002 and np.max(abs(state['lapse'])) < .0001),
        no_spatial_neck_bound=state['base']['growth_over_expansion'] < .025 and state['q'].max() < 10
                              and analytic_mode_bound < 14.2 and radial_derivative_bound < 1)
    summary = dict(initial_central_null=float(state['null_plus'][0, 0]),
        final_central_null=float(state['null_plus'][-1, 0]),
        initial_central_density=float(state['density'][0, 0]),
        largest_absolute_compactness=float(np.max(abs(state['compactness_excess']))),
        largest_absolute_lapse=float(np.max(abs(state['lapse']))),
        largest_absolute_zeta=float(np.max(abs(state['zeta']))),
        largest_sampled_radial_derivative=float(np.max(abs(state['radial_gradient_times_r']))),
        all_time_mode_bound=analytic_mode_bound, all_time_radial_derivative_bound=radial_derivative_bound,
        initial_width_times_cutoff=state['radius'], width_growth_factor=math.exp(4),
        initial_width_over_background_horizon=h*state['radius'],
        width_reaches_background_horizon_at_Ht=math.log(1/(h*state['radius'])))
    arrays = {name:value for name,value in state.items() if isinstance(value,np.ndarray)}
    state_path = project_path('localized_expanding_packet_state.npz')
    np.savez_compressed(state_path, **arrays)
    result = dict(complete=True, passed=all(checks.values()), checks=checks,
        controls=controls, resolution_changes=changes, summary=summary,
        parameters=dict(epsilon=state['epsilon'], cutoff_over_planck=.01, H_over_cutoff=h,
                        spectral_window=[8.,10.], time_range_Ht=[0.,4.], quadrature_orders=[96,128]),
        source_sha256=sha(__file__), input_sha256={name:sha(project_path(name)) for name in ('expanding_field_completion.py','project_paths.py')},
        baseline_sha256=sha(baseline_path), state_sha256=sha(state_path),
        scope=['Smooth band-limited spherical perturbation of the derived homogeneous expanding background, at first order only.',
               'Total stress and flux follow from full linear Einstein projections; nonlinear stress and backreaction are not enclosed.',
               'Compensated energy and weak excess compactness do not create a spatial throat. The analytic derivative bound excludes a neck within this linear packet ansatz.',
               'The packet has rapidly decaying tails, not compact spatial support or physical walls. Its physical width expands; a cosmological horizon is not a wormhole mouth.',
               'Material production, actual cutoff, nonlinear stability, bridge support and chronology change remain unproved.'])
    project_path('localized_expanding_packet_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Localized packet checks failed; evidence retained.')
    plot(result, state)
    print(json.dumps(dict(passed=True, checks=checks, summary=summary, controls=controls, resolution_changes=changes), indent=2))
    return result


if __name__ == '__main__':
    run()
