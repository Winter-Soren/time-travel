"""Nonlinear spherical initial constraints for the illustrative action144.

This solves one initial slice, not the evolution or material realizability.
Units: Lambda=1, a=1, beta=H*r*w(r/R), prescribed zeta=-epsilon*exp(-x*x/2).
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_bvp

from expanding_field_completion import coefficients
from project_paths import project_path


def background():
    prior_path = project_path('expanding_field_completion_results.json')
    prior = json.loads(prior_path.read_text(encoding='utf-8'))
    source = project_path('expanding_field_completion.py')
    if not prior['passed'] or hashlib.sha256(source.read_bytes()).hexdigest() != prior['source_sha256']:
        raise ArithmeticError('Background source provenance failed.')
    h = prior['optima'][1]['H_over_cutoff']
    base = coefficients(.01, h)
    P, c, X = 1e4, base['X'], base['X']
    Z = c*X*X+h*h*(P+3*c)*(P+12*c)
    return dict(P=P, c=c, X=X, H=h, radius=1/base['jeans_wavenumber_over_cutoff'],
                rho=base['added_vacuum_energy_over_cutoff4'],
                J=P*P*h*(P+9*c)/Z, W=base['fourth_derivative'], G=base['gradient_over_cutoff2'])


def profile(x, epsilon, radius):
    z = -epsilon*np.exp(-x*x/2)
    zx, zxx = -x*z, (x*x-1)*z
    u, ux = x*zx, x*(x*x-2)*z
    lap = (x*x-3)*z
    lapx = (5*x-x**3)*z
    curvature = np.exp(-2*z)*(-4*lap-2*zx*zx)/radius**2
    curvature_x = np.exp(-2*z)*(-4*lapx-4*zx*zxx+8*zx*lap+4*zx**3)/radius**2
    return z, u, ux, curvature, curvature_x


def fields(x, w, wx, epsilon, base):
    P, c, X, h = [base[key] for key in ('P', 'c', 'X', 'H')]
    z, u, ux, curvature, curvature_x = profile(x, epsilon, base['radius'])
    t, d = 1-w*(1+u), x*wx
    e = 3*t-d
    A, B = 3*c*h*h*e*e-.75*X*X, P*h*h*(3*t*t-2*t*d)+.5*X
    C = .5*P*curvature+.25-base['rho']
    discriminant = B*B-4*A*C
    if np.any(discriminant <= 0) or np.any(A >= 0):
        raise ArithmeticError('The lapse branch is not the continued background branch.')
    y = (-B-np.sqrt(discriminant))/(2*A)
    if np.any(y <= 0):
        raise ArithmeticError('Nonpositive inverse lapse squared.')
    N = 1/np.sqrt(y)
    Fy = 2*A*y+B
    Ft = 18*c*h*h*e*y*y+P*h*h*(6*t-2*d)*y
    Fd = -6*c*h*h*e*y*y-2*P*h*h*t*y
    D, K = P*t*N+3*c*e/N, N**3/(2*Fy)
    denominator = c+D*K*Fd
    if np.any(denominator <= 0):
        raise ArithmeticError('Momentum reduction became singular or changed branch.')
    numerator = (-(P*N*N+3*c)*w*ux-c*(4+3*u)*wx
                 -D*K*(Ft*(-(1+u)*wx-w*ux)+Fd*wx+.5*P*curvature_x))
    wxx = numerator/(x*denominator)
    Nx = K*(Ft*(-(1+u)*wx-w*ux)+Fd*(wx+x*wxx)+.5*P*curvature_x)
    return dict(zeta=z, u=u, ux=ux, curvature=curvature, t=t, d=d, e=e,
                lapse=N, lapse_x=Nx, wxx=wxx, discriminant=discriminant,
                implicit_denominator=denominator, lapse_polynomial=A*y*y+B*y+C)


def solve(epsilon=.001, *, xmax=12., xmin=1e-5, tolerance=2e-10):
    if not (0 <= epsilon <= .001 and 0 < xmin < .01 and xmax >= 8 and tolerance > 0):
        raise ValueError('Require 0<=epsilon<=.001, 0<xmin<.01, xmax>=8 and positive tolerance.')
    base = background()
    mesh = np.unique(np.r_[np.geomspace(xmin, .1, 60), np.linspace(.1, xmax, 550)])
    amplitude = -base['J']*epsilon/(base['P']*base['H']*base['radius']**2)
    guess = np.vstack((amplitude*np.exp(-mesh*mesh/2), -mesh*amplitude*np.exp(-mesh*mesh/2)))
    def rhs(x, values):
        w, wx = values
        return np.vstack((wx, fields(x, w, wx, epsilon, base)['wxx']))
    # For an even center w=w0+b*x^2/2+..., the singular equation gives
    # b=lim(rhs_wxx(x,w0,0))/5. Finite-x evaluation changes b by O(x^2).
    # The weak exterior is constant + x^-3; exclude the constant at infinity.
    def boundary(left, right):
        center_slope = fields(xmin, left[0], 0., epsilon, base)['wxx']/5
        return np.array([left[1]-xmin*center_slope, xmax*right[1]+3*right[0]])
    solution = solve_bvp(rhs, boundary, mesh, guess, tol=tolerance, max_nodes=20000)
    if not solution.success:
        raise ArithmeticError(solution.message)
    return solution, base


def evaluate(solution, base, epsilon, x):
    w, wx = solution.sol(x)
    f = fields(x, w, wx, epsilon, base)
    P, h, R = [base[key] for key in ('P', 'H', 'radius')]
    z, N, u, t, d = [f[key] for key in ('zeta', 'lapse', 'u', 't', 'd')]
    radius = R*x*np.exp(z)
    F = -h*(w*(1+u)+(N-1))/N
    compactness = radius**2*(2*h*F+F*F)-2*u-u*u
    delta_rho = .5*P*f['curvature']+P*h*h*((3*t*t-2*t*d)/N**2-3)
    flux = 2*P*h*np.exp(-z)/R*(w*f['ux']/N+t*f['lapse_x']/N**2)
    lap = epsilon/R**2*(3-x*x)*np.exp(-x*x/2)
    linear = -2*base['G']/P*u-2*h*h*base['W']/P*(R*x)**2*lap
    return dict(x=x, w=w, w_x=wx, zeta=z, lapse=N, density_excess=delta_rho,
                flux=flux, compactness_excess=compactness, linear_compactness=linear,
                areal_radius=radius, areal_derivative_factor=1+u, **{
                    key:f[key] for key in ('lapse_x','curvature','t','d','e','discriminant','implicit_denominator')})


def residuals(solution, base, epsilon):
    """Check the original equations and a separate geometric mass identity."""
    P, c, X, h = [base[key] for key in ('P', 'c', 'X', 'H')]
    x = np.linspace(.01, 10., 1801)
    f = evaluate(solution, base, epsilon, x)
    N, t, d, e = [f[key] for key in ('lapse','t','d','e')]
    local_X = X/N**2
    QE = h*h*(-6*t*t+4*t*d)
    lapse_terms = np.array([.5*P*(f['curvature']-QE/N**2),
        (local_X-1)**2/4-local_X*(local_X-1)-base['rho'], 3*c*h*h*e*e/N**4])
    lapse_relative = float(np.max(abs(lapse_terms.sum(axis=0)))/np.max(np.sum(abs(lapse_terms),axis=0)))
    u = profile(x, epsilon, base['radius'])[1]
    shift_term = P*(1+u)*f['w_x']/N
    def canonical(xx):
        w, wx = solution.sol(xx)
        v = fields(xx, w, wx, epsilon, base)
        return P*(v['t']/v['lapse']-1)+c*(v['e']/v['lapse']**3-3)
    momentum = []
    matter_flux = []
    for step in (.002, .001):
        derivative = (canonical(x-2*step)-8*canonical(x-step)
                      +8*canonical(x+step)-canonical(x+2*step))/(12*step)
        momentum.append(float(np.max(abs(derivative+shift_term))/max(np.max(abs(shift_term)), 1e-25)))
        def clock_trace(xx):
            ww, vv = solution.sol(xx)
            ff = fields(xx, ww, vv, epsilon, base)
            return ff['e']/ff['lapse']**3
        trace_x = (clock_trace(x-2*step)-8*clock_trace(x-step)
                   +8*clock_trace(x+step)-clock_trace(x+2*step))/(12*step)
        flux_from_action = 2*c*h*np.exp(-f['zeta'])/base['radius']*trace_x
        matter_flux.append(float(np.max(abs(flux_from_action-f['flux']))/np.max(abs(f['flux']))))
    density_from_action = (base['rho']+local_X*(local_X-1)-(local_X-1)**2/4
                           -3*c*h*h*e*e/N**4-3*P*h*h)
    matter_density = float(np.max(abs(density_from_action-f['density_excess']))/np.max(abs(f['density_excess'])))
    xx = np.linspace(solution.x[0], 10., 8001)
    v = evaluate(solution, base, epsilon, xx)
    radius_x = base['radius']*np.exp(v['zeta'])*v['areal_derivative_factor']
    integrand = v['areal_radius']**2*(v['density_excess']*radius_x
        -v['flux']*base['radius']*np.exp(v['zeta'])*v['areal_radius']*h*v['t']/v['lapse'])
    boundary_mass = P*(v['areal_radius'][-1]*v['compactness_excess'][-1]
                       -v['areal_radius'][0]*v['compactness_excess'][0])
    mass_relative = float(abs(simpson(integrand,x=xx)-boundary_mass)/simpson(abs(integrand),x=xx))
    return dict(original_lapse_relative=lapse_relative,
                original_momentum_relative_five_point=momentum,
                action_density_relative=matter_density, action_flux_relative_five_point=matter_flux,
                geometric_mass_integral_relative=mass_relative)


def plot(state):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12., 4.3), layout='constrained')
    x = state['x']
    axes[0].plot(x, state['density_excess']*1e6, color='#205e80')
    axes[0].axhline(0, color='black', lw=.7)
    axes[0].set(xlim=(0,5), xlabel='Radius / initial packet width',
                ylabel='Density excess / cutoff^4 (x 1,000,000)', title='Constraint-compatible density')
    axes[1].plot(x, state['compactness_excess']*1e6, label='Nonlinear initial constraints')
    axes[1].plot(x, state['linear_compactness']*1e6, '--', label='Earlier linear estimate')
    axes[1].set(xlim=(0,6), xlabel='Radius / initial packet width',
                ylabel='Excess compactness (x 1,000,000)', title='The gravitational estimate changes')
    axes[1].legend(fontsize=7)
    axes[2].plot(x, state['areal_derivative_factor'], color='#b85432')
    axes[2].axhline(0, color='black', ls='--', lw=.8, label='Necessary for a spatial neck')
    axes[2].set(xlim=(0,5), ylim=(-.05,1.1), xlabel='Radius / initial packet width',
                ylabel='Areal-radius derivative / exp(zeta)', title='No spatial neck on this slice')
    axes[2].legend(fontsize=7)
    fig.suptitle('A consistent initial slice for a hypothetical field', fontsize=13)
    fig.supxlabel('Initial constraints only: nonlinear evolution, null stress, stability, material realization and backward passage remain unproved.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('nonlinear_packet_constraints_progress.'+extension),dpi=160)
    plt.close(fig)


def run():
    epsilon = .001
    coarse, base = solve(epsilon)
    fine, _ = solve(epsilon, tolerance=1e-11)
    outer, _ = solve(epsilon, xmax=18., tolerance=1e-11)
    center, _ = solve(epsilon, xmin=1e-6, tolerance=1e-11)
    x = np.linspace(1e-5, 10., 2001)
    state = evaluate(fine, base, epsilon, x)
    changes = {}
    for name, solution in [('resolution',coarse),('outer_radius',outer),('center_cutoff',center)]:
        other = evaluate(solution,base,epsilon,x)
        changes[name] = {key:float(np.max(abs(other[key]-state[key])))
                        for key in ('w','lapse','density_excess','compactness_excess')}
    controls = residuals(fine,base,epsilon)
    controls['maximum_bvp_residual'] = float(np.max(fine.rms_residuals))
    amplitudes = []
    for amplitude in (epsilon, .0001, .00001):
        solution = fine if amplitude == epsilon else solve(amplitude,tolerance=1e-11)[0]
        v = evaluate(solution,base,amplitude,x)
        error = float(np.max(abs(v['compactness_excess']-v['linear_compactness'])))
        amplitudes.append(dict(epsilon=amplitude, largest_absolute_compactness=float(np.max(abs(v['compactness_excess']))),
            largest_linear_compactness=float(np.max(abs(v['linear_compactness']))),
            largest_difference=error, difference_over_maximum_linear=error/float(np.max(abs(v['linear_compactness']))),
            density_near_center=float(v['density_excess'][0])))
    zero, _ = solve(0.)
    zero_fields = evaluate(zero,base,0.,x)
    controls['zero_amplitude_maximum_compactness'] = float(np.max(abs(zero_fields['compactness_excess'])))
    controls['zero_amplitude_maximum_density_excess'] = float(np.max(abs(zero_fields['density_excess'])))
    correction_ratios = [amplitudes[i]['largest_difference']/amplitudes[i+1]['largest_difference'] for i in range(2)]
    summary = dict(density_near_center=float(state['density_excess'][0]),
        maximum_density_excess=float(np.max(state['density_excess'])),
        maximum_absolute_compactness=float(np.max(abs(state['compactness_excess']))),
        difference_from_linear_over_maximum_linear=amplitudes[0]['difference_over_maximum_linear'],
        maximum_lapse_deviation=float(np.max(abs(state['lapse']-1))),
        minimum_lapse=float(np.min(state['lapse'])),
        minimum_discriminant=float(np.min(state['discriminant'])),
        minimum_implicit_denominator=float(np.min(state['implicit_denominator'])),
        minimum_areal_derivative_factor=float(np.min(state['areal_derivative_factor'])),
        near_center_shift_over_Hr=float(state['w'][0]),
        exterior_shift_tail_x3w=float(state['w'][-1]*x[-1]**3))
    checks = dict(original_lapse=controls['original_lapse_relative'] < 3e-14,
        independent_momentum=max(controls['original_momentum_relative_five_point']) < 2e-7,
        action_density=controls['action_density_relative'] < 2e-9,
        action_flux=max(controls['action_flux_relative_five_point']) < 2e-6,
        independent_mass_identity=controls['geometric_mass_integral_relative'] < 2e-7,
        bvp_converged=controls['maximum_bvp_residual'] < 1.01e-11,
        regular_branch=summary['minimum_lapse'] > .99 and summary['minimum_discriminant'] > .2
                       and summary['minimum_implicit_denominator'] > 3000,
        resolution=changes['resolution']['compactness_excess'] < 3e-11,
        exterior_refinement=changes['outer_radius']['compactness_excess'] < 3e-11,
        center_refinement=changes['center_cutoff']['compactness_excess'] < 3e-11,
        homogeneous_control=controls['zero_amplitude_maximum_compactness'] < 1e-12
                            and controls['zero_amplitude_maximum_density_excess'] < 1e-13,
        nonlinear_correction_is_quadratic=all(95 < r < 105 for r in correction_ratios),
        linear_limit=amplitudes[-1]['difference_over_maximum_linear'] < .004,
        no_initial_neck=summary['minimum_areal_derivative_factor'] >= 1)
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    state_path = project_path('nonlinear_packet_constraints_state.npz')
    np.savez_compressed(state_path, **state, bvp_x=fine.x, bvp_y=fine.y)
    result = dict(complete=True, passed=all(checks.values()),checks=checks,summary=summary,
        controls=controls, absolute_refinement_changes=changes, amplitude_controls=amplitudes,
        quadratic_correction_ratios=correction_ratios, background=base,
        parameters=dict(epsilon=epsilon, zetadot=0, scale_factor=1, profile='Exact Gaussian, no spectral truncation',
                        outer_radii=[12.,18.],inner_cutoffs=[1e-5,1e-6],tolerances=[2e-10,1e-11]),
        source_sha256=sha(__file__),state_sha256=sha(state_path),
        inputs_sha256={name:sha(project_path(name)) for name in
            ('expanding_field_completion.py','expanding_field_completion_results.json','project_paths.py')},
        scope=['Numerically consistent nonlinear initial constraints of the stated hypothetical action, not a proof of continuum existence or uniqueness.',
               'A regular center and a decaying asymptotic shift are imposed through tested finite-boundary approximations.',
               'Negative density relative to the expanding background is not negative total density, null-energy violation, or evidence for a throat.',
               'Prescribed zeta=-epsilon exp(-x^2/2) implies 1+r*zeta_r>=1 exactly for epsilon>=0; this slice cannot have a spatial neck.',
               'No nonlinear null stress, full evolution, long-term stability, UV completion, material recipe or chronology change has been established.'])
    project_path('nonlinear_packet_constraints_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Nonlinear initial-data checks failed; saved evidence retained.')
    plot(state)
    print(json.dumps({key:result[key] for key in ('passed','checks','summary','controls','absolute_refinement_changes','amplitude_controls')},indent=2))
    return result


if __name__ == '__main__':
    run()
