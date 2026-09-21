"""Coupled scalar perturbations of the explicitly vacuum-supported action.

Natural-unit output uses Lambda=1. This is a homogeneous expanding exterior,
not a localized NEC source, manufactured material, or chronology-changing medium.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from project_paths import project_path


def symbolic_checks():
    """Start from a geometric ADM ansatz, not the claimed dispersion law."""
    import sympy as s
    e, z, v, g, d, b, shift, n = s.symbols('e z v g d b shift n', real=True)
    P, c, H, X, L, a = s.symbols('P c H X Lambda a', positive=True)
    rho = s.symbols('rho', real=True)
    lapse = 1+e*n
    # One-direction scalar mode; rotational symmetry supplies arbitrary k.
    A = H+e*v-e**2*b*g
    trace = 3*A-e*shift
    extrinsic = -6*A**2+4*A*e*shift
    curvature = s.exp(-2*e*z)/a**2*(-4*e*d-2*e**2*g**2)
    density = s.exp(3*e*z)*(P/2*(lapse*curvature+extrinsic/lapse)
        +lapse*(L**4*(X/lapse**2-1)**2/4-rho)-c*trace**2/lapse**3)
    quadratic = s.expand(s.diff(density, e, 2).subs(e, 0)/2)
    # Background: c=alpha*Lambda^2*X, with X held constant in time.
    background_rho = L**4*(X-1)**2/4+(3*P+9*c)*H**2
    quadratic = s.expand(quadratic.subs(rho, background_rho))
    # Spatial IBP: b*z_x -> -z*b_x, z*z_xx -> -z_x^2.
    # Time IBP with measure a^3: z*z_t -> -3H*z^2/2.
    integrated = s.expand(quadratic.subs(b*g, -z*shift).subs(z*d, -g**2)
                          .subs(z*v, -s.Rational(3, 2)*H*z**2))
    # The shift-current equation is L^4 X(X-1)=18 c H^2.
    E, B = P+3*c, P+9*c
    U = L**4*X*(3*X-1)/2
    raw = (-3*E*v**2+2*E*v*shift-c*shift**2+6*H*B*n*v
           -2*H*B*n*shift+(U-3*P*H**2-54*c*H**2)*n**2
           +P*g**2/a**2-2*P*n*d/a**2)
    residual = s.factor(integrated-raw)
    # Divisibility by the actual background relation avoids a circular rewrite.
    relation = L**4*X*(X-1)-18*c*H**2
    _, remainder = s.div(residual, relation, X)
    V = s.symbols('V', positive=True)
    # Substitute the current relation in the expanded polynomial.
    on_background = s.expand(raw).subs(L**4*X**2, V)
    on_background = s.expand(on_background.subs(L**4*X, V-18*c*H**2))
    constraints = s.solve([s.diff(on_background, n), s.diff(on_background, shift)], (n, shift))
    reduced = s.factor(on_background.subs(constraints))
    Z = c*V+H**2*(P+3*c)*(P+12*c)
    K = P*((P+3*c)*V+108*c**2*H**2)/Z
    J, W = P**2*H*(P+9*c)/Z, P**2*c/Z
    target = K*v**2-2*J*v*d/a**2+P*g**2/a**2-W*d**2/a**4
    G = P*c*(V+6*H**2*(P+6*c))/Z
    static_K = 3*P+P**2/c
    # Minisuperspace lapse/current/scale equations of the changed action.
    alpha = s.symbols('alpha', positive=True)
    x0 = 1+18*alpha*H**2/L**2
    vacuum = (3*P+9*alpha*L**2)*H**2+243*alpha**2*H**4
    px = (x0-1)/2
    f = (x0-1)**2/4-vacuum/L**4
    checks = dict(
        geometric_ADM_expansion=s.simplify(remainder) == 0,
        actual_background_lapse=s.simplify((3*P+27*alpha*L**2*x0)*H**2
                              -L**4*(2*x0*px-f)) == 0,
        actual_background_current=s.simplify(px-9*alpha*H**2/L**2) == 0,
        actual_background_scale=s.simplify((3*P+9*alpha*L**2*x0)*H**2+L**4*f) == 0,
        full_constraint_elimination=s.simplify(reduced-target) == 0,
        expanding_time_IBP=s.simplify(P-H*J-G) == 0,
        fixed_metric_decoupling=s.simplify(s.limit(K, P, s.oo)-V/H**2) == 0
                     and s.simplify(s.limit(G, P, s.oo)-6*c) == 0
                     and s.simplify(s.limit(W, P, s.oo)-c/H**2) == 0,
        static_limit=s.simplify(K.subs(H, 0)-static_K) == 0
                     and s.simplify(W.subs(H, 0)-P**2/V) == 0
                     and s.simplify(G.subs(H, 0)-P) == 0,
        zero_higher_derivative_limit=s.simplify(K.subs(c, 0)-V/H**2) == 0
                     and s.simplify(W.subs(c, 0)) == 0 and s.simplify(G.subs(c, 0)) == 0)
    return checks


def coefficients(cutoff_over_planck, h, alpha=1.):
    """All dimensionful coefficients are divided by their Lambda powers."""
    if (not all(math.isfinite(x) and x > 0 for x in (cutoff_over_planck, alpha))
            or not math.isfinite(h) or h < 0):
        raise ValueError('Require finite positive cutoff/Planck, alpha and finite H/Lambda>=0.')
    P, X = cutoff_over_planck**-2, 1+18*alpha*h*h
    c, V = alpha*X, X*X
    Z = c*V+h*h*(P+3*c)*(P+12*c)
    K = P*((P+3*c)*V+108*c*c*h*h)/Z
    W = P*P*c/Z
    G = P*c*(V+6*h*h*(P+6*c))/Z
    A, kJ2 = W/K, G/W
    gamma = math.sqrt(A)*kJ2/2
    vacuum = (3*P+9*alpha)*h*h+243*alpha*alpha*h**4
    return dict(cutoff_over_planck=cutoff_over_planck, H_over_cutoff=h, alpha=alpha,
        X=X, kinetic_over_cutoff2=K, gradient_over_cutoff2=G, fourth_derivative=W,
        dispersion_A_times_cutoff2=A, jeans_wavenumber_over_cutoff=math.sqrt(kJ2),
        instantaneous_peak_growth_over_cutoff=gamma,
        growth_over_expansion=gamma/h if h else None,
        added_vacuum_energy_over_cutoff4=vacuum,
        two_derivative_NEC_over_cutoff4=X*(X-1))


def evolve(r, *, rtol=2e-12, atol=2e-13, max_step=.12):
    """Curvature mode at band entry: z(0)=1,z'(0)=0, tau=Ht."""
    if not math.isfinite(r) or r < 0:
        raise ValueError('Require finite nonnegative growth/expansion ratio.')
    tau = np.linspace(0., 16., 161)
    def rhs(t, y):
        q = 4*r*r*(math.exp(-2*t)-math.exp(-4*t))
        return y[1], -3*y[1]+q*y[0]
    solution = solve_ivp(rhs, (0., tau[-1]), (1., 0.), method='DOP853',
                         t_eval=tau, rtol=rtol, atol=atol, max_step=max_step)
    if not solution.success:
        raise ArithmeticError(solution.message)
    z, velocity = solution.y[:, -1]
    estimate = float(z+velocity/3)
    return dict(growth_over_expansion=r, tau=tau.tolist(), mode=solution.y[0].tolist(),
        final_velocity=float(velocity), late_mode_estimate=estimate,
        analytic_lower=1+r*r/3, analytic_upper=math.exp(r*r/3),
        analytic_log_tail_factor=r*r/3*(2*math.exp(-32)-math.exp(-64)),
        scope='Numerical ODE values are not interval enclosures; analytic bounds apply to the stated band-entry data.')


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), layout='constrained')
    rows = result['scan']
    axes[0].loglog([x['H_over_cutoff'] for x in rows],
                  [math.expm1(x['growth_over_expansion']**2/3) for x in rows], color='#205e80')
    optimum = result['optima'][1]
    axes[0].axvline(optimum['H_over_cutoff'], ls='--', color='#b85432', label='Best tested expansion rate')
    axes[0].set(xlabel='Expansion rate H / cutoff', ylabel='Bound on fractional mode amplification',
                title='Derived coupled theory; cutoff / Planck = 0.01')
    axes[0].legend(fontsize=8)
    for row in result['evolutions']:
        axes[1].plot(row['tau'], row['mode'], label=f"H / cutoff = {row['H_over_cutoff']:.4g}")
    axes[1].set(xlim=(0., 6.), xlabel='Time in expansion units H t', ylabel='Curvature mode / initial amplitude',
                title='Growth saturates for the specified initial mode')
    axes[1].legend(fontsize=8)
    fig.suptitle('An expanding exterior passes a conditional finite-growth test', fontsize=13)
    fig.supxlabel('Added vacuum energy is required. No localized negative-energy source, nonlinear stability, material or time passage established.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('expanding_field_completion_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    checks = symbolic_checks()
    scan = [coefficients(.01, float(h)) for h in np.geomspace(1e-5, .03, 81)]
    optima = []
    for ratio in (.001, .01, .05):
        fit = minimize_scalar(lambda log_h: coefficients(ratio, math.exp(log_h))['growth_over_expansion'],
                              bounds=(math.log(1e-6), math.log(.05)), method='bounded',
                              options=dict(xatol=1e-12))
        if not fit.success:
            raise ArithmeticError(fit.message)
        optima.append(coefficients(ratio, math.exp(fit.x)))
    evolutions, differences = [], []
    for h in (1e-4, optima[1]['H_over_cutoff'], .03):
        row = coefficients(.01, h)
        first = evolve(row['growth_over_expansion'])
        fine = evolve(row['growth_over_expansion'], rtol=8e-14, atol=1e-14, max_step=.06)
        differences.append(abs(first['late_mode_estimate']-fine['late_mode_estimate']))
        evolutions.append(dict(first, H_over_cutoff=h))
    checks.update(positive_coupled_kinetic=all(r['kinetic_over_cutoff2'] > 0 for r in scan),
        residual_instantaneous_jeans_band=all(r['gradient_over_cutoff2'] > 0 and r['fourth_derivative'] > 0 for r in scan),
        mode_refinement=max(differences) < 2e-10,
        analytic_growth_bounds=all(r['analytic_lower'] <= r['late_mode_estimate'] <= r['analytic_upper'] for r in evolutions),
        no_instability_control=max(abs(x-1) for x in evolve(0.)['mode']) < 1e-14,
        interior_optima=all(1e-6 < r['H_over_cutoff'] < .05 for r in optima))
    result = dict(complete=True, passed=all(checks.values()), checks=checks, scan=scan,
        optima=optima, evolutions=evolutions, refinement_differences=differences,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={'project_paths.py': hashlib.sha256(project_path('project_paths.py').read_bytes()).hexdigest()},
        sources=['https://arxiv.org/abs/hep-th/0312099'],
        scope=['Exact background and quadratic scalar coefficients for the displayed truncated action with an explicit constant vacuum-energy term.',
               'The upper bound concerns a normalized curvature mode initially at band entry with zero time derivative; arbitrary initial data require their own normalization.',
               'Positive homogeneous expansion does not create a localized NEC-violating region: total background rho+pressure is zero.',
               'Cutoff-scale ratios are illustrative; strong-coupling scale, nonlinear stability, finite region joins, matter production and chronology change remain unproved.'])
    project_path('expanding_field_completion_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Expanding-completion checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True, checks=checks, optima=optima, refinement_differences=differences), indent=2))
    return result


if __name__ == '__main__':
    run()
