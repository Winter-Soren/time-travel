"""A conditional localized-field hypothesis, not a material or new law.

The prescribed quadratic action has positive K,a, negative gradient
coefficient -sigma in a finite slab, and positive p_out outside. Its exact
linear-stability criterion extends to smoother compact blobs by comparison.
No covariant background, gravitational stability or physical source is supplied.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from project_paths import project_path


def critical_width(sigma, p_out, a):
    if not all(math.isfinite(x) and x > 0 for x in (sigma, p_out, a)):
        raise ValueError('sigma, p_out and a must be finite and positive.')
    return math.sqrt(a/sigma)*math.atan(math.sqrt(p_out/sigma))


def lowest_auxiliary_eigenvalue(width, sigma, p_out, a):
    """Ground level of T=-a*d_x^2+p, NOT a fourth-order growth rate."""
    critical_width(sigma, p_out, a)
    if not math.isfinite(width) or width <= 0:
        raise ValueError('width must be finite and positive.')
    def matching(energy):
        inside = math.sqrt(max(0., sigma+energy))
        outside = math.sqrt(max(0., p_out-energy))
        return width*inside/math.sqrt(a)-math.atan2(outside, inside)
    return brentq(matching, -sigma, p_out, xtol=1e-15, rtol=1e-14)


def transverse_identity(width=.5, sigma=1., p_out=1., a=1., q=.3):
    """Independent integrals of the full slab form and its T decomposition."""
    def parts(x):
        v = math.exp(-x*x/2)
        first, second = -x*v, (x*x-1)*v
        p = -sigma if abs(x) < width else p_out
        full = a*(second-q*q*v)**2+p*(first*first+q*q*v*v)
        decomposed = (a*second*second+p*first*first
                      +q*q*(2*a*first*first+p*v*v)+a*q**4*v*v)
        return full, decomposed
    bounds = (-np.inf, -width, width, np.inf)
    values = [sum(quad(lambda x: parts(x)[i], lo, hi, epsabs=1e-12, epsrel=1e-12)[0]
                  for lo, hi in zip(bounds[:-1], bounds[1:])) for i in (0, 1)]
    return dict(q=q, full=values[0], decomposed=values[1], error=abs(values[0]-values[1]))


def hardy_control(radius=1.):
    """Independent radial integrals for pi=exp(-r^2/2) in three dimensions."""
    def integral(function, hi=np.inf):
        return 4*math.pi*quad(function, 0., hi, epsabs=1e-11, epsrel=1e-12)[0]
    inside = integral(lambda r: r**4*math.exp(-r*r), radius)
    weighted = integral(lambda r: r*r*math.exp(-r*r))
    hessian = integral(lambda r: r*r*(r**4-2*r*r+3)*math.exp(-r*r))
    laplacian = integral(lambda r: r*r*(r*r-3)**2*math.exp(-r*r))
    return dict(inside=inside, weighted=weighted, hessian=hessian, laplacian=laplacian,
                satisfies_chain=inside <= radius**2*weighted <= 4*radius**2*hessian,
                hessian_laplacian_error=abs(hessian-laplacian))


def gravity_completion():
    """Eliminate the actual lapse/shift constraints of the proposed action."""
    import sympy as sp
    planck2, coefficient, cutoff, alpha, k = sp.symbols('P c Lambda alpha k', positive=True)
    velocity, laplacian, gradient2, shift, lapse = sp.symbols('v d g2 b n', real=True)
    action = (-3*planck2*velocity**2+2*planck2*velocity*shift+planck2*gradient2
              -2*planck2*lapse*laplacian+cutoff**4*lapse**2
              -coefficient*(3*velocity-shift)**2)
    constraints = sp.solve([sp.diff(action, lapse), sp.diff(action, shift)], (lapse, shift))
    reduced = sp.expand(action.subs(constraints))
    kinetic = 3*planck2+planck2**2/coefficient
    expected = kinetic*velocity**2+planck2*gradient2-planck2**2/cutoff**4*laplacian**2
    frequency = ((planck2**2/cutoff**4*k**4-planck2*k**2)/kinetic).subs(coefficient, alpha*cutoff**2)
    target = alpha/cutoff**2/(1+3*alpha*cutoff**2/planck2)*k**2*(k**2-cutoff**4/planck2)
    ratio = .01  # Illustrative Lambda/Mbar_Pl; not a measured material parameter.
    prefactor = 1/(1+3*ratio**2)
    def omega2_over_cutoff2(k_over_cutoff):
        return prefactor*k_over_cutoff**2*(k_over_cutoff**2-ratio**2)
    gamma = ratio**2/(2*math.sqrt(1+3*ratio**2))
    optimum = ratio/math.sqrt(2)
    return dict(constraint_identity=sp.simplify(reduced-expected) == 0,
        dispersion_identity=sp.simplify(frequency-target) == 0,
        unstable_band_checked=all(omega2_over_cutoff2(f*ratio) < 0 for f in (.1, .5, .9)),
        outside_band_positive=all(omega2_over_cutoff2(f*ratio) > 0 for f in (1.1, 2.)),
        growth_identity_error=abs(gamma**2+omega2_over_cutoff2(optimum)),
        illustrative_cutoff_over_reduced_planck=ratio, jeans_wavenumber_over_cutoff=ratio,
        fastest_growth_over_cutoff=gamma,
        decision='Reject this minimal action with an infinite static ghost-vacuum exterior as a globally linearly stable medium.',
        scope='Exact quadratic ADM result for the displayed truncated action; not a general exclusion of other theories or a prediction of immediate collapse.')


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), layout='constrained')
    rows = result['width_scan']
    axes[0].plot([r['width_over_critical'] for r in rows], [r['auxiliary_energy_over_sigma'] for r in rows])
    axes[0].axhline(0, color='black', lw=.8)
    axes[0].axvline(1, color='#b85432', ls='--', label='Exact onset in this model')
    axes[0].set(xlabel='Half-width / critical half-width', ylabel='Lowest auxiliary level / sigma',
                title='Localizing the negative coefficient can remove instability')
    axes[0].legend(fontsize=8)
    k = np.linspace(0., 1.35, 150)
    frequency = 4*k*k*(k*k-1)
    axes[1].plot(k, frequency, color='#b85432')
    axes[1].fill_between(k, frequency, 0., where=frequency < 0, color='#f0b8a3', label='Exponential growth')
    axes[1].axhline(0, color='black', lw=.8)
    axes[1].set(xlabel='Wavenumber / gravitational threshold', ylabel='Squared frequency / peak growth rate squared',
                title='Minimal completion: its exterior is unstable')
    axes[1].legend(fontsize=8)
    fig.suptitle('Reduced field model and its minimal gravitational completion', fontsize=13)
    fig.supxlabel('Left: prescribed coefficients. Right: gravity creates an exterior Jeans instability. No stable material or time passage established.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('localized_field_hypothesis_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    sigma, p_out, a = 1e-4, 1., 1.
    critical = critical_width(sigma, p_out, a)
    factors = np.linspace(.5, 2., 61)
    scan = [dict(width_over_critical=float(f),
                 auxiliary_energy_over_sigma=lowest_auxiliary_eigenvalue(f*critical, sigma, p_out, a)/sigma)
            for f in factors]
    at = lowest_auxiliary_eigenvalue(critical, sigma, p_out, a)
    below = lowest_auxiliary_eigenvalue(.99*critical, sigma, p_out, a)
    above = lowest_auxiliary_eigenvalue(1.01*critical, sigma, p_out, a)
    identities = [transverse_identity(q=q) for q in (0., .3, 1.2)]
    hardy = hardy_control()
    gravity = gravity_completion()
    zero_exterior = dict(sigma=sigma, alpha=1., K=1., outside_gradient_coefficient=0.,
        sufficient_radius_times_cutoff=.5/math.sqrt(sigma), example_radius_times_cutoff=40.,
        core_radius_times_cutoff=30., smoothing_length_times_cutoff=10.,
        retained_bending_margin=1-4*sigma*40**2,
        scope='Three-dimensional compact blob: sufficient Hardy comparison, not a necessary radius limit or a gravitational stability result.')
    example = dict(sigma=sigma, p_out=p_out, alpha=1., K=1., X=1.,
        critical_L_times_cutoff=critical, outer_L_times_cutoff=100., core_L_times_cutoff=80.,
        smoothing_length_times_cutoff=20.,
        auxiliary_energy=lowest_auxiliary_eigenvalue(100., sigma, p_out, a),
        conditional_NEC_magnitude_over_cutoff4=2*sigma,
        homogeneous_fastest_growth_over_cutoff=sigma/2,
        homogeneous_fastest_wavenumber_over_cutoff=math.sqrt(sigma/2))
    checks = dict(exact_symmetric_threshold=abs(critical_width(1., 1., 1.)-math.pi/4) < 1e-15,
        zero_at_critical=abs(at) < 1e-13, stable_side=below > 0, unstable_side=above < 0,
        decreasing_auxiliary_ground_level=all(b['auxiliary_energy_over_sigma'] < x['auxiliary_energy_over_sigma']
                                              for x, b in zip(scan[:-1], scan[1:])),
        all_transverse_identity=max(row['error'] for row in identities) < 2e-12,
        smooth_blob_comparison_margin=example['outer_L_times_cutoff'] < critical and example['auxiliary_energy'] > 0,
        three_dimensional_hardy_control=hardy['satisfies_chain'] and hardy['hessian_laplacian_error'] < 2e-11,
        zero_exterior_compact_window=zero_exterior['retained_bending_margin'] > 0,
        gravity_constraint_elimination=gravity['constraint_identity'],
        gravity_dispersion_identity=gravity['dispersion_identity'],
        gravity_jeans_band=gravity['unstable_band_checked'] and gravity['outside_band_positive'],
        gravity_growth_rate=gravity['growth_identity_error'] < 1e-20)
    result = dict(complete=True, passed=all(checks.values()), checks=checks, example=example,
        critical_controls=dict(at=at, below=below, above=above), transverse_controls=identities, width_scan=scan,
        zero_exterior_blob=zero_exterior, hardy_control=hardy, gravity_completion=gravity,
        conditional_cutoff_over_reduced_planck_for_B_alpha_X_1=math.sqrt(2)/math.pi,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=['https://arxiv.org/abs/hep-th/0312099', 'https://arxiv.org/abs/hep-th/0606090'],
        scope=['Analytic criterion for a prescribed decoupled quadratic action; no claim of a new physical law or novelty.',
               'The auxiliary T eigenvalue tests the quadratic form sign; it is not a physical fluctuation frequency or growth rate.',
               'Smooth compact blobs inherit nonnegativity only when their p is pointwise above the comparison slab and K,a remain positive with the stated values.',
               'The separate zero-exterior 3D Hardy criterion is sufficient only; it does not prove instability above that radius.',
               'The NEC and Planck comparison require an actual P(X) background, negligible omitted background stress, and a specified geometric source coefficient B.',
               'The displayed minimal covariant completion is rejected for strict global linear stability because its infinite static exterior has a Jeans band.',
               'No covariant completion, nonlinear/gravitational stability, causal UV completion, material generation, or backward passage established.'])
    project_path('localized_field_hypothesis_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Hypothesis checks failed; numerical evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'], checks=checks, example=example), indent=2))
    return result


if __name__ == '__main__':
    run()
