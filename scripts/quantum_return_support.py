"""Conditional return-force bound and a phase-sensitive support decision.

This reuses the existing loop/reservoir normalization. It does not compute
the balanced-dihole return operator or identify classical flux with channels.
The matrix inequality is analytic; the controls and sensitivity values here
are floating diagnostics, not new interval certificates or a physical throat.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from project_paths import project_path
from scattering_vacuum import LoopReservoir


def support_weight(modulus, phase=0.):
    """l^2 E_return' for one zero-delay Dirac loop of length pi*l."""
    model = LoopReservoir(modulus, 1., phase)
    with mp.workdps(50):
        return float(-model.vacuum_energy(quadrature=False)/mp.pi)


def force_kernel(return_matrix, x):
    """Actual phase-dependent integrand and its singular-value upper bound.

    x=pi*l*kappa. The matrix includes the spin/Wilson phase; NS is -I.
    The exterior is fixed under the length variation. Matrix dimension
    counts paired Dirac channels, not two independently counted chiralities.
    """
    R = np.asarray(return_matrix, complex)
    x = float(x)
    if (R.ndim != 2 or R.shape[0] != R.shape[1] or len(R) == 0
            or not np.isfinite(R).all() or not math.isfinite(x) or x <= 0):
        raise ValueError('Require a finite nonempty square matrix and x>0.')
    singular = np.linalg.svd(R, compute_uv=False)
    if singular[0] > 1.:
        raise ValueError('The claimed Euclidean contraction hypothesis fails.')
    U = math.exp(-x)*R
    actual = -np.trace(np.linalg.solve(np.eye(len(R))-U, U)).real
    scaled = math.exp(-x)*singular
    bound = np.sum(scaled/(1.+scaled))
    return float(2*x*actual/math.pi**2), float(2*x*bound/math.pi**2)


def controls():
    rng = np.random.default_rng(129)
    cases = 0
    worst = 0.
    for n in (2, 4, 8):
        for _ in range(4):
            R = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            R *= .95/np.linalg.svd(R, compute_uv=False)[0]
            for x in (.01, .2, 1., 5.):
                actual, bound = force_kernel(R, x)
                worst = max(worst, actual-bound)
                assert actual <= bound+3e-14
                cases += 1
    # A nonzero nilpotent return has no closed winding; singular values
    # give an upper bound, not the actual force or a return probability.
    nilpotent = np.array([[0., .8], [0., 0.]])
    value, bound = force_kernel(nilpotent, .7)
    assert abs(value) < 1e-15 and bound > .01
    # Test the new length derivative against an independently differentiated
    # energy using the already established scalar loop model.
    length, delay, modulus, phase = 1.4, .6, .4, .37
    predicted = support_weight(modulus, phase)*(math.pi*length/(math.pi*length+delay))**2
    energy = lambda l: float(LoopReservoir(modulus, math.pi*l, phase, delay).vacuum_energy(quadrature=False))
    h = 2e-5
    derivative = length**2*(energy(length-2*h)-8*energy(length-h)
                           +8*energy(length+h)-energy(length+2*h))/(12*h)
    assert abs(derivative-predicted) < 2e-11
    direct = quad(lambda x: force_kernel(np.array([[-.4]]), x)[0],
                  0., 80., epsabs=2e-13, epsrel=2e-13)[0]
    assert abs(direct-support_weight(.4)) < 2e-13
    assert abs(support_weight(1.)-1/6) < 1e-15
    assert support_weight(.04) > 0 > support_weight(.04, math.pi)
    return dict(passed=True, nonnormal_pointwise_cases=cases,
                maximum_upper_bound_violation=worst,
                nilpotent_actual=value, nilpotent_upper=bound,
                delayed_derivative_error=abs(derivative-predicted),
                scalar_force_integral_error=abs(direct-support_weight(.4)))


def run():
    started = time.perf_counter()
    path = project_path('dihole_return_geometry_results.json')
    assert hashlib.sha256(path.read_bytes()).hexdigest() == '7f39c4cd67628961a1c8622fb67638b1247e59941bbae594feae64c0496b34bc'
    previous = json.loads(path.read_text())
    assert previous['passed']
    f = previous['geometry']['returning_flux_fraction']
    threshold = 1/24
    total = lambda r, phase=0.: f/6+(1-f)*support_weight(r, phase)
    critical = brentq(lambda r: total(r)-threshold, .02, .05, xtol=2e-15)
    checks = controls()
    rows = []
    for modulus, phase in ((0.,0.),(.01,0.),(.03,0.),(.04,0.),(.04,math.pi)):
        capacity = total(modulus,phase)
        margin = capacity-threshold
        rows.append(dict(residual_amplitude=modulus, phase=phase,
                         constant_vertex_return_probability=modulus**2,
                         scaled_return_force_per_q=capacity,
                         force_margin_per_q=margin,
                         optimistic_stationary_length_for_C_over_q_1=2/margin if margin>0 else None))
    result = dict(passed=checks['passed'], controls=checks,
        classical_fraction_used_only_as_unproved_proxy=f,
        ideal_all_channel_anomaly_threshold_per_q=threshold,
        critical_optimally_phased_residual_amplitude=critical,
        corresponding_constant_vertex_probability=critical**2,
        examples=rows,
        matrix_bound='-Re Tr[U(I-U)^-1] <= sum sigma_j(U)/(1+sigma_j(U)), ||U||<1; nonnormal matrices allowed.',
        force='l^2 E_return_prime = -(2/pi^2) integral_0^infinity x Re Tr[U(I-U)^-1] dx, U=exp(-x) R(i x/(pi*l)), with R fixed under l variation.',
        varying_exterior_term='If R=R_l, add (2*l^2/pi) integral exp(-pi*l*kappa) Re Tr[(I-U)^-1 partial_l R] d kappa. Contractivity alone does not bound this term.',
        envelope='If sigma_j(R(i*kappa))<=r_j exp(-kappa*L_j), then l^2 E_return_prime<=sum_j g(r_j)*(pi*l/(pi*l+L_j))^2; g(r)=-2 Li2(-r)/pi^2.',
        rejection='Within the separately justified E=C/l^2+q/(24*l)+E_return family, an all-l upper bound <=q/24 excludes stationarity for C>0. All other l-dependent energy terms must be included.',
        normalization='Each dimension represents one Dirac channel: two conjugate complex chiral symmetric-sea contributions. The matrix includes spin/Wilson holonomy.',
        assumptions=['A valid relative determinant and contour rotation for a specified common self-adjoint system/state/subtraction.',
            'Canonical Euclidean return contraction, fixed exterior and channel identification under the length variation, and enough integrability to differentiate.',
            'For infinite channel systems: trace-class or appropriate regulated relative determinant plus uniform summable bounds. Finite matrix checks do not establish this.',
            'The q/(24*l) anomaly and C/l^2 classical term require their own common-geometry derivation; a fixed-boundary cylinder does not automatically supply them.'],
        limitations=['No actual balanced-dihole quantum return operator or uniform singular-value envelope has been calculated.',
            'The critical amplitude uses an optimistic zero-delay, perfectly aligned phase model and an unproved classical-flux channel proxy.',
            'The displayed probability is only r^2 for that specified constant vertex. General real-axis probabilities do not determine imaginary-axis return or vacuum force.',
            'Exceeding the upper bound threshold does not prove a minimum; an illustrative minimum is confined to its one-coordinate toy energy.',
            'No absolute four-dimensional stress, common throat/exterior solution, gravitational stability, material recipe, or backward passage follows.'],
        sources=['https://doi.org/10.1073/pnas.35.7.408', 'https://arxiv.org/abs/1606.07834',
                 'https://arxiv.org/abs/0711.2395', 'https://arxiv.org/abs/1807.04726'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in
            ('scattering_vacuum.py','dihole_return_geometry_results.json','project_paths.py')},
        runtime_seconds=time.perf_counter()-started)
    project_path('quantum_return_support_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    draw(result)
    print(json.dumps(dict(passed=result['passed'],critical_amplitude=critical,
                         examples=rows,controls=checks),indent=2))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f = result['classical_fraction_used_only_as_unproved_proxy']
    amplitude = np.linspace(0.,.07,100)
    fig,ax = plt.subplots(figsize=(8.7,4.8))
    for phase,label,color in ((0.,'All weak returns optimally phased','#007e75'),
                              (math.pi,'Same strengths, opposite phase','#9d4263')):
        value = [f/6+(1-f)*support_weight(float(r),phase) for r in amplitude]
        ax.plot(100*amplitude,value,label=label,color=color,lw=2)
    ax.axhline(1/24,color='#353b45',ls='--',label='Necessary anomaly threshold')
    root=result['critical_optimally_phased_residual_amplitude']
    ax.plot(100*root,1/24,'o',color='#007e75')
    ax.annotate(f'Optimistic crossing: {100*root:.3f}% amplitude',
                (100*root,1/24),xytext=(.3,.0465),arrowprops=dict(arrowstyle='->'),fontsize=9)
    ax.set(xlabel='Return amplitude in nominally open channels (%)',
           ylabel=r'$\ell^2 E_{\mathrm{return}}^\prime/q$',
           title='The quantum support test needs phase, not just a channel count')
    ax.grid(alpha=.2);ax.legend(loc='lower left',fontsize=9)
    fig.text(.5,.035,'Conditional zero-delay model. Classical flux is an unproved channel proxy.',ha='center',fontsize=9)
    fig.text(.5,.007,'No physical return amplitude, supported throat or backward travel has been established.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,1))
    for ext in ('png','pdf'):fig.savefig(project_path('quantum_return_support_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
