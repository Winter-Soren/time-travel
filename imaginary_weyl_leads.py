"""Decaying Weyl boundary maps at E=i*kappa, kappa>0.

The real-frequency characteristic bases are unchanged.  These maps specify
upper-half-plane resolvent boundary conditions, not scattering probabilities.
The product mouth leads remain the declared constant-radius model; the tail
retains the actual flat point-pair potential only between R and the far L.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp
from scipy.special import kve

from magnetic_tail_map import _AngularCoupling
from spherical_weyl_leads import _channels, _positive


def outer_decaying_ratios(kappa, R, mu=.5, nmax=12):
    """Free spherical b=gamma*a at E=i*kappa in the existing outer basis.

    With x=kappa*R and r_n=K_(n-1/2)(x)/K_(n+1/2)(x),
    gamma_n=(r_n-1)/(r_n+1).  The recurrence is positive in r_n;
    for x>1 it is evaluated directly in gamma to avoid subtracting 1.
    """
    kappa, R = _positive(kappa, 'kappa'), _positive(R, 'R')
    channels = _channels(mu, nmax)
    x = kappa*R
    if not np.isfinite(x) or x <= 0:
        raise ValueError('kappa*R must be finite and positive')
    first = int(channels[0])
    result = np.empty(len(channels), dtype=complex)
    if x <= 1:
        ratio = x/(1+x)
        for n in range(1, int(channels[-1])+1):
            if n > 1:
                ratio = x/(2*n-1+x*ratio)
            if n >= first:
                result[n-first] = (ratio-1)/(ratio+1)
    else:
        half_inverse = .5/x
        gamma = -half_inverse/(1+half_inverse)
        for n in range(1, int(channels[-1])+1):
            if n > 1:
                term = (2*n-1)*half_inverse*(1-gamma)
                gamma = -(term+gamma)/(1+term)
            if n >= first:
                result[n-first] = gamma
    return result


def mouth_decaying_ratios(kappa, b, nmax=12):
    """Paired product-mouth modes n=1,...,nmax; no zero mode is included.

    m=sqrt(n(n+1))/b, p=i*sqrt(kappa^2+m^2), beta=m/(i*kappa+p).
    The lead coordinate points out of the exterior into the mouth, and
    incoming_characteristic=beta*outgoing_characteristic.  The unchanged
    upper zero mode enters the exterior and has no homogeneous decaying
    continuation; the lower zero mode leaves it and decays into its lead.
    """
    kappa, b = _positive(kappa, 'kappa'), _positive(b, 'b')
    channels = _channels(.5, nmax).astype(float)
    masses = np.sqrt(channels)*np.sqrt(channels+1)/b
    if not np.all(np.isfinite(masses)):
        raise ValueError('Product-lead masses must be finite')
    scale = np.maximum(kappa, masses)
    kk, mm = kappa/scale, masses/scale
    return -1j*mm/(kk+np.hypot(kk, mm))


def tail_map(kappa, R, L, q=1., d=1., mu=.5, nmax=12,
             angular_order=64, rtol=1e-10, atol=1e-12):
    """Return (Z, diagnostics) for the actual pair on R<=r<=L.

    Z is in the continuous spherical coefficient basis, before any facet
    Gram conversion.  At L, free decaying data replace the remaining tail.
    Tolerance and radius refinements are numerical controls, not error proofs.
    """
    kappa, R, L, rtol, atol = map(float, (kappa, R, L, rtol, atol))
    angular = _AngularCoupling(q, d, mu, nmax, angular_order)
    if (not np.all(np.isfinite([kappa, R, L, rtol, atol])) or kappa <= 0
            or not L > R > angular.d/2 or not 0 < rtol < 1 or not 0 < atol < 1):
        raise ValueError('Require kappa>0, L>R>d/2 and positive ODE tolerances<1')
    count = len(angular.channels)
    initial = np.diag(outer_decaying_ratios(kappa, L, mu, nmax))

    def rhs(radius, packed):
        Z = packed.reshape(count, count)
        C = angular.matrix(radius)
        return (C.conj().T+2*kappa*Z-Z@C@Z).reshape(-1)

    solved = solve_ivp(rhs, (L, R), initial.reshape(-1), method='DOP853',
                       rtol=rtol, atol=atol)
    if not solved.success:
        raise RuntimeError('Imaginary magnetic tail failed: '+solved.message)
    Z = solved.y[:, -1].reshape(count, count)
    if not np.all(np.isfinite(Z)):
        raise ArithmeticError('Nonfinite imaginary magnetic tail map')
    metric = np.eye(count)-Z.conj().T@Z
    eigenvalues = np.linalg.eigvalsh((metric+metric.conj().T)/2)
    diagnostics = dict(kappa=kappa, R=R, L=L, q=float(q), d=float(d),
        mu=float(mu), nmax=angular.nmax, channels=angular.channels.tolist(),
        angular_order=angular.order, rtol=rtol, atol=atol,
        nfev=int(solved.nfev), accepted_steps=len(solved.t)-1,
        angular_gram_error=angular.gram_error,
        radial_current_metric_eigenvalues=eigenvalues.tolist(),
        minimum_radial_current_metric_eigenvalue=float(eigenvalues[0]),
        maximum_singular_value=float(np.linalg.svd(Z, compute_uv=False)[0]),
        eigenvalues_clipped=False,
        interpretation='Resolvent current metric, not a transmission probability')
    return Z, diagnostics


def _complex_record(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def _free_and_mouth_controls():
    bessel_error = n1_error = mouth_error = mouth_integral_error = 0.
    for x in (.003, .3, 3., 30., 300.):
        channels = _channels(.5, 12)
        ratios = kve(channels-.5, x)/kve(channels+.5, x)
        exact = (ratios-1)/(ratios+1)
        result = outer_decaying_ratios(x, 1.)
        bessel_error = max(bessel_error, float(np.max(abs(result-exact))))
        n1_error = max(n1_error, float(abs(result[0]+1/(1+2*x))))
    for kappa in (.01, .4, 4.):
        b = .2
        masses = np.sqrt(np.arange(1, 9)*np.arange(2, 10))/b
        omega = np.hypot(kappa, masses)
        beta = mouth_decaying_ratios(kappa, b, 8)
        # Independent constant-lead Hamiltonian: (p sigma_z+m sigma_x)v=E v.
        first = 1j*omega+masses*beta-1j*kappa
        second = masses-1j*omega*beta-1j*kappa*beta
        mouth_error = max(mouth_error,
            float(np.max(np.maximum(abs(first), abs(second))/(omega+masses))))
        # Its decaying profile exp(-omega*y) satisfies the integrated current law.
        mouth_integral_error = max(mouth_integral_error, float(np.max(abs(
            1-abs(beta)**2-kappa*(1+abs(beta)**2)/omega))))
    # Check the recurrence retains a nonzero small gamma without Hankel overflow.
    large = outer_decaying_ratios(1e200, 1., nmax=3)
    return dict(modified_Bessel_maximum_error=bessel_error,
        free_n1_analytic_error=n1_error, product_lead_eigenvector_error=mouth_error,
        product_lead_integrated_current_error=mouth_integral_error,
        large_argument_scaled_gamma=(1e200*large.real).tolist())


def _linear_current_control():
    """Independent linear propagation plus a dense-solution current integral."""
    kappa, R, L, count = .4, 2., 5., 3
    angular = _AngularCoupling(nmax=count)
    initial_Z = np.diag(outer_decaying_ratios(kappa, L, .5, count))
    initial = np.vstack((np.eye(count), initial_Z))

    def rhs(radius, packed):
        a, b = packed.reshape(2, count, count)
        C = angular.matrix(radius)
        return np.vstack((-kappa*a+C@b, C.conj().T@a+kappa*b)).reshape(-1)

    solved = solve_ivp(rhs, (L, R), initial.reshape(-1), method='DOP853',
        rtol=8e-13, atol=8e-15, dense_output=True)
    if not solved.success:
        raise RuntimeError(solved.message)
    a, b = solved.y[:, -1].reshape(2, count, count)
    direct = np.linalg.solve(a.T, b.T).T
    riccati, _ = tail_map(kappa, R, L, nmax=count, rtol=8e-13, atol=8e-15)

    def density(radius):
        aa, bb = solved.sol(radius).reshape(2, count, count)
        return 2*kappa*(aa.conj().T@aa+bb.conj().T@bb)

    integrated, quadrature_error = quad_vec(density, R, L, epsabs=2e-12,
                                           epsrel=2e-12)
    initial_current = np.eye(count)-initial_Z.conj().T@initial_Z
    final_current = a.conj().T@a-b.conj().T@b
    metric = np.diag([1.]*count+[-1.]*count)
    C = angular.matrix(R)
    generator = np.block([[-kappa*np.eye(count), C],
                          [C.conj().T, kappa*np.eye(count)]])
    return dict(kappa=kappa, R=R, L=L, nmax=count,
        maximum_map_difference=float(np.max(abs(direct-riccati))),
        integrated_current_identity_absolute_error=float(np.max(abs(
            final_current-initial_current-integrated))),
        integrated_current_identity_relative_error=float(np.max(abs(
            final_current-initial_current-integrated))/np.max(abs(final_current))),
        density_quadrature_error_estimate=float(quadrature_error),
        generator_current_identity_error=float(np.max(abs(
            generator.conj().T@metric+metric@generator+2*kappa*np.eye(2*count)))),
        final_plus_block_condition=float(np.linalg.cond(a)), nfev=int(solved.nfev))


def checks():
    """Bounded isolated resolvent checks; no real-energy or historical reruns."""
    fine = dict(rtol=3e-12, atol=3e-14)
    cases, maps = {}, {}

    def evaluate(name, kappa=.2, R=2., L=32., **options):
        Z, diagnostics = tail_map(kappa, R, L, **options)
        cases[name], maps[name] = diagnostics, Z
        return Z

    free = evaluate('free', q=0., **fine)
    free_error = float(np.max(abs(free-np.diag(outer_decaying_ratios(.2, 2.)))))
    far_rows, previous = [], None
    for L in (8., 16., 32., 64.):
        Z = evaluate('L'+str(int(L)), L=L, **fine)
        row = dict(L=L)
        if previous is not None:
            row['maximum_entry_change'] = float(np.max(abs(Z-previous)))
        far_rows.append(row)
        previous = Z
    coarse = evaluate('ode_default')
    angular = evaluate('angular32', angular_order=32, **fine)
    n8 = evaluate('n8', nmax=8, **fine)
    n16 = evaluate('n16', nmax=16, **fine)
    tolerance_error = float(np.max(abs(coarse-maps['L32'])))
    angular_error = float(np.max(abs(angular-maps['L32'])))
    channel_rows = []
    for lower, higher, lo, hi in ((8, 12, n8, maps['L32']),
                                 (12, 16, maps['L32'], n16)):
        difference = lo-hi[:lower, :lower]
        channel_rows.append(dict(lower_nmax=lower, higher_nmax=higher,
            low4_maximum_change=float(np.max(abs(difference[:4, :4]))),
            full_common_block_maximum_change=float(np.max(abs(difference)))))
    evaluate('kappa_004_R4', kappa=.04, R=4., L=64., **fine)
    evaluate('kappa_08_R2', kappa=.8, **fine)
    special = _free_and_mouth_controls()
    linear = _linear_current_control()
    minimum_metric = min(c['minimum_radial_current_metric_eigenvalue']
                         for c in cases.values())
    tests = dict(
        modified_Bessel_reference=special['modified_Bessel_maximum_error'] < 2e-13,
        exact_free_n1=special['free_n1_analytic_error'] < 2e-15,
        large_argument_recurrence=np.max(abs(np.asarray(
            special['large_argument_scaled_gamma'])+np.arange(1, 4)/2)) < 2e-15,
        mouth_decaying_eigenvector=special['product_lead_eigenvector_error'] < 2e-14,
        mouth_integrated_current=special['product_lead_integrated_current_error'] < 2e-14,
        free_tail_map=free_error < 2e-10,
        independent_linear_propagation=linear['maximum_map_difference'] < 2e-10,
        integrated_radial_current=linear['integrated_current_identity_relative_error'] < 2e-10,
        generator_current_identity=linear['generator_current_identity_error'] < 2e-14,
        angular_quadrature=angular_error < 2e-10,
        ode_tolerance=tolerance_error < 2e-9,
        far_radius_refinement=(far_rows[3]['maximum_entry_change']
            < far_rows[2]['maximum_entry_change'] < far_rows[1]['maximum_entry_change']),
        final_far_radius_change=far_rows[-1]['maximum_entry_change'] < 2e-8,
        low_channel_refinement=channel_rows[-1]['low4_maximum_change'] < 2e-10,
        positive_resolvent_current=minimum_metric > -2e-10)
    tests = {name: bool(value) for name, value in tests.items()}
    if not all(tests.values()):
        raise AssertionError({name: value for name, value in tests.items() if not value})
    here = Path(__file__).resolve()
    return dict(passed=True, tests=tests, special_function_controls=special,
        independent_linear_system=linear, free_map_maximum_error=free_error,
        refinement=dict(kappa=.2, base_R=2., far_radius=far_rows,
            angular_order_32_to_64_maximum_change=angular_error,
            ode_default_to_refined_maximum_change=tolerance_error,
            channel_cutoff=channel_rows),
        minimum_observed_current_metric_eigenvalue=minimum_metric,
        cases=cases, reference_map=_complex_record(maps['L64']),
        source_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256((here.parent/name).read_bytes()).hexdigest()
            for name in ('magnetic_tail_map.py', 'spherical_weyl_leads.py',
                         'exterior_boundary_map.py')},
        equations=dict(energy='E=i*kappa, kappa>0; time convention exp(-iEt)',
            free_gamma='gamma=(K_(n-1/2)(kR)-K_(n+1/2)(kR))/(K_(n-1/2)(kR)+K_(n+1/2)(kR))',
            mouth_beta='-i*m/(kappa+sqrt(kappa^2+m^2)), m=sqrt(n(n+1))/b',
            radial='a_prime=-kappa*a+C*b; b_prime=C_dagger*a+kappa*b',
            riccati='Z_prime=C_dagger+2*kappa*Z-Z*C*Z',
            current='d(a_dagger*a-b_dagger*b)/dr=-2*kappa*(a_dagger*a+b_dagger*b)'),
        scope=[
            'Upper-half-plane resolvent data in the unchanged spherical characteristic basis. Imaginary-frequency coefficients and I-Z_dagger*Z are not transmission probabilities.',
            'Flat, static point-pair exterior, fixed half-integer mu and finite retained angular channels. The product mouth lead is still a declared isolated constant-radius model.',
            'The actual magnetic potential is retained on R<=r<=L. Free decay beyond L remains a truncation; numerical L refinement is not an analytic infinity-error bound.',
            'Channel refinements apply to the reported low block; the highest common channels are reported separately. No whole-operator, all-energy, ultraviolet, vacuum-stress, stability or backward-time conclusion follows.',
            'No facet Gram conversion or eigenvalue clipping is performed. ODE tolerances and refinement differences are not certified exact-error bounds.'
        ])


if __name__ == '__main__':
    result = checks()
    target = Path(__file__).with_name('imaginary_weyl_leads_results.json')
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in ('passed', 'tests',
        'special_function_controls', 'independent_linear_system', 'refinement',
        'minimum_observed_current_metric_eigenvalue', 'source_sha256')}, indent=2))
