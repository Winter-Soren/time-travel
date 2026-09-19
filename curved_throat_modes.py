"""Finite-profile Weyl channel scattering, without growing transfer matrices.

H=-i*sigma_3*d/dx+V*sigma_1, with time convention exp(-i*E*t).
The default V=lambda_n*sec(x/ell)/ell is a prescribed product throat.
Its finite cuts and any later sewing to an exterior are separate assumptions.
"""
import hashlib
import json
import operator
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(name+' must be finite and positive')
    return value


def _energy(value):
    value = complex(value)
    if (not np.isfinite(value.real) or not np.isfinite(value.imag)
            or value.imag < 0 or (value.real != 0 and value.imag != 0)):
        raise ValueError('Energy must be real or purely imaginary with Im(E)>0')
    return value


def _complex_record(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def throat_scattering(energy, n, ell=1., rho_c=10., q=1,
                      rtol=1e-10, atol=1e-12, profile='secant', b=None):
    """Return (S, diagnostics), incoming [f_plus(left),f_minus(right)].

    Outgoing is [f_minus(left),f_plus(right)].  The mass-pair convention is
    sigma_1, as for the declared product mouth leads.  The secant interval
    is x/ell in [-atan(rho_c),atan(rho_c)]; the constant option uses the
    same interval and V=sqrt(n(n+q))/b.  All energies use one common clock.

    For n=0 the returned paired matrix is a convenient formal control.
    A single Weyl field has only its allowed chiral zero mode, not both
    entries as independent physical channels.
    """
    energy = _energy(energy)
    ell, rho_c = _positive(ell, 'ell'), _positive(rho_c, 'rho_c')
    try:
        n, q = operator.index(n), operator.index(q)
    except TypeError as exc:
        raise ValueError('n and q must be integers') from exc
    if n < 0 or q <= 0:
        raise ValueError('Require n>=0 and q>0')
    rtol, atol = _positive(rtol, 'rtol'), _positive(atol, 'atol')
    if rtol >= 1 or atol >= 1:
        raise ValueError('ODE tolerances must be less than one')
    if profile not in ('secant', 'constant'):
        raise ValueError("profile must be 'secant' or 'constant'")
    if profile == 'constant':
        if b is None:
            raise ValueError('The constant profile requires its radius b')
        b = _positive(b, 'b')
    elif b is not None:
        raise ValueError('b is used only by the constant profile')
    cut = float(np.arctan(rho_c))
    length = 2*ell*cut
    lam = np.sqrt(float(n))*np.sqrt(float(n+q))
    scaled_energy = energy*ell
    if (not cut < np.pi/2 or not np.isfinite(length)
            or not np.isfinite(lam) or not np.isfinite(scaled_energy)):
        raise ValueError('Parameters exceed the finite floating-point interval')

    # Extend the right end of the interval.  t=exp(log_t) never comes from
    # division by an exponentially large transfer-matrix element.
    if n == 0:
        log_t = 1j*energy*length
        right = left = 0j
        nfev, steps = 0, 0
    else:
        def rhs(angle, state):
            right, log_t, left = state
            mass = lam/np.cos(angle) if profile == 'secant' else ell*lam/b
            return np.array((2j*scaled_energy*right-1j*mass*(1+right*right),
                1j*scaled_energy-1j*mass*right,
                -1j*mass*np.exp(2*log_t)), dtype=complex)

        solved = solve_ivp(rhs, (-cut, cut), np.zeros(3, dtype=complex),
                           method='DOP853', rtol=rtol, atol=atol)
        if not solved.success:
            raise RuntimeError('Throat scattering integration failed: '+solved.message)
        right, log_t, left = solved.y[:, -1]
        nfev, steps = int(solved.nfev), len(solved.t)-1
    transmission = np.exp(log_t)
    S = np.array(((left, transmission), (transmission, right)), dtype=complex)
    if not np.all(np.isfinite(S)) or not np.isfinite(log_t):
        raise ArithmeticError('Nonfinite throat scattering result')
    defect = np.eye(2)-S.conj().T@S
    eigenvalues = np.linalg.eigvalsh((defect+defect.conj().T)/2)
    diagnostics = dict(energy=_complex_record(energy), n=n, q=q, ell=ell,
        rho_c=rho_c, angular_eigenvalue=float(lam), optical_length=length,
        profile=profile, constant_radius=b if profile == 'constant' else None,
        rtol=rtol, atol=atol, nfev=nfev, accepted_steps=steps,
        log_transmission=_complex_record(log_t),
        transmission_underflowed=bool(transmission == 0),
        endpoint_reflection_symmetry_error=float(abs(left-right)),
        unitarity_defect_maximum=float(np.max(abs(defect))),
        passivity_defect_eigenvalues=eigenvalues.tolist(),
        minimum_passivity_defect_eigenvalue=float(eigenvalues[0]),
        singular_values_clipped=False,
        zero_mode_scope=('Formal paired zero control; retain only the physical Weyl chirality'
                         if n == 0 else 'Nonzero angular pair'),
        normalization='Canonical current per common coordinate time, f_dagger*sigma_3*f')
    return S, diagnostics


def _constant_exact(energy, mass, length):
    """Independent closed transfer exponential, divided before it can overflow."""
    energy = complex(energy)
    square = float(np.real(mass*mass-energy*energy))
    if square > 0:
        decay = np.sqrt(square)
        e2 = np.exp(-2*decay*length)
        denominator = 1+e2-1j*energy/decay*(1-e2)
        transmission = 2*np.exp(-decay*length)/denominator
        reflection = -1j*mass/decay*(1-e2)/denominator
    elif square < 0:
        momentum = np.sqrt(-square)
        integral = np.sin(momentum*length)/momentum
        denominator = np.cos(momentum*length)-1j*energy*integral
        transmission = 1/denominator
        reflection = -1j*mass*integral/denominator
    else:
        denominator = 1-1j*energy*length
        transmission = 1/denominator
        reflection = -1j*mass*length/denominator
    return np.array(((reflection, transmission),
                     (transmission, reflection)), dtype=complex)


def _direct_transfer_control(energy):
    """Low-action full transfer propagation, independent of the embedding RHS."""
    ell, rho_c, n, q = .7, .8, 1, 1
    cut = np.arctan(rho_c)
    lam = np.sqrt(n*(n+q))
    sigma3 = np.diag([1., -1.])
    sigma2 = np.array(((0., -1j), (1j, 0.)))

    def rhs(angle, flattened):
        matrix = flattened.reshape(2, 2)
        generator = 1j*energy*ell*sigma3+lam/np.cos(angle)*sigma2
        return (generator@matrix).reshape(-1)

    solved = solve_ivp(rhs, (-cut, cut), np.eye(2, dtype=complex).reshape(-1),
        method='DOP853', rtol=5e-13, atol=5e-15, dense_output=True)
    if not solved.success:
        raise RuntimeError(solved.message)
    matrix = solved.y[:, -1].reshape(2, 2)
    a, bb, c, d = matrix.reshape(-1)
    direct = np.array(((-c/d, 1/d), (np.linalg.det(matrix)/d, bb/d)))
    S, _ = throat_scattering(energy, n, ell, rho_c, q,
                             rtol=5e-13, atol=5e-15)
    start = np.array(((1., 0.), (S[0, 0], S[0, 1])), dtype=complex)

    def density(angle):
        columns = solved.sol(angle).reshape(2, 2)@start
        return ell*columns.conj().T@columns

    norm, quadrature_error = quad_vec(density, -cut, cut, epsabs=2e-12,
                                     epsrel=2e-12)
    absorption = np.eye(2)-S.conj().T@S
    identity_error = float(np.max(abs(absorption-2*complex(energy).imag*norm)))
    return dict(energy=_complex_record(energy), ell=ell, rho_c=rho_c, n=n,
        scattering_maximum_difference=float(np.max(abs(direct-S))),
        transfer_determinant_error=float(abs(np.linalg.det(matrix)-1)),
        integrated_norm_identity_error=identity_error,
        norm_quadrature_error_estimate=float(quadrature_error),
        norm_matrix=_complex_record(norm), nfev=int(solved.nfev))


def checks():
    """Small profile controls and synthetic finite-boundary comparison rows."""
    fine = dict(rtol=3e-12, atol=3e-14)
    zero_error = constant_error = refinement_error = symmetry_error = 0.
    real_unitarity = 0.
    minimum_imaginary = np.inf
    rows = []
    b, rho_c = .2, 4.
    ell = b*np.hypot(1., rho_c)
    for energy in (.4, .4j):
        for n in (0, 1, 2, 4, 12):
            S, info = throat_scattering(energy, n, ell, rho_c, **fine)
            product, product_info = throat_scattering(energy, n, ell, rho_c,
                profile='constant', b=b, **fine)
            rows.append(dict(energy=_complex_record(energy), n=n,
                curved_scattering=_complex_record(S),
                constant_scattering=_complex_record(product),
                diagnostics=info, constant_diagnostics=product_info))
            symmetry_error = max(symmetry_error, info['endpoint_reflection_symmetry_error'])
            if complex(energy).imag == 0:
                real_unitarity = max(real_unitarity, info['unitarity_defect_maximum'])
            else:
                minimum_imaginary = min(minimum_imaginary,
                                        info['minimum_passivity_defect_eigenvalue'])
            if n == 0:
                exact = np.array(((0., np.exp(1j*energy*info['optical_length'])),
                                  (np.exp(1j*energy*info['optical_length']), 0.)))
                zero_error = max(zero_error, float(np.max(abs(S-exact))))
            else:
                closed = _constant_exact(energy, np.sqrt(n*(n+1))/b,
                                         info['optical_length'])
                constant_error = max(constant_error, float(np.max(abs(product-closed))))
    # Include propagating and threshold constant channels, not only evanescence.
    for energy in (0., np.sqrt(2), 3., -.4):
        S, info = throat_scattering(energy, 1, ell=.7, rho_c=.8,
            profile='constant', b=1., **fine)
        exact = _constant_exact(energy, np.sqrt(2), info['optical_length'])
        constant_error = max(constant_error, float(np.max(abs(S-exact))))
    refinement_rows = []
    for energy, n in ((.4, 1), (.4j, 1), (4., 4), (.4j, 12)):
        coarse, _ = throat_scattering(energy, n, ell, rho_c)
        precise, _ = throat_scattering(energy, n, ell, rho_c, **fine)
        error = float(np.max(abs(coarse-precise)))
        refinement_error = max(refinement_error, error)
        refinement_rows.append(dict(energy=_complex_record(energy), n=n,
                                    maximum_change=error))
    direct = [_direct_transfer_control(.4), _direct_transfer_control(.3j)]
    high_S, high = throat_scattering(.4j, 160, ell=1., rho_c=10., **fine)
    tests = dict(
        exact_zero_channel=zero_error < 2e-15,
        analytic_constant_profile=constant_error < 2e-10,
        independent_short_transfer=max(r['scattering_maximum_difference'] for r in direct) < 2e-10,
        independent_transfer_determinant=max(r['transfer_determinant_error'] for r in direct) < 2e-11,
        integrated_norm_identity=max(r['integrated_norm_identity_error'] for r in direct) < 2e-10,
        real_energy_unitarity=real_unitarity < 2e-9,
        imaginary_energy_passivity=minimum_imaginary > -2e-10,
        symmetric_profile_reflections=symmetry_error < 2e-9,
        tolerance_refinement=refinement_error < 2e-9,
        large_barrier_log_transmission=(high['log_transmission']['real'] < -800
            and np.all(np.isfinite(high_S))
            and high['minimum_passivity_defect_eigenvalue'] > -2e-10))
    tests = {name: bool(value) for name, value in tests.items()}
    if not all(tests.values()):
        raise AssertionError({name: value for name, value in tests.items() if not value})
    here = Path(__file__).resolve()
    note = here.with_name('CURVED_THROAT_MATCHING.md')
    return dict(passed=True, tests=tests,
        controls=dict(zero_mode_error=zero_error, constant_profile_error=constant_error,
            real_energy_unitarity_error=real_unitarity,
            minimum_imaginary_passivity_eigenvalue=float(minimum_imaginary),
            reflection_symmetry_error=symmetry_error,
            independent_short_transfer=direct, tolerance_refinement=refinement_rows,
            high_barrier=high, high_barrier_scattering=_complex_record(high_S)),
        synthetic_calibration=dict(q=1, sphere_radius=b, rho_c=rho_c,
            ell=ell, optical_length=2*ell*np.arctan(rho_c), boundary_lapse=1.,
            meaning='ell=b*sqrt(1+rho_c^2) is a synthetic N_cut=1 calibration, not a matched MMP collar'),
        samples=rows,
        source_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
        input_sha256={'CURVED_THROAT_MATCHING.md': hashlib.sha256(note.read_bytes()).hexdigest()},
        equations=dict(potential='V_n=sqrt(n(n+q))*sec(x/ell)/ell',
            generator='f_prime=(i*E*sigma_3+V*sigma_2)*f',
            embedding='r_right_prime=2iE*r_right-iV*(1+r_right^2); log_t_prime=iE-iV*r_right; r_left_prime=-iV*exp(2log_t)',
            current='j_prime=-2*Im(E)*f_dagger*f',
            scattering='[f_minus(left),f_plus(right)]=S*[f_plus(left),f_minus(right)]',
            norm_identity='I-S_dagger*S=2*Im(E)*integral(F_dagger*F dx), F for the two unit incoming characteristics'),
        scope=[
            'A prescribed, static, spherical optical potential with an unchanged angular bundle. This does not solve Einstein-Maxwell equations or determine the quantum source.',
            'Nonzero angular pairs use the sigma_1 mass convention. A single Weyl zero sector has one chirality; its formal paired control matrix must not create an extra physical channel.',
            'Canonical current is per common coordinate time. A local physical spinor requires the radius, lapse, spin-frame and gauge transitions before sewing to an exterior.',
            'The constant-profile comparison and N_cut=1 calibration are synthetic finite-boundary tests. No physical collar, endpoint phase, spin structure, or transparent junction is inferred.',
            'At imaginary energy these are resolvent data, not propagation probabilities. Numerical passivity defects are reported without clipping; tolerances are not certified error bounds.',
            'Log transmission is retained even if its floating-point exponential underflows. No all-frequency, full angular, ultraviolet, vacuum-stress, stability or backward-time conclusion follows.'
        ])


if __name__ == '__main__':
    result = checks()
    Path(__file__).with_name('curved_throat_modes_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'controls',
                                                 'synthetic_calibration', 'source_sha256')}, indent=2))
