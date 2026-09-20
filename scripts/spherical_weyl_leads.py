"""Free outgoing spherical leads for H = sigma . (-i grad), E > 0.

The time convention is exp(-i E t).  Spinor harmonics follow Holder, Huang
and Ostrovsky, Phys. Rev. B 96, 174205 (2017), Appendix A, (A2)--(A12):
https://arxiv.org/pdf/1704.05481 .  In particular sigma_r Omega_{+/-n}
= -Omega_{-/+n}.  The outgoing spinor is proportional to
h^(1)_{n-1} Omega_{-n} - i h^(1)_n Omega_{+n}.  This minus-i phase follows
the first-order equations (A8)--(A12), and is independently checked here;
the opposite sign displayed in the extracted arXiv v1 (2.6) must not be
used with these same harmonic conventions.

The characteristic map is derived here from those free solutions.  It is
exact only if the region beyond the matching sphere is free.  Vanishing
net magnetic flux does not eliminate a monopole pair's exterior dipole
field: that tail needs radius refinement or explicit continuation.

At phi=0, u is the Cartesian-spin-frame meridional spinor with boundary
measure rho ds = R**2 sin(theta) dtheta.  The returned characteristic
bases are sqrt(2*pi) eta_+/- / R and have unit norm in that measure.
They are flux-orthonormal boundary bases, not individually radial waves.
No physical mouth boundary condition or throat model is supplied here.
"""

import operator

import numpy as np
from scipy.special import eval_jacobi, eval_legendre, sph_harm_y, spherical_jn, spherical_yn


def _channels(mu, nmax):
    mu = float(mu)
    if not np.isfinite(mu):
        raise ValueError("mu must be a finite half-odd-integer")
    twice = int(round(2 * mu))
    if abs(2 * mu - twice) > 1e-12 or twice % 2 != 1:
        raise ValueError("mu must be a half-odd-integer")
    try:
        nmax = operator.index(nmax)
    except TypeError as exc:
        raise ValueError("nmax must be an integer") from exc
    first = (abs(twice) + 1) // 2
    if nmax < first:
        raise ValueError("nmax must include at least one channel")
    return np.arange(first, nmax + 1, dtype=int)


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")
    return value


def _spinor_harmonics(mu, theta, channels):
    """Return Omega_-n, Omega_+n at phi=0; shape (2, points, channels)."""
    lower = np.empty((2, theta.size, len(channels)), dtype=complex)
    upper = np.empty_like(lower)
    m0, m1 = int(round(mu - .5)), int(round(mu + .5))

    def harmonic(ell, order):
        if abs(order) > ell:
            return np.zeros_like(theta, dtype=complex)
        return sph_harm_y(ell, order, theta, 0.)

    for column, n in enumerate(channels):
        j = n - .5
        lower[0, :, column] = np.sqrt((j + mu) / (2*j)) * harmonic(n-1, m0)
        lower[1, :, column] = np.sqrt((j - mu) / (2*j)) * harmonic(n-1, m1)
        upper[0, :, column] = -np.sqrt((j-mu+1) / (2*j+2)) * harmonic(n, m0)
        upper[1, :, column] = np.sqrt((j+mu+1) / (2*j+2)) * harmonic(n, m1)
    return lower, upper


def outer_characteristics(mu, theta, radius, nmax):
    """Return (plus, minus), each (2, npoints, nchannels).

    Channels are n = |mu| + 1/2, ..., nmax.  sigma_r plus=plus and
    sigma_r minus=-minus.  For an outer trace u, project its two
    coefficients with the rho ds weighted inner product; outgoing
    matching sets minus coefficients to outgoing_ratios * plus ones.
    """
    channels = _channels(mu, nmax)
    radius = _positive(radius, "radius")
    theta = np.asarray(theta, dtype=float)
    if (theta.ndim != 1 or not np.all(np.isfinite(theta))
            or np.any(theta < 0) or np.any(theta > np.pi)):
        raise ValueError("theta must be a one-dimensional array in [0, pi]")
    lo, hi = _spinor_harmonics(float(mu), theta, channels)
    factor = np.sqrt(np.pi) / radius
    return factor * (lo-hi), factor * (lo+hi)


def outgoing_ratios(energy, radius, mu, nmax):
    """Exact free outgoing characteristic ratios at real positive E R.

    gamma_n=(h_(n-1)-i h_n)/(h_(n-1)+i h_n), all h of the first kind.
    Recursing t_n=h_(n-1)/h_n avoids factorial-sized Hankel functions:
    t_1=z/(1-i*z), t_n=z/(2*n-1-z*t_(n-1)).
    """
    energy = _positive(energy, "energy")
    radius = _positive(radius, "radius")
    channels = _channels(mu, nmax)
    z = energy * radius
    if not np.isfinite(z) or z <= 0:
        raise ValueError("energy * radius must be finite and positive")
    inverse_ratio = z/(1-1j*z) if z <= 1 else 1/(1/z-1j)
    result = np.empty(len(channels), dtype=complex)
    first = int(channels[0])
    for n in range(1, int(channels[-1])+1):
        if n > 1:
            inverse_ratio = z/(2*n-1-z*inverse_ratio)
        if n >= first:
            result[n-first] = (inverse_ratio-1j)/(inverse_ratio+1j)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError("nonfinite outgoing characteristic ratio")
    return result


def _mouth_arguments(theta, radius, port):
    radius = _positive(radius, "radius")
    theta = np.asarray(theta, dtype=float)
    if (theta.ndim != 1 or not np.all(np.isfinite(theta))
            or np.any(theta < 0) or np.any(theta > np.pi)):
        raise ValueError("theta must be a one-dimensional array in [0, pi]")
    if port not in ("upper", "lower"):
        raise ValueError("port must be 'upper' or 'lower'")
    radial_plus = np.stack((np.cos(theta/2), np.sin(theta/2)))
    radial_minus = np.stack((-np.sin(theta/2), np.cos(theta/2)))
    return theta, radius, radial_plus, radial_minus


def mouth_characteristics(theta, radius, port, nmax):
    """q=1, global mu=1/2 product-lead pairs, n=1,...,nmax.

    Return plus/minus arrays (2, points, channels), normalized in rho ds.
    Their signs refer to the exterior inner normal Gamma=-sigma_r.
    Theta is measured from global +z at BOTH mouths.  The isolated lead
    gauges are Aphi_upper=sin(theta/2)**2 and Aphi_lower=cos(theta/2)**2.
    Pair phases give the transverse matrix sqrt(n*(n+1))*sigma_x/radius.

    This is a declared constant-radius cylindrical lead with an isolated
    monopole angular operator, not a derivation of the MMP throat or of
    the partner-field correction at a finite mouth.  No zero mode is
    included in these paired arrays; use zero_mode separately.
    """
    theta, radius, er_plus, er_minus = _mouth_arguments(theta, radius, port)
    channels = _channels(.5, nmax)
    legendre = np.empty((theta.size, len(channels)))
    jacobi = np.empty_like(legendre)
    cosine, sine = np.cos(theta), np.sin(theta)
    for column, n in enumerate(channels):
        factor = np.sqrt((2*n+1)/2)/radius
        legendre[:, column] = factor*eval_legendre(n, cosine)
        derivative = (n+1)/2*eval_jacobi(n-1, 1, 1, cosine)
        jacobi[:, column] = factor*sine*derivative/np.sqrt(n*(n+1))
    if port == "upper":
        return (1j*er_minus[:, :, None]*jacobi,
                er_plus[:, :, None]*legendre)
    return (er_minus[:, :, None]*legendre,
            1j*er_plus[:, :, None]*jacobi)


def zero_mode(theta, radius, port):
    """One normalized chiral mouth mode, shape (2, points).

    Upper is pure Gamma-minus (enters the exterior); lower is pure
    Gamma-plus (leaves it).  The lower overall phase is the e_r-minus
    convention (-sin(theta/2), cos(theta/2)), with no extra sign.
    """
    _, radius, er_plus, er_minus = _mouth_arguments(theta, radius, port)
    return (er_plus if port == "upper" else er_minus)/(np.sqrt(2)*radius)


def mouth_outgoing_ratios(energy, radius, nmax):
    """Product-lead beta=m/(E+p), p=sqrt((E+i0)**2-m**2), m=sqrt(n(n+1))/b.

    The lead coordinate increases out of the computational domain into
    the mouth.  Positive real p gives outward propagation, positive
    imaginary p gives decay.  In this pair-phase convention the no-input
    condition is incoming_characteristic=beta*outgoing_characteristic.
    For 0<E<m this is an evanescent relation, not a reflection probability.
    At E=m the returned continuous limit has zero flux.  If E>m, a unit-
    flux incident lead amplitude I instead adds sqrt(2*p/(E+p))*I to
    incoming_characteristic-beta*outgoing_characteristic.
    """
    energy = _positive(energy, "energy")
    radius = _positive(radius, "radius")
    channels = _channels(.5, nmax).astype(float)
    masses = np.sqrt(channels)*np.sqrt(channels+1)/radius
    if not np.all(np.isfinite(masses)):
        raise ValueError("lead masses must be finite")
    result = np.empty(len(channels), dtype=complex)
    propagating = energy >= masses
    t = masses[propagating]/energy
    result[propagating] = t/(1+np.sqrt(1-t)*np.sqrt(1+t))
    t = energy/masses[~propagating]
    result[~propagating] = 1/(t+1j*np.sqrt(1-t)*np.sqrt(1+t))
    return result


def _mouth_self_check():
    nodes, weights = np.polynomial.legendre.leggauss(96)
    theta = np.arccos(nodes)
    radius, nmax = .7, 8

    def gamma(values, angles):
        c, s = np.cos(angles)[:, None], np.sin(angles)[:, None]
        return -np.stack((c*values[0]+s*values[1],
                         s*values[0]-c*values[1]))

    gram_error = sign_error = zero_error = 0.
    for port in ("upper", "lower"):
        plus, minus = mouth_characteristics(theta, radius, port, nmax)
        zero = zero_mode(theta, radius, port)
        both = np.concatenate((zero[:, :, None], plus, minus), axis=2)
        gram = radius**2*np.einsum("spn,p,spm->nm", both.conj(), weights, both)
        gram_error = max(gram_error, float(np.max(abs(gram-np.eye(gram.shape[0])))))
        sign_error = max(sign_error, float(np.max(abs(gamma(plus, theta)-plus))),
                         float(np.max(abs(gamma(minus, theta)+minus))))
        expected = -zero if port == "upper" else zero
        zero_error = max(zero_error,
            float(np.max(abs(gamma(zero[:, :, None], theta)[:, :, 0]-expected))))

    # Apply K=-i sigma_x d/dtheta + sigma_y V to the radial spin-frame
    # functions, rather than checking only normalization or phase labels.
    angles = np.linspace(.12, np.pi-.12, 41)
    step = 1e-4
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    lambdas = np.sqrt(np.arange(1, nmax+1)*np.arange(2, nmax+2))
    angular_error = 0.
    for port in ("upper", "lower"):
        def frame_functions(at):
            plus, minus = mouth_characteristics(at, radius, port, nmax)
            c, s = np.cos(at/2)[:, None], np.sin(at/2)[:, None]
            def transform(values):
                return radius*np.sqrt(np.sin(at))[None, :, None]*np.stack(
                    (c*values[0]+s*values[1], -s*values[0]+c*values[1]))
            return np.stack((transform(plus), transform(minus)))

        functions = frame_functions(angles)
        derivative = (-frame_functions(angles+2*step)+8*frame_functions(angles+step)
                      -8*frame_functions(angles-step)+frame_functions(angles-2*step))/(12*step)
        potential = (.5 if port == "upper" else -.5)/np.tan(angles)
        applied = (-1j*np.einsum("ab,sbpn->sapn", sigma_x, derivative)
                   + np.einsum("ab,sbpn,p->sapn", sigma_y, functions, potential))
        expected = functions[::-1]*lambdas
        angular_error = max(angular_error,
            float(np.max(abs(applied-expected))/np.max(abs(expected))))

    dispersion_error = propagating_flux_error = evanescent_flux_error = 0.
    minimum_decay_rate = np.inf
    for energy in (.3, 1., 3., 20.):
        masses = np.sqrt(np.arange(1, nmax+1)*np.arange(2, nmax+2))/radius
        momentum = np.sqrt((energy+0j)**2-masses**2)
        beta = mouth_outgoing_ratios(energy, radius, nmax)
        outgoing = np.stack((np.ones_like(beta), beta))
        applied = np.stack((momentum+ masses*beta, masses-momentum*beta))
        scale = energy+np.max(masses)
        dispersion_error = max(dispersion_error,
            float(np.max(abs(applied-energy*outgoing))/scale))
        for mass, p, bb in zip(masses, momentum, beta):
            if p.imag > 0:
                minimum_decay_rate = min(minimum_decay_rate, float(p.imag))
                evanescent_flux_error = max(evanescent_flux_error, float(abs(1-abs(bb)**2)))
            else:
                p = p.real
                norm = np.sqrt(2*p*(energy+p))
                basis = np.array([[energy+p, mass], [mass, energy+p]])/norm
                metric = basis.T@np.diag([1., -1.])@basis
                propagating_flux_error = max(propagating_flux_error,
                    float(np.max(abs(metric-np.diag([1., -1.])))))
    example_error = float(abs(mouth_outgoing_ratios(1., 1., 1)[0]-(1-1j)/np.sqrt(2)))
    assert gram_error < 2e-12 and sign_error < 2e-13 and zero_error < 2e-13
    assert angular_error < 2e-8
    assert dispersion_error < 2e-14 and propagating_flux_error < 2e-13
    assert evanescent_flux_error < 2e-14 and minimum_decay_rate > 0
    assert example_error < 2e-14
    return {
        "weighted_basis_and_zero_mode_gram_max_error": gram_error,
        "characteristic_eigenvalue_max_error": sign_error,
        "zero_mode_chirality_max_error": zero_error,
        "direct_angular_pair_relative_residual": angular_error,
        "lead_dispersion_scaled_residual": dispersion_error,
        "propagating_flux_metric_max_error": propagating_flux_error,
        "evanescent_flux_max_error": evanescent_flux_error,
        "smallest_checked_positive_decay_rate": minimum_decay_rate,
        "q1_n1_evanescent_example_error": example_error,
        "scope": "declared q=1, mu=1/2 constant-radius monopole leads; no physical throat match",
    }


def self_check():
    """Focused independent normalization, first-order PDE and flux checks."""
    nodes, weights = np.polynomial.legendre.leggauss(96)
    theta = np.arccos(nodes)
    radius = 1.7

    def sigma_r(values):
        c, s = np.cos(theta)[:, None], np.sin(theta)[:, None]
        return np.stack((c*values[0]+s*values[1],
                         s*values[0]-c*values[1]))

    gram_error = spin_error = characteristic_error = 0.
    for mu in (.5, -.5, 1.5, -2.5):
        channels = _channels(mu, 7)
        lo, hi = _spinor_harmonics(mu, theta, channels)
        spin_error = max(spin_error, float(np.max(abs(sigma_r(lo)+hi))),
                         float(np.max(abs(sigma_r(hi)+lo))))
        plus, minus = outer_characteristics(mu, theta, radius, 7)
        both = np.concatenate((plus, minus), axis=2)
        gram = radius**2 * np.einsum("spn,p,spm->nm", both.conj(), weights, both)
        gram_error = max(gram_error, float(np.max(abs(gram-np.eye(gram.shape[0])))))
        characteristic_error = max(characteristic_error,
            float(np.max(abs(sigma_r(plus)-plus))),
            float(np.max(abs(sigma_r(minus)+minus))))

    ratio_error = n1_error = radial_error = wronskian_error = flux_error = 0.
    largest_ratio = 0.
    for z in (.02, .2, 1., 2.5, 9.3, 80.):
        gamma = outgoing_ratios(z/radius, radius, .5, 12)
        n = np.arange(1, 13)
        h0 = spherical_jn(n-1, z)+1j*spherical_yn(n-1, z)
        h1 = spherical_jn(n, z)+1j*spherical_yn(n, z)
        reference = (h0-1j*h1)/(h0+1j*h1)
        ratio_error = max(ratio_error, float(np.max(abs(gamma-reference))))
        n1_error = max(n1_error, float(abs(gamma[0]+1/(1-2j*z))))
        largest_ratio = max(largest_ratio, float(np.max(abs(gamma))))
        d0 = spherical_jn(n-1, z, True)+1j*spherical_yn(n-1, z, True)
        d1 = spherical_jn(n, z, True)+1j*spherical_yn(n, z, True)
        scale = abs(d0)+abs(d1)+(n+1)*(abs(h0)+abs(h1))/z
        radial_error = max(radial_error,
            float(np.max(abs(d0-(n-1)*h0/z+h1)/scale)),
            float(np.max(abs(d1+(n+1)*h1/z-h0)/scale)))
        wronskian = z*z*(spherical_jn(n, z)*spherical_yn(n-1, z)
                        - spherical_jn(n-1, z)*spherical_yn(n, z))
        wronskian_error = max(wronskian_error, float(np.max(abs(wronskian-1))))
        if z in (2.5, 9.3):
            lo, hi = _spinor_harmonics(.5, theta, np.arange(1, 6))
            wave = (z/radius)/np.sqrt(2) * (lo*h0[:5]-1j*hi*h1[:5])
            flux = 2*np.pi*radius**2*np.einsum(
                "spn,p,spn->n", wave.conj(), weights, sigma_r(wave)).real
            flux_error = max(flux_error, float(np.max(abs(flux-1))))

    # Independently differentiate full Cartesian spinors, including their
    # azimuthal phases; this catches a wrong radial sign or spin frame.
    pauli = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                      [[1, 0], [0, -1]]], dtype=complex)
    point = np.array([.52, .31, .79])
    energy, step = 1.7, 1e-4
    cartesian_error = 0.
    for mu, n in ((.5, 1), (.5, 3), (-1.5, 3)):
        def wave_at(position):
            rr = np.linalg.norm(position)
            angle = np.arccos(position[2]/rr)
            phi = np.arctan2(position[1], position[0])
            lo, hi = _spinor_harmonics(mu, np.array([angle]), [n])
            h0 = spherical_jn(n-1, energy*rr)+1j*spherical_yn(n-1, energy*rr)
            h1 = spherical_jn(n, energy*rr)+1j*spherical_yn(n, energy*rr)
            phases = np.exp(1j*np.array([mu-.5, mu+.5])*phi)
            return phases*(lo[:, 0, 0]*h0-1j*hi[:, 0, 0]*h1)

        applied = np.zeros(2, dtype=complex)
        for axis in range(3):
            shift = np.eye(3)[axis]*step
            derivative = (-wave_at(point+2*shift)+8*wave_at(point+shift)
                          -8*wave_at(point-shift)+wave_at(point-2*shift))/(12*step)
            applied += -1j*pauli[axis]@derivative
        expected = energy*wave_at(point)
        cartesian_error = max(cartesian_error,
            float(np.linalg.norm(applied-expected)/np.linalg.norm(expected)))

    for z in (1e-250, 1e-12, 1e12, 1e250):
        gamma = outgoing_ratios(z, 1., .5, 64)
        largest_ratio = max(largest_ratio, float(np.max(abs(gamma))))
    rejected = 0
    for energy in (0., -1., np.inf, np.nan):
        try:
            outgoing_ratios(energy, 1., .5, 1)
        except ValueError:
            rejected += 1

    assert gram_error < 2e-12 and spin_error < 2e-13
    assert characteristic_error < 2e-13
    assert ratio_error < 3e-14 and n1_error < 3e-14
    assert radial_error < 2e-14 and wronskian_error < 3e-13
    assert flux_error < 2e-10 and cartesian_error < 2e-8
    assert largest_ratio <= 1+2e-14 and rejected == 4
    return {
        "weighted_basis_gram_max_error": gram_error,
        "spin_flip_max_error": spin_error,
        "characteristic_eigenvalue_max_error": characteristic_error,
        "hankel_ratio_reference_max_error": ratio_error,
        "n1_closed_form_max_error": n1_error,
        "radial_weyl_scaled_residual": radial_error,
        "cartesian_weyl_relative_residual": cartesian_error,
        "spherical_wronskian_max_error": wronskian_error,
        "unit_outgoing_flux_max_error": flux_error,
        "largest_ratio_modulus_including_extreme_arguments": largest_ratio,
        "invalid_energies_rejected": rejected,
        "scope": "exact free exterior map; magnetic dipole tail not included",
        "mouth_leads": _mouth_self_check(),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2))
