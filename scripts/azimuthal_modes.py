"""q=1 isolated-mouth modes for global mu=-1/2, 1/2, 3/2.

The monopole Dirac spectrum and zero-mode counting are standard; see
Borokhov, Kapustin and Wu, hep-th/0206054, Section 4 and Appendix:
https://arxiv.org/pdf/hep-th/0206054 .  The explicit phases and endpoint
reflection below are our derivation for the existing meridional solver.

Theta is measured from global +z at BOTH mouths.  Aphi is (1-cos(theta))/2
at the upper and (1+cos(theta))/2 at the lower isolated mouth.  Set
m=mu-1/2.  In the sqrt(sin(theta))-rescaled radial spin frame the angular
operator is K=-i*sigma_x*d_theta+sigma_y*(m +/- cos(theta)/2)/sin(theta).
Without that rescaling, d_theta must include +cot(theta)/2.

The returned Cartesian-frame spinors are evaluated at phi=0.  Their
boundary measure is rho ds=b**2 sin(theta)dtheta, and plus/minus means
the exterior-domain inner normal Gamma=-sigma_r.  Nonzero pairs have
K_pair=sqrt(n*(n+1))*sigma_x, n=1,...,nmax, in all three sectors.  The
physical transverse mass is K_pair/b.  Only mu=1/2 has a chiral zero mode.

For |m|=1 we use the same overall minus of the Condon--Shortley
order-|m| associated-Legendre convention at both ends and for both partners.
The existing mu=1/2 functions are delegated without any numerical change.

These are isolated monopole mouth bases.  The partner's field at a finite
mouth and a physical curved junction are not solved here.  An optional
reflected product-frame sewing has paired cross-transmission signs
(-1)**(n+abs(m)); it does NOT modify these local bases.  The inherited
same-label sewing remains a separate, explicit model choice.
"""

import json
import operator

import numpy as np
from scipy.special import eval_jacobi

from spherical_weyl_leads import _mouth_arguments
from spherical_weyl_leads import mouth_characteristics as _legacy_mouth
from spherical_weyl_leads import zero_mode as _legacy_zero


def _sector(mu):
    mu = float(mu)
    if mu not in (-.5, .5, 1.5):
        raise ValueError("mu must be one of -0.5, 0.5, 1.5")
    return int(mu-.5)


def _count(nmax):
    try:
        nmax = operator.index(nmax)
    except TypeError as exc:
        raise ValueError("nmax must be a positive integer") from exc
    if nmax < 1:
        raise ValueError("nmax must include at least one paired mode")
    return nmax


def _scalar_pair(m, theta, nmax):
    """Unweighted f and g_m; pole-safe Jacobi forms for m=+/-1 only."""
    n = np.arange(1, nmax+1)[None, :]
    c, s = np.cos(theta)[:, None], np.sin(theta)[:, None]
    norm = np.sqrt((2*n+1)/2)
    lam = np.sqrt(n*(n+1))
    # f=c_n sin(theta) P_n'(cos(theta))/lambda_n.
    f = norm*s*(n+1)/2*eval_jacobi(n-1, 1, 1, c)/lam
    # g_m=-(f'-m*csc(theta)*f)/lambda_n.  The factored forms
    # avoid subtracting nearly equal Legendre values at a pole.
    if m == 1:
        g = norm*np.sin(theta[:, None]/2)**2*eval_jacobi(n-1, 2, 0, c)
    elif m == -1:
        g = -norm*np.cos(theta[:, None]/2)**2*eval_jacobi(n-1, 0, 2, c)
    else:
        raise ValueError("this scalar helper is only for m=+/-1")
    return f, g


def mouth_characteristics(mu, theta, radius, port, nmax):
    """Return (plus, minus), each shape (2, npoints, nmax), n=1..nmax.

    Gamma=-sigma_r is the exterior inner normal.  No zero mode or
    reflected-sewing parity is included.  The two added sectors have the
    same mouth n=1 onset despite their different FREE OUTER mode onset.
    """
    m, nmax = _sector(mu), _count(nmax)
    if m == 0:
        return _legacy_mouth(theta, radius, port, nmax)
    theta, radius, er_plus, er_minus = _mouth_arguments(theta, radius, port)
    f, g = _scalar_pair(m if port == "upper" else -m, theta, nmax)
    f, g = f/radius, g/radius
    if port == "upper":
        return 1j*er_minus[:, :, None]*g, er_plus[:, :, None]*f
    return er_minus[:, :, None]*f, 1j*er_plus[:, :, None]*g


def zero_mode(mu, theta, radius, port):
    """Return the frozen mu=1/2 zero spinor (2, npoints), otherwise None.

    The upper zero is Gamma-minus and the lower zero is Gamma-plus.
    Input validation is retained even in a sector without a zero mode.
    """
    m = _sector(mu)
    if m == 0:
        return _legacy_zero(theta, radius, port)
    _mouth_arguments(theta, radius, port)
    return None


def reflection_parities(mu, nmax):
    """Explicit reflected-frame cross-transmission signs, n=1..nmax.

    Theta_lower=pi-theta_upper; the spin lift -exp(i*theta*sigma_y)
    fixes the zero-channel sign to +1 and gives this sign on BOTH paired
    cross transmissions.  Reflection blocks are unchanged.  This is a
    declared isolated-sphere product-frame gluing, not a physical junction.
    """
    m, nmax = _sector(mu), _count(nmax)
    return 1-2*((np.arange(1, nmax+1)+abs(m)) % 2)


def basis_checks():
    """Focused angular, polar and spin-frame controls; no exterior solve."""
    nodes, weights = np.polynomial.legendre.leggauss(96)
    theta = np.arccos(nodes)
    radius, nmax = .7, 8
    lam = np.sqrt(np.arange(1, nmax+1)*np.arange(2, nmax+2))
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)

    def gamma(values, angles):
        c, s = np.cos(angles)[:, None], np.sin(angles)[:, None]
        return -np.stack((c*values[0]+s*values[1], s*values[0]-c*values[1]))

    gram_error = sign_error = angular_error = explicit_error = 0.
    reflection_error = polar_scaling_error = 0.
    legacy_exact = True
    zero_counts = {}
    pole_rows = []
    angles, step = np.linspace(.12, np.pi-.12, 41), 1e-4
    for mu in (-.5, .5, 1.5):
        m = _sector(mu)
        parities = reflection_parities(mu, nmax)
        for port in ("upper", "lower"):
            plus, minus = mouth_characteristics(mu, theta, radius, port, nmax)
            zero = zero_mode(mu, theta, radius, port)
            both = np.concatenate((plus, minus), axis=2)
            if zero is not None:
                both = np.concatenate((zero[:, :, None], both), axis=2)
                expected = -zero if port == "upper" else zero
                sign_error = max(sign_error, float(np.max(abs(
                    gamma(zero[:, :, None], theta)[:, :, 0]-expected))))
            gram = radius**2*np.einsum("spn,p,spm->nm", both.conj(), weights, both)
            gram_error = max(gram_error, float(np.max(abs(gram-np.eye(len(gram))))))
            sign_error = max(sign_error, float(np.max(abs(gamma(plus, theta)-plus))),
                             float(np.max(abs(gamma(minus, theta)+minus))))

            def frame(at):
                pp, mm = mouth_characteristics(mu, at, radius, port, nmax)
                c, s = np.cos(at/2)[:, None], np.sin(at/2)[:, None]
                def rotate(values):
                    return radius*np.sqrt(np.sin(at))[None, :, None]*np.stack(
                        (c*values[0]+s*values[1], -s*values[0]+c*values[1]))
                return np.stack((rotate(pp), rotate(mm)))

            f0 = frame(angles)
            df = (-frame(angles+2*step)+8*frame(angles+step)
                  -8*frame(angles-step)+frame(angles-2*step))/(12*step)
            sign = 1 if port == "upper" else -1
            potential = (m+sign*.5*np.cos(angles))/np.sin(angles)
            action = (-1j*np.einsum("ab,sbpn->sapn", sx, df)
                      + np.einsum("ab,sbpn,p->sapn", sy, f0, potential))
            target = f0[::-1]*lam
            angular_error = max(angular_error, float(np.max(abs(action-target))/np.max(abs(target))))

            # n=1 elementary functions are independent of Jacobi evaluation.
            _, _, ep, em = _mouth_arguments(angles, radius, port)
            if m:
                f = np.sqrt(3)/2*np.sin(angles)
                g = np.sqrt(3/8)*((m if port == "upper" else -m)-np.cos(angles))
            else:
                f = np.sqrt(3/2)*np.cos(angles)
                g = np.sqrt(3)/2*np.sin(angles)
            exact = ((1j*em*g/radius, ep*f/radius) if port == "upper"
                     else (em*f/radius, 1j*ep*g/radius))
            observed = mouth_characteristics(mu, angles, radius, port, 1)
            explicit_error = max(explicit_error, *(float(np.max(abs(a[:, :, 0]-b)))
                                                  for a, b in zip(observed, exact)))

            endpoints = np.array([0., np.pi])
            ends = mouth_characteristics(mu, endpoints, radius, port, nmax)
            assert all(np.all(np.isfinite(a)) for a in ends)
            if not m:
                for at in (theta, endpoints):
                    legacy_exact &= all(np.array_equal(a, b) for a, b in zip(
                        mouth_characteristics(mu, at, radius, port, nmax),
                        _legacy_mouth(at, radius, port, nmax)))
                    legacy_exact &= np.array_equal(zero_mode(mu, at, radius, port),
                                                   _legacy_zero(at, radius, port))

            # All channel columns obey the component powers dictated by the
            # local gauge patch.  Use their maximum so isolated nodes do not
            # masquerade as a different leading pole power.
            for pole, a_axis in ((0, 0 if port == "upper" else 1),
                                 (1, 1 if port == "upper" else 0)):
                powers = [int(abs(mu-.5-a_axis)), int(abs(mu+.5-a_axis))]
                eps = 1e-4
                at = np.array([eps, eps/2]) if pole == 0 else np.pi-np.array([eps, eps/2])
                near = np.concatenate(mouth_characteristics(mu, at, radius, port, nmax), axis=2)
                for component, power in enumerate(powers):
                    peak = np.max(abs(near[component]), axis=1)
                    observed_ratio = peak[1]/peak[0]
                    polar_scaling_error = max(polar_scaling_error,
                        float(abs(observed_ratio/2**(-power)-1)))
                    if power:
                        assert max(abs(a[component, pole]).max() for a in ends) < 2e-12
                pole_rows.append({"mu": mu, "port": port, "pole": pole,
                                  "component_powers": powers})

        # Direct geometry: -S sends upper +/- to reflected lower -/+
        # with the same parity.  This does not modify the raw returned bases.
        up = mouth_characteristics(mu, angles, radius, "upper", nmax)
        low = mouth_characteristics(mu, np.pi-angles, radius, "lower", nmax)
        c, s = np.cos(angles)[:, None], np.sin(angles)[:, None]
        for a, b in zip(up, low[::-1]):
            mapped = -np.stack((c*a[0]+s*a[1], -s*a[0]+c*a[1]))
            reflection_error = max(reflection_error, float(np.max(abs(mapped-b*parities))))
        if not m:
            zup = zero_mode(mu, angles, radius, "upper")
            zlo = zero_mode(mu, np.pi-angles, radius, "lower")
            mapped = -np.stack((c[:, 0]*zup[0]+s[:, 0]*zup[1],
                               -s[:, 0]*zup[0]+c[:, 0]*zup[1]))
            reflection_error = max(reflection_error, float(np.max(abs(mapped-zlo))))

        # The two formal chiral solutions have these norm exponents at the
        # poles.  Integrability requires BOTH exponents > -1.
        f_exponents, g_exponents = (2*m+1, 1-2*m), (-1-2*m, -1+2*m)
        count = int(min(f_exponents) > -1)+int(min(g_exponents) > -1)
        assert count == int(m == 0)
        zero_counts[str(mu)] = count

    assert legacy_exact
    assert gram_error < 2e-12 and sign_error < 2e-13
    assert angular_error < 2e-8 and explicit_error < 2e-13
    assert reflection_error < 2e-13 and polar_scaling_error < 2e-6
    return {
        "passed": True,
        "sectors": [-.5, .5, 1.5], "paired_channels": list(range(1, nmax+1)),
        "legacy_mu_half_bitwise_equal": bool(legacy_exact),
        "weighted_gram_max_error": gram_error, "characteristic_max_error": sign_error,
        "angular_action_relative_error": angular_error, "explicit_n1_max_error": explicit_error,
        "reflected_spin_lift_max_error": reflection_error,
        "polar_scaling_relative_error": polar_scaling_error,
        "zero_mode_counts": zero_counts, "pole_component_powers": pole_rows,
        "reflection_parities": {str(mu): reflection_parities(mu, nmax).tolist()
                                for mu in (-.5, .5, 1.5)},
        "conventions": "q=1; m=mu-1/2; Gamma=-sigma_r; rho ds=b^2 sin(theta)dtheta; paired mass=sqrt(n(n+1))/b",
        "source": "https://arxiv.org/pdf/hep-th/0206054 (monopole Dirac spectrum); explicit gauge/phase/reflection calculation is ours",
        "scope": "Three azimuthal sectors complete the n=1 triplet but not higher angular multiplets. The exterior mixes n within each sector. Reflected gluing is a declared isolated-sphere product-frame model, not physical junction matching."
    }


if __name__ == "__main__":
    print(json.dumps(basis_checks(), indent=2))
