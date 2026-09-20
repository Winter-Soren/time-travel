"""Normalized q=1 isolated-mouth modes for every half-integer global mu.

Set m=mu-1/2. Paired mouth labels are absolute n=max(1,abs(m)),...,N,
with angular eigenvalue sqrt(n(n+1)); only m=0 has a chiral zero mode.
The m=0,+/-1 implementations are delegated unchanged to azimuthal_modes.
The monopole spectrum is given in Borokhov, Kapustin and Wu, Section 4.1:
https://arxiv.org/pdf/hep-th/0206054 . Explicit phases below are our fixed
gauge convention, not an inference of a physical throat junction.

Theta is measured from global +z at BOTH mouths. Returned Cartesian-frame
spinors are at phi=0 and have unit boundary norm in rho ds=b^2 sin(theta)
dtheta. The exterior inner characteristic is Gamma=-sigma_r. The angular
operator in the unrescaled radial frame includes +cot(theta)/2; removing
that term requires multiplying the spinor by sqrt(sin(theta)).

Complete mouth multiplets n<=N require m=-N,...,N. The free outer labels
are different: j=abs(mu)+1/2,...,J. A common J>=N+1 is necessary to keep
an outer channel in every retained sector; it is not a convergence bound.
"""

import json
import operator

import numpy as np
from scipy.special import eval_jacobi, gammaln, lpmv

import azimuthal_modes as _legacy
from spherical_weyl_leads import _mouth_arguments


def _sector(mu):
    mu = float(mu)
    if not np.isfinite(mu):
        raise ValueError("mu must be a finite half-odd-integer")
    twice = int(round(2*mu))
    if abs(2*mu-twice) > 1e-12 or twice % 2 != 1:
        raise ValueError("mu must be a half-odd-integer")
    return (twice-1)//2


def mouth_channels(mu, nmax):
    """Absolute paired labels n=max(1,abs(mu-1/2)),...,nmax."""
    first = max(1, abs(_sector(mu)))
    try:
        nmax = operator.index(nmax)
    except TypeError as exc:
        raise ValueError("nmax must be an integer") from exc
    if nmax < first:
        raise ValueError("nmax must include at least one paired mouth mode")
    return np.arange(first, nmax+1, dtype=int)


def _scalar_pair(m, theta, channels):
    """Pole-safe f and g_m for abs(m)>=2, with sin(theta) scalar norm."""
    s = abs(m)
    n = channels[None, :]
    j = n-s
    c, sn = np.cos(theta)[:, None], np.sin(theta)[:, None]
    log_c = .5*(np.log(2*n+1)+gammaln(n-s+1)+gammaln(n+s+1)
                -(2*s+1)*np.log(2.)-2*gammaln(n+1))
    norm = np.exp(log_c)
    f = norm*sn**s*eval_jacobi(j, s, s, c)
    sign = 1 if m > 0 else -1
    half = np.sin(theta[:, None]/2) if sign > 0 else np.cos(theta[:, None]/2)
    g = (sign*2*n*norm/np.sqrt(n*(n+1))*sn**(s-1)*half**2
         *eval_jacobi(j, s+sign, s-sign, c))
    # Both exact scalar pole limits vanish for abs(m)>=2. Avoid sin(pi)
    # roundoff in the public endpoint values without altering nearby values.
    poles = (theta == 0) | (theta == np.pi)
    f[poles], g[poles] = 0., 0.
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(g)):
        raise ArithmeticError("nonfinite normalized mouth basis")
    return f, g


def mouth_characteristics(mu, theta, radius, port, nmax):
    """Return (plus, minus), each (2, points, len(mouth_channels(mu,N))).

    Gamma=-sigma_r; K_pair=sqrt(n(n+1))*sigma_x/radius. No zero mode or
    reflection-sewing parity is included in these arrays.
    """
    m = _sector(mu)
    channels = mouth_channels(mu, nmax)
    if abs(m) <= 1:
        return _legacy.mouth_characteristics(m+.5, theta, radius, port, nmax)
    theta, radius, ep, em = _mouth_arguments(theta, radius, port)
    f, g = _scalar_pair(m if port == "upper" else -m, theta, channels)
    f, g = f/radius, g/radius
    if port == "upper":
        return 1j*em[:, :, None]*g, ep[:, :, None]*f
    return em[:, :, None]*f, 1j*ep[:, :, None]*g


def zero_mode(mu, theta, radius, port):
    """The frozen m=0 zero spinor, or None; validate every sector's inputs."""
    m = _sector(mu)
    if abs(m) <= 1:
        return _legacy.zero_mode(m+.5, theta, radius, port)
    _mouth_arguments(theta, radius, port)
    return None


def reflection_parities(mu, nmax):
    """Declared reflected-frame signs (-1)^(n+abs(m)) on both transmissions.

    Reflections are unchanged. This separate endpoint-sewing choice leaves
    the local bases untouched and does not establish a physical junction.
    """
    m = _sector(mu)
    channels = mouth_channels(mu, nmax)
    if abs(m) <= 1:
        return _legacy.reflection_parities(m+.5, nmax)
    return 1-2*((channels+abs(m)) % 2)


def basis_checks():
    """Focused independent angular controls; no exterior or source solves."""
    c, weights = np.polynomial.legendre.leggauss(160)
    theta = np.arccos(c)
    angles = np.linspace(.12, np.pi-.12, 65)
    radius, nmax, step = .7, 12, 1e-4
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    errors = dict(associated_legendre=0., weighted_gram=0., characteristic=0.,
                  angular_action_relative=0., reflection=0., polar_limit=0.)
    legacy_exact = True
    rows = []

    def characteristic(values, at):
        cc, ss = np.cos(at)[:, None], np.sin(at)[:, None]
        return -np.stack((cc*values[0]+ss*values[1],
                          ss*values[0]-cc*values[1]))

    def legendre_pair(m, at, labels):
        s, cc, ss = abs(m), np.cos(at), np.sin(at)
        ff, gg = [], []
        for n in labels:
            # Independent associated Legendre evaluation/derivative, including
            # the conventional Condon--Shortley phase removed at BOTH ends.
            norm = (-1.)**s*np.exp(.5*(np.log((2*n+1)/2)
                         +gammaln(n-s+1)-gammaln(n+s+1)))
            p = lpmv(s, n, cc)
            prev = lpmv(s, n-1, cc) if n > s else np.zeros_like(cc)
            dp = (n*cc*p-(n+s)*prev)/(cc*cc-1)
            f, ft = norm*p, -ss*norm*dp
            ff.append(f)
            gg.append(-(ft-m*f/ss)/np.sqrt(n*(n+1)))
        return np.array(ff).T, np.array(gg).T

    for m in range(-8, 9):
        mu = m+.5
        labels = mouth_channels(mu, nmax)
        lam = np.sqrt(labels*(labels+1))
        parity = reflection_parities(mu, nmax)
        assert len(labels) == nmax-max(1, abs(m))+1
        for port in ("upper", "lower"):
            pp, mm = mouth_characteristics(mu, theta, radius, port, nmax)
            zero = zero_mode(mu, theta, radius, port)
            both = np.concatenate((pp, mm), axis=2)
            if zero is not None:
                both = np.concatenate((zero[:, :, None], both), axis=2)
            assert (zero is not None) == (m == 0)
            gram = radius**2*np.einsum("spn,p,spm->nm", both.conj(), weights, both)
            errors["weighted_gram"] = max(errors["weighted_gram"],
                                           float(np.max(abs(gram-np.eye(len(gram))))))
            errors["characteristic"] = max(errors["characteristic"],
                float(np.max(abs(characteristic(pp, theta)-pp))),
                float(np.max(abs(characteristic(mm, theta)+mm))))
            if zero is not None:
                expected = -zero if port == "upper" else zero
                errors["characteristic"] = max(errors["characteristic"], float(np.max(
                    abs(characteristic(zero[:, :, None], theta)[:, :, 0]-expected))))

            ff, gg = legendre_pair(m if port == "upper" else -m, angles, labels)
            _, _, ep, em = _mouth_arguments(angles, radius, port)
            exact = ((1j*em[:, :, None]*gg/radius, ep[:, :, None]*ff/radius)
                     if port == "upper" else
                     (em[:, :, None]*ff/radius, 1j*ep[:, :, None]*gg/radius))
            observed = mouth_characteristics(mu, angles, radius, port, nmax)
            errors["associated_legendre"] = max(errors["associated_legendre"],
                *(float(np.max(abs(a-b))) for a, b in zip(observed, exact)))

            def frame(at):
                pairs = mouth_characteristics(mu, at, radius, port, nmax)
                ch, sh = np.cos(at/2)[:, None], np.sin(at/2)[:, None]
                return np.stack([radius*np.sqrt(np.sin(at))[None, :, None]*np.stack(
                    (ch*p[0]+sh*p[1], -sh*p[0]+ch*p[1])) for p in pairs])

            framed = frame(angles)
            derivative = (-frame(angles+2*step)+8*frame(angles+step)
                          -8*frame(angles-step)+frame(angles-2*step))/(12*step)
            sign = 1 if port == "upper" else -1
            potential = (m+sign*.5*np.cos(angles))/np.sin(angles)
            action = (-1j*np.einsum("ab,sbpn->sapn", sx, derivative)
                      +np.einsum("ab,sbpn,p->sapn", sy, framed, potential))
            target = framed[::-1]*lam
            errors["angular_action_relative"] = max(errors["angular_action_relative"],
                float(np.max(abs(action-target))/np.max(abs(target))))

            # Analytic endpoint limits, independent of evaluation at a small
            # angle. Only abs(m)=0,1 can have nonzero scalar pole values.
            at = np.array([0., np.pi])
            ff, gg = np.zeros((2, len(labels))), np.zeros((2, len(labels)))
            cn = np.sqrt((2*labels+1)/2)
            mp = m if port == "upper" else -m
            if m == 0:
                ff[0], ff[1] = cn, (-1.)**labels*cn
            elif mp == 1:
                gg[1] = (-1.)**(labels-1)*cn
            elif mp == -1:
                gg[0] = -cn
            ep, em = np.array([[1., 0.], [0., 1.]]), np.array([[0., -1.], [1., 0.]])
            poles = ((1j*em[:, :, None]*gg/radius, ep[:, :, None]*ff/radius)
                     if port == "upper" else
                     (em[:, :, None]*ff/radius, 1j*ep[:, :, None]*gg/radius))
            actual_poles = mouth_characteristics(mu, at, radius, port, nmax)
            errors["polar_limit"] = max(errors["polar_limit"],
                *(float(np.max(abs(a-b))) for a, b in zip(actual_poles, poles)))
            if abs(m) >= 2:
                assert all(np.array_equal(a, np.zeros_like(a)) for a in actual_poles)
            else:
                for test_angles in (theta, at):
                    legacy_exact &= all(np.array_equal(a, b) for a, b in zip(
                        mouth_characteristics(mu, test_angles, radius, port, nmax),
                        _legacy.mouth_characteristics(mu, test_angles, radius, port, nmax)))
                    old_zero = _legacy.zero_mode(mu, test_angles, radius, port)
                    new_zero = zero_mode(mu, test_angles, radius, port)
                    legacy_exact &= (old_zero is None and new_zero is None
                                     or np.array_equal(old_zero, new_zero))
                legacy_exact &= np.array_equal(parity, _legacy.reflection_parities(mu, nmax))

        upper = mouth_characteristics(mu, angles, radius, "upper", nmax)
        lower = mouth_characteristics(mu, np.pi-angles, radius, "lower", nmax)
        cc, ss = np.cos(angles)[:, None], np.sin(angles)[:, None]
        for a, b in zip(upper, lower[::-1]):
            mapped = -np.stack((cc*a[0]+ss*a[1], -ss*a[0]+cc*a[1]))
            errors["reflection"] = max(errors["reflection"],
                                        float(np.max(abs(mapped-b*parity))))
        rows.append({"m": m, "mu": mu, "mouth_first": int(labels[0]),
                     "outer_first": int(abs(mu)+.5), "zero": m == 0})

    # All finite mouth multiplets through N are retained without an empty
    # extreme outer port, using a separate outer absolute cutoff J=N+1.
    complete_n, outer_j = 8, 9
    counts = [len(mouth_channels(m+.5, complete_n)) for m in range(-complete_n, complete_n+1)]
    assert sum(counts) == complete_n*(complete_n+2)
    assert all(int(abs(m+.5)+.5) <= outer_j for m in range(-complete_n, complete_n+1))
    assert legacy_exact
    assert errors["associated_legendre"] < 5e-11
    assert errors["weighted_gram"] < 3e-12
    assert errors["characteristic"] < 3e-13
    assert errors["angular_action_relative"] < 2e-8
    assert errors["reflection"] < 3e-12
    assert errors["polar_limit"] < 3e-12
    return {"passed": True, "m_range": [-8, 8], "mouth_nmax": nmax,
            "legacy_three_sectors_bitwise_equal": bool(legacy_exact),
            "max_errors": errors, "channels": rows,
            "complete_cut_control": {"mouth_N": complete_n, "outer_J": outer_j,
                                     "paired_labels_across_sectors": sum(counts)},
            "scope": "Isolated q=1 angular bases and declared reflected-frame signs only. No exterior solve, complete source, physical junction or infinite-cutoff convergence is claimed."}


if __name__ == "__main__":
    print(json.dumps(basis_checks(), separators=(",", ":"), allow_nan=False))
