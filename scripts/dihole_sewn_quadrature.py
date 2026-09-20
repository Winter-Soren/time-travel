"""Tensor Gauss quadrature on a rectangle outside two actual round mouth cuts.

No polygon, mesh ownership, or point-rejection mask enters this rule.
Below the circle radius use rho=radius*sin(t), so circle half-height is
radius*cos(t). Above it use ordinary radial panels. Returned weights are
drho*dz: callers supply rho*P for the Euclidean norm or rho/P for the
differential-action form, and no extra 2*pi for a normalized fixed-mu mode.
Geometry and weights are binary64 numerical values, not interval enclosures.
"""
import json
import math
import operator

import numpy as np


def exterior_nodes(order, rho_bounds=(.005, .5), z_bounds=(-.85, .85),
                   radius=.05, k=.5):
    """Return (points[n,2], positive drho*dz weights[n]).

    The centers are (rho,z)=(0,+/-k). The rectangle must contain both
    circles in z and extend beyond their radius in rho; its positive
    rho lower bound truncates their caps. Radial panels include R,2R,4R.
    Axial panels include +/-k and +/-k+/-3R and zero. Extra radial splits
    handle any axial panel crossing a circular endpoint. Thus each tensor
    patch has smooth bounds even for nondefault disjoint mouth sizes.
    """
    try:
        if isinstance(order, bool):
            raise TypeError
        order = operator.index(order)
    except TypeError as exc:
        raise ValueError('order must be an integer >=2.') from exc
    if order < 2:
        raise ValueError('order must be an integer >=2.')
    rho_bounds, z_bounds = np.asarray(rho_bounds, float), np.asarray(z_bounds, float)
    if rho_bounds.shape != (2,) or z_bounds.shape != (2,):
        raise ValueError('rho_bounds and z_bounds must each contain two endpoints.')
    radius, k = float(radius), float(k)
    rlo, rhi = rho_bounds
    zlo, zhi = z_bounds
    if (not np.all(np.isfinite(np.r_[rho_bounds, z_bounds, radius, k]))
            or not (0 < rlo < rhi and 0 < radius < k and rhi > radius
                    and zlo <= -k-radius and zhi >= k+radius)):
        raise ValueError('Require positive rho, disjoint mouths, and both full z extents inside the rectangle.')

    gauss_x, gauss_w = np.polynomial.legendre.leggauss(order)
    unit_x, unit_w = (gauss_x+1)/2, gauss_w/2
    zcuts = sorted({float(z) for z in (-k-3*radius, -k, -k+3*radius,
                                      0., k-3*radius, k, k+3*radius)
                    if zlo < z < zhi})
    rcuts = {float(rlo), float(rhi)}
    rcuts.update(float(r) for r in (radius, 2*radius, 4*radius) if rlo < r < rhi)
    # A fixed z split can change which side of a circular endpoint it lies
    # on. Split rho at every such crossing before choosing tensor patches.
    for center in (-k, k):
        for z in zcuts:
            distance = abs(z-center)
            if distance < radius:
                crossing = math.sqrt((radius-distance)*(radius+distance))
                if rlo < crossing < rhi:
                    rcuts.add(crossing)
    rcuts = sorted(rcuts)

    points, weights = [], []
    for lo, hi in zip(rcuts[:-1], rcuts[1:]):
        below = hi <= radius
        if below:
            tlo, thi = math.asin(lo/radius), math.asin(hi/radius)
            angles = tlo+(thi-tlo)*unit_x
            rr = radius*np.sin(angles)
            hh = radius*np.cos(angles)
            rw = (thi-tlo)*unit_w*hh
        else:
            rr = lo+(hi-lo)*unit_x
            rw = (hi-lo)*unit_w
            hh = np.zeros(order)
        for rho, halfheight, radial_weight in zip(rr, hh, rw):
            if below:
                segments = ((zlo, -k-halfheight), (-k+halfheight, k-halfheight),
                            (k+halfheight, zhi))
            else:
                segments = ((zlo, zhi),)
            for bottom, top in segments:
                cuts = [bottom]+[z for z in zcuts if bottom < z < top]+[top]
                for left, right in zip(cuts[:-1], cuts[1:]):
                    if right <= left:
                        continue
                    zz = left+(right-left)*unit_x
                    points.append(np.column_stack((np.full(order, rho), zz)))
                    weights.append(radial_weight*(right-left)*unit_w)
    if not points:
        raise ValueError('Empty integration region.')
    return np.vstack(points), np.concatenate(weights)


def self_check():
    """Independent circle-antiderivative moments; no coefficient or PDE solve."""
    rlo, rhi, zlo, zhi, radius, k = .005, .5, -.85, .85, .05, .5
    halfheight = math.sqrt((radius-rlo)*(radius+rlo))
    cap_area = radius**2*math.acos(rlo/radius)-rlo*halfheight
    cap_rho_moment = 2*halfheight**3/3
    exact = np.array([(rhi-rlo)*(zhi-zlo)-2*cap_area,
                      (rhi*rhi-rlo*rlo)*(zhi-zlo)/2-2*cap_rho_moment,
                      (rhi-rlo)*(zhi*zhi-zlo*zlo)/2])
    rows = []
    for order in (8, 16, 24):
        points, weights = exterior_nodes(order)
        actual = np.array([np.sum(weights), weights@points[:, 0], weights@points[:, 1]])
        distance = np.minimum(np.hypot(points[:, 0], points[:, 1]-k),
                              np.hypot(points[:, 0], points[:, 1]+k))
        rows.append(dict(order=order, node_count=len(weights), moments=actual.tolist(),
            area_error=float(abs(actual[0]-exact[0])),
            rho_moment_error=float(abs(actual[1]-exact[1])),
            z_moment_error=float(abs(actual[2]-exact[2])),
            centroid=(actual[1:]/actual[0]).tolist(),
            minimum_weight=float(np.min(weights)), minimum_mouth_distance=float(np.min(distance))))
    tests = dict(exact_area=max(r['area_error'] for r in rows) < 3e-12,
                 exact_rho_moment=max(r['rho_moment_error'] for r in rows) < 3e-12,
                 exact_z_moment=max(r['z_moment_error'] for r in rows) < 3e-12,
                 positive_weights=min(r['minimum_weight'] for r in rows) > 0,
                 actual_round_exclusions=min(r['minimum_mouth_distance'] for r in rows) >= radius-2e-15,
                 bounded_default_size=rows[-1]['node_count'] == 18432)
    return dict(passed=bool(all(tests.values())), tests={key: bool(value) for key, value in tests.items()},
                exact_moments=exact.tolist(), exact_centroid=(exact[1:]/exact[0]).tolist(), controls=rows,
                weights='drho*dz; no rho, metric factor, or azimuthal multiplier included',
                scope='Exact circular-domain parameterization with numerical Gauss integration. Moment checks do not certify a field/form quadrature error.')


if __name__ == '__main__':
    report = self_check()
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report['passed']:
        raise SystemExit('Round-mouth quadrature controls failed.')
