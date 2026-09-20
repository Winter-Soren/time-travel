"""Small binary64 interval primitives for the compact dihole quadrature.

An interval is a tuple (lower, upper) of finite NumPy float64 arrays.  A
point(x) means the EXACT binary64 value supplied, not the mathematical decimal
which may have produced it.  Operations do not mutate their operands.

Arithmetic contract: NumPy's separate add/subtract/multiply/divide ufuncs use
IEEE754 binary64 round-to-nearest, gradual underflow, and no reassociation or
fused operations across ufunc calls.  Each result is enlarged by nextafter in
both directions, including subnormal results.  Overflow is rejected.  sqrt
uses NumPy only for an initial approximation: the returned endpoints are
verified by outward-rounded DIVISION, so no unverified sqrt rounding claim or
underflow-prone squared-endpoint comparison is needed.

matmul uses explicit interval outer products and additions, never BLAS.
sum_interval uses a balanced pairwise tree, never an opaque floating reduction.
These primitives do not certify an integrand, quadrature remainder, physical
parameters, or source model; the caller must supply those separate arguments.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

import numpy as np


_NEG_INF = np.float64(-np.inf)
_POS_INF = np.float64(np.inf)
_ARITHMETIC = dict(over='raise', divide='raise', invalid='raise', under='ignore')


def _array(value):
    a = np.asarray(value)
    if a.dtype.kind in 'iu':
        if np.any(a > 2**53) or (a.dtype.kind == 'i' and np.any(a < -2**53)):
            raise ValueError('Integer inputs must lie in [-2**53,2**53]; supply intended binary64 explicitly.')
    elif a.dtype.kind != 'f' or a.dtype.itemsize > 8:
        raise TypeError('Expected real binary64-compatible numeric input.')
    a = np.asarray(a, dtype=np.float64)
    if not np.isfinite(a).all():
        raise ValueError('Interval endpoints must be finite.')
    return a


def point(x):
    """Exact dyadic input; a float32 is promoted exactly, integers are checked."""
    a = _array(x)
    return a, a


def _ends(a):
    if not isinstance(a, tuple) or len(a) != 2:
        raise TypeError('Expected the interval tuple (lower, upper); use point for constants.')
    lo, hi = np.broadcast_arrays(_array(a[0]), _array(a[1]))
    if np.any(lo > hi):
        raise ValueError('Interval lower endpoint exceeds upper endpoint.')
    return lo, hi


def _finish(lo, hi):
    if not (np.isfinite(lo).all() and np.isfinite(hi).all()):
        raise FloatingPointError('A finite outward enclosure overflowed.')
    return lo, hi


def add(a, b):
    al, ah = _ends(a); bl, bh = _ends(b)
    with np.errstate(**_ARITHMETIC):
        return _finish(np.nextafter(al+bl, _NEG_INF), np.nextafter(ah+bh, _POS_INF))


def sub(a, b):
    al, ah = _ends(a); bl, bh = _ends(b)
    with np.errstate(**_ARITHMETIC):
        return _finish(np.nextafter(al-bh, _NEG_INF), np.nextafter(ah-bl, _POS_INF))


def _four_corners(al, ah, bl, bh, operation):
    p0, p1 = operation(al, bl), operation(al, bh)
    p2, p3 = operation(ah, bl), operation(ah, bh)
    lo = np.minimum(np.minimum(p0, p1), np.minimum(p2, p3))
    hi = np.maximum(np.maximum(p0, p1), np.maximum(p2, p3))
    return np.nextafter(lo, _NEG_INF), np.nextafter(hi, _POS_INF)


def mul(a, b):
    al, ah = _ends(a); bl, bh = _ends(b)
    with np.errstate(**_ARITHMETIC):
        return _finish(*_four_corners(al, ah, bl, bh, np.multiply))


def div(a, b):
    al, ah = _ends(a); bl, bh = _ends(b)
    if np.any((bl <= 0) & (bh >= 0)):
        raise ZeroDivisionError('Denominator interval contains zero.')
    with np.errstate(**_ARITHMETIC):
        return _finish(*_four_corners(al, ah, bl, bh, np.divide))


def square(a):
    lo, hi = _ends(a)
    with np.errstate(**_ARITHMETIC):
        x, y = lo*lo, hi*hi
        lower = np.maximum(0., np.nextafter(np.minimum(x, y), _NEG_INF))
        lower = np.where((lo <= 0) & (hi >= 0), 0., lower)
        return _finish(lower, np.nextafter(np.maximum(x, y), _POS_INF))


def sqrt(a):
    """Enclose sqrt by certified division tests, including subnormal inputs."""
    lo, hi = _ends(a)
    if np.any(lo < 0):
        raise ValueError('sqrt requires a nonnegative interval.')
    with np.errstate(**_ARITHMETIC):
        lower = np.maximum(0., np.nextafter(np.sqrt(lo), _NEG_INF))
        upper = np.where(hi == 0, 0., np.nextafter(np.sqrt(hi), _POS_INF))
        for _ in range(8):
            # lower <= down(lo/lower) => lower**2 <= lo exactly; likewise
            # upper >= up(hi/upper).  Division avoids square underflow.
            dl = np.divide(lo, lower, out=np.zeros_like(lower), where=lower > 0)
            du = np.divide(hi, upper, out=np.zeros_like(upper), where=upper > 0)
            bad_lo = (lower > 0) & (lower > np.nextafter(dl, _NEG_INF))
            bad_hi = (upper > 0) & (upper < np.nextafter(du, _POS_INF))
            if not (np.any(bad_lo) or np.any(bad_hi)):
                return _finish(lower, upper)
            lower = np.where(bad_lo, np.maximum(0., np.nextafter(lower, _NEG_INF)), lower)
            upper = np.where(bad_hi, np.nextafter(upper, _POS_INF), upper)
    raise ArithmeticError('sqrt approximation could not be certified by division.')


def _sum_axis(lo, hi, axis):
    lo, hi = np.moveaxis(lo, axis, -1), np.moveaxis(hi, axis, -1)
    if lo.shape[-1] == 0:
        zero = np.zeros(lo.shape[:-1], dtype=np.float64)
        return zero, zero.copy()
    while lo.shape[-1] > 1:
        n = lo.shape[-1] // 2
        new_lo = np.nextafter(lo[..., :2*n:2]+lo[..., 1:2*n:2], _NEG_INF)
        new_hi = np.nextafter(hi[..., :2*n:2]+hi[..., 1:2*n:2], _POS_INF)
        if lo.shape[-1] % 2:
            new_lo = np.concatenate((new_lo, lo[..., -1:]), axis=-1)
            new_hi = np.concatenate((new_hi, hi[..., -1:]), axis=-1)
        lo, hi = new_lo, new_hi
    return lo[..., 0], hi[..., 0]


def sum_interval(a, axis=None):
    """Balanced sums over all entries, one axis, or a tuple of distinct axes."""
    lo, hi = _ends(a)
    if axis is None:
        lo, hi = lo.reshape(-1), hi.reshape(-1)
        axes = (0,)
    else:
        axes = axis if isinstance(axis, tuple) else (axis,)
        if any(isinstance(i, bool) or not isinstance(i, (int, np.integer)) for i in axes):
            raise TypeError('axis must be an integer, tuple of integers, or None.')
        if any(not -lo.ndim <= i < lo.ndim for i in axes):
            raise ValueError('axis is outside the input dimensions.')
        axes = tuple(int(i) % lo.ndim for i in axes)
        if len(set(axes)) != len(axes):
            raise ValueError('Repeated reduction axis.')
    with np.errstate(**_ARITHMETIC):
        for ax in sorted(axes, reverse=True):
            lo, hi = _sum_axis(lo, hi, ax)
        return _finish(lo, hi)


def matmul(a, b):
    """Verified two-dimensional product by k interval outer products; no BLAS."""
    al, ah = _ends(a); bl, bh = _ends(b)
    if al.ndim != 2 or bl.ndim != 2 or al.shape[1] != bl.shape[0]:
        raise ValueError('matmul requires compatible two-dimensional arrays.')
    shape = al.shape[0], bl.shape[1]
    lower, upper = np.zeros(shape), np.zeros(shape)
    if al.shape[1] == 0:
        return lower, upper
    p0, p1, p2, p3, term_lo, term_hi = [np.empty(shape) for _ in range(6)]
    with np.errstate(**_ARITHMETIC):
        for j in range(al.shape[1]):
            np.multiply(al[:, j, None], bl[None, j, :], out=p0)
            np.multiply(al[:, j, None], bh[None, j, :], out=p1)
            np.multiply(ah[:, j, None], bl[None, j, :], out=p2)
            np.multiply(ah[:, j, None], bh[None, j, :], out=p3)
            np.minimum(p0, p1, out=term_lo)
            np.minimum(term_lo, p2, out=term_lo)
            np.minimum(term_lo, p3, out=term_lo)
            np.maximum(p0, p1, out=term_hi)
            np.maximum(term_hi, p2, out=term_hi)
            np.maximum(term_hi, p3, out=term_hi)
            np.nextafter(term_lo, _NEG_INF, out=term_lo)
            np.nextafter(term_hi, _POS_INF, out=term_hi)
            if j == 0:
                lower[...] = term_lo; upper[...] = term_hi
            else:
                np.add(lower, term_lo, out=lower)
                np.add(upper, term_hi, out=upper)
                np.nextafter(lower, _NEG_INF, out=lower)
                np.nextafter(upper, _POS_INF, out=upper)
        return _finish(lower, upper)


def run_checks(benchmark=False):
    """Exact-rational enclosure checks; optional one representative contraction."""
    started = time.perf_counter()
    from mpmath import mp

    def exact(x):
        return Fraction.from_float(float(x))

    def encloses(result, lower, upper):
        rl, rh = result
        assert exact(rl) <= lower <= upper <= exact(rh)

    tiny = float.fromhex('0x0.0000000000001p-1022')
    cases = [(-1.25, 1.75, -.5, 2.), (1e-200, 2e-200, -2e-200, -1e-200),
             (tiny, 2*tiny, .25, .5), (1e150, 1.5e150, -.75, -.5),
             (2.**52, 2.**52+1, -2.**52, -2.**52+1)]
    for al, ah, bl, bh in cases:
        a, b = (np.asarray(al), np.asarray(ah)), (np.asarray(bl), np.asarray(bh))
        xl, xh, yl, yh = map(exact, (al, ah, bl, bh))
        encloses(add(a,b), xl+yl, xh+yh)
        encloses(sub(a,b), xl-yh, xh-yl)
        products = [x*y for x in (xl,xh) for y in (yl,yh)]
        encloses(mul(a,b), min(products), max(products))
        if not bl <= 0 <= bh:
            quotients = [x/y for x in (xl,xh) for y in (yl,yh)]
            encloses(div(a,b), min(quotients), max(quotients))
        squares = [xl*xl,xh*xh]
        encloses(square(a), 0 if al <= 0 <= ah else min(squares), max(squares))

    square_inputs = [0., tiny, 2*tiny, np.finfo(float).tiny, 1e-200,
                     np.nextafter(1.,0.), 1., np.nextafter(1.,2.), 2., np.finfo(float).max]
    sl, sh = sqrt(point(square_inputs))
    for x, lower, upper in zip(square_inputs, sl, sh):
        assert 0 <= exact(lower)**2 <= exact(x) <= exact(upper)**2

    values = np.array([[2.**53, 1., -2.**53, tiny], [3., -1., -2., -tiny]])
    for ax in (None, 0, 1, (0,1)):
        got = sum_interval(point(values), axis=ax)
        groups = [values.ravel()] if ax is None or isinstance(ax,tuple) else list(np.moveaxis(values, ax, -1).reshape(-1,values.shape[ax]))
        for i, row in enumerate(groups):
            target = sum(map(exact,row),Fraction())
            encloses((np.ravel(got[0])[i],np.ravel(got[1])[i]),target,target)

    rng = np.random.default_rng(119)
    a, b = rng.normal(size=(5,7)), rng.normal(size=(7,4))
    a[0] = [2.**50,1.,-2.**50,1e-200,tiny,0.,-tiny]
    ai = (np.nextafter(a,-np.inf),np.nextafter(a,np.inf))
    bi = (np.nextafter(b,-np.inf),np.nextafter(b,np.inf))
    product = matmul(ai,bi)
    with mp.workdps(90):
        for i in range(5):
            for j in range(4):
                low, high = Fraction(), Fraction()
                direct = mp.mpf(0)
                for k in range(7):
                    corners = [exact(x)*exact(y) for x in (ai[0][i,k],ai[1][i,k]) for y in (bi[0][k,j],bi[1][k,j])]
                    low += min(corners); high += max(corners)
                    direct += mp.mpf(float(a[i,k]))*mp.mpf(float(b[k,j]))
                encloses((product[0][i,j],product[1][i,j]),low,high)
                assert mp.mpf(float(product[0][i,j])) <= direct <= mp.mpf(float(product[1][i,j]))

    bad_calls = [lambda: point(np.nan), lambda: point(2**53+1),
                 lambda: div(point(1.),(-1.,1.)), lambda: sqrt((-1.,2.)),
                 lambda: add((2.,1.),point(0.)),
                 lambda: matmul(point(np.ones((2,3))),point(np.ones((2,3))))]
    rejected = 0
    for call in bad_calls:
        try:
            call()
        except (ValueError,TypeError,ZeroDivisionError,FloatingPointError):
            rejected += 1
    assert rejected == len(bad_calls)
    with np.errstate(**_ARITHMETIC):
        assert np.add(np.float64(1),np.float64(2**-53)) == 1
        assert np.add(np.float64(tiny),np.float64(tiny)) == 2*tiny
        assert np.multiply(np.float64(tiny),np.float64(1)) == tiny

    pilot = None
    if benchmark:
        a, b = rng.normal(size=(648,128)), rng.normal(size=(128,648))
        t0 = time.perf_counter()
        product = matmul(point(a),point(b))
        pilot = dict(shape=[648,128,648], seconds=time.perf_counter()-t0,
                     largest_interval_width=float(np.max(product[1]-product[0])),
                     proof='Explicit outward interval outer products, not empirical agreement with BLAS.')
    tests = dict(exact_rational_elementary_enclosures=True,
                 subnormal_and_cancellation_enclosures=True, division_certified_square_roots=True,
                 pairwise_sum_enclosures=True, exact_rational_matrix_enclosures=True,
                 independent_mpmath_matrix_values=True, invalid_inputs_rejected=True,
                _round_to_nearest_and_gradual_underflow_smoke=True)
    return dict(passed=True,tests=tests,benchmark=pilot,runtime_seconds=time.perf_counter()-started,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Arithmetic building block under the stated IEEE754 contract; no quadrature or physical certificate by itself.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark',action='store_true')
    args = parser.parse_args()
    print(json.dumps(run_checks(args.benchmark),indent=2))
