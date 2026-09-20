"""Validated GL/Hermite tables for the compact curved-residual quadrature.

Only one-dimensional tables are constructed.  No PDE, residual integral or
historical driver is run.  Each returned numerical quantity is a tuple of
outward binary64 arrays (lower, upper).  Internally, isolated mpmath ordinary
and interval contexts are used; no global precision is modified.

Approximate Legendre roots merely propose brackets.  Interval evaluation
proves opposite nonzero signs at every bracket's endpoints.  There are n
disjoint brackets inside (-1,1) for the degree-n polynomial, so these contain
all n roots, exactly one each.  Interval differentiation then encloses the
positive weights 2/((1-x^2)*P_n'(x)^2).  Hermite functions and derivatives
are evaluated at the high-precision root intervals, BEFORE binary64 export.

The normalization is phi_j(x)=H_j(x)*exp(-x^2/2)/sqrt(2^j*j!*sqrt(pi));
phi_j'=sqrt(j/2)*phi_{j-1}-sqrt((j+1)/2)*phi_{j+1}.
The nine width-two panels are the rule bounded in dihole_quadrature_remainder.
Analytic quadrature truncation, physical coefficients and tensor contractions
are outside this helper's scope.  The enclosures assume mpmath's documented
interval arithmetic; no binary64 BLAS assumption is used here.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext


PANEL_EDGES = tuple(range(-9, 10, 2))


def _positive_integer(name, value, minimum=1):
    if isinstance(value, bool) or int(value) != value or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}.')
    return int(value)


def _outward(value):
    lower, upper = float(value.a), float(value.b)
    if not math.isfinite(lower) or not math.isfinite(upper):
        raise OverflowError('An interval endpoint exceeds finite binary64 range.')
    return math.nextafter(lower, -math.inf), math.nextafter(upper, math.inf)


def _exact_mpf(iv, value):
    """Enclose the exact dyadic stored by an ordinary high-precision mpf."""
    sign, mantissa, exponent, _ = value._mpf_
    return (-1 if sign else 1)*iv.mpf(mantissa)*iv.mpf(2)**exponent


def _legendre(x, order):
    """P_n and P_{n-1}, with exact integer recurrence coefficients."""
    previous, current = 1+x*0, x
    if order == 1:
        return current, previous
    for degree in range(2, order+1):
        previous, current = current, ((2*degree-1)*x*current-(degree-1)*previous)/degree
    return current, previous


def _certified_gauss(order, digits):
    # Interval three-term recurrence loses correlation even at very narrow
    # point enclosures.  Explicit guard digits handle its finite-order
    # cancellation; the actual signs must still pass, independently of this
    # heuristic precision choice.  Approximate endpoints are exact dyadics
    # within this higher-precision context.
    iv = MPIntervalContext(); iv.dps = digits+max(40,order)
    mp = MPContext(); mp.dps = digits+30
    tolerance = mp.mpf(10)**(-digits-20)
    guesses = []
    for j in range(1, order+1):
        x = mp.cos(mp.pi*(mp.mpf(j)-mp.mpf(1)/4)/(mp.mpf(order)+mp.mpf(1)/2))
        for _ in range(60):
            pn, previous = _legendre(x, order)
            derivative = order*(x*pn-previous)/(x*x-1)
            step = pn/derivative
            x -= step
            if abs(step) < tolerance:
                break
        else:
            raise ArithmeticError('Approximate Gauss-root iteration did not converge.')
        guesses.append(x)
    guesses.sort()
    roots, weights, sign_pairs, attempts = [], [], [], []
    for guess in guesses:
        # Wide compared with interval roundoff, narrow compared with binary64.
        radius = mp.mpf(10)**(-digits+8)
        for attempt in range(12):
            left = _exact_mpf(iv, guess-radius)
            right = _exact_mpf(iv, guess+radius)
            pl = _legendre(left, order)[0]
            pr = _legendre(right, order)[0]
            left_sign = 1 if pl > 0 else (-1 if pl < 0 else 0)
            right_sign = 1 if pr > 0 else (-1 if pr < 0 else 0)
            if left_sign*right_sign == -1:
                break
            radius *= 10
        else:
            raise ArithmeticError('Could not validate a Gauss-root sign bracket.')
        root = iv.mpf([left.a, right.b])
        if not (-1 < root < 1) or (roots and not roots[-1].b < root.a):
            raise ArithmeticError('Gauss-root brackets are not disjoint inside (-1,1).')
        pn, previous = _legendre(root, order)
        derivative = order*(root*pn-previous)/(root*root-1)
        if not (derivative > 0 or derivative < 0):
            raise ArithmeticError('Weight derivative interval includes zero.')
        weight = 2/((1-root*root)*derivative**2)
        if not weight > 0:
            raise ArithmeticError('Gauss weight is not provably positive.')
        roots.append(root); weights.append(weight)
        sign_pairs.append([left_sign, right_sign]); attempts.append(attempt+1)
    # A completeness/count assertion is essential: individual brackets alone
    # would not show that no root was omitted or counted twice.
    assert len(roots) == order
    weight_sum = sum(weights, iv.mpf(0))
    assert weight_sum.a <= 2 <= weight_sum.b
    return iv, roots, weights, dict(
        normalized_root_sign_pairs=sign_pairs,
        root_bracket_attempts=attempts,
        normalized_weight_sum_interval=list(_outward(weight_sum)),
        root_count=order, all_brackets_disjoint=True,
        all_brackets_inside_unit_interval=True, all_weights_positive=True)


def _hermites(iv, x, count):
    """Values 0..count, giving the extra function required by derivatives."""
    first = iv.exp(-x*x/2)/iv.sqrt(iv.sqrt(iv.pi))
    values = [first, iv.sqrt(2)*x*first]
    for j in range(1, count):
        values.append(iv.sqrt(iv.mpf(2)/(j+1))*x*values[-1]
                      -iv.sqrt(iv.mpf(j)/(j+1))*values[-2])
    derivatives = []
    for j in range(count):
        value = -iv.sqrt(iv.mpf(j+1)/2)*values[j+1]
        if j:
            value += iv.sqrt(iv.mpf(j)/2)*values[j-1]
        derivatives.append(value)
    return values[:count], derivatives


def interval_grid(order=72, count=128, digits=70):
    """Return outward (lo,hi) arrays and metadata for the fixed nine panels.

    x, weights: (9*order,); basis, derivative: (9*order,count).
    The nodes are sorted in increasing order.  ``weights`` are the actual
    physical-panel weights: every panel has half-width one, so no scale is
    missing.  ``panel_ids`` is an integer array, and ``panel_edges`` is exact.
    """
    started = time.perf_counter()
    order = _positive_integer('order', order)
    count = _positive_integer('count', count)
    digits = _positive_integer('digits', digits, 30)
    iv, roots, weights, diagnostics = _certified_gauss(order, digits)
    nx = 9*order
    xlo, xhi, wlo, whi = (np.empty(nx, dtype=np.float64) for _ in range(4))
    blo, bhi, dlo, dhi = (np.empty((nx,count), dtype=np.float64) for _ in range(4))
    panel_ids = np.repeat(np.arange(9, dtype=np.int64), order)
    for panel, (left, right) in enumerate(zip(PANEL_EDGES, PANEL_EDGES[1:])):
        assert right-left == 2
        midpoint = (left+right)//2
        for local, (root, weight) in enumerate(zip(roots, weights)):
            row = panel*order+local
            x = root+midpoint
            xlo[row], xhi[row] = _outward(x)
            wlo[row], whi[row] = _outward(weight)
            basis, derivative = _hermites(iv, x, count)
            for column, value in enumerate(basis):
                blo[row,column], bhi[row,column] = _outward(value)
            for column, value in enumerate(derivative):
                dlo[row,column], dhi[row,column] = _outward(value)
    if not (np.all(xhi[:-1] < xlo[1:]) and np.all(wlo > 0)):
        raise ArithmeticError('Outward binary64 tables lost node ordering or positive weights.')
    if not all(np.isfinite(a).all() for a in (xlo,xhi,wlo,whi,blo,bhi,dlo,dhi)):
        raise OverflowError('Nonfinite table output.')
    metadata = dict(order_per_panel=order, count=count, requested_digits=digits,
        interval_digits=iv.dps,
        approximate_root_digits=digits+30, width=1, total_nodes=nx,
        node_weight_method='Disjoint interval Legendre sign brackets; exact derivative weight formula.',
        Hermite_method='Normalized interval recurrence through index count, then exact ladder derivative.',
        rounding='Every exported endpoint is moved one binary64 nextafter outwards.',
        tables_cover='Exact GL roots/weights and exact width1 Hermite values there.',
        excludes='Tensor contractions, curved coefficients, quadrature truncation, PDE and physical parameter errors.',
        max_basis_interval_width=float(np.max(bhi-blo)),
        max_derivative_interval_width=float(np.max(dhi-dlo)),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime_seconds=time.perf_counter()-started, **diagnostics)
    return dict(x=(xlo,xhi), weights=(wlo,whi), basis=(blo,bhi),
        derivative=(dlo,dhi), panel_ids=panel_ids,
        panel_edges=np.array(PANEL_EDGES,dtype=np.int64), metadata=metadata)


def self_check():
    """Tiny exact GL2 and normalized Hermite checks; no full grid or solve."""
    iv, roots, weights, metadata = _certified_gauss(2, 60)
    exact_root = 1/iv.sqrt(3)
    assert roots[0].a <= (-exact_root).a and roots[0].b >= (-exact_root).b
    assert roots[1].a <= exact_root.a and roots[1].b >= exact_root.b
    moment_widths = []
    for power in range(4):
        moment = sum((w*x**power for x,w in zip(roots,weights)), iv.mpf(0))
        target = iv.mpf(0) if power % 2 else iv.mpf(2)/(power+1)
        assert moment.a <= target.a and moment.b >= target.b
        moment_widths.append(float(moment.delta))
    h, d = _hermites(iv, iv.mpf(0), 4)
    phi0 = 1/iv.sqrt(iv.sqrt(iv.pi))
    targets = [phi0, iv.mpf(0), -phi0/iv.sqrt(2), iv.mpf(0)]
    derivative_targets = [iv.mpf(0), iv.sqrt(2)*phi0, iv.mpf(0), -iv.sqrt(3)*phi0]
    for value, target in zip(h+d, targets+derivative_targets):
        assert value.a <= target.b and value.b >= target.a
    # A nonzero point independently checks the Gaussian and polynomial
    # factors, rather than relying only on zeros/parity at the origin.
    x = iv.mpf(2)/3
    h, d = _hermites(iv, x, 4)
    gaussian = iv.exp(-iv.mpf(2)/9)/iv.sqrt(iv.sqrt(iv.pi))
    polynomials = [iv.mpf(1), 2*x, 4*x*x-2, 8*x**3-12*x]
    polynomial_derivatives = [iv.mpf(0), iv.mpf(2), 8*x, 24*x*x-12]
    for j, (polynomial, dp) in enumerate(zip(polynomials,polynomial_derivatives)):
        normalization = iv.sqrt((2**j)*math.factorial(j))
        expected = gaussian*polynomial/normalization
        expected_derivative = gaussian*(dp-x*polynomial)/normalization
        assert h[j].a <= expected.b and h[j].b >= expected.a
        assert d[j].a <= expected_derivative.b and d[j].b >= expected_derivative.a
    # Gaussian moments give exact normalization/orthogonality for phi0,phi1:
    # integral exp(-x^2)=sqrt(pi), integral x^2 exp(-x^2)=sqrt(pi)/2.
    norm0 = phi0**2*iv.sqrt(iv.pi)
    norm1 = 2*phi0**2*iv.sqrt(iv.pi)/2
    assert norm0.a <= 1 <= norm0.b and norm1.a <= 1 <= norm1.b
    grid = interval_grid(order=2,count=4,digits=50)
    assert grid['basis'][0].shape == (18,4)
    assert set(grid['panel_ids']) == set(range(9))
    return dict(passed=True, exact_GL2_roots=True, polynomial_moments_through_degree3=True,
        normalized_Hermite_and_derivative_at_zero=True,
        independent_nonzero_Hermite_polynomials=True,
        Gaussian_normalization=True, outward_grid_shapes_and_ordering=True,
        moment_interval_widths=moment_widths,
        small_grid_metadata=grid['metadata'])


if __name__ == '__main__':
    print(json.dumps(self_check(),indent=2))
