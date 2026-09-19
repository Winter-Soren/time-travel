"""Origin-regular Dirichlet polynomial basis on the unit disk.

Modes are unnormalized, ordered by n=m+2*j, then increasing m, with real
before imaginary parts. Thus the first three are (1-r^2), x*(1-r^2),
y*(1-r^2). All derivatives are Cartesian; the third coordinate is z.
"""
from operator import index

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_jacobi


def _nonnegative_integer(value):
    value = index(value)
    if value < 0:
        raise ValueError("degree and padding must be nonnegative integers")
    return value


def _modes(degree, max_angular=None):
    for n in range(degree + 1):
        for m in range(n % 2, n + 1, 2):
            if max_angular is not None and m > max_angular:
                continue
            j = (n - m) // 2
            yield m, j, False
            if m:
                yield m, j, True


def quadrature(degree, padding=4, max_angular=None):
    """Return flat x, y, weights arrays for disk area integration.

    Gauss-Legendre integrates s=r^2 in [0,1]; uniform angles integrate theta.
    dA=ds*dtheta/2. Padding adds radial points and angular resolution.
    Nonpolynomial coefficient integrals still require refinement checks.
    """
    degree = _nonnegative_integer(degree)
    padding = _nonnegative_integer(padding)
    angular_degree = degree if max_angular is None else min(degree, _nonnegative_integer(max_angular))
    nodes, radial_weights = leggauss((degree + 4) // 2 + padding)
    s = (nodes + 1) / 2
    angles = 2 * np.pi * np.arange(2 * angular_degree + 5 + 2 * padding)
    angles /= len(angles)
    r = np.sqrt(s)[:, None]
    x = r * np.cos(angles)[None, :]
    y = r * np.sin(angles)[None, :]
    weights = np.broadcast_to(
        (radial_weights * np.pi / (2 * len(angles)))[:, None], x.shape
    )
    return x.ravel(), y.ravel(), weights.ravel()


def basis(degree, x, y, max_angular=None, compact=False):
    """Return B[P,N], grad[P,3,N], hess[P,3,3,N], with z derivatives zero.

    Coordinates broadcast, then flatten into P points. Mode (m,j) is the
    real/imaginary part of (1-s)*(x+i*y)^m*P_j^(2,m)(2*s-1).
    No radius division occurs, so the formulas are regular at the origin.
    An optional angular cap retains all radial modes with m <= max_angular;
    angular convergence must be checked separately when this is used.
    Compact storage returns grad[P,2,N] and hess[P,3,N], the latter in
    (xx,xy,yy) order. It omits only identically zero or duplicate derivatives.
    """
    degree = _nonnegative_integer(degree)
    if max_angular is not None:
        max_angular = _nonnegative_integer(max_angular)
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    x, y = x.ravel(), y.ravel()
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("coordinates must be finite")
    s = x * x + y * y
    w = x + 1j * y
    modes = list(_modes(degree, max_angular))
    count = len(modes)
    B = np.empty((x.size, count))
    grad = np.empty((x.size, 2, count)) if compact else np.zeros((x.size, 3, count))
    hess = np.empty((x.size, 3, count)) if compact else np.zeros((x.size, 3, 3, count))
    for column, (m, j, imaginary) in enumerate(modes):
        R = eval_jacobi(j, 2, m, 2 * s - 1)
        Rp = (j + m + 3) * eval_jacobi(j - 1, 3, m + 1, 2 * s - 1) if j else 0
        Rpp = ((j + m + 3) * (j + m + 4)
               * eval_jacobi(j - 2, 4, m + 2, 2 * s - 1)) if j > 1 else 0
        Q, Qp, Qpp = (1 - s) * R, -R + (1 - s) * Rp, -2 * Rp + (1 - s) * Rpp
        A = w ** m
        Ax = m * w ** (m - 1) if m else np.zeros_like(w)
        Axx = m * (m - 1) * w ** (m - 2) if m > 1 else np.zeros_like(w)
        value = A * Q
        dx = Ax * Q + 2 * x * A * Qp
        dy = 1j * Ax * Q + 2 * y * A * Qp
        dxx = Axx * Q + 4 * x * Ax * Qp + A * (2 * Qp + 4 * x * x * Qpp)
        dxy = 1j * Axx * Q + 2 * (y + 1j * x) * Ax * Qp + 4 * x * y * A * Qpp
        dyy = -Axx * Q + 4j * y * Ax * Qp + A * (2 * Qp + 4 * y * y * Qpp)
        part = np.imag if imaginary else np.real
        B[:, column] = part(value)
        grad[:, 0, column], grad[:, 1, column] = part(dx), part(dy)
        if compact:
            hess[:, 0, column], hess[:, 1, column], hess[:, 2, column] = part(dxx), part(dxy), part(dyy)
        else:
            hess[:, 0, 0, column] = part(dxx)
            hess[:, 0, 1, column] = hess[:, 1, 0, column] = part(dxy)
            hess[:, 1, 1, column] = part(dyy)
    return B, grad, hess


def _check():
    import sympy as sp

    # Independent polynomial differentiation includes radial and angular modes.
    xx, yy = sp.symbols("x y", real=True)
    radius2 = xx * xx + yy * yy
    px = np.array([0., .13, -.27, .4, -.61])
    py = np.array([0., -.31, .19, .1, -.2])
    degree = 8
    B, grad, hess = basis(degree, px, py)
    expected_first = np.column_stack((1-px*px-py*py,
                                     px*(1-px*px-py*py), py*(1-px*px-py*py)))
    np.testing.assert_allclose(B[:, :3], expected_first, atol=1e-15)
    selected = {(0, 0), (1, 0), (0, 3), (1, 2), (2, 2), (3, 1), (6, 1), (8, 0)}
    checked = 0
    for column, (m, j, imaginary) in enumerate(_modes(degree)):
        if (m, j) not in selected:
            continue
        complex_poly = (xx + sp.I*yy)**m * sp.jacobi(j, 2, m, 2*radius2-1)
        poly = sp.expand((1-radius2)*(sp.im(complex_poly) if imaginary else sp.re(complex_poly)))
        for variables, actual in (((), B[:, column]),
                                  ((xx,), grad[:, 0, column]), ((yy,), grad[:, 1, column]),
                                  ((xx, xx), hess[:, 0, 0, column]),
                                  ((xx, yy), hess[:, 0, 1, column]),
                                  ((yy, yy), hess[:, 1, 1, column])):
            reference = sp.lambdify((xx, yy), sp.diff(poly, *variables) if variables else poly, "numpy")
            np.testing.assert_allclose(actual, reference(px, py), rtol=3e-12, atol=3e-12)
        checked += 1
    assert np.all(np.isfinite(B)) and np.all(np.isfinite(grad)) and np.all(np.isfinite(hess))
    assert np.all(grad[:, 2] == 0) and np.all(hess[:, 2] == 0) and np.all(hess[:, :, 2] == 0)
    np.testing.assert_array_equal(hess, np.swapaxes(hess, 1, 2))
    theta = np.linspace(0, 2*np.pi, 97)
    np.testing.assert_allclose(basis(degree, np.cos(theta), np.sin(theta))[0], 0, atol=2e-14)
    qx, qy, weights = quadrature(degree, padding=0)
    np.testing.assert_allclose(weights.sum(), np.pi, atol=2e-14)
    np.testing.assert_allclose(np.dot(weights, qx*qx + qy*qy), np.pi/2, atol=2e-14)
    values = basis(degree, qx, qy)[0]
    gram = values.T @ (weights[:, None] * values)
    # Analytic disk norms include pi for angular cos/sin, 2*pi for m=0.
    norms = [np.pi * (j+1)*(j+2) / ((2*j+m+3)*(j+m+1)*(j+m+2))
             * (1 if m == 0 else .5) for m, j, _ in _modes(degree)]
    np.testing.assert_allclose(gram, np.diag(norms), rtol=2e-12, atol=2e-14)
    print(f"PASS: {checked} modes versus SymPy through second derivatives; "
          f"{B.shape[1]} modes finite at origin, zero at boundary, and disk-orthogonal.")


if __name__ == "__main__":
    _check()
