"""Curved Ori Einstein-scalar initial data using smooth disk polynomials.

The basis enforces Cartesian regularity at the origin and Dirichlet data on
the circle. Solves nonlinear constraints only; no exterior/evolution claim.
Units 8*pi*G=c=1, z period 2*pi. One complex field denotes two real scalars.
"""
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import newton_krylov
from scipy.sparse.linalg import LinearOperator

from coupled_ori import exponent
from disk_basis import basis, quadrature


def scaled_basis(degree, x, y, radius=1., max_angular=None, compact=False):
    """Disk basis in physical coordinates, including Cartesian derivative scaling."""
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('radius must be finite and positive')
    values, first, second = basis(degree, np.asarray(x)/radius, np.asarray(y)/radius,
                                  max_angular=max_angular, compact=compact)
    first /= radius
    second /= radius**2
    return values, first, second


def scaled_quadrature(degree, padding=4, radius=1., source_radius=None, max_angular=None):
    """Area quadrature; split r^2 at an interior compact-source cutoff.

    Each radial panel uses the original number of Gauss points. This resolves
    the nonsmooth cutoff as a panel endpoint, without changing the trial basis.
    """
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('radius must be finite and positive')
    if source_radius is not None and (not np.isfinite(source_radius) or source_radius <= 0):
        raise ValueError('source_radius must be finite and positive')
    x, y, weights = quadrature(degree, padding=padding, max_angular=max_angular)
    if source_radius is None or source_radius >= radius:
        return radius*x, radius*y, radius**2*weights
    s = x*x+y*y
    cuts = (0., source_radius**2/radius**2, 1.)
    xs, ys, ws = [], [], []
    for left, right in zip(cuts[:-1], cuts[1:]):
        factor = radius*np.sqrt((left+(right-left)*s)/s)
        xs.append(factor*x)
        ys.append(factor*y)
        ws.append(radius**2*(right-left)*weights)
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)


def metric_geometry(g, dg, ddg):
    """Connection and Ricci from metric derivatives; derivative index first."""
    gi = np.linalg.inv(g)
    dgi = -np.einsum('...ia,...lab,...bj->...lij', gi, dg, gi)
    connection = np.empty(g.shape[:-2]+(3, 3, 3))
    dc = np.empty(g.shape[:-2]+(3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                connection[..., i, j, k] = sum(gi[..., i, m]*(dg[..., j, m, k]+
                    dg[..., k, m, j]-dg[..., m, j, k])/2 for m in range(3))
                for l in range(3):
                    dc[..., l, i, j, k] = sum((dgi[..., l, i, m]*(dg[..., j, m, k]+
                        dg[..., k, m, j]-dg[..., m, j, k])+gi[..., i, m]*(
                        ddg[..., l, j, m, k]+ddg[..., l, k, m, j]-ddg[..., l, m, j, k]))/2
                        for m in range(3))
    ricci = np.empty_like(g)
    for i in range(3):
        for j in range(3):
            ricci[..., i, j] = sum(dc[..., k, k, i, j]-dc[..., j, k, i, k] + sum(
                connection[..., k, k, l]*connection[..., l, i, j]-
                connection[..., k, j, l]*connection[..., l, i, k] for l in range(3))
                for k in range(3))
    return gi, dgi, connection, dc, np.einsum('...ij,...ij->...', gi, ricci)


def scalar_window(x, y, source_radius=1., taper_kind='polynomial'):
    """Compact scalar window and its Cartesian first derivatives.

    The smooth option is a distinct seed: it matches the polynomial window's
    central value and second spatial derivatives, but not its whole profile.
    """
    if not np.isfinite(source_radius) or source_radius <= 0:
        raise ValueError('source_radius must be finite and positive')
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    s = (x*x+y*y)/source_radius**2
    if taper_kind == 'polynomial':
        gap = np.maximum(0., 1-s)
        window, derivative_s = gap**4, -4*gap**3
    elif taper_kind == 'smooth':
        window, derivative_s = np.zeros_like(s), np.zeros_like(s)
        inside = s < 1
        gap = 1-s[inside]
        window[inside] = np.exp(-4*s[inside]/gap)
        derivative_s[inside] = -4*window[inside]/gap**2
    else:
        raise ValueError("taper_kind must be 'polynomial' or 'smooth'")
    return window, 2*x*derivative_s/source_radius**2, 2*y*derivative_s/source_radius**2


def seed_data(x, y, epsilon, a=.02, e=1/16, t=-1., source_radius=1., taper_kind='polynomial'):
    x, y = np.asarray(x), np.asarray(y)
    if not np.isfinite(source_radius) or source_radius <= 0:
        raise ValueError('source_radius must be finite and positive')
    P, Q = 2*e-a, 2*e+a
    D = -t+(e-P*P)*x*x+(e-Q*Q)*y*y
    assert np.min(D) > 0
    g = np.zeros((x.size, 3, 3))
    dg, ddg = np.zeros((x.size, 3, 3, 3)), np.zeros((x.size, 3, 3, 3, 3))
    g[:, 0, 0] = g[:, 1, 1] = 1
    g[:, 0, 2] = g[:, 2, 0] = P*x
    g[:, 1, 2] = g[:, 2, 1] = Q*y
    g[:, 2, 2] = e*(x*x+y*y)-t
    dg[:, 0, 0, 2] = dg[:, 0, 2, 0] = P
    dg[:, 1, 1, 2] = dg[:, 1, 2, 1] = Q
    dg[:, 0, 2, 2], dg[:, 1, 2, 2] = 2*e*x, 2*e*y
    ddg[:, 0, 0, 2, 2] = ddg[:, 1, 1, 2, 2] = 2*e
    gi, dgi, conn, dc, R = metric_geometry(g, dg, ddg)
    N = np.zeros_like(g)
    N[:, 0, 0], N[:, 1, 1] = P, Q
    N[:, 0, 2] = N[:, 2, 0] = e*x
    N[:, 1, 2] = N[:, 2, 1] = e*y
    N[:, 2, 2] = e*(P*x*x+Q*y*y)
    dN = np.zeros_like(dg)
    dN[:, 0, 0, 2] = dN[:, 0, 2, 0] = e
    dN[:, 1, 1, 2] = dN[:, 1, 2, 1] = e
    dN[:, 0, 2, 2], dN[:, 1, 2, 2] = 2*e*P*x, 2*e*Q*y
    dD = np.stack((2*(e-P*P)*x, 2*(e-Q*Q)*y, np.zeros_like(x)), axis=-1)
    Z = np.diag([0., 0., 1.])
    Kcov = (N+D[:, None, None]*Z/2)/np.sqrt(D)[:, None, None]
    dKcov = ((dN+dD[:, :, None, None]*Z/2)/np.sqrt(D)[:, None, None, None]-
              Kcov[:, None]*dD[:, :, None, None]/(2*D[:, None, None, None]))
    K = np.einsum('pij,pij->p', gi, Kcov)
    dK = np.einsum('plij,pij->pl', dgi, Kcov)+np.einsum('pij,plij->pl', gi, dKcov)
    A = gi@Kcov@gi-K[:, None, None]*gi/3
    dA = np.empty_like(dg)
    for l in range(3):
        dA[:, l] = (dgi[:, l]@Kcov@gi+gi@dKcov[:, l]@gi+gi@Kcov@dgi[:, l]-
                     dK[:, l, None, None]*gi/3-K[:, None, None]*dgi[:, l]/3)
    gradK = np.einsum('pij,pj->pi', gi, dK)
    expected_R = -2*e/D+2*t*(e+P*Q)/D**2
    assert np.max(abs(R-expected_R)) < 1e-12
    b, d = -a/(1+np.sqrt(1+8*a)), a/(1+np.sqrt(1-8*a))
    S = t+(a/2-e+b)*x*x+(-a/2-e+d)*y*y
    Sx, Sy = 2*(a/2-e+b)*x, 2*(-a/2-e+d)*y
    window, wx, wy = scalar_window(x, y, source_radius, taper_kind)
    p = exponent(a, 0)+2j
    phase = epsilon*np.exp(p*np.log(-S))
    field = phase*window
    C = np.stack((phase*(wx+window*p*Sx/S), phase*(wy+window*p*Sy/S), 1j*field), axis=-1)
    Pi = (D*field*p/S-P*x*C[:, 0]-Q*y*C[:, 1]+C[:, 2])/np.sqrt(D)
    G2 = np.einsum('pij,pi,pj->p', gi, C.conj(), C).real
    source = np.einsum('pij,pj->pi', gi, (Pi.conj()[:, None]*C).real)
    return dict(x=x, y=y, g=g, dg=dg, ddg=ddg, gi=gi, dgi=dgi, conn=conn, dc=dc,
        R=R, K=K, dK=dK, Kcov=Kcov, A=A, dA=dA, gradK=gradK,
        field=field, C=C, Pi=Pi, G2=G2, source=source, det=D)


def longitudinal(data, W, dW, ddW):
    DW = dW+np.einsum('pijk,pk->pji', data['conn'], W)
    dDW = ddW+np.einsum('plijk,pk->plji', data['dc'], W)+np.einsum('pijk,plk->plji', data['conn'], dW)
    raised = np.einsum('pik,pkj->pij', data['gi'], DW)
    divW = np.trace(DW, axis1=-2, axis2=-1)
    ddivW = np.trace(dDW, axis1=-2, axis2=-1)
    L = raised+raised.transpose(0, 2, 1)-(2/3)*data['gi']*divW[:, None, None]
    d_raised = np.einsum('plik,pkj->plij', data['dgi'], DW)+np.einsum('pik,plkj->plij', data['gi'], dDW)
    dL = d_raised+d_raised.transpose(0, 1, 3, 2)-(2/3)*(data['dgi']*divW[:, None, None, None]+
            data['gi'][:, None]*ddivW[:, :, None, None])
    return L, dL


def divergence(conn, tensor, derivative):
    out = np.empty((tensor.shape[0], 3))
    for i in range(3):
        out[:, i] = sum(derivative[:, j, i, j]+sum(conn[:, i, j, k]*tensor[:, k, j]+
            conn[:, j, j, k]*tensor[:, i, k] for k in range(3)) for j in range(3))
    return out


def fields(matrices, coefficients):
    B, dB, ddB = matrices
    c = coefficients.reshape(B.shape[1], 4)
    values = B@c
    if ddB.ndim == 3:
        # Contract before expanding to full three-dimensional tensor jets.
        # The z derivatives vanish; mixed Cartesian derivatives are symmetric.
        derivative = np.zeros((B.shape[0], 3, 4))
        second = np.zeros((B.shape[0], 3, 3, 4))
        derivative[:, :2] = np.einsum('plb,bf->plf', dB, c)
        compact_second = np.einsum('plb,bf->plf', ddB, c)
        second[:, 0, 0] = compact_second[:, 0]
        second[:, 0, 1] = second[:, 1, 0] = compact_second[:, 1]
        second[:, 1, 1] = compact_second[:, 2]
    else:
        derivative = np.einsum('plb,bf->plf', dB, c)
        second = np.einsum('plmb,bf->plmf', ddB, c)
    return values, derivative, second


def residual(data, matrices, coefficients, physical=False, pointwise=False):
    values, derivative, second = fields(matrices, coefficients)
    lam, W = values[:, 0], values[:, 1:]
    dl, ddl = derivative[:, :, 0], second[:, :, :, 0]
    psi = np.exp(lam)
    L, dL = longitudinal(data, W, derivative[:, :, 1:], second[:, :, :, 1:])
    Ahat, dAhat = data['A']+L, data['dA']+dL
    if not physical:
        dpsi = psi[:, None]*dl
        ddpsi = psi[:, None, None]*(ddl+np.einsum('pi,pj->pij', dl, dl))
        lap = np.einsum('pij,pij->p', data['gi'], ddpsi-np.einsum('pkij,pk->pij', data['conn'], dpsi))
        A2 = np.einsum('pij,pij->p', data['g']@Ahat@data['g'], Ahat)
        H = (lap-data['R']*psi/8+(A2+abs(data['Pi'])**2)*psi**-7/8-
             data['K']**2*psi**5/12+data['G2']*psi/8)
        M = divergence(data['conn'], L, dL)-(2/3)*(psi**6-1)[:, None]*data['gradK']+data['source']
        return np.column_stack((H, np.einsum('pij,pj->pi', data['g'], M)))
    # Reconstruct physical tensors and compute ADM constraints directly.
    g = psi[:, None, None]**4*data['g']
    dg = psi[:, None, None, None]**4*(data['dg']+4*dl[:, :, None, None]*data['g'][:, None])
    ddg = psi[:, None, None, None, None]**4*(data['ddg']+
        4*dl[:, :, None, None, None]*data['dg'][:, None]+
        4*dl[:, None, :, None, None]*data['dg'][:, :, None]+
        (4*ddl+16*np.einsum('pi,pj->pij', dl, dl))[:, :, :, None, None]*data['g'][:, None, None])
    gi, dgi, conn, _, R = metric_geometry(g, dg, ddg)
    Bcov = data['g']@Ahat@data['g']
    Kcov = psi[:, None, None]**-2*Bcov+data['K'][:, None, None]*g/3
    dKcov = np.empty_like(dg)
    for l in range(3):
        dBcov = data['dg'][:, l]@Ahat@data['g']+data['g']@dAhat[:, l]@data['g']+data['g']@Ahat@data['dg'][:, l]
        dKcov[:, l] = (psi[:, None, None]**-2*(dBcov-2*dl[:, l, None, None]*Bcov)+
                       (data['dK'][:, l, None, None]*g+data['K'][:, None, None]*dg[:, l])/3)
    Kup = gi@Kcov@gi
    dKup = np.stack([dgi[:, l]@Kcov@gi+gi@dKcov[:, l]@gi+gi@Kcov@dgi[:, l] for l in range(3)], axis=1)
    K = np.einsum('pij,pij->p', gi, Kcov)
    dK = np.einsum('plij,pij->pl', dgi, Kcov)+np.einsum('pij,plij->pl', gi, dKcov)
    Pi = psi**-6*data['Pi']
    rho2 = abs(Pi)**2+np.einsum('pij,pi,pj->p', gi, data['C'].conj(), data['C']).real
    H = R+K*K-np.einsum('pij,pij->p', Kcov, Kup)-rho2
    M = (divergence(conn, Kup, dKup)-np.einsum('pij,pj->pi', gi, dK)+
         np.einsum('pij,pj->pi', gi, (Pi.conj()[:, None]*data['C']).real))
    Mnorm = np.sqrt(np.maximum(0, np.einsum('pij,pi,pj->p', g, M, M)))
    if not all(np.all(np.isfinite(value)) for value in (H, Mnorm, psi, W, rho2)):
        raise ValueError('Nonfinite physical initial-data constraints or fields')
    result = dict(max_H=float(np.max(abs(H))), max_M=float(np.max(Mnorm)),
        min_psi=float(psi.min()), max_psi_minus_one=float(np.max(abs(psi-1))),
        max_W=float(np.max(abs(W))), max_energy=float(np.max(rho2/2)))
    return (result, H, Mnorm) if pointwise else result


def validation_points(radius=1., source_radius=1.):
    """Common off-grid cloud, including both sides of a compact-source cutoff."""
    large = source_radius < radius
    xf, yf, _ = scaled_quadrature(56 if large else 40, padding=10, radius=radius,
                                  source_radius=source_radius)
    theta = (np.arange(73)+.37)*2*np.pi/73
    radii = radius*np.array([0., 1e-8, 1e-5, .001, .1, .5, .99, 1.])
    if large:
        source_rings = source_radius*np.r_[.1, .25, .5, .75, .9, .95, .99,
            1-1e-5, 1-1e-8, 1., 1+1e-8, 1+1e-5, 1.01, 1.05, 1.1, 1.25]
        radii = np.unique(np.r_[radii, source_rings[source_rings <= radius]])
    xf = np.r_[xf, (radii[:, None]*np.cos(theta)).ravel()]
    yf = np.r_[yf, (radii[:, None]*np.sin(theta)).ravel()]
    return xf, yf


def validate_solution(degree, epsilon, coefficients, radius=1., source_radius=1., batch_size=512,
                      max_angular=None, taper_kind='polynomial'):
    """Reconstruct physical ADM tensors in bounded-memory off-grid batches."""
    xf, yf = validation_points(radius, source_radius)
    physical, uncorrected = None, None
    locations = {}
    for start in range(0, len(xf), batch_size):
        x, y = xf[start:start+batch_size], yf[start:start+batch_size]
        data = seed_data(x, y, epsilon, source_radius=source_radius, taper_kind=taper_kind)
        matrices = scaled_basis(degree, x, y, radius, max_angular=max_angular, compact=True)
        actual, H, Mnorm = residual(data, matrices, coefficients, physical=True, pointwise=True)
        vacuum_metric = residual(data, matrices, np.zeros_like(coefficients), physical=True)
        for key, samples in (('max_H', abs(H)), ('max_M', Mnorm)):
            if physical is None or actual[key] > physical[key]:
                i = int(np.argmax(samples))
                locations[key+'_location'] = [float(x[i]), float(y[i])]
        if physical is None:
            physical, uncorrected = actual, vacuum_metric
        else:
            for aggregate, values in ((physical, actual), (uncorrected, vacuum_metric)):
                for key, value in values.items():
                    aggregate[key] = min(aggregate[key], value) if key.startswith('min_') else max(aggregate[key], value)
    physical.update(locations)
    return physical, uncorrected, len(xf)


def projected_constraint_check(degree, epsilon, coefficients, padding=2, radius=1., source_radius=1.,
                               max_angular=None, taper_kind='polynomial', batch_size=256):
    """Check fixed coefficients on a chosen quadrature without a dense factorization.

    This is a quadrature check, not a replacement for direct off-grid ADM
    validation and not a new solve when the projected residual is unresolved.
    """
    x, y, weights = scaled_quadrature(degree, padding, radius, source_radius, max_angular)
    projected = np.zeros_like(np.asarray(coefficients, dtype=float)).reshape(-1, 4)
    for start in range(0, len(x), batch_size):
        xx, yy = x[start:start+batch_size], y[start:start+batch_size]
        data = seed_data(xx, yy, epsilon, source_radius=source_radius, taper_kind=taper_kind)
        matrices = scaled_basis(degree, xx, yy, radius, max_angular=max_angular, compact=True)
        weight = weights[start:start+batch_size]*np.sqrt(data['det'])
        projected += matrices[0].T@(weight[:, None]*residual(data, matrices, coefficients))
    if not np.all(np.isfinite(projected)):
        raise ValueError('Nonfinite projected initial-data constraints')
    return dict(max_projected_residual=float(np.max(abs(projected))), quadrature_points=len(x),
                quadrature_padding=padding, max_angular=max_angular, degree=degree)


def solve_case(degree, epsilon, padding=2, radius=1., source_radius=1., initial_coefficients=None,
               max_angular=None, taper_kind='polynomial', on_converged=None):
    """Solve one finite disk; a warm start must use the same angular cap/order.

    A shorter coefficient vector is zero-padded in total-degree ordering.
    Source radius and taper are independent of the outer Dirichlet boundary.
    An optional callback(coefficients, projected_error) runs after nonlinear
    convergence, before physical validation. Such a checkpoint is provisional.
    """
    x, y, weights = scaled_quadrature(degree, padding=padding, radius=radius,
                                     source_radius=source_radius, max_angular=max_angular)
    data = seed_data(x, y, epsilon, source_radius=source_radius, taper_kind=taper_kind)
    matrices = scaled_basis(degree, x, y, radius, max_angular=max_angular, compact=True)
    B, _, ddB = matrices
    weight = weights*np.sqrt(data['det'])
    calls = [0]
    def projected(coefficients):
        calls[0] += 1
        return (B.T@(weight[:, None]*residual(data, matrices, coefficients))).ravel()
    start = time.perf_counter()
    initial = np.zeros(4*B.shape[1])
    if initial_coefficients is not None:
        supplied = np.asarray(initial_coefficients, dtype=float).ravel()
        if supplied.size > initial.size or supplied.size % 4 or not np.all(np.isfinite(supplied)):
            raise ValueError('initial_coefficients must contain complete modes no larger than the requested basis')
        initial[:supplied.size] = supplied
    if np.max(abs(projected(initial))) < 1e-14:
        coefficients = initial
    else:
        lap0 = ddB[:, 0]+ddB[:, 2]
        lap0 *= weight[:, None]
        factor = lu_factor(B.T@lap0)
        del lap0
        pre = LinearOperator((4*B.shape[1], 4*B.shape[1]), matvec=lambda v:
            lu_solve(factor, v.reshape(B.shape[1], 4)).ravel(), dtype=float)
        coefficients = newton_krylov(projected, initial, inner_M=pre, f_tol=1e-12, maxiter=30)
    error = float(np.max(abs(projected(coefficients))))
    assert error < 2e-12
    if on_converged is not None:
        on_converged(coefficients.copy(), error)
    physical, uncorrected, points = validate_solution(degree, epsilon, coefficients, radius, source_radius,
                                                      max_angular=max_angular, taper_kind=taper_kind)
    result = dict(degree=degree, epsilon=epsilon, quadrature_padding=padding,
        radius=radius, source_radius=source_radius, max_angular=max_angular, taper_kind=taper_kind,
        scalar_taper=('max(1-s,0)^4; scalar C3 and normal momentum C2 at s=1' if taper_kind == 'polynomial'
                      else 'exp(-4*s/(1-s)) for s<1, else0; scalar and momentum C-infinity; distinct seed'),
        quadrature_points=len(x), unknowns=coefficients.size, projected_residual=error,
        residual_evaluations=calls[0], elapsed_seconds=time.perf_counter()-start,
        independent_points=points, physical_constraints=physical, uncorrected=uncorrected)
    return result, coefficients


def check_reconstruction():
    x, y = np.array([0., .2, -.6, .91]), np.array([0., -.3, .1, .2])
    data, matrices = seed_data(x, y, .05), basis(4, x, y)
    coefficients = 1e-4*np.cos(np.arange(4*matrices[0].shape[1])+.37)
    raw = residual(data, matrices, coefficients)
    psi = np.exp(fields(matrices, coefficients)[0][:, 0])
    direct = residual(data, matrices, coefficients, physical=True)
    expected_H = np.max(abs(-8*psi**-5*raw[:, 0]))
    expected_M = np.max(psi**-8*np.sqrt(np.einsum('pij,pi,pj->p', data['gi'], raw[:, 1:], raw[:, 1:])))
    assert abs(direct['max_H']-expected_H) < 1e-12
    assert abs(direct['max_M']-expected_M) < 1e-12


def main():
    check_reconstruction()
    runs, coefficients = [], None
    for degree in (6, 10, 14, 18, 22):
        result, coefficients = solve_case(degree, .03)
        runs.append(result)
        print(json.dumps(result), flush=True)
    errors = [max(r['physical_constraints'][key] for key in ('max_H', 'max_M')) for r in runs]
    checks = dict(refinement=errors[-1] < errors[0]/100, accuracy=errors[-1] < 1e-7,
                  positive_metric=all(r['physical_constraints']['min_psi'] > 0 for r in runs),
                  physical_reconstruction_identity=True)
    quadrature_runs = []
    for padding in (6, 10):
        result, padded = solve_case(22, .03, padding=padding)
        result['max_coefficient_change_from_padding2'] = float(np.max(abs(padded-coefficients)))
        quadrature_runs.append(result)
        print(json.dumps(result), flush=True)
    checks['quadrature_insensitive'] = all(r['max_coefficient_change_from_padding2'] < 1e-10
                                           for r in quadrature_runs)
    amplitude_runs, strong_coefficients, strong_degree = [], None, None
    for epsilon in (0., .003, .01, .1, .3, 1.):
        attempts = []
        for degree in (22, 26, 30, 34, 38):
            result, candidate = solve_case(degree, epsilon, padding=6)
            attempts.append(result)
            print(json.dumps(result), flush=True)
            if max(result['physical_constraints'][key] for key in ('max_H', 'max_M')) < 1e-7:
                break
        passed = max(result['physical_constraints'][key] for key in ('max_H', 'max_M')) < 1e-7
        amplitude_runs.append(dict(epsilon=epsilon, passed=passed, attempts=attempts))
        if epsilon == 1.:
            strong_coefficients, strong_degree = candidate, degree
    checks['finite_amplitudes'] = all(run['passed'] for run in amplitude_runs)
    vacuum = amplitude_runs[0]['attempts'][-1]['physical_constraints']
    checks['vacuum_recovery'] = max(vacuum[key] for key in ('max_H', 'max_M', 'max_psi_minus_one')) < 1e-12
    small = [run['attempts'][-1]['physical_constraints']['max_psi_minus_one']/run['epsilon']**2
             for run in amplitude_runs if run['epsilon'] in (.003, .01)]
    checks['small_amplitude_quadratic_response'] = abs(small[0]/small[1]-1) < .005
    strong_quad, strong_padded = solve_case(strong_degree, 1., padding=10)
    strong_quad['max_coefficient_change_from_padding6'] = float(np.max(abs(strong_padded-strong_coefficients)))
    checks['strong_quadrature_insensitive'] = (strong_quad['max_coefficient_change_from_padding6'] < 1e-8 and
        max(strong_quad['physical_constraints'][key] for key in ('max_H', 'max_M')) < 1e-7)
    print(json.dumps(strong_quad), flush=True)
    output = Path(__file__).resolve().parent
    np.savez_compressed(output/'ori_constraint_coefficients.npz', coefficients=coefficients,
                        degree=22, epsilon=.03, a=.02, e=1/16, initial_t=-1., radius=1., wave_number=1,
                        source_radius=1., max_angular=-1, taper_kind='polynomial',
                        scalar_taper='max(1-r^2/source_radius^2,0)^4',
                        strong_coefficients=strong_padded, strong_degree=strong_degree, strong_epsilon=1.)
    report = dict(scope='Finite-disk nonlinear curved Ori initial constraints; NO TIME EVOLUTION',
        convention='8*pi*G=c=1; Kij=-1/2 Lie_n gamma; Pi=n(phi); complex scalar denotes two real fields',
        method='Origin-regular Cartesian disk-polynomial Galerkin, independently differentiated physical ADM residual',
        boundary='psi=1 and W=0 on radius1; neither exterior matching nor a vacuum collar is established',
        source='Choquet-Bruhat, Isenberg, Pollack (2007), arXiv:gr-qc/0610045 equations4-13',
        refinement=runs, constraint_errors=errors, checks=checks,
        quadrature_refinement=quadrature_runs, amplitudes=amplitude_runs,
        strong_quadrature_check=strong_quad,
        validation='Common 3944-point cloud, including origin, shrinking rings and circular boundary; sampled maxima, not certified continuum sup bounds',
        previous_attempts='Square collocation failed independent corner checks; parity-only disk collocation failed origin checks; both superseded, not counted as verified data',
        remaining=['Well-posed curved-core evolution with justified boundary data',
                   'Exterior completion and physical preparation', 'Generic stability, quantum effects and human traversal'])
    (output/'ori_constraint_results.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    assert all(checks.values()), 'Independent constraints not yet resolved; inspect saved results.'
    print(json.dumps({'all_checks': checks, 'saved': str(output/'ori_constraint_coefficients.npz')}), flush=True)


if __name__ == '__main__':
    main()
