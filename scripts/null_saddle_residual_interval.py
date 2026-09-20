"""Interval residuals for fixed dyadic whole-plane Hermite trial spinors.

This module is independent of the Galerkin solver.  Its operator is exactly
H=-i AU*d_u-i AV*d_v+SY*(eta-2*u*v), with normalized Hermite width ``width``.
For real binary64 arrays a,b, the actual trial fields are psi=i*a and chi=i*b.
An isolated mpmath interval context encloses their FULL continuum coefficient
residuals, including the first omitted rows and their corner.  No simulation,
historical module, result file, or numerical solution is loaded or written.

The mathematical interpretation uses the self-adjoint whole-plane operator.
The exported enclosures assume the documented outward rounding of mpmath's
interval arithmetic; input floats define exact dyadic numbers, not uncertain
measurements or intended decimal values.  There is no angular/lead truncation
in the action on each finite Hermite trial.
"""
from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np
from mpmath.ctx_iv import MPIntervalContext


def _outward_bounds(value):
    """Export binary64 endpoints outwards, including conversion rounding."""
    low, high = float(value.a), float(value.b)
    if not math.isfinite(low) or not math.isfinite(high):
        raise OverflowError('Interval endpoint is outside finite binary64 range.')
    low = 0. if value.a == 0 else math.nextafter(low, -math.inf)
    high = 0. if value.b == 0 else math.nextafter(high, math.inf)
    return [low, high]


def _inputs(a, b, kappa, eta, width, digits):
    arrays = []
    for name, value in (('a', a), ('b', b)):
        raw = np.asarray(value)
        if np.iscomplexobj(raw):
            raise ValueError(f'{name} must be real; the trial spinor is i*{name}.')
        array = np.asarray(raw, dtype=np.float64)
        if (array.ndim != 3 or array.shape[0] != array.shape[1]
                or array.shape[0] < 1 or array.shape[2] != 2
                or not np.isfinite(array).all()):
            raise ValueError(f'{name} must be a finite real array of shape (N,N,2).')
        arrays.append(array)
    if arrays[0].shape != arrays[1].shape:
        raise ValueError('The two trial arrays must have identical shapes.')
    kappa, eta, width = map(float, (kappa, eta, width))
    if (not np.isfinite([kappa, eta, width]).all() or kappa <= 0 or width <= 0):
        raise ValueError('Require finite kappa>0, width>0, and eta.')
    if isinstance(digits, bool) or int(digits) != digits or digits < 20:
        raise ValueError('digits must be an integer >=20, sufficient for exact binary64 inputs.')
    return *arrays, kappa, eta, width, int(digits)


def _dyadic_array(array, context):
    """mpmath converts each finite binary64 number to its exact dyadic value."""
    return [[[context.mpf(float(value)) for value in spin] for spin in row]
            for row in array]


def _residual(coefficients, signed_kappa, eta, width, roots, context):
    """Return f+K*c+signed_kappa*c on indices 0..N, where H=i*K.

    The two passes apply the Hermite raising/lowering identities separately.
    Combining the Pauli 1/sqrt(2) with the derivative 1/sqrt(2) leaves the
    exact kinetic factor 1/(2*width); no rounded Pauli matrix is imported.
    """
    n = len(coefficients)
    zero = context.mpf(0)
    pair_zero = (zero, zero)
    xsum = [[pair_zero for _ in range(n)] for _ in range(n+1)]
    xdiff = [[pair_zero for _ in range(n)] for _ in range(n+1)]
    for p in range(n+1):
        for q in range(n):
            upper = ([roots[p+1]*v for v in coefficients[p+1][q]]
                     if p+1 < n else pair_zero)
            lower = ([roots[p]*v for v in coefficients[p-1][q]]
                     if p else pair_zero)
            xsum[p][q] = (upper[0]+lower[0], upper[1]+lower[1])
            xdiff[p][q] = (upper[0]-lower[0], upper[1]-lower[1])

    inverse_twice_width = 1/(2*width)
    width_squared = width**2
    result = []
    retained_squared, shell_squared = zero, zero
    for p in range(n+1):
        row = []
        for q in range(n+1):
            here = coefficients[p][q] if p < n and q < n else pair_zero
            du = xdiff[p][q] if q < n else pair_zero
            high_v = ([roots[q+1]*v for v in coefficients[p][q+1]]
                      if p < n and q+1 < n else pair_zero)
            low_v = ([roots[q]*v for v in coefficients[p][q-1]]
                     if p < n and q else pair_zero)
            dv = (high_v[0]-low_v[0], high_v[1]-low_v[1])
            high_uv = ([roots[q+1]*v for v in xsum[p][q+1]]
                       if q+1 < n else pair_zero)
            low_uv = ([roots[q]*v for v in xsum[p][q-1]]
                      if q else pair_zero)
            uv = (high_uv[0]+low_uv[0], high_uv[1]+low_uv[1])

            # K=-AU*d_u-AV*d_v+J*(eta-2uv), J=SY/i=[[0,-1],[1,0]].
            r0 = ((-du[0]-du[1]+dv[0]-dv[1])*inverse_twice_width
                  -eta*here[1]+width_squared*uv[1]+signed_kappa*here[0])
            r1 = ((-du[0]+du[1]-dv[0]-dv[1])*inverse_twice_width
                  +eta*here[0]-width_squared*uv[0]+signed_kappa*here[1])
            if p == 0 and q == 0:
                r0 += 1
            row.append((r0, r1))
            squared = r0**2+r1**2
            if p < n and q < n:
                retained_squared += squared
            else:
                shell_squared += squared
        result.append(row)
    return result, retained_squared, shell_squared


def enclose_residual(a, b, kappa, eta=1., width=1., digits=35):
    """Enclose residuals and a dual-corrected imaginary Green observable.

    rho=f-(H-i*kappa)psi=f+K*a-kappa*a;
    sigma=f-(H+i*kappa)chi=f+K*b+kappa*b.
    With <x,y>=x^dagger*y, Im(G_corrected)=a[0,0,0]-sum(b*rho).
    The exact observable lies in a complex disk about that corrected value
    with radius <=||rho||*||sigma||/kappa; its imaginary-part enclosure is
    also returned.  Retained solve defects are never discarded.
    """
    started = time.perf_counter()
    a, b, kappa, eta, width, digits = _inputs(a, b, kappa, eta, width, digits)
    context = MPIntervalContext()
    context.dps = digits
    aa, bb = _dyadic_array(a, context), _dyadic_array(b, context)
    kk, ee, ww = map(context.mpf, (kappa, eta, width))
    roots = [context.sqrt(j) for j in range(len(a)+1)]
    rho, rho_in, rho_out = _residual(aa, -kk, ee, ww, roots, context)
    _, sigma_in, sigma_out = _residual(bb, kk, ee, ww, roots, context)
    rho_norm = context.sqrt(rho_in+rho_out)
    sigma_norm = context.sqrt(sigma_in+sigma_out)
    correction = context.mpf(0)
    for p in range(len(a)):
        for q in range(len(a)):
            correction += bb[p][q][0]*rho[p][q][0]+bb[p][q][1]*rho[p][q][1]
    center = aa[0][0][0]-correction
    radius = rho_norm*sigma_norm/kk
    rho_bounds, sigma_bounds = _outward_bounds(rho_norm), _outward_bounds(sigma_norm)
    center_bounds, radius_bounds = _outward_bounds(center), _outward_bounds(radius)
    imaginary_bounds = [_outward_bounds(center-radius)[0], _outward_bounds(center+radius)[1]]
    return dict(
        N=len(a), kappa=kappa, eta=eta, width=width, digits=digits,
        residual_norm_upper=rho_bounds[1], adjoint_residual_norm_upper=sigma_bounds[1],
        solution_error_upper=_outward_bounds(rho_norm/kk)[1],
        G_disk_radius_upper=radius_bounds[1],
        corrected_imaginary_G_interval=center_bounds, imaginary_G_interval=imaginary_bounds,
        residual_norm_interval=rho_bounds, adjoint_residual_norm_interval=sigma_bounds,
        residual_retained_norm_interval=_outward_bounds(context.sqrt(rho_in)),
        residual_shell_norm_interval=_outward_bounds(context.sqrt(rho_out)),
        adjoint_retained_norm_interval=_outward_bounds(context.sqrt(sigma_in)),
        adjoint_shell_norm_interval=_outward_bounds(context.sqrt(sigma_out)),
        dual_correction_interval=_outward_bounds(correction),
        corrected_ImG_interval=center_bounds, dual_disk_radius_upper=radius_bounds[1],
        dual_disk_radius_interval=radius_bounds,
        ImG_interval=imaginary_bounds,
        input_interpretation='Finite binary64 trial arrays/parameters are exact dyadic numbers.',
        residual_definition='rho=f+Ka-kappa*a; sigma=f+Kb+kappa*b; H=iK; f000=1.',
        backend='Isolated mpmath.ctx_iv.MPIntervalContext; outward binary64 export.',
        scope='Full continuum Hermite action on fixed Schwartz trials; disk bound uses whole-plane self-adjointness.',
        seconds=time.perf_counter()-started)


def _direct_real_residual(c, signed_kappa, eta, width):
    """Independent binary64 reference: explicit shifts then Pauli matrices."""
    n = len(c)
    shape = (n+1, n+1, 2)
    du, dv, product = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    for p in range(n):
        for q in range(n):
            for step, sign in ((-1, 1.), (1, -1.)):
                if p+step >= 0:
                    du[p+step, q] += sign*math.sqrt(p+(step+1)//2)/(math.sqrt(2)*width)*c[p, q]
                if q+step >= 0:
                    dv[p, q+step] += sign*math.sqrt(q+(step+1)//2)/(math.sqrt(2)*width)*c[p, q]
            for sp in (-1, 1):
                for sq in (-1, 1):
                    if p+sp >= 0 and q+sq >= 0:
                        product[p+sp, q+sq] += (width**2/2*math.sqrt(p+(sp+1)//2)
                                               *math.sqrt(q+(sq+1)//2)*c[p, q])
    au = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    av = np.array([[-1., 1.], [1., 1.]])/math.sqrt(2)
    j = np.array([[0., -1.], [1., 0.]])
    padded = np.zeros(shape)
    padded[:n, :n] = c
    residual = -du@au.T-dv@av.T+(eta*padded-2*product)@j.T+signed_kappa*padded
    residual[0, 0, 0] += 1
    return residual


def self_check(benchmark=False):
    """Tiny fixed-vector/action checks; optional N32/64 interval timing only."""
    from mpmath import mp
    with mp.workdps(80):
        a, b = np.zeros((1, 1, 2)), np.zeros((1, 1, 2))
        a[0, 0, 0], b[0, 0, 0] = .125, -.25
        fixed = enclose_residual(a, b, .5)
        t, u, k = mp.mpf(1)/8, -mp.mpf(1)/4, mp.mpf(1)/2
        exact_rho = mp.sqrt((1-k*t)**2+3*t*t)
        exact_sigma = mp.sqrt((1+k*u)**2+3*u*u)
        exact_center = t-u*(1-k*t)
        for key, target in (('residual_norm_interval', exact_rho),
                            ('adjoint_residual_norm_interval', exact_sigma),
                            ('corrected_ImG_interval', exact_center),
                            ('dual_disk_radius_interval', exact_rho*exact_sigma/k)):
            assert mp.mpf(fixed[key][0]) <= target <= mp.mpf(fixed[key][1]), key

    rng = np.random.default_rng(116082)
    a = rng.integers(-64, 65, (4, 4, 2)).astype(float)/1024
    b = rng.integers(-64, 65, (4, 4, 2)).astype(float)/1024
    enclosed = enclose_residual(a, b, .75, eta=.5, width=1.25)
    rho = _direct_real_residual(a, -.75, .5, 1.25)
    sigma = _direct_real_residual(b, .75, .5, 1.25)
    numerical_center = a[0, 0, 0]-float(np.sum(b*rho[:4, :4]))
    errors = dict(primal_norm=abs(enclosed['residual_norm_upper']-np.linalg.norm(rho)),
                  adjoint_norm=abs(enclosed['adjoint_residual_norm_upper']-np.linalg.norm(sigma)),
                  corrected_center=abs(np.mean(enclosed['corrected_ImG_interval'])-numerical_center),
                  corner=abs(np.linalg.norm(rho[-1, -1])-1.25**2*4*np.linalg.norm(a[-1, -1])))
    assert max(errors.values()) < 2e-14, errors
    zero = enclose_residual(np.zeros((2, 2, 2)), np.zeros((2, 2, 2)), .5)
    assert zero['residual_shell_norm_interval'] == [0., 0.]
    assert zero['adjoint_shell_norm_interval'] == [0., 0.]
    assert zero['corrected_ImG_interval'] == [0., 0.]
    timings = []
    if benchmark:
        for n in (32, 64):
            aa = rng.integers(-64, 65, (n, n, 2)).astype(float)/1024
            bb = rng.integers(-64, 65, (n, n, 2)).astype(float)/1024
            row = enclose_residual(aa, bb, .5)
            timings.append(dict(N=n, seconds=row['seconds'],
                                projected_N256_seconds=row['seconds']*(256/n)**2))
    return dict(passed=True, fixed_N1=fixed, direct_action_errors=errors,
                zero_trial_control=True, timings=timings,
                note='No Galerkin solve; benchmarks use fixed dyadic vectors.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark', action='store_true', help='Time only N32/64 interval checks.')
    args = parser.parse_args()
    print(json.dumps(self_check(args.benchmark), indent=2, allow_nan=False))
