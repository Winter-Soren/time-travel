"""Frequency tails of exact binary64 Hermitian matrix polynomials only.

Inputs are already in one common orthonormal basis, including the fixed-space
variation. Assembly/whitening, continuum, exterior and regulator-limit errors
are NOT enclosed. Interval arithmetic inherits dihole_interval_arithmetic's
IEEE754 contract. Returned tails include both signs with measure dκ/(2π).
"""
import json
import time
from fractions import Fraction

import numpy as np
from scipy.linalg import eigh

import dihole_interval_arithmetic as ia


# Adjacent dyadics enclosing pi; the self-check proves this using Machin's series.
PI = tuple(np.asarray(float.fromhex(s)) for s in
           ('0x1.921fb54442d18p+1', '0x1.921fb54442d19p+1'))
p = ia.point


def _gersh(matrix):
    lo, hi = matrix
    radius = np.maximum(abs(lo), abs(hi)).copy()
    np.fill_diagonal(radius, 0.)
    radius = ia.sum_interval(p(radius), axis=1)
    return (float(np.min(ia.sub(p(np.diag(lo)), radius)[0])),
            float(np.max(ia.add(p(np.diag(hi)), radius)[1])))


def _lower(matrix):
    real = (matrix.real if not np.any(matrix.imag) else
            np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]]))
    _, Q = eigh(real, check_finite=False)  # A trial basis, never an eigenvalue proof.
    gram = _gersh(ia.matmul(p(Q.T), p(Q)))
    if gram[0] <= 0:
        raise ArithmeticError('Trial basis invertibility was not certified.')
    L = _gersh(ia.matmul(p(Q.T), ia.matmul(p(real), p(Q))))[0]
    # Rayleigh quotient: negative L requires the LOWER Gram endpoint.
    bound = float(ia.div(p(L), p(gram[1] if L >= 0 else gram[0]))[0])
    return bound, dict(congruence_lower=L, gram_lower=gram[0], gram_upper=gram[1],
                       real_dimension=len(real))


def _frobenius(matrix):
    if not np.any(matrix):
        return 0.
    return float(ia.sqrt(ia.sum_interval(ia.add(ia.square(p(matrix.real)),
                                               ia.square(p(matrix.imag)))))[1])


def coefficient_bounds(polynomial):
    """q=[C,B,A], derivative=[C1,B1,A1]; no silent Hermitian symmetrization."""
    matrices = [np.asarray(x) for key in ('q', 'derivative') for x in polynomial[key]]
    if len(polynomial['q']) != 3 or len(polynomial['derivative']) != 3:
        raise ValueError('Each polynomial must have three ascending coefficients.')
    shape = matrices[0].shape
    for x in matrices:
        if (x.dtype.kind not in 'fc' or x.dtype.itemsize > (16 if x.dtype.kind == 'c' else 8)
                or x.ndim != 2 or x.shape != shape or shape[0] != shape[1] or not shape[0]
                or not np.isfinite(x).all() or not np.array_equal(x, x.conj().T)):
            raise ValueError('Require finite, exactly Hermitian binary64-compatible square arrays.')
    C, B, A, C1, B1, A1 = [x.astype(np.complex128) for x in matrices]
    a, pa = _lower(A); c0, pc = _lower(C)
    return dict(a=a, c0=c0, b=_frobenius(B), c=max(0., -c0),
                d0=_frobenius(C1), d1=_frobenius(B1), d2=_frobenius(A1), n=len(A),
                eigenvalue_proofs=dict(A=pa, C=pc),
                scope='Exact supplied dyadic matrices; no assembly or whitening enclosure.')


def _power(x, n):
    result = p(1.)
    for _ in range(n):
        result = ia.mul(result, x)
    return result


def _inputs(bounds, K, M):
    if not (np.isfinite(K) and np.isfinite(M) and K > 0 and M > 0 and bounds['a'] > 0):
        raise ValueError('Require finite K,M>0 and a certified positive a.')
    if any(not np.isfinite(bounds[k]) or bounds[k] < 0 for k in ('b','c','d0','d1','d2')):
        raise ValueError('Norm and negative-part bounds must be finite and nonnegative.')
    return p(float(K)), p(float(M))


def high_tail(bounds, K, M=1):
    """Upper bound on ∫_|κ|>=K |Fκ| dκ/(2π), Fκ=-3M^6 Tr(deltaP/product)."""
    k, mass = _inputs(bounds, K, M)
    chi = ia.sub(ia.sub(p(bounds['a']), ia.div(p(bounds['b']), k)),
                 ia.div(p(bounds['c']), ia.square(k)))
    if chi[0] <= 0:
        raise ValueError('The chosen cutoff does not certify a positive quadratic gap.')
    terms = [ia.div(p(bounds['d'+str(i)]), ia.mul(p(float(7-i)), _power(k, 7-i)))
             for i in range(3)]
    total = ia.add(ia.add(terms[0], terms[1]), terms[2])
    prefactor = ia.div(ia.mul(p(3.*bounds['n']), _power(mass, 6)), ia.mul(PI, _power(chi, 4)))
    return float(ia.mul(prefactor, total)[1])


def low_tail(bounds, K, M=1):
    """Upper bound on ∫_|κ|<=K |Fκ| dκ/(2π) if the coefficient gap is positive."""
    k, mass = _inputs(bounds, K, M)
    ratio = ia.div(p(bounds['b']), ia.mul(p(2.), p(bounds['a'])))
    x = tuple(np.maximum(0., np.minimum(float(K), endpoint)) for endpoint in ratio)
    ell = ia.sub(ia.add(p(bounds['c0']), ia.mul(p(bounds['a']), ia.square(x))), ia.mul(p(bounds['b']), x))
    if ell[0] <= 0:
        raise ValueError('No positive low-frequency gap was certified; this is inconclusive.')
    denom = p(1.)
    for j in range(4):
        denom = ia.mul(denom, ia.add(ell, ia.mul(p(float(j)), ia.square(mass))))
    terms = [ia.div(ia.mul(p(bounds['d'+str(i)]), _power(k, i+1)), p(float(i+1))) for i in range(3)]
    prefactor = ia.div(ia.mul(p(3.*bounds['n']), _power(mass, 6)), ia.mul(PI, denom))
    return float(ia.mul(prefactor, ia.add(ia.add(terms[0], terms[1]), terms[2]))[1])


def self_check():
    from mpmath import mp
    started = time.perf_counter()
    def atan_bracket(q):
        s = sum((Fraction((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(20)), Fraction())
        return s, s+Fraction(1, 41*q**41)
    a, b = atan_bracket(5), atan_bracket(239)
    assert Fraction(float(PI[0])) < 16*a[0]-4*b[1] < 16*a[1]-4*b[0] < Fraction(float(PI[1]))
    with mp.workdps(75):
        one, zero = np.ones((1,1)), np.zeros((1,1))
        scalar = coefficient_bounds(dict(q=[one,zero,one], derivative=[zero,zero,one]))
        f = lambda k: 3*k*k/mp.fprod(k*k+j for j in (1,2,3,4))/mp.pi
        actual_high, actual_low = mp.quad(f,[2,mp.inf]), mp.quad(f,[0,mp.mpf('.25')])
        assert actual_high <= high_tail(scalar,2.) and actual_low <= low_tail(scalar,.25)
        A=np.array([[2.,.25+.125j],[.25-.125j,1.5]])
        B=np.array([[.25,.25j],[-.25j,-.125]])
        C=np.array([[1.,.25+.125j],[.25-.125j,-.25]])
        dense=coefficient_bounds(dict(q=[C,B,A],derivative=[B,C,A/8]))
        conv=lambda X:mp.matrix([[mp.mpc(float(z.real),float(z.imag)) for z in row] for row in X])
        aa,bb,cc=map(conv,(A,B,C)); residual=mp.mpf(0); ratios=[]
        assert dense['a'] <= mp.eighe(aa,eigvals_only=True)[0] and dense['c0'] <= mp.eighe(cc,eigvals_only=True)[0] < 0
        assert mp.norm(aa*bb-bb*aa)>0
        for key, matrix in (('b',bb),('d0',bb),('d1',cc),('d2',aa/8)):
            assert max(abs(v) for v in mp.eighe(matrix,eigvals_only=True)) <= dense[key]
        for k in (4,8,16):
            P=aa*k*k+bb*k+cc; inverse=P**-1
            residual=max(residual,mp.norm(P*inverse-mp.eye(2)))
            gap=dense['a']*k*k-dense['b']*k-dense['c']
            assert gap>0 and mp.mpf(gap) <= mp.eighe(P,eigvals_only=True)[0]
            kernel=mp.eye(2)
            for j in range(4):kernel=kernel*(P+j*mp.eye(2))**-1
            response=3*abs(mp.fsum(((aa*k*k/8+cc*k+bb)*kernel)[i,i] for i in range(2)))
            estimate=6*(dense['d2']*k*k+dense['d1']*k+dense['d0'])/mp.mpf(gap)**4
            ratios.append(float(response/estimate)); assert response<=estimate
        assert residual<mp.mpf('1e-70') and high_tail(dense,4.)>0
    return dict(passed=True,scalar_high=float(actual_high),scalar_high_upper=high_tail(scalar,2.),
                scalar_low=float(actual_low),scalar_low_upper=low_tail(scalar,.25),
                dense_bounds=dense,dense_pointwise_bound_ratios=ratios,
                dense_inverse_residual=str(residual),seconds=time.perf_counter()-started)


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2))
