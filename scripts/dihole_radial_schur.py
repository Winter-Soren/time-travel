"""Exact finite added-mode accounting, evaluated in binary64 arithmetic.

The raw fixed-space Gram matrix G has delta G=0. For each PV shift j*M**2,
partition Q+j*M**2*G into A,V; V^*,H before any whitening. With
S=A-V H^{-1}V^*, delta S=a-v H^{-1}V^*-V H^{-1}v^*
+V H^{-1}d H^{-1}V^*, the added response is the sum of a direct H trace
and the retained Schur-feedback trace. All j*M**2*G_LH terms are included.

The normalization is fibre_response's -1/2 sum_j (1,-3,3,-1)_j trace;
there is NO frequency-integration factor 1/pi here. This is an algebraic
identity for these finite matrices, not an omitted-infinite-mode bound,
an assembly-error certificate, or a renormalized physical stress.
"""
import json
import math
from fractions import Fraction

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh

from dihole_bridge_enrichment import canonical
from dihole_sewn_form import fibre_response


def _hermitian(matrix):
    return (matrix + matrix.conj().T) / 2


def source_change(raw, kappa, keep, mass=1.):
    """Return direct-high + feedback = full minus retained at one frequency.

    raw={'norm':G, 'q':[C,B,A], 'derivative':[C1,B1,A1]}; keep is a
    nonempty proper subset of column indices, in any order. Inputs must
    be finite and exactly Hermitian as stored. Generalized Schur gaps
    use G_LL; the reported identity errors are numerical, not intervals.
    """
    if not np.isfinite(kappa) or not np.isfinite(mass) or mass <= 0:
        raise ValueError('Use finite real kappa and positive finite mass.')
    gram = np.asarray(raw['norm'], dtype=complex)
    n = len(gram)
    if gram.shape != (n, n) or n < 2:
        raise ValueError('norm must be square with at least two columns.')
    arrays = {'norm': [gram]}
    for key in ('q', 'derivative'):
        arrays[key] = [np.asarray(a, dtype=complex) for a in raw[key]]
        if len(arrays[key]) != 3:
            raise ValueError('q and derivative each need three coefficients.')
    for values in arrays.values():
        for matrix in values:
            if (matrix.shape != (n, n) or not np.isfinite(matrix).all()
                    or not np.array_equal(matrix, matrix.conj().T)):
                raise ValueError('All stored matrices must be finite and exactly Hermitian.')
    retained = np.asarray(keep)
    if (retained.ndim != 1 or retained.dtype.kind not in 'iu'
            or not 0 < len(retained) < n or len(np.unique(retained)) != len(retained)
            or np.any(retained < 0) or np.any(retained >= n)):
        raise ValueError('keep must contain distinct valid integer indices, with a nonempty complement.')
    retained = retained.astype(int)
    added = np.setdiff1d(np.arange(n), retained)
    checked_raw = {'norm': gram, 'q': arrays['q'], 'derivative': arrays['derivative']}
    q, dq = [sum(kappa**i * a for i, a in enumerate(arrays[key]))
             for key in ('q', 'derivative')]
    ix_l, ix_h, ix_lh = [np.ix_(*indices) for indices in
                         ((retained, retained), (added, added), (retained, added))]
    a, v, d = dq[ix_l], dq[ix_lh], dq[ix_h]
    records, direct_terms, feedback_terms = [], [], []
    for j, coefficient in enumerate((1., -3., 3., -1.)):
        shifted = q + (j * mass**2) * gram
        low, coupling, high = shifted[ix_l], shifted[ix_lh], shifted[ix_h]
        factor_h = cho_factor(high, lower=True, check_finite=False)
        x = cho_solve(factor_h, coupling.conj().T, check_finite=False)
        correction = coupling @ x
        schur = _hermitian(low - correction)
        factor_s = cho_factor(schur, lower=True, check_finite=False)
        factor_l = cho_factor(low, lower=True, check_finite=False)
        vx = v @ x
        delta_correction = -vx - vx.conj().T + x.conj().T @ d @ x
        # S^-1-A^-1=S^-1(V H^-1 V*) A^-1 avoids subtracting two traces.
        feedback_matrix = correction @ cho_solve(factor_l, a, check_finite=False) + delta_correction
        feedback_trace = np.trace(cho_solve(factor_s, feedback_matrix, check_finite=False)).real
        direct_trace = np.trace(cho_solve(factor_h, d, check_finite=False)).real
        direct_terms.append(-.5 * coefficient * float(direct_trace))
        feedback_terms.append(-.5 * coefficient * float(feedback_trace))
        spectrum = eigh(schur, gram[ix_l], eigvals_only=True, check_finite=False)
        records.append({'j': j, 'mass_squared': j * mass**2,
                        'direct_trace': float(direct_trace), 'feedback_trace': float(feedback_trace),
                        'schur_gap': float(spectrum[0]),
                        'schur_condition': float(spectrum[-1] / spectrum[0])})
    # Keep common-P cancellation inside the positive rational function for
    # the direct-high term; retain the alternating-shift check separately.
    direct_high = fibre_response(canonical(checked_raw, added), kappa, mass)
    feedback = math.fsum(feedback_terms)
    total = direct_high + feedback
    full = fibre_response(canonical(checked_raw, np.arange(n)), kappa, mass)
    low = fibre_response(canonical(checked_raw, retained), kappa, mass)
    independent = full - low
    return {'direct_high': direct_high, 'feedback': feedback, 'total_added': total,
            'independent_full_minus_retained': independent, 'error': abs(total-independent),
            'direct_high_shift_sum': math.fsum(direct_terms),
            'direct_high_shift_error': abs(math.fsum(direct_terms)-direct_high),
            'full_response': full, 'retained_response': low,
            'minimum_schur_gap': min(r['schur_gap'] for r in records),
            'maximum_schur_condition': max(r['schur_condition'] for r in records),
            'dimension': n, 'retained_dimension': len(retained), 'added_dimension': len(added),
            'gram_cross_frobenius': float(np.linalg.norm(gram[ix_lh])), 'shifts': records,
            'scope': 'Numerical identity for added finite columns; delta G=0; no 1/pi, infinite-tail, assembly, or physical-error bound.'}


def self_check():
    zero = np.zeros((2, 2))
    raw = {'norm': np.eye(2), 'q': [np.array([[2., .5], [.5, 3.]]), zero, np.eye(2)],
           'derivative': [np.array([[.25, .125], [.125, -.375]]), zero, np.diag([.125, .25])]}
    scalar = source_change(raw, .5, [0])
    exact_full, exact_low, exact_high = Fraction(0), Fraction(0), Fraction(0)
    for j, coefficient in enumerate((1, -3, 3, -1)):
        a, h, c = Fraction(9, 4)+j, Fraction(13, 4)+j, Fraction(1, 2)
        da, dh, dc = Fraction(9, 32), Fraction(-5, 16), Fraction(1, 8)
        weight = Fraction(-coefficient, 2)
        exact_full += weight*(h*da+a*dh-2*c*dc)/(a*h-c*c)
        exact_low += weight*da/a
        exact_high += weight*dh/h
    scalar_error = max(abs(scalar['total_added']-float(exact_full-exact_low)),
                       abs(scalar['direct_high']-float(exact_high)))

    rng = np.random.default_rng(140)
    z = rng.normal(size=(6, 6)) + 1j*rng.normal(size=(6, 6))
    transform = np.eye(6) + .12*z
    def pull(a):
        return _hermitian(transform.conj().T @ a @ transform)
    dense_raw = {'norm': pull(np.eye(6)),
                 'q': [pull(z.conj().T@z+2*np.eye(6)), pull(.1*_hermitian(z)), pull(np.eye(6))],
                 'derivative': [pull(.2*_hermitian(z)), pull(.03*np.eye(6)), pull(.01*_hermitian(z))]}
    keep, kappa, mass = np.array([4, 0, 2]), .7, 1.3
    dense = source_change(dense_raw, kappa, keep, mass)
    # Independently reconstruct the entire inverse from its block factors.
    added = np.setdiff1d(np.arange(6), keep)
    order = np.r_[keep, added]
    pencil = sum(kappa**j*a for j, a in enumerate(dense_raw['q'])) + 2*mass**2*dense_raw['norm']
    pencil = pencil[np.ix_(order, order)]
    a, v, h = pencil[:3, :3], pencil[:3, 3:], pencil[3:, 3:]
    x = np.linalg.solve(h, v.conj().T)
    sinv = np.linalg.inv(a-v@x)
    inverse = np.block([[sinv, -sinv@x.conj().T],
                        [-x@sinv, np.linalg.inv(h)+x@sinv@x.conj().T]])
    inverse_residual = float(np.linalg.norm(pencil@inverse-np.eye(6)))
    checks = {'exact_scalar_decomposition': scalar_error < 2e-16,
              'dense_nonorthogonal_full_minus_low': dense['error'] < 2e-14,
              'nonzero_mass_gram_cross_block': dense['gram_cross_frobenius'] > .1,
              'independent_block_inverse': inverse_residual < 2e-13,
              'positive_generalized_schur_gaps': min(scalar['minimum_schur_gap'], dense['minimum_schur_gap']) > 0}
    if not all(checks.values()):
        raise AssertionError({'checks': checks, 'scalar_error': scalar_error,
                              'dense_error': dense['error'], 'inverse_residual': inverse_residual})
    return {'checks': checks, 'scalar_exact_total': float(exact_full-exact_low),
            'scalar_error': scalar_error, 'dense_error': dense['error'],
            'block_inverse_residual': inverse_residual, 'dense_decomposition': dense}


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2, allow_nan=False))
