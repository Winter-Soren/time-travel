"""Finite-Hermite continuation in one actual bilinear null-saddle arm.

The normal coordinate r increases outwards; tau follows the CCW tangent.
The mass is eta+2*sign*r*tau, with sign=-1 on u arms and +1 on v arms.
The moving Hermite basis retains its full projected width/center connection.
At Rfar the arm is closed by a frozen linear-mass lead, then its impedance,
incoming source, and outgoing zero-channel functional are evolved inwards.

This is an E=0, finite-channel, finite-Rfar model on the whole transverse
line. The four separate arms overlap outside the square near its corners;
their independent use is not an exact decomposition of the infinite plane.
No result files or figures are written by this module.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_hermite

from null_saddle_leads import NORMALS, FrozenLead, frozen_lead


def hermite_connection(count, r, eta, slope_sign):
    """Return K_ij=<phi_i,d_r phi_j>, i,j=0,...,count-1, in L2(d tau).

    omega=2*r, tau0=-eta/(2*sign*r). The first band is the width change;
    the second is translation of tau0. Truncation preserves skew symmetry.
    """
    if (isinstance(count, bool) or int(count) != count or count < 1
            or not np.isfinite([r, eta]).all() or r <= 0
            or slope_sign not in (-1, 1)):
        raise ValueError('Require positive integer count, r>0, finite eta and sign +/-1.')
    count = int(count)
    connection = np.zeros((count, count))
    shift = -eta/(2*slope_sign*r**1.5)
    for j in range(count):
        if j >= 2:
            connection[j-2, j] = math.sqrt(j*(j-1))/(4*r)
        if j+2 < count:
            connection[j+2, j] = -math.sqrt((j+1)*(j+2))/(4*r)
        if j:
            connection[j-1, j] = shift*math.sqrt(j)
        if j+1 < count:
            connection[j+1, j] = -shift*math.sqrt(j+1)
    return connection


def projected_matrices(lead, r):
    """Return Kplus, Kminus, Lambda in the existing FrozenLead ordering."""
    connection = hermite_connection(lead.pairs+1, r, lead.eta, int(np.sign(lead.slope)))
    kp = connection[np.ix_(lead.plus_indices, lead.plus_indices)]
    km = connection[np.ix_(lead.minus_indices, lead.minus_indices)]
    coupling = np.zeros((len(lead.plus_indices), len(lead.minus_indices)))
    coupling[lead.plus_pair_indices, lead.minus_pair_indices] = lead.masses*math.sqrt(r/lead.R)
    return kp, km, coupling


def coefficient_generator(lead, r):
    """Complex a/b generator at E=0; preserves Sigma=diag(Iplus,-Iminus)."""
    kp, km, coupling = projected_matrices(lead, r)
    return np.block([[-kp, -1j*coupling], [1j*coupling.T, -km]])


def _unpack(state, np_, nm):
    nmatrix = np_*nm
    return (state[:nmatrix].reshape(nm, np_),
            state[nmatrix:nmatrix+nm], state[nmatrix+nm:])


def _cap_state(lead):
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    y = np.zeros((nm, np_))
    y[lead.minus_pair_indices, lead.plus_pair_indices] = 1.
    source, output = np.zeros(nm), np.zeros(np_)
    if lead.incoming_zero_index is not None:
        source[lead.incoming_zero_index] = 1.
    if lead.outgoing_zero_index is not None:
        output[lead.outgoing_zero_index] = 1.
    return np.concatenate([y.ravel(), source, output])


def _derivative(state, np_, nm, kp, km, coupling):
    y, source, output = _unpack(state, np_, nm)
    # Z=-iY. Source is the physical incoming-zero amplitude at the far cap;
    # output @ a(r) is the physical outgoing-zero amplitude at that cap.
    dy = -coupling.T-km@y+y@kp+y@coupling@y
    dg = (-km+y@coupling)@source
    dt = output@(kp+coupling@y)
    return np.concatenate([dy.ravel(), dg, dt])


def _current_errors(lead, state):
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    y, source, output = _unpack(state, np_, nm)
    if lead.outgoing_zero_index is not None:
        return dict(outgoing_current=float(np.linalg.norm(np.eye(np_)-y.T@y-np.outer(output, output), 2)),
                    massive_coisometry=float(np.linalg.norm(y@y.T-np.eye(nm), 2)),
                    outgoing_zero_norm=abs(float(output@output)-1.))
    return dict(massive_isometry=float(np.linalg.norm(y.T@y-np.eye(np_), 2)),
                incoming_orthogonality=float(np.linalg.norm(y.T@source)),
                incoming_zero_norm=abs(float(source@source)-1.))


def tail_port(face, R, eta, pairs=3, Rfar=12., rtol=2e-10, atol=None):
    """Return a matched analytic-coordinate port for one bilinear arm at E=0.

    Keys: ``lead`` (FrozenLead at R), ``Z`` (Nm,Np), ``source`` (Nm,),
    ``output`` (Np,), and JSON-safe ``metadata``. On v arms,
    b=Z@a+source*incoming_zero. On u arms, b=Z@a and the far outgoing-zero
    amplitude is output@a. The unused source/output is an explicit zero array.
    No same-arm zero input-to-output direct term exists in this chiral model.

    DG callers must apply their distinct finite-face Gram transforms:
    Znorm=Rminus@Z@inv(Rplus), driven source=Rminus@source, and
    far output=output@solve(Rplus, outgoing_normalized_coefficients).
    """
    lead = frozen_lead(face, R, eta, pairs=pairs, energy=0.)
    if (not math.isfinite(float(Rfar)) or Rfar < R or not math.isfinite(float(rtol))
            or not 1e-14 <= rtol < 1):
        raise ValueError('Require finite Rfar>=R and 1e-14<=rtol<1.')
    atol = max(1e-14, .02*rtol) if atol is None else float(atol)
    if not math.isfinite(atol) or not 0 < atol < 1:
        raise ValueError('atol must be finite and between0 and1.')
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    cap = _cap_state(lead)

    def rhs(r, state):
        return _derivative(state, np_, nm, *projected_matrices(lead, r))

    if Rfar == R:
        samples = cap[:, None]
        near = cap.copy()
        nfev, steps = 0, 0
    else:
        solution = solve_ivp(rhs, (float(Rfar), float(R)), cap, method='DOP853',
                             rtol=float(rtol), atol=atol)
        if not solution.success or not np.isfinite(solution.y).all():
            raise RuntimeError(f'Bilinear-arm continuation failed: {solution.message}')
        samples, near = solution.y, solution.y[:, -1]
        nfev, steps = solution.nfev, len(solution.t)-1
    y, source, output = _unpack(near, np_, nm)
    errors = _current_errors(lead, near)
    along_error = max(max(_current_errors(lead, samples[:, j]).values())
                      for j in range(samples.shape[1]))
    # This is an integrator acceptance diagnostic, not a mathematical enclosure.
    tolerance = max(2e-8, 200*rtol)
    if along_error > tolerance:
        raise RuntimeError(f'Finite-channel current check failed: {along_error} > {tolerance}')
    yf, gf, tf = _unpack(cap, np_, nm)
    metadata = dict(face=face, R=float(R), Rfar=float(Rfar), eta=float(eta),
                    pairs=int(pairs), energy=0., method='DOP853', rtol=float(rtol), atol=atol,
                    function_evaluations=int(nfev), accepted_steps=int(steps),
                    incoming_zero_index=lead.incoming_zero_index,
                    outgoing_zero_index=lead.outgoing_zero_index,
                    near_current_errors=errors, maximum_sampled_current_error=float(along_error),
                    current_acceptance_tolerance=tolerance,
                    Z_change_from_frozen=float(np.linalg.norm(y-yf)),
                    source_change_from_bare_zero=float(np.linalg.norm(source-gf)),
                    output_change_from_bare_zero=float(np.linalg.norm(output-tf)),
                    cap='Frozen linear-mass lead at Rfar; E=0 outgoing prescription.',
                    norm='Full transverse-line L2(d tau), analytic FrozenLead phases/order.',
                    direct_zero_input_output_term=False,
                    scope=['Actual bilinear arm coefficients and exact projected moving-basis connection.',
                           'Finite Hermite cutoff and finite Rfar; no omitted-channel or infinite-tail bound.',
                           'Separate full-line arms overlap near exterior corners; no global plane decomposition certificate.',
                           'Current identities apply before restriction to finite faces; preserve both DG Gram transforms.'])
    return dict(lead=lead, Z=-1j*y.copy(), source=source.astype(complex),
                output=output.astype(complex), metadata=metadata)


def self_check():
    """Small independent differential/current/frozen-limit controls; no DG."""
    errors = dict(connection_finite_difference=0., full_derivative_reconstruction=0.,
                  connection_skew=0., generator_current=0., frozen_rhs=0.,
                  zero_length_port=0., continued_current=0., tolerance_refinement=0.)
    nodes, weights = roots_hermite(36)
    for face in NORMALS:
        for eta in (-.7, 0., .7):
            r, maximum_column, count = 2.5, 3, 6
            lead = frozen_lead(face, r, eta, pairs=count, energy=0.)
            part = 0 if lead.slope < 0 else 1
            tau = lead.center+nodes/math.sqrt(lead.frequency)
            points = r*lead.normal+tau[:, None]*lead.tangent
            f = lead.basis_function(points)[part][:, :, :count]
            quadrature = weights*np.exp(nodes**2)/math.sqrt(lead.frequency)
            step = r*1e-4
            shifted = [frozen_lead(face, r+j*step, eta, pairs=count).basis_function(points)[part]
                       for j in (-2, -1, 1, 2)]
            df = (shifted[0]-8*shifted[1]+8*shifted[2]-shifted[3])/(12*step)
            df = df[:, :, :maximum_column+1]
            numerical = np.einsum('sni,n,snj->ij', f.conj(), quadrature, df)
            exact = hermite_connection(count, r, eta, int(np.sign(lead.slope)))
            errors['connection_finite_difference'] = max(errors['connection_finite_difference'],
                float(np.max(abs(numerical-exact[:, :maximum_column+1]))))
            reconstruction = np.einsum('sni,ij->snj', f, exact[:, :maximum_column+1])
            errors['full_derivative_reconstruction'] = max(errors['full_derivative_reconstruction'],
                float(np.sqrt(np.sum(quadrature[None, :, None]*abs(df-reconstruction)**2))))
            errors['connection_skew'] = max(errors['connection_skew'], float(np.max(abs(exact+exact.T))))
            generator = coefficient_generator(lead, r)
            np_, nm = len(lead.plus_indices), len(lead.minus_indices)
            signature = np.diag(np.r_[np.ones(np_), -np.ones(nm)])
            errors['generator_current'] = max(errors['generator_current'],
                float(np.max(abs(generator.conj().T@signature+signature@generator))))
            _, _, coupling = projected_matrices(lead, r)
            frozen = _derivative(_cap_state(lead), np_, nm, np.zeros((np_, np_)),
                                  np.zeros((nm, nm)), coupling)
            errors['frozen_rhs'] = max(errors['frozen_rhs'], float(np.max(abs(frozen))))
            near_cap = tail_port(face, r, eta, pairs=3, Rfar=r)
            errors['zero_length_port'] = max(errors['zero_length_port'],
                                             float(np.max(abs(near_cap['Z']-near_cap['lead'].Z))))
    examples = []
    for face in ('u+', 'v+'):
        base = tail_port(face, 2.5, .7, pairs=3, Rfar=8.)
        refined = tail_port(face, 2.5, .7, pairs=3, Rfar=8., rtol=2e-12)
        errors['continued_current'] = max(errors['continued_current'],
                                          base['metadata']['maximum_sampled_current_error'])
        for key in ('Z', 'source', 'output'):
            errors['tolerance_refinement'] = max(errors['tolerance_refinement'],
                                                  float(np.linalg.norm(base[key]-refined[key])))
        examples.append(base['metadata'])
    for key, value in errors.items():
        tolerance = 2e-9 if key in ('connection_finite_difference', 'full_derivative_reconstruction') else 2e-8
        if value > tolerance:
            raise AssertionError(f'{key}: {value} > {tolerance}')
    return dict(passed=True, maximum_errors=errors, examples=examples,
                scope='Finite projected arm controls only; no full-plane scattering certificate.')


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2, allow_nan=False))
