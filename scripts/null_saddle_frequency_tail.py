"""Frequency-dependent finite-Hermite continuation in a bilinear saddle arm.

The moving-basis matrices are reused from null_saddle_tail. The returned
zero-channel source and readout are anchored at the near plane R: both raw
far-cap amplitudes have exp(i*E*(Rfar-R)) divided out. The anchored ODEs are
integrated directly, without exponentially small/large numerical cancellation.

Real energies must lie below the first massive threshold at R. Upper-half-
plane energies use the analytic decaying cap branch; their anchored amplitudes
are not probabilities or flux-normalized scattering amplitudes. Finite Hermite
cutoff, finite cap, and overlapping whole-line arms remain model limitations.
This module writes no artifacts and runs no DG calculation.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.integrate import solve_ivp

from null_saddle_leads import frozen_lead
from null_saddle_tail import projected_matrices, tail_port


def cap_impedance(masses, energy):
    """Outgoing gamma=lambda/(E+i*sqrt(lambda**2-E**2)), Re(sqrt)>0.

    This equals (E-i*sqrt(lambda**2-E**2))/lambda while avoiding subtraction
    of nearly equal imaginary terms when E=i*kappa and kappa is large.
    """
    masses = np.asarray(masses, float)
    energy = complex(energy)
    if (not np.isfinite(masses).all() or np.any(masses <= 0)
            or not np.isfinite([energy.real, energy.imag]).all() or energy.imag < 0):
        raise ValueError('Require positive finite masses and finite Im(E)>=0.')
    if energy.imag == 0 and np.any(abs(energy.real) >= masses):
        raise ValueError('Real cap energy must be strictly below every massive threshold.')
    decay = np.sqrt(masses.astype(complex)**2-energy**2)
    if np.any(decay.real <= 0):
        raise ValueError('The decaying square-root branch has no positive real part.')
    return masses/(energy+1j*decay), decay


def coefficient_generator(lead, r, energy):
    """G for [a,b]'=G[a,b]; G†Sigma+SigmaG=-2*Im(E)*I."""
    kp, km, coupling = projected_matrices(lead, r)
    return np.block([[1j*energy*np.eye(len(kp))-kp, -1j*coupling],
                     [1j*coupling.T, -1j*energy*np.eye(len(km))-km]])


def _unpack(state, np_, nm):
    count = np_*nm
    return state[:count].reshape(nm, np_), state[count:count+nm], state[count+nm:]


def _cap_state(lead, Rfar, energy):
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    masses = lead.masses*math.sqrt(Rfar/lead.R)
    gamma, _ = cap_impedance(masses, energy)
    z = np.zeros((nm, np_), complex)
    z[lead.minus_pair_indices, lead.plus_pair_indices] = gamma
    source, output = np.zeros(nm, complex), np.zeros(np_, complex)
    if lead.incoming_zero_index is not None:
        source[lead.incoming_zero_index] = 1.
    if lead.outgoing_zero_index is not None:
        output[lead.outgoing_zero_index] = 1.
    return np.concatenate([z.ravel(), source, output])


def _derivative(state, np_, nm, kp, km, coupling, energy):
    z, source, output = _unpack(state, np_, nm)
    b = -1j*coupling
    dz = (1j*coupling.T-2j*energy*z-km@z+z@kp-z@b@z)
    # The scalar -iE terms in the raw source/readout ODEs are removed by
    # near-plane anchoring. Z itself is invariant under this normalization.
    dg = (-km-z@b)@source
    dt = output@(kp-b@z)
    return np.concatenate([dz.ravel(), dg, dt])


def _real_current_errors(lead, state):
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    z, source, output = _unpack(state, np_, nm)
    if lead.outgoing_zero_index is not None:
        return dict(outgoing_current=float(np.linalg.norm(np.eye(np_)-z.conj().T@z
                                                         -np.outer(output.conj(), output), 2)),
                    massive_coisometry=float(np.linalg.norm(z@z.conj().T-np.eye(nm), 2)),
                    outgoing_zero_norm=abs(float(np.vdot(output, output).real)-1.))
    return dict(massive_isometry=float(np.linalg.norm(z.conj().T@z-np.eye(np_), 2)),
                incoming_orthogonality=float(np.linalg.norm(z.conj().T@source)),
                incoming_zero_norm=abs(float(np.vdot(source, source).real)-1.))


def frequency_port(face, R, eta, energy, pairs=4, Rfar=12., rtol=2e-10, atol=None):
    """Return lead,Z,source,output,metadata in analytic near-face coordinates.

    ``lead`` is FrozenLead(...,energy=0) used ONLY for its energy-independent
    Hermite/spin basis. The port's actual energy is recorded in metadata.
    On v arms b=Z@a+source*incoming_zero; on u arms b=Z@a and output@a is
    the outgoing-zero readout. Unused source/output arrays are zero. Both
    zero amplitudes are referenced to R, not the arbitrary far cap.

    Apply the existing separate DG Gram transforms at R; no new pointwise
    normalization is implicit. Real-current identities are checked only for
    real subthreshold energies. Im(E)>0 checks homogeneous Z contraction.
    """
    lead = frozen_lead(face, R, eta, pairs=pairs, energy=0.)
    energy, Rfar, rtol = complex(energy), float(Rfar), float(rtol)
    if (not np.isfinite([energy.real, energy.imag, Rfar, rtol]).all()
            or energy.imag < 0 or Rfar < lead.R or not 1e-14 <= rtol < 1):
        raise ValueError('Require finite Im(E)>=0, Rfar>=R and 1e-14<=rtol<1.')
    if energy.imag == 0 and abs(energy.real) >= lead.masses[0]:
        raise ValueError('Real energy must lie below the first massive threshold at R.')
    atol = max(1e-14, .02*rtol) if atol is None else float(atol)
    if not math.isfinite(atol) or not 0 < atol < 1:
        raise ValueError('atol must be finite and between 0 and 1.')
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    cap = _cap_state(lead, Rfar, energy)

    def rhs(r, state):
        return _derivative(state, np_, nm, *projected_matrices(lead, r), energy)

    if Rfar == lead.R:
        samples, near, nfev, steps = cap[:, None], cap.copy(), 0, 0
    else:
        solution = solve_ivp(rhs, (Rfar, lead.R), cap, method='DOP853', rtol=rtol, atol=atol)
        if not solution.success or not np.isfinite(solution.y).all():
            raise RuntimeError(f'Frequency-arm continuation failed: {solution.message}')
        samples, near, nfev, steps = solution.y, solution.y[:, -1], solution.nfev, len(solution.t)-1
    z, source, output = _unpack(near, np_, nm)
    maximum_norm = max(float(np.linalg.norm(_unpack(samples[:, j], np_, nm)[0], 2))
                       for j in range(samples.shape[1]))
    contraction_excess = max(0., maximum_norm-1.)
    tolerance = max(2e-8, 200*rtol)
    real_errors, maximum_real_error = None, None
    if energy.imag == 0:
        real_errors = _real_current_errors(lead, near)
        maximum_real_error = max(max(_real_current_errors(lead, samples[:, j]).values())
                                 for j in range(samples.shape[1]))
        if maximum_real_error > tolerance:
            raise RuntimeError(f'Real-energy current check failed: {maximum_real_error} > {tolerance}')
    if contraction_excess > tolerance:
        raise RuntimeError(f'Homogeneous impedance contraction failed: {contraction_excess} > {tolerance}')
    metadata = dict(face=face, R=lead.R, Rfar=Rfar, eta=lead.eta, pairs=lead.pairs,
                    energy=dict(real=energy.real, imag=energy.imag), method='DOP853', rtol=rtol, atol=atol,
                    function_evaluations=int(nfev), accepted_steps=int(steps),
                    incoming_zero_index=lead.incoming_zero_index,
                    outgoing_zero_index=lead.outgoing_zero_index,
                    near_real_current_errors=real_errors, maximum_sampled_real_current_error=maximum_real_error,
                    maximum_sampled_Z_norm=maximum_norm, contraction_excess=contraction_excess,
                    current_acceptance_tolerance=tolerance,
                    anchored_source_norm=float(np.linalg.norm(source)),
                    anchored_output_norm=float(np.linalg.norm(output)),
                    anchor='Both raw cap zero amplitudes divided by exp(i*E*(Rfar-R)); anchored ODEs integrated directly.',
                    norm='Full transverse-line analytic basis; retain separate finite-face DG Gram transforms.',
                    cap='Frozen massive outgoing branch at Rfar; scalar sqrt has positive real part.',
                    direct_zero_input_output_term=False,
                    scope=['Finite Hermite cutoff and finite far cap; no omitted-mode or infinite-tail certificate.',
                           'Separate whole-line arms overlap near corners; not an exact decomposition of the bilinear plane.',
                           'Anchored upper-half-plane source/readout are not probabilities or flux-normalized amplitudes.',
                           'Only real subthreshold energies have the stated one-chiral-channel current identities.'])
    return dict(lead=lead, Z=z.copy(), source=source.copy(), output=output.copy(), metadata=metadata)


def _transfer_control(face, energy):
    """Short direct coefficient transport; does not use the Riccati equations."""
    near, far = 2.5, 2.7
    port = frequency_port(face, near, .7, energy, pairs=3, Rfar=far, rtol=2e-12)
    lead = port['lead']
    np_, nm = len(lead.plus_indices), len(lead.minus_indices)
    size = np_+nm
    def rhs(r, state):
        return (coefficient_generator(lead, r, energy)@state.reshape(size, size)).ravel()
    solution = solve_ivp(rhs, (near, far), np.eye(size, dtype=complex).ravel(),
                         method='DOP853', rtol=2e-12, atol=4e-14)
    if not solution.success:
        raise RuntimeError(solution.message)
    transfer = solution.y[:, -1].reshape(size, size)
    aa, ab = transfer[:np_, :np_], transfer[:np_, np_:]
    ba, bb = transfer[np_:, :np_], transfer[np_:, np_:]
    zcap, gcap, tcap = _unpack(_cap_state(lead, far, energy), np_, nm)
    denominator = bb-zcap@ab
    z = np.linalg.solve(denominator, zcap@aa-ba)
    # Independent reconstruction divides the raw amplitudes afterwards. The
    # short interval keeps this control well conditioned at complex energy.
    anchor = np.exp(-1j*energy*(far-near))
    source = anchor*np.linalg.solve(denominator, gcap)
    output = anchor*tcap@(aa+ab@z)
    return max(float(np.linalg.norm(port[key]-value))
               for key, value in (('Z', z), ('source', source), ('output', output)))


def self_check():
    """Bounded finite-matrix checks, including absorption and phase anchoring."""
    errors = dict(generator_absorption=0., cap_decay=0., frozen_rhs=0.,
                  zero_energy_reproduction=0., direct_transfer=0., real_current=0.,
                  contraction_excess=0., conjugation=0., imaginary_axis_structure=0.)
    energies = (0., .35, -.35, .4j, .3+.2j)
    for face in ('u+', 'v+'):
        lead = frozen_lead(face, 2.5, .7, pairs=3)
        np_, nm = len(lead.plus_indices), len(lead.minus_indices)
        signature = np.diag(np.r_[np.ones(np_), -np.ones(nm)])
        for energy in energies:
            generator = coefficient_generator(lead, 2.7, energy)
            identity = generator.conj().T@signature+signature@generator+2*complex(energy).imag*np.eye(np_+nm)
            errors['generator_absorption'] = max(errors['generator_absorption'], float(np.max(abs(identity))))
            gamma, decay = cap_impedance(lead.masses, energy)
            for mass, g, d in zip(lead.masses, gamma, decay):
                mode = np.array([1., g])
                residual = np.array([[1j*d, mass], [mass, -1j*d]])@mode-energy*mode
                errors['cap_decay'] = max(errors['cap_decay'], float(np.max(abs(residual))))
            _, _, coupling = projected_matrices(lead, lead.R)
            frozen = _derivative(_cap_state(lead, lead.R, energy), np_, nm,
                                  np.zeros((np_, np_)), np.zeros((nm, nm)), coupling, energy)
            errors['frozen_rhs'] = max(errors['frozen_rhs'], float(np.max(abs(frozen))))
            port = frequency_port(face, 2.5, .7, energy, pairs=3, Rfar=4.)
            errors['contraction_excess'] = max(errors['contraction_excess'], port['metadata']['contraction_excess'])
            if complex(energy).imag == 0:
                errors['real_current'] = max(errors['real_current'], port['metadata']['maximum_sampled_real_current_error'])
            errors['direct_transfer'] = max(errors['direct_transfer'], _transfer_control(face, energy))
        zero, old = frequency_port(face, 2.5, .7, 0., pairs=3, Rfar=8.), tail_port(face, 2.5, .7, pairs=3, Rfar=8.)
        for key in ('Z', 'source', 'output'):
            errors['zero_energy_reproduction'] = max(errors['zero_energy_reproduction'], float(np.linalg.norm(zero[key]-old[key])))
        positive = frequency_port(face, 2.5, .7, .3+.2j, pairs=3, Rfar=4.)
        partner = frequency_port(face, 2.5, .7, -.3+.2j, pairs=3, Rfar=4.)
        for key in ('Z', 'source', 'output'):
            sign = -1 if key == 'Z' else 1
            errors['conjugation'] = max(errors['conjugation'],
                float(np.linalg.norm(partner[key]-sign*positive[key].conj())))
        imaginary = frequency_port(face, 2.5, .7, .4j, pairs=3, Rfar=4.)
        errors['imaginary_axis_structure'] = max(errors['imaginary_axis_structure'],
            float(np.max(abs(imaginary['Z'].real))), float(np.max(abs(imaginary['source'].imag))),
            float(np.max(abs(imaginary['output'].imag))))
    for name, value in errors.items():
        if value > 2e-8:
            raise AssertionError(f'{name}: {value} > 2e-8')
    return dict(passed=True, maximum_errors=errors,
                scope='Finite projected arms; real current and complex absorption are distinct controls.')


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2, allow_nan=False))
