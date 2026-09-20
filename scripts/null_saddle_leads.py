"""Exact frozen-face Hermite channels for the local null-saddle Dirac model.

The interior model is H=-i AU*d_u-i AV*d_v+SY*(eta-2*u*v). Each
axis-aligned box face is continued by freezing its normal coordinate. This
produces a translation-invariant, linear-mass lead, not the true bilinear
exterior. Channels are orthonormal on the WHOLE transverse line. Restricting
them to a finite face changes their two characteristic Gram matrices; callers
must retain those matrices and must not presume evanescent flux cancellation.

No result files are written. Run this module for small independent checks.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.special import roots_hermite


SX = np.array([[0., 1.], [1., 0.]], complex)
SY = np.array([[0., -1j], [1j, 0.]], complex)
SZ = np.diag([1., -1.]).astype(complex)
AU, AV = (SX+SZ)/math.sqrt(2), (SX-SZ)/math.sqrt(2)
NORMALS = {'u+': (1., 0.), 'u-': (-1., 0.),
           'v+': (0., 1.), 'v-': (0., -1.)}


def _hermite_functions(x, maximum):
    """Normalized oscillator functions in dimensionless x, by recurrence."""
    values = np.empty((len(x), maximum+1))
    values[:, 0] = np.exp(-x*x/2)/math.pi**.25
    if maximum:
        values[:, 1] = math.sqrt(2)*x*values[:, 0]
    for n in range(2, maximum+1):
        values[:, n] = (math.sqrt(2/n)*x*values[:, n-1]
                        -math.sqrt((n-1)/n)*values[:, n-2])
    return values


class FrozenLead:
    """Analytic full-line channels and an outward-decaying massive lead map.

    ``basis_function(points)`` takes Cartesian points of shape (npoints, 2)
    and returns (plus, minus), each shaped (2, npoints, nchannels). Plus/minus
    are the +1/-1 eigenspaces of the OUTWARD normal Pauli matrix.

    In analytic channel coordinates, incoming = Z @ outgoing. On v faces,
    add ``incoming_zero_amplitude * e_zero`` to this relation to prescribe the
    one incoming propagating channel. On u faces the outgoing zero is free.
    At E=0 the zero lead is understood in the outgoing/limiting-absorption
    sense; its longitudinal profile is constant, not square-integrable.
    """

    def __init__(self, face, R, eta, pairs=2, energy=0.):
        if face not in NORMALS:
            raise ValueError('face must be u+, u-, v+, or v-')
        if isinstance(pairs, bool) or int(pairs) != pairs or pairs < 1:
            raise ValueError('pairs must be a positive integer')
        if any(not math.isfinite(float(x)) for x in (R, eta, energy)) or R <= 0:
            raise ValueError('R must be positive; R, eta and real energy finite')
        self.face, self.R, self.eta = face, float(R), float(eta)
        self.pairs, self.energy = int(pairs), float(energy)
        self.normal = np.array(NORMALS[face])
        self.tangent = np.array([-self.normal[1], self.normal[0]])
        self.sigma_n = self.normal[0]*AU+self.normal[1]*AV
        self.sigma_t = self.tangent[0]*AU+self.tangent[1]*AV
        self.slope = -2*self.R*(self.normal[0]**2-self.normal[1]**2)
        self.center = -self.eta/self.slope
        self.frequency = abs(self.slope)
        self.masses = np.sqrt(2*self.frequency*np.arange(1, self.pairs+1))
        if abs(self.energy) >= self.masses[0]:
            raise ValueError('real energy must be strictly below the first massive threshold')

        # Fix an absolute phase, then use the oriented transverse frame to fix
        # the relative phase so H_transverse swaps each pair with +lambda_j.
        _, spinors = np.linalg.eigh(self.sigma_n.real)
        self.plus_spinor = spinors[:, 1].astype(complex)
        pivot = np.argmax(abs(self.plus_spinor))
        if self.plus_spinor[pivot].real < 0:
            self.plus_spinor *= -1
        self.minus_spinor = 1j*np.sign(self.slope)*self.sigma_t@self.plus_spinor
        if self.slope < 0:
            self.plus_indices = np.arange(self.pairs+1)
            self.minus_indices = np.arange(self.pairs)
            self.plus_pair_indices = np.arange(1, self.pairs+1)
            self.minus_pair_indices = np.arange(self.pairs)
            self.outgoing_zero_index, self.incoming_zero_index = 0, None
        else:
            self.plus_indices = np.arange(self.pairs)
            self.minus_indices = np.arange(self.pairs+1)
            self.plus_pair_indices = np.arange(self.pairs)
            self.minus_pair_indices = np.arange(1, self.pairs+1)
            self.outgoing_zero_index, self.incoming_zero_index = None, 0
        self.decay_rates = np.sqrt(self.masses**2-self.energy**2)
        self.Z = np.zeros((len(self.minus_indices), len(self.plus_indices)), complex)
        self.Z[self.minus_pair_indices, self.plus_pair_indices] = (
            self.energy-1j*self.decay_rates)/self.masses

    def _scalar_functions(self, points, derivative=False):
        points = np.asarray(points, float)
        if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
            raise ValueError('points must be a finite (npoints, 2) array')
        x = math.sqrt(self.frequency)*(points@self.tangent-self.center)
        values = _hermite_functions(x, self.pairs+int(derivative))
        if not derivative:
            return self.frequency**.25*values
        result = np.empty((len(points), self.pairs+1))
        for n in range(self.pairs+1):
            result[:, n] = -math.sqrt((n+1)/2)*values[:, n+1]
            if n:
                result[:, n] += math.sqrt(n/2)*values[:, n-1]
        return self.frequency**.75*result

    def basis_function(self, points):
        """Full-line normalized characteristic traces, evaluated at points."""
        values = self._scalar_functions(points)
        return (self.plus_spinor[:, None, None]*values[None, :, self.plus_indices],
                self.minus_spinor[:, None, None]*values[None, :, self.minus_indices])

    def transverse_derivative(self, points):
        """Analytic derivative with respect to the oriented tangent coordinate."""
        values = self._scalar_functions(points, derivative=True)
        return (self.plus_spinor[:, None, None]*values[None, :, self.plus_indices],
                self.minus_spinor[:, None, None]*values[None, :, self.minus_indices])

    def facet_gram(self, order=96):
        """Numerical Grams on [-R,R]; no normalization or tail certification."""
        if isinstance(order, bool) or int(order) != order or order < 2:
            raise ValueError('quadrature order must be an integer >=2')
        nodes, weights = np.polynomial.legendre.leggauss(int(order))
        points = self.R*self.normal+(self.R*nodes[:, None])*self.tangent
        return tuple(np.einsum('sni,n,snj->ij', f.conj(), self.R*weights, f)
                     for f in self.basis_function(points))

    def metadata(self):
        """JSON-safe analytic labels; Z itself is available as a complex array."""
        pairs = [f'pair{j}' for j in range(1, self.pairs+1)]
        return dict(face=self.face, R=self.R, eta=self.eta, energy=self.energy,
                    normal=self.normal.tolist(), tangent=self.tangent.tolist(),
                    angular_mass_slope=self.slope, transverse_center=self.center,
                    masses=self.masses.tolist(), decay_rates=self.decay_rates.tolist(),
                    outgoing_zero_index=self.outgoing_zero_index,
                    incoming_zero_index=self.incoming_zero_index,
                    plus_labels=(['zero']+pairs if self.slope < 0 else pairs),
                    minus_labels=(pairs if self.slope < 0 else ['zero']+pairs),
                    plus_hermite_indices=self.plus_indices.tolist(),
                    minus_hermite_indices=self.minus_indices.tolist(),
                    plus_pair_indices=self.plus_pair_indices.tolist(),
                    minus_pair_indices=self.minus_pair_indices.tolist(),
                    Z_real=self.Z.real.tolist(), Z_imag=self.Z.imag.tolist(),
                    norm_convention='L2(d tau) on the whole transverse real line',
                    frozen_lead_only=True,
                    finite_face_warning='Finite-face characteristic Grams differ; retain both transforms.')


def frozen_lead(face, R, eta, pairs=2, energy=0.):
    """Build one face's analytic channels; see FrozenLead for shapes/order."""
    return FrozenLead(face, R, eta, pairs, energy)


def self_check():
    """Small full-line quadrature and independent differential controls."""
    errors = dict(pauli_orientation=0., full_line_gram=0., normal_eigenvectors=0.,
                  actual_face_mass=0., transverse_finite_difference=0.,
                  decaying_lead_equation=0., zero_transverse_equation=0.)
    errors['pauli_orientation'] = float(np.max(abs(AU@AV-1j*SY)))
    nodes, weights = roots_hermite(32)
    facet_rows = []
    for face in NORMALS:
        for eta in (-.8, 0., .8):
            for energy in (0., .3):
                lead = frozen_lead(face, 2.5, eta, pairs=3, energy=energy)
                tau = lead.center+nodes/math.sqrt(lead.frequency)
                points = lead.R*lead.normal+tau[:, None]*lead.tangent
                full_weights = weights*np.exp(nodes**2)/math.sqrt(lead.frequency)
                plus, minus = lead.basis_function(points)
                for sign, f in ((1., plus), (-1., minus)):
                    gram = np.einsum('sni,n,snj->ij', f.conj(), full_weights, f)
                    errors['full_line_gram'] = max(errors['full_line_gram'],
                                                  float(np.max(abs(gram-np.eye(gram.shape[0])))))
                    errors['normal_eigenvectors'] = max(errors['normal_eigenvectors'],
                        float(np.max(abs(np.einsum('ab,bnj->anj', lead.sigma_n, f)-sign*f))))

                tau = lead.center+np.array([-1.2, -.45, .17, .7, 1.4])/math.sqrt(lead.frequency)
                points = lead.R*lead.normal+tau[:, None]*lead.tangent
                mass = lead.eta-2*points[:, 0]*points[:, 1]
                errors['actual_face_mass'] = max(errors['actual_face_mass'],
                    float(np.max(abs(mass-(lead.eta+lead.slope*tau)))))
                step = 2e-5/math.sqrt(lead.frequency)
                shifted = [lead.basis_function(points+j*step*lead.tangent) for j in (-2, -1, 1, 2)]
                basis = lead.basis_function(points)
                acted = []
                for part in (0, 1):
                    derivative = (shifted[0][part]-8*shifted[1][part]
                                  +8*shifted[2][part]-shifted[3][part])/(12*step)
                    acted.append(-1j*np.einsum('ab,bnj->anj', lead.sigma_t, derivative)
                                 +np.einsum('ab,bnj,n->anj', SY, basis[part], mass))
                for j, (ip, im) in enumerate(zip(lead.plus_pair_indices, lead.minus_pair_indices)):
                    error = max(np.max(abs(acted[0][:, :, ip]-lead.masses[j]*basis[1][:, :, im])),
                                np.max(abs(acted[1][:, :, im]-lead.masses[j]*basis[0][:, :, ip])))
                    errors['transverse_finite_difference'] = max(errors['transverse_finite_difference'], float(error))
                    z = lead.Z[im, ip]
                    field = basis[0][:, :, ip]+z*basis[1][:, :, im]
                    hfield = (1j*lead.decay_rates[j]*(lead.sigma_n@field)
                              +acted[0][:, :, ip]+z*acted[1][:, :, im])
                    errors['decaying_lead_equation'] = max(errors['decaying_lead_equation'],
                                                          float(np.max(abs(hfield-energy*field))))
                part = 0 if lead.outgoing_zero_index is not None else 1
                errors['zero_transverse_equation'] = max(errors['zero_transverse_equation'],
                                                         float(np.max(abs(acted[part][:, :, 0]))))
                if energy == 0.:
                    gp, gm = lead.facet_gram()
                    # These are measured finite-face omitted norms, not estimates
                    # on omitted channels or on the bilinear exterior solution.
                    facet_rows.append(dict(face=face, eta=eta,
                        plus_missing_norms=(1-np.diag(gp).real).tolist(),
                        minus_missing_norms=(1-np.diag(gm).real).tolist()))
    for key, value in errors.items():
        tolerance = 2e-9 if 'equation' in key or 'finite_difference' in key else 2e-13
        if value > tolerance:
            raise AssertionError(f'{key}: {value} exceeds {tolerance}')
    return dict(passed=True, cases=24, maximum_errors=errors,
                finite_face_tail_estimates=facet_rows,
                scope='Exact frozen full-line lead basis; no bilinear-tail or angular-completeness certificate.')


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2, allow_nan=False))
