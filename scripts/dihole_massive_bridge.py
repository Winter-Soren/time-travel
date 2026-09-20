"""Physical Dirac mass on the existing prescribed balanced-dihole bridge.

In the original R/L ordering H_m=[[H_0,m*N],[m*N,-H_0]].  The canonical
current variable therefore has a transverse coupling m*P, not m*N or m*W.
Reorder the characteristic components as

  plus  = [R_up(n+1), L_down(n)]
  minus = [R_down(n), L_up(n+1)].

The two characteristic spaces both have dimension 2*n+1.  This helper
includes both physical Weyls explicitly; it supplies no extra determinant
multiplicity.  It reuses the existing finite angular basis, actual metric
and stable embedding.  Neither an exterior state, regulator subtraction,
renormalized stress nor a self-consistent throat is computed here.
"""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import block_diag

from balanced_dihole import symbolic_fields
from dihole_bridge_scattering import BridgeMatrices, collar_mask, complex_record, embed
from dihole_direct_source import direct_response
from dihole_stretched_collar import exact_target_parameters, field_values, gauss, radius_profile
from project_paths import project_path
from spherical_weyl_leads import _mouth_arguments, mouth_characteristics, zero_mode


class MassiveBridgeMatrices:
    """Fixed physical mass; coefficients(s,epsilon=0) is embed-compatible.

    The same compact paired variation as BridgeMatrices is used.  Radius,
    angular metric and magnetic potential stay fixed: all clock blocks gain
    exp(-2*epsilon*f), and the entire transverse block, including m*P,
    gains exp(-epsilon*f).  Cached unvaried blocks do not depend on epsilon.
    """

    def __init__(self, function, parameters, size, order, mass=0.):
        if isinstance(size, bool) or int(size) != size or size < 1:
            raise ValueError('size must be a positive integer.')
        if isinstance(order, bool) or int(order) != order or order <= size+1:
            raise ValueError('order must be an integer greater than size+1.')
        if isinstance(mass, bool) or not math.isfinite(float(mass)) or mass < 0:
            raise ValueError('Physical mass must be finite and nonnegative.')
        if parameters.get('q') != 1.:
            raise ValueError('This angular basis is the q=1, mu=1/2 sector.')
        self.function, self.parameters = function, dict(parameters)
        self.size, self.order, self.mass = int(size), int(order), float(mass)
        self.weyl = BridgeMatrices(function, self.parameters, self.size, self.order)
        self.plus = self.minus = 2*self.size+1
        self.cosine, self.weights = gauss(self.order)
        theta = np.arccos(self.cosine)
        _, _, er_plus, er_minus = _mouth_arguments(theta, 1., 'upper')
        plus, minus = mouth_characteristics(theta, 1., 'upper', self.size)
        upper = np.sum(er_plus[:, :, None]*minus, axis=0).real
        self.lower = np.sum(er_minus[:, :, None]*plus, axis=0).imag
        zeroth = np.sum(er_plus*zero_mode(theta, 1., 'upper'), axis=0)
        self.upper = np.column_stack((zeroth, upper))
        # Scalar profiles omit sqrt(sin theta); weights integrate dcos(theta).
        self.basis_gram_error = float(max(
            np.max(abs(self.upper.conj().T@(self.weights[:, None]*self.upper)-np.eye(self.size+1))),
            np.max(abs(self.lower.conj().T@(self.weights[:, None]*self.lower)-np.eye(self.size)))))

    @property
    def gram_error(self):
        return max(self.basis_gram_error, self.weyl.gram_error)

    @lru_cache(maxsize=4096)
    def mass_blocks_at_radius(self, radius):
        """Actual positive P compressions in the unchanged up/down basis."""
        data = field_values(self.function, self.parameters, float(radius), self.cosine)
        meridian = np.asarray(data['meridian_factor'], float)
        if not np.all(np.isfinite(meridian)) or np.any(meridian <= 0):
            raise ValueError('The bridge meridian factor must be finite and positive.')
        weighted = (self.weights*meridian)[:, None]
        return (self.upper.conj().T@(weighted*self.upper),
                self.lower.conj().T@(weighted*self.lower))

    @lru_cache(maxsize=4096)
    def at_radius(self, radius):
        wp, wm, cross = self.weyl.at_radius(float(radius))
        pp, pm = self.mass_blocks_at_radius(float(radius))
        full = np.block([[cross, self.mass*pp],
                         [self.mass*pm, -cross.conj().T]])
        return block_diag(wp, wm), block_diag(wm, wp), full

    def coefficients(self, s, epsilon=0.):
        s, epsilon = float(s), float(epsilon)
        if not math.isfinite(s) or not math.isfinite(epsilon):
            raise ValueError('Coordinate and variation must be finite.')
        u = min(.05, max(0., abs(s)-.05))
        radius = float(radius_profile(u)[0])
        wp, wm, cross = self.at_radius(radius)
        factor = math.exp(-epsilon*collar_mask(s))
        return factor*factor*wp, factor*factor*wm, factor*cross


def massless_permutations(size):
    """Index a blockdiag(S_R,S_L) into the massive helper's port ordering."""
    p, m = size+1, size
    total = p+m
    incoming = np.r_[np.arange(p), np.arange(total, total+m),
                     np.arange(p, total), np.arange(total+m, 2*total)]
    outgoing = np.r_[np.arange(m), np.arange(total, total+p),
                     np.arange(m, total), np.arange(total+p, 2*total)]
    return outgoing, incoming


def _coefficient_controls(model):
    """Reorder the physical R/L Hamiltonian independently of the block recipe."""
    p, m = model.size+1, model.size
    n = p+m
    original_j = np.diag(np.r_[np.ones(p), -np.ones(m)])
    order = np.r_[np.arange(p), np.arange(n+p, 2*n),
                  np.arange(p, n), np.arange(n, n+p)]
    rows = []
    for radius in (.025, .037, .05):
        wp, wm, x = model.weyl.at_radius(radius)
        pp, pm = model.mass_blocks_at_radius(radius)
        transverse = np.block([[np.zeros((p, p)), x],
                               [x.conj().T, np.zeros((m, m))]])
        metric = block_diag(pp, pm)
        original = np.block([[transverse, model.mass*metric],
                             [model.mass*metric, -transverse]])
        clockwise = block_diag(block_diag(wp, wm), block_diag(wp, wm))
        wplus, wminus, full = model.at_radius(radius)
        reordered = np.block([[np.zeros((n, n)), full],
                              [full.conj().T, np.zeros((n, n))]])
        rows.append(dict(radius=radius,
            transverse_permutation_error=float(np.max(abs(original[np.ix_(order, order)]-reordered))),
            clock_permutation_error=float(np.max(abs(clockwise[np.ix_(order, order)]-block_diag(wplus, wminus)))),
            current_permutation_error=float(np.max(abs(block_diag(original_j, -original_j)[np.ix_(order, order)]
                                                       -np.diag(np.r_[np.ones(n), -np.ones(n)])))),
            minimum_mass_compression_eigenvalue=float(min(np.linalg.eigvalsh(pp)[0], np.linalg.eigvalsh(pm)[0])),
            minimum_clock_eigenvalue=float(min(np.linalg.eigvalsh(wplus)[0], np.linalg.eigvalsh(wminus)[0]))))
    base = model.coefficients(-.08)
    scaling = []
    for epsilon in (-.025, .025):
        varied = model.coefficients(-.08, epsilon)
        f = math.exp(-epsilon*collar_mask(-.08))
        scaling.append(float(max(np.max(abs(a-scale*b)) for a, b, scale in
                                 zip(varied, base, (f*f, f*f, f)))))
    return dict(radial_rows=rows, variation_scaling_errors=scaling,
                angular_gram_error=model.gram_error)


def _constant_mass_control():
    """Independent flat constant-lapse dispersion: E^2=E0^2+(mass*N)^2."""
    p, m = 3, 2
    lapse, meridian, mass, momentum = .7, 2.3, 1.3, .4
    weight = meridian/lapse
    x = np.vstack((np.zeros((1, m)), np.diag([1.7, 2.9])))
    full = np.block([[x, mass*meridian*np.eye(p)],
                     [mass*meridian*np.eye(m), -x.T]])
    count = p+m
    h = (np.diag(np.r_[np.full(count, momentum), np.full(count, -momentum)])
         +np.block([[np.zeros((count, count)), full], [full.T, np.zeros((count, count))]]))/weight
    h0 = (np.diag(np.r_[np.full(p, momentum), np.full(m, -momentum)])
          +np.block([[np.zeros((p, p)), x], [x.T, np.zeros((m, m))]]))/weight
    positive = np.sqrt(np.linalg.eigvalsh(h0)**2+(mass*lapse)**2)
    expected = np.sort(np.r_[-positive, positive])
    actual = np.linalg.eigvalsh(h)
    return dict(lapse=lapse, meridian=meridian, mass=mass, momentum=momentum,
                spectrum=actual.tolist(), expected=expected.tolist(),
                maximum_spectrum_error=float(np.max(abs(actual-expected))))


def checks():
    started = perf_counter()
    function, _ = symbolic_fields()
    parameters = exact_target_parameters()
    size, order, physical_mass = 2, 32, 1.3
    zero = MassiveBridgeMatrices(function, parameters, size, order, mass=0.)
    massive = MassiveBridgeMatrices(function, parameters, size, order, mass=physical_mass)
    coefficient_controls = _coefficient_controls(massive)
    permutation = massless_permutations(size)
    p, m, count = size+1, size, 2*size+1

    def left_coefficients(s):
        wp, wm, cross = zero.weyl.coefficients(s)
        return wm, wp, -cross.conj().T

    rows = []
    for energy in (.2, .1j):
        right, _ = embed(zero.weyl.coefficients, p, m, energy)
        left, _ = embed(left_coefficients, m, p, energy)
        old = block_diag(right, left)[np.ix_(*permutation)]
        joined, zero_info = embed(zero.coefficients, count, count, energy)
        result, info = embed(massive.coefficients, count, count, energy)
        rows.append(dict(energy=complex_record(energy),
                         massless_permutation_error=float(np.max(abs(joined-old))),
                         massless_diagnostics=zero_info, massive_diagnostics=info,
                         massive_scattering=complex_record(result)))

    kappa = .1
    direct_zero, _ = direct_response(zero.coefficients, count, count, kappa)
    direct_right, _ = direct_response(zero.weyl.coefficients, p, m, kappa)
    direct_left, _ = direct_response(left_coefficients, m, p, kappa)
    sum_keys = ('direct_logdet', 'direct_dlog', 'power_subtracted_dlog', 'volume_factor')
    direct_additivity = {key: abs(direct_zero[key]-direct_right[key]-direct_left[key]) for key in sum_keys}
    direct, reflection = direct_response(massive.coefficients, count, count, kappa)
    imaginary = rows[1]['massive_scattering']
    imaginary = np.array(imaginary['real'])+1j*np.array(imaginary['imag'])
    reflection_error = float(np.max(abs(reflection-imaginary[count:, count:])))
    eta = 1e-4
    varied = [direct_response(lambda s, e=e: massive.coefficients(s, e), count, count, kappa)[0]['direct_logdet']
              for e in (-eta, eta)]
    derivative_error = abs(direct['direct_dlog']-(varied[1]-varied[0])/(2*eta))
    constant = _constant_mass_control()
    tests = dict(
        physical_block_permutation=max(r['transverse_permutation_error'] for r in coefficient_controls['radial_rows']) < 2e-12,
        clock_and_current_permutation=max(max(r['clock_permutation_error'], r['current_permutation_error'])
                                         for r in coefficient_controls['radial_rows']) < 2e-12,
        positive_mass_and_clock=min(min(r['minimum_mass_compression_eigenvalue'], r['minimum_clock_eigenvalue'])
                                    for r in coefficient_controls['radial_rows']) > 0,
        angular_normalization=coefficient_controls['angular_gram_error'] < 2e-12,
        paired_variation_scaling=max(coefficient_controls['variation_scaling_errors']) < 2e-12,
        massless_two_weyl_scattering=max(r['massless_permutation_error'] for r in rows) < 2e-8,
        real_energy_unitarity=rows[0]['massive_diagnostics']['unitarity_error'] < 2e-8,
        imaginary_energy_passivity=rows[1]['massive_diagnostics']['minimum_passivity_defect'] > -2e-9,
        constant_mass_dispersion=constant['maximum_spectrum_error'] < 2e-12,
        direct_massless_additivity=max(direct_additivity.values()) < 2e-8,
        direct_reflection_agreement=reflection_error < 2e-8,
        direct_metric_derivative=derivative_error < 3e-6)
    inputs = ('balanced_dihole.py', 'dihole_stretched_collar.py', 'dihole_cylinder_transfer.py',
              'dihole_bridge_scattering.py', 'dihole_direct_source.py', 'spherical_weyl_leads.py', 'project_paths.py')
    return dict(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                parameters=dict(q=1, mu=.5, size=size, angular_order=order, physical_mass=physical_mass),
                ordering=dict(plus='R_up(n+1), L_down(n)', minus='R_down(n), L_up(n+1)',
                              scattering_input='plus_left, minus_right', scattering_output='minus_left, plus_right'),
                coefficient_controls=coefficient_controls, scattering_controls=rows,
                constant_mass_control=constant,
                direct_controls=dict(kappa=kappa, result=direct, massless_additivity_errors=direct_additivity,
                                     reflection_error=reflection_error, metric_derivative_fd_error=derivative_error),
                scope=['Actual prescribed q1, mu1/2 bridge coefficients, finite angular compression and numerical ODE checks.',
                       'Both physical Weyls are present explicitly; do not multiply determinant responses by a second Weyl factor.',
                       'Imaginary-frequency scattering is resolvent data, not a propagation probability.',
                       'No massive exterior, Pauli-Villars subtraction, covariant normalization, complete stress or Einstein solution is computed.'],
                input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in inputs},
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                runtime_seconds=perf_counter()-started)


if __name__ == '__main__':
    report = checks()
    project_path('dihole_massive_bridge_results.json').write_text(
        json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=report['passed'], tests=report['tests'],
                         direct_controls=report['direct_controls'], runtime_seconds=report['runtime_seconds']), indent=2))
    if not report['passed']:
        raise SystemExit('Massive bridge checks failed; saved diagnostic evidence retained.')
