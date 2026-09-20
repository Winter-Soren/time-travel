"""Pointwise bridge data combining sewn traces with independent zero-trace modes.

The global columns are exactly the 138 linear trace lifts. Compact columns
are the 137 regular monopole/sine family and vanish at both cuts. All columns
share split angular and radial quadrature before their differential actions
are squared. The two polar panels are needed by the compact family even
though the exterior trace vanishes there. No exterior assembly, projected
first-order square, source integral or regulator limit is implemented here.
"""
from functools import lru_cache
import json
import math
import operator

import numpy as np

from balanced_dihole import symbolic_fields
from dihole_azimuthal_bridge import AzimuthalBridgeMatrices
from dihole_bridge_scattering import collar_mask
from dihole_euclidean_form import mask_derivative
from dihole_sewn_form import RHO_BOUNDS, SX4, SY4, SZ4, matrix_action, trace_values
from dihole_stretched_collar import exact_target_parameters, field_values, gauss, radius_profile


def _positive_integer(value, name):
    try:
        if isinstance(value, bool):
            raise TypeError
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(name+' must be a positive integer.') from exc
    if result < 1:
        raise ValueError(name+' must be a positive integer.')
    return result


def _arguments(angular_modes, radial_modes, degrees):
    angular_modes = _positive_integer(angular_modes, 'angular_modes')
    radial_modes = _positive_integer(radial_modes, 'radial_modes')
    if len(degrees) != 2:
        raise ValueError('degrees must contain the rho and z polynomial counts.')
    degrees = tuple(_positive_integer(d, 'degrees entry') for d in degrees)
    return angular_modes, radial_modes, degrees


@lru_cache(maxsize=1)
def _geometry():
    return symbolic_fields()[0], exact_target_parameters()


def column_metadata(mu, angular_modes=2, radial_modes=2, degrees=(2, 3)):
    """Explicit column labels, suitable for nested subspace selectors.

    Global order: component, rho degree, z degree. Compact order: R/L,
    upper/lower, absolute angular l, radial sine number (starting at one).
    """
    angular_modes, radial_modes, degrees = _arguments(angular_modes, radial_modes, degrees)
    # Reuse the same half-odd-integer validation and label convention as the
    # actual profile helper; no physical angular degeneracy is inserted.
    from complete_angular_modes import _sector
    a = _sector(mu)
    first = max(1, abs(a))
    lower = list(range(first, first+angular_modes))
    upper = ([0]+lower) if a == 0 else lower.copy()
    columns = []
    for component in range(4):
        for rho_degree in range(degrees[0]):
            for z_degree in range(degrees[1]):
                columns.append(dict(index=len(columns), family='global', component=component,
                                    rho_degree=rho_degree, z_degree=z_degree))
    global_count = len(columns)
    for chirality in ('R', 'L'):
        for spin, labels in (('upper', upper), ('lower', lower)):
            for ell in labels:
                for radial_mode in range(1, radial_modes+1):
                    columns.append(dict(index=len(columns), family='compact', chirality=chirality,
                                        spin=spin, l=ell, radial_mode=radial_mode))
    cmax = math.sqrt(1-(RHO_BOUNDS[0]/.05)**2)
    return dict(mu=a+.5, a=a, angular_modes=angular_modes, radial_modes=radial_modes,
                degrees=list(degrees), p=len(upper), m=len(lower), plus_labels=upper, minus_labels=lower,
                global_columns=global_count, compact_columns=len(columns)-global_count,
                total_columns=len(columns), columns=columns,
                cosine_panels=[[-1., -cmax], [-cmax, cmax], [cmax, 1.]],
                radial_panels=[[-.1, -.05], [-.05, .05], [.05, .1]],
                trace_phase='right global trace multiplied by exp(-i*phase); compact traces zero')


def _two_weyls(data, action=False):
    """137 ordering: R columns first, L columns second; T_L=-T_R."""
    zero = np.zeros_like(data)
    left = -data if action else data
    return np.concatenate((np.concatenate((data, zero), axis=2),
                           np.concatenate((zero, left), axis=2)), axis=1)


def bridge_data(mu, phase, angular_modes=2, radial_modes=2, order=48, degrees=(2, 3)):
    """Return the nine 137 assemble keys, with global columns before compact.

    value/radial/angular have shape (9*order**2,4,total_columns).
    Scalar fields and weights have shape (9*order**2,). Weights integrate
    ds*dcos(theta); no extra sin(theta), radius factor or 2*pi is needed.
    The caller must append zero compact columns to its exterior data.
    """
    angular_modes, radial_modes, degrees = _arguments(angular_modes, radial_modes, degrees)
    order = _positive_integer(order, 'order')
    if order <= radial_modes:
        raise ValueError('order must exceed radial_modes.')
    phase = float(phase)
    if not math.isfinite(phase):
        raise ValueError('phase must be finite.')
    function, parameters = _geometry()
    bridge = AzimuthalBridgeMatrices(function, parameters, angular_modes, order, mu=mu)
    nodes, weights = gauss(order)
    cmax = math.sqrt(1-(RHO_BOUNDS[0]/.05)**2)
    intervals = ((-1., -cmax), (-cmax, cmax), (cmax, 1.))
    cosine = np.concatenate([(hi+lo)/2+(hi-lo)*nodes/2 for lo, hi in intervals])
    angular_weights = np.concatenate([(hi-lo)*weights/2 for lo, hi in intervals])
    theta = np.arccos(cosine)
    sine = np.sin(theta)
    left, left_theta = trace_values(theta, 'left', degrees)
    right, right_theta = [np.exp(-1j*phase)*v for v in trace_values(theta, 'right', degrees)]
    global_s = (right-left)/.2

    upper, lower = bridge.profiles(theta)
    p, m, n = bridge.p, bridge.m, bridge.n
    profiles = np.zeros((len(theta), 2, n), complex)
    profiles[:, 0, :p], profiles[:, 1, p:] = upper, 1j*lower
    reference = np.zeros_like(profiles)
    lam = np.sqrt(bridge.minus_labels*(bridge.minus_labels+1))
    reference[:, 0, p:] = upper[:, int(bridge.a == 0):]*lam
    reference[:, 1, int(bridge.a == 0):p] = 1j*lower*lam
    sy_profiles = np.stack((-1j*profiles[:, 1], 1j*profiles[:, 0]), axis=1)
    modes = np.arange(1, radial_modes+1)
    rows = {key: [] for key in ('value', 'radial', 'angular', 'P', 'W', 'f', 'fp', 'weight')}

    for lo, hi in ((-.1, -.05), (-.05, .05), (.05, .1)):
        for s, weight in zip((hi+lo)/2+(hi-lo)*nodes/2, (hi-lo)*weights/2):
            t = (s+.1)/.2
            global_value = (1-t)*left+t*right
            global_theta = (1-t)*left_theta+t*right_theta
            radius = float(radius_profile(min(.05, max(0., abs(s)-.05)))[0])
            fields = field_values(function, parameters, radius, cosine)
            lapse, meridian, potential = [np.asarray(fields[key], float) for key in
                                         ('lapse', 'meridian_factor', 'unit_A_phi')]
            angular_coefficient = lapse*meridian*(bridge.mu-potential)/sine
            global_angular = (-1j*matrix_action(SX4, global_theta
                              +cosine[:, None, None]*global_value/(2*sine[:, None, None]))
                              +angular_coefficient[:, None, None]*matrix_action(SY4, global_value))/radius
            global_radial = -1j*matrix_action(SZ4, global_s)

            delta = (lapse*meridian*(bridge.mu-potential)
                     -(bridge.a+.5*np.cos(theta)))/sine
            compact_angular = (reference+delta[:, None, None]*sy_profiles)/radius
            radial = np.sqrt(10.)*np.sin(modes*np.pi*(s+.1)/.2)
            derivative = np.sqrt(10.)*modes*np.pi/.2*np.cos(modes*np.pi*(s+.1)/.2)
            compact_value = (profiles[..., None]*radial).reshape(len(theta), 2, -1)
            compact_radial = (-1j*np.array([1., -1.])[None, :, None, None]
                              *profiles[..., None]*derivative).reshape(len(theta), 2, -1)
            compact_angular = (compact_angular[..., None]*radial).reshape(len(theta), 2, -1)
            for key, global_data, compact in (
                    ('value', global_value, _two_weyls(compact_value)),
                    ('radial', global_radial, _two_weyls(compact_radial, action=True)),
                    ('angular', global_angular, _two_weyls(compact_angular, action=True))):
                rows[key].append(np.concatenate((global_data, compact), axis=2))
            for key, data in (('P', meridian), ('W', meridian/lapse),
                              ('f', np.full(len(theta), collar_mask(s))),
                              ('fp', np.full(len(theta), mask_derivative(s))),
                              ('weight', weight*angular_weights)):
                rows[key].append(data)
    return {key: np.concatenate(value) for key, value in rows.items()}


def self_check():
    """Small separate-block and nested-column controls; no exterior/source solve."""
    from dihole_euclidean_form import assemble, relative_error, trial_data
    from dihole_sewn_form import bridge_data as old_global_data

    function, parameters = _geometry()
    records = []
    for mu, phase in ((-.5, 0.), (.5, .7), (1.5, np.pi)):
        data = bridge_data(mu, phase, angular_modes=1, radial_modes=1, order=24, degrees=(1, 2))
        metadata = column_metadata(mu, 1, 1, (1, 2))
        cut = metadata['global_columns']
        full = assemble(data, .7, mass=1.3)
        global_form = assemble(old_global_data(mu, phase, order=24, degrees=(1, 2)), .7, mass=1.3)
        old_bridge = AzimuthalBridgeMatrices(function, parameters, 1, 48, mu=mu)
        compact_form = assemble(trial_data(old_bridge, radial_modes=1, radial_order=24), .7, mass=1.3)
        keys = ('q', 'norm', 'derivative', 'weak_dirac', 'direct_massive')
        records.append(dict(mu=mu, phase=phase, dimension=metadata['total_columns'],
            global_block_errors={key: relative_error(full[key][:cut, :cut], global_form[key]) for key in keys},
            compact_block_quadrature_changes={key: relative_error(full[key][cut:, cut:], compact_form[key]) for key in keys},
            mixed_norm=float(np.linalg.norm(full['norm'][:cut, cut:])),
            scalar_weight_sum=float(np.sum(data['weight']))))

    small = bridge_data(.5, .7, angular_modes=1, radial_modes=1, order=12, degrees=(1, 2))
    large = bridge_data(.5, .7, angular_modes=2, radial_modes=2, order=12, degrees=(1, 2))
    low_meta, high_meta = column_metadata(.5, 1, 1, (1, 2)), column_metadata(.5, 2, 2, (1, 2))
    def key(record):
        return tuple((name, value) for name, value in record.items() if name != 'index')
    lookup = {key(record): record['index'] for record in high_meta['columns']}
    indices = [lookup[key(record)] for record in low_meta['columns']]
    nested = {name: float(np.max(abs(small[name]-large[name][:, :, indices])))
              for name in ('value', 'radial', 'angular')}
    tests = dict(
        unchanged_global_block=max(value for row in records for value in row['global_block_errors'].values()) < 2e-10,
        compact_recovers_137=max(value for row in records for value in row['compact_block_quadrature_changes'].values()) < 2e-6,
        mixed_terms_present=min(row['mixed_norm'] for row in records) > 1e-8,
        full_sphere_and_radial_weight=max(abs(row['scalar_weight_sum']-.4) for row in records) < 2e-14,
        nested_compact_labels=max(nested.values()) < 2e-11,
        column_dimensions=all(row['dimension'] == 8+2*(2+int(row['mu'] == .5)) for row in records))
    return dict(passed=bool(all(tests.values())), tests={name: bool(value) for name, value in tests.items()},
                controls=records, nested_pointwise_errors=nested,
                scope='Pointwise shared-grid and finite form controls only; no sewn exterior integral, spectral convergence, or physical source.')


if __name__ == '__main__':
    report = self_check()
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report['passed']:
        raise SystemExit('Enriched bridge controls failed.')
