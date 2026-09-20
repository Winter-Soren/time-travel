"""Memory-bounded radial enrichment of the frozen 139 pointwise bridge data.

Only the Nr=1 fields are materialized by the existing angular/geometry helper.
Compact values and angular actions factor as profile times the first positive
interior sine. Recover that profile per small group of radial slices, then
expand the requested sines. The full-Dirac radial action uses
Sigma=diag(1,-1,-1,1). Global trace-lift columns and all scalar data are unchanged.
No angular operator, geometry, exterior, or projected square is reimplemented.
"""
import json
import operator

import numpy as np

from dihole_enriched_bridge import bridge_data, column_metadata
from dihole_stretched_collar import gauss


SIGMA = np.array([1., -1., -1., 1.])
RADIAL_PANELS = ((-.1, -.05), (-.05, .05), (.05, .1))
FIELD_KEYS = ('value', 'radial', 'angular')


def _integer(value, name):
    try:
        if isinstance(value, bool):
            raise TypeError
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError(name+' must be a positive integer.') from exc
    if value < 1:
        raise ValueError(name+' must be a positive integer.')
    return value


def radial_recovery_bounds(order=64, radial_modes=12):
    """Numerical scale diagnostics and analytic sine inequalities at these nodes.

    sin(x)>=2*min(x,pi-x)/pi and |sin(n*x)/sin(x)|<=n hold for 0<x<pi.
    These formula checks are not interval enclosures of the Gauss nodes.
    """
    order, radial_modes = _integer(order, 'order'), _integer(radial_modes, 'radial_modes')
    if order <= radial_modes:
        raise ValueError('order must exceed radial_modes, as in frozen bridge_data.')
    nodes = gauss(order)[0]
    s = np.concatenate([(hi+lo)/2+(hi-lo)*nodes/2 for lo, hi in RADIAL_PANELS])
    first = np.sqrt(10.)*np.sin(np.pi*(s+.1)/.2)
    if np.any(~np.isfinite(first)) or np.any(first <= np.finfo(float).tiny):
        raise ArithmeticError('The first sine is not positive and normal at every interior radial node.')
    distances = np.minimum(s+.1, .1-s)
    labels = np.arange(1, radial_modes+1)
    sines = np.sqrt(10.)*np.sin(labels[None, :]*np.pi*(s[:, None]+.1)/.2)
    relative = abs(sines/first[:, None])/labels[None, :]
    return dict(order=order, radial_modes=radial_modes, radial_slices=len(s),
                minimum_endpoint_distance=float(np.min(distances)),
                minimum_sine1=float(np.min(first)),
                maximum_profile_recovery_amplification=float(1/np.min(first)),
                analytic_sine1_lower_bound_at_float_nodes=float(np.sqrt(10.)*10*np.min(distances)),
                maximum_sine_ratio_divided_by_mode=float(np.max(relative)),
                maximum_exact_sine_derivative=np.sqrt(10.)*radial_modes*np.pi/.2,
                scope='Numerical node diagnostics; sine inequalities are analytic, but endpoints/roundoff are not enclosed.')


def bridge_chunks(mu, phase, angular_modes=3, radial_modes=12, order=64,
                  degrees=(3, 4), slices_per_chunk=4):
    """Yield the nine assemble-compatible data arrays, global columns first.

    Column order is exactly column_metadata(mu,angular_modes,radial_modes,degrees).
    A chunk has at most slices_per_chunk*3*order angular/radial points.
    The baseline Nr=1 arrays remain resident; high-Nr fields do not. Scalar
    arrays are read-only views by convention: callers must not mutate them.
    """
    order, radial_modes = _integer(order, 'order'), _integer(radial_modes, 'radial_modes')
    slices_per_chunk = _integer(slices_per_chunk, 'slices_per_chunk')
    if order <= radial_modes:
        raise ValueError('order must exceed radial_modes, as in frozen bridge_data.')
    metadata = column_metadata(mu, angular_modes, radial_modes, degrees)
    global_count = metadata['global_columns']
    baseline_count = metadata['compact_columns']//radial_modes
    nodes = gauss(order)[0]
    s = np.concatenate([(hi+lo)/2+(hi-lo)*nodes/2 for lo, hi in RADIAL_PANELS])
    first = np.sqrt(10.)*np.sin(np.pi*(s+.1)/.2)
    if np.any(~np.isfinite(first)) or np.any(first <= np.finfo(float).tiny):
        raise ArithmeticError('Cannot safely recover profiles from a nonpositive/subnormal first sine.')
    baseline = bridge_data(mu, phase, angular_modes=angular_modes, radial_modes=1,
                           order=order, degrees=degrees)
    angular_points = 3*order
    expected_shape = (len(s)*angular_points, 4, global_count+baseline_count)
    if any(baseline[key].shape != expected_shape for key in FIELD_KEYS):
        raise ValueError('Frozen bridge_data shape/order no longer matches this streaming adapter.')
    modes = np.arange(1, radial_modes+1)
    for lo in range(0, len(s), slices_per_chunk):
        hi = min(len(s), lo+slices_per_chunk)
        rows = slice(lo*angular_points, hi*angular_points)
        ns = hi-lo
        angle = modes[None, :]*np.pi*(s[lo:hi, None]+.1)/.2
        sines = np.sqrt(10.)*np.sin(angle)
        derivatives = np.sqrt(10.)*modes[None, :]*np.pi/.2*np.cos(angle)
        denominator = first[lo:hi, None, None, None]
        profiles = (baseline['value'][rows, :, global_count:]
                    .reshape(ns, angular_points, 4, baseline_count)/denominator)
        angular_profiles = (baseline['angular'][rows, :, global_count:]
                            .reshape(ns, angular_points, 4, baseline_count)/denominator)
        compact = dict(
            value=profiles[..., None]*sines[:, None, None, None, :],
            radial=-1j*SIGMA[None, None, :, None, None]*profiles[..., None]
                   *derivatives[:, None, None, None, :],
            angular=angular_profiles[..., None]*sines[:, None, None, None, :])
        result = {key: value[rows] for key, value in baseline.items() if key not in FIELD_KEYS}
        for key in FIELD_KEYS:
            expanded = compact[key].reshape(ns*angular_points, 4, baseline_count*radial_modes)
            result[key] = np.concatenate((baseline[key][rows, :, :global_count], expanded), axis=2)
        yield result


def self_check():
    """Small direct-139 comparisons and independent full-form sums; no exterior."""
    from dihole_euclidean_form import assemble, relative_error

    records = []
    for mu, nr, order in ((.5, 3, 16), (-.5, 8, 32)):
        kwargs = dict(mu=mu, phase=.7, angular_modes=1, radial_modes=nr,
                      order=order, degrees=(1, 2))
        direct = bridge_data(**kwargs)
        expected = assemble(direct, .7, epsilon=.013, mass=1.3)
        assembled, offset, largest_chunk = None, 0, 0
        field_errors = {key: 0. for key in FIELD_KEYS}
        relative_fields = {key: 0. for key in FIELD_KEYS}
        scalar_equal, global_equal = True, True
        global_count = column_metadata(mu, 1, nr, (1, 2))['global_columns']
        for chunk in bridge_chunks(**kwargs, slices_per_chunk=4):
            count = len(chunk['P'])
            at = slice(offset, offset+count)
            for key in FIELD_KEYS:
                difference = float(np.max(abs(chunk[key]-direct[key][at])))
                field_errors[key] = max(field_errors[key], difference)
                relative_fields[key] = max(relative_fields[key], difference/max(float(np.max(abs(direct[key][at]))), 1e-300))
                global_equal &= np.array_equal(chunk[key][:, :, :global_count], direct[key][at, :, :global_count])
            scalar_equal &= all(np.array_equal(chunk[key], direct[key][at])
                                for key in direct if key not in FIELD_KEYS)
            form = assemble(chunk, .7, epsilon=.013, mass=1.3)
            if assembled is None:
                assembled = form
            else:
                for key in assembled:
                    assembled[key] += form[key]
            largest_chunk = max(largest_chunk, sum(value.nbytes for value in chunk.values()))
            offset += count
        records.append(dict(mu=mu, radial_modes=nr, order=order, points=offset,
            field_absolute_errors=field_errors, field_relative_errors=relative_fields,
            form_errors={key: relative_error(assembled[key], expected[key]) for key in expected},
            global_columns_bitwise_equal=bool(global_equal), scalar_data_bitwise_equal=bool(scalar_equal),
            all_points_consumed=offset == len(direct['P']), largest_chunk_array_bytes=largest_chunk))
    recovery = radial_recovery_bounds(64, 12)
    modes = np.arange(1, 13)
    numerical_right_trace = float(np.max(abs(np.sqrt(10.)*np.sin(modes*np.pi))))
    tests = dict(
        unchanged_global_columns=all(row['global_columns_bitwise_equal'] for row in records),
        unchanged_scalar_data=all(row['scalar_data_bitwise_equal'] for row in records),
        all_points=all(row['all_points_consumed'] for row in records),
        direct_pointwise_agreement=max(v for row in records for v in row['field_relative_errors'].values()) < 3e-13,
        direct_form_agreement=max(v for row in records for v in row['form_errors'].values()) < 3e-12,
        positive_normal_recovery=recovery['minimum_sine1'] >= recovery['analytic_sine1_lower_bound_at_float_nodes'] > 0,
        sine_ratio_bound=recovery['maximum_sine_ratio_divided_by_mode'] < 1+2e-12,
        analytic_zero_trace_roundoff=numerical_right_trace < 2e-14)
    return dict(passed=bool(all(tests.values())), tests={key: bool(value) for key, value in tests.items()},
                controls=records, production_node_recovery=recovery,
                analytic_zero_trace_numerical_residue=numerical_right_trace,
                ordering=column_metadata(.5, 1, 3, (1, 2))['columns'],
                scope='Streamed pointwise assembly reproduces frozen139 at fixed quadrature; no radial, exterior, spectral, or source convergence is inferred.')


if __name__ == '__main__':
    report = self_check()
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report['passed']:
        raise SystemExit('Radial streaming controls failed.')
