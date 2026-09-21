"""Stream zero-trace residual candidates on the unchanged checkpoint157 bridge.

Local columns0:324 reproduce mouth_fields exactly. The next176 columns extend
its22 compact angular profiles to radial sines13:20; the final96 carry angular
l=6,7 in both spins/Weyls and radial sines1:12. Only an Nr=1 seed is retained.
Geometry, angular actions, polar panels and spin/gauge conventions are reused.
These are finite form-domain trial fields, not a source or resolvent solve.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from dihole_enriched_bridge import bridge_data, column_metadata as base_metadata
from dihole_mouth_fields import bridge_chunks as old_bridge_chunks
from dihole_radial_stream import FIELD_KEYS, RADIAL_PANELS, SIGMA, _integer
from dihole_stretched_collar import gauss
from project_paths import project_path


OLD_COLUMNS, CANDIDATE_COLUMNS = 324, 272


def column_metadata():
    """New-column labels only; index is relative to this272-column pool."""
    rows = []
    for group, angular, radial in (('radial', 5, range(13, 21)),
                                    ('angular', 7, range(1, 13))):
        for row in base_metadata(.5, angular, 1, (1, 1))['columns']:
            if row['family'] != 'compact' or (group == 'angular' and row['l'] < 6):
                continue
            for mode in radial:
                index = len(rows)
                rows.append(dict(index=index, group=group, chirality=row['chirality'],
                    spin=row['spin'], l=row['l'], radial_mode=mode,
                    local_index=OLD_COLUMNS+index, global_index=356+index))
    return rows


def local_to_global():
    """The32 exterior-only old columns312:344 do not occur on this bridge."""
    return np.r_[np.arange(312), np.arange(344, 356), np.arange(356, 628)]


def _radial_factors(s, modes):
    s, modes = np.asarray(s), np.asarray(modes)
    angle = modes[None, :]*np.pi*(s[:, None]+.1)/.2
    value = np.sqrt(10.)*np.sin(angle)
    # The endpoint trace is analytically zero, not a floating sin(n*pi).
    value[(s == -.1) | (s == .1), :] = 0.
    derivative = np.sqrt(10.)*modes[None, :]*np.pi/.2*np.cos(angle)
    return value, derivative


def _seed(order):
    """Recover only compact Nr=1 profiles/actions, discarding global fields."""
    metadata = base_metadata(.5, 7, 1, (1, 1))
    baseline = bridge_data(.5, 0., 7, 1, order, (1, 1))
    nodes = gauss(order)[0]
    s = np.concatenate([(hi+lo)/2+(hi-lo)*nodes/2 for lo, hi in RADIAL_PANELS])
    first = _radial_factors(s, [1])[0][:, 0]
    if np.any(first <= np.finfo(float).tiny):
        raise ArithmeticError('The interior first sine cannot safely recover profiles.')
    count, start, angles = metadata['compact_columns'], metadata['global_columns'], 3*order
    expected = (len(s)*angles, 4, start+count)
    if any(baseline[key].shape != expected for key in FIELD_KEYS):
        raise ValueError('Frozen Nr=1 ordering changed.')
    # Angular values are independent of s. Actions retain the actual collar
    # coefficients at every slice; no angular differential formula is copied.
    profiles = baseline['value'][:angles, :, start:].copy()/first[0]
    angular = baseline['angular'][:, :, start:].copy()
    lookup = {(r['chirality'], r['spin'], r['l']): r['index']-start
              for r in metadata['columns'] if r['family'] == 'compact'}
    return s, first, profiles, angular, lookup


def bridge_chunks(order=64, slices_per_chunk=4):
    """Yield standard assemble data, old324 then272 candidates, fixed phase0.

    Each yielded field has shape(points,4,596). Scalar data are unchanged
    views from mouth_fields and must not be mutated. No exterior is added.
    """
    order, slices_per_chunk = _integer(order, 'order'), _integer(slices_per_chunk, 'slices_per_chunk')
    if order <= 20:
        raise ValueError('order must exceed the highest candidate radial mode20.')
    metadata = column_metadata()
    s, first, profiles, angular, lookup = _seed(order)
    indices = np.array([lookup[(r['chirality'], r['spin'], r['l'])] for r in metadata])
    modes = np.array([r['radial_mode'] for r in metadata])
    angles, consumed = 3*order, 0
    selected_profiles = profiles[:, :, indices]
    for old in old_bridge_chunks(phase=0., order=order, slices_per_chunk=slices_per_chunk):
        ns, remainder = divmod(len(old['P']), angles)
        if remainder or old['value'].shape[-1] != OLD_COLUMNS:
            raise ValueError('Frozen157 chunk shape/order changed.')
        hi = consumed+ns
        rows = slice(consumed*angles, hi*angles)
        sine, derivative = _radial_factors(s[consumed:hi], modes)
        action = angular[rows].reshape(ns, angles, 4, -1)/first[consumed:hi, None, None, None]
        candidates = dict(
            value=selected_profiles[None, ...]*sine[:, None, None, :],
            radial=-1j*SIGMA[None, None, :, None]*selected_profiles[None, ...]
                   *derivative[:, None, None, :],
            angular=action[:, :, :, indices]*sine[:, None, None, :])
        result = dict(old)
        for key in FIELD_KEYS:
            result[key] = np.concatenate((old[key], candidates[key].reshape(ns*angles, 4, -1)), axis=2)
        consumed = hi
        yield result
    if consumed != len(s):
        raise ArithmeticError('Incomplete radial stream.')


def self_check():
    """Pointwise/nesting controls only: no Gram, source or exterior solve."""
    from dihole_radial_stream import bridge_chunks as reference_stream, radial_recovery_bounds

    start = perf_counter()
    metadata = column_metadata()
    reference_labels = base_metadata(.5, 7, 20, (1, 1))['columns']
    lookup = {(r['chirality'], r['spin'], r['l'], r['radial_mode']): r['index']
              for r in reference_labels if r['family'] == 'compact'}
    selected = np.array([lookup[(r['chirality'], r['spin'], r['l'], r['radial_mode'])] for r in metadata])
    paired = {(r['spin'], r['l'], r['radial_mode']): r['index']
              for r in metadata if r['chirality'] == 'L'}
    right = [r['index'] for r in metadata if r['chirality'] == 'R']
    left = [paired[(r['spin'], r['l'], r['radial_mode'])]
            for r in metadata if r['chirality'] == 'R']
    old_exact, scalar_exact, pairing, points, largest = True, True, 0., 0, 0
    errors = {key: 0. for key in FIELD_KEYS}
    order = 24
    streams = (bridge_chunks(order, 3), old_bridge_chunks(0., order, slices_per_chunk=3),
               reference_stream(.5, 0., 7, 20, order, (1, 1), slices_per_chunk=3))
    from itertools import zip_longest
    for new, old, reference in zip_longest(*streams):
        if new is None or old is None or reference is None:
            raise ArithmeticError('Stream lengths differ.')
        points += len(new['P'])
        largest = max(largest, sum(v.nbytes for v in new.values()))
        scalar_exact &= all(np.array_equal(old[k], new[k]) for k in old if k not in FIELD_KEYS)
        for key in FIELD_KEYS:
            old_exact &= np.array_equal(old[key], new[key][:, :, :OLD_COLUMNS])
            candidate = new[key][:, :, OLD_COLUMNS:]
            expected = reference[key][:, :, selected]
            errors[key] = max(errors[key], float(np.max(abs(candidate-expected))/(1+np.max(abs(expected)))))
            sign = 1 if key == 'value' else -1
            pairing = max(pairing, float(np.max(abs(candidate[:, :2, right]-sign*candidate[:, 2:, left]))))
    modes = np.arange(1, 21)
    trace, endpoint_derivative = _radial_factors(np.array([-.1, .1]), modes)
    endpoint_expected = np.sqrt(10.)*modes*np.pi/.2*np.array([np.ones(20), (-1.)**modes])
    samples = np.array([-.079, -.013, .046, .093])
    h = 2e-6
    a, b, c, d = [_radial_factors(samples+j*h, modes)[0] for j in (-2, -1, 1, 2)]
    difference = (a-8*b+8*c-d)/(12*h)
    exact = _radial_factors(samples, modes)[1]
    derivative_error = float(np.max(abs(difference-exact))/(1+np.max(abs(exact))))
    mapping = local_to_global()
    recovery = radial_recovery_bounds(64, 20)
    rejected = 0
    for bad in (True, 20, 0, 24.5):
        try:
            next(bridge_chunks(bad))
        except ValueError:
            rejected += 1
    tests = dict(old324_bitwise=old_exact, scalar_data_bitwise=scalar_exact,
        same_streamed_candidates=max(errors.values()) < 3e-13,
        exact_chiral_pairing=pairing == 0., exact_zero_trace=np.all(trace == 0.),
        endpoint_derivatives=np.max(abs(endpoint_derivative-endpoint_expected)) < 3e-12,
        radial_derivative=derivative_error < 2e-9,
        candidate_counts=len(metadata) == 272 and sum(r['group'] == 'radial' for r in metadata) == 176,
        unique_labels=len(set(lookup_key(r) for r in metadata)) == 272,
        old_global_mapping=np.array_equal(mapping[:324], np.r_[np.arange(312), np.arange(344, 356)]),
        candidate_mapping=np.array_equal(mapping[324:], np.arange(356, 628)),
        all_points=points == 9*order**2,
        production_recovery=recovery['minimum_sine1'] > 0 and recovery['maximum_sine_ratio_divided_by_mode'] < 1+2e-12,
        invalid_orders=rejected == 4)
    return dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        seconds=perf_counter()-start, order=order, points=points, largest_yielded_chunk_bytes=largest,
        candidate_relative_errors=errors, chiral_error=pairing, radial_derivative_error=derivative_error,
        endpoint_derivative_error=float(np.max(abs(endpoint_derivative-endpoint_expected))),
        production_recovery=recovery, columns=metadata, local_to_global=mapping.tolist(),
        scope=['Fixed157 geometry, mu=1/2 and phase0; no exterior added.',
               'Sine endpoints are exactly zero; normal derivatives need not vanish for the first-order form domain.',
               'All profiles/actions reuse the existing regular monopole and actual angular-coefficient helpers.',
               'Floating pointwise controls; no form, source, residual reduction, or continuum convergence claimed.'])


def lookup_key(row):
    return tuple(row[k] for k in ('chirality', 'spin', 'l', 'radial_mode'))


if __name__ == '__main__':
    result = self_check()
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['input_sha256'] = {name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
        for name in ('dihole_mouth_fields.py', 'dihole_radial_stream.py', 'dihole_enriched_bridge.py',
                     'dihole_azimuthal_bridge.py', 'dihole_stretched_collar.py', 'project_paths.py')}
    project_path('dihole_residual_candidates_checks.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('columns', 'local_to_global')}, indent=2))
    if not result['passed']:
        raise SystemExit('Residual candidate field checks failed.')
