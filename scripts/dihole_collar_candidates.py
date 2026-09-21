"""Localized collar and high radial candidates on the frozen158 field stream.

The old596 local columns and scalar arrays are retained exactly. New fields
are132 collar profiles followed by132 global sines21:26, all using the same22
regular angular profiles already present at radial mode1. Their values and
actual angular actions are recovered from each old chunk; no new geometry or
angular operator is evaluated. All new mouth traces vanish. This is a finite
first-order form-domain provider, not a source or residual-reduction result.
"""
import hashlib
import json
from itertools import zip_longest
from pathlib import Path
from time import perf_counter

import numpy as np

from dihole_enriched_bridge import column_metadata as base_metadata
from dihole_radial_stream import FIELD_KEYS, RADIAL_PANELS, SIGMA, _integer
from dihole_residual_candidates import bridge_chunks as old_bridge_chunks, _radial_factors
from dihole_stretched_collar import gauss, smooth_step
from project_paths import project_path


OLD_COLUMNS, CANDIDATE_COLUMNS = 596, 264
WIDTHS = (.0025, .005, .01)


def _profiles():
    return [r for r in base_metadata(.5, 5, 12, (3, 4))['columns']
            if r['family'] == 'compact' and r['radial_mode'] == 1]


def column_metadata():
    """New-column labels only; index is relative to the264-column addition."""
    rows = []
    for side in ('left', 'right'):
        for width in WIDTHS:
            for p in _profiles():
                i = len(rows)
                rows.append(dict(index=i, group='collar', side=side, width=width,
                    chirality=p['chirality'], spin=p['spin'], l=p['l'], radial_mode=None,
                    local_index=OLD_COLUMNS+i, global_index=628+i))
    for p in _profiles():
        for mode in range(21, 27):
            i = len(rows)
            rows.append(dict(index=i, group='global', side=None, width=None,
                chirality=p['chirality'], spin=p['spin'], l=p['l'], radial_mode=mode,
                local_index=OLD_COLUMNS+i, global_index=628+i))
    return rows


def local_to_global():
    """Skip only the32 old exterior-only columns312:344."""
    return np.r_[np.arange(312), np.arange(344, 628), np.arange(628, 892)]


def collar_factor(s, side, width):
    """h(s) and dh/ds, with exact zero trace and flat inner cutoff.

    h=2/sqrt(width)*(d/width)*exp(-d/width)*(1-S(d/.05)).
    The nonzero normal derivative at the mouth is allowed by the form domain.
    """
    if side not in ('left', 'right') or not np.isfinite(width) or width <= 0:
        raise ValueError('Require left/right and a finite positive width.')
    s = np.asarray(s, float)
    if s.ndim != 1 or not np.isfinite(s).all():
        raise ValueError('s must be a finite one-dimensional array.')
    d, sign = (s+.1, 1.) if side == 'left' else (.1-s, -1.)
    value, derivative = np.zeros_like(d), np.zeros_like(d)
    active = (d >= 0) & (d <= .05)
    distance = d[active]
    step, step_prime = smooth_step(distance/.05)
    chi, chip = 1-step, -step_prime/.05
    v = distance/width
    exponential = np.exp(-v)
    scale = 2/np.sqrt(width)
    value[active] = scale*v*exponential*chi
    derivative[active] = sign*scale/width*exponential*((1-v)*chi+distance*chip)
    return value, derivative


def radial_factors(s):
    """All264 radial values/derivatives in column_metadata order."""
    values, derivatives = [], []
    for side in ('left', 'right'):
        for width in WIDTHS:
            h, hp = collar_factor(s, side, width)
            values.append(np.repeat(h[:, None], 22, axis=1))
            derivatives.append(np.repeat(hp[:, None], 22, axis=1))
    global_value, global_derivative = _radial_factors(np.asarray(s), np.arange(21, 27))
    values.append(np.tile(global_value, (1, 22)))
    derivatives.append(np.tile(global_derivative, (1, 22)))
    return np.concatenate(values, axis=1), np.concatenate(derivatives, axis=1)


def bridge_chunks(order=64, slices_per_chunk=4):
    """Standard data with860 local columns; old596 are bitwise unchanged."""
    order, slices_per_chunk = _integer(order, 'order'), _integer(slices_per_chunk, 'slices_per_chunk')
    if order <= 26:
        raise ValueError('order must exceed the highest radial candidate mode26.')
    metadata, seed = column_metadata(), _profiles()
    seed_columns = np.array([r['index'] for r in seed])
    lookup = {(r['chirality'], r['spin'], r['l']): i for i, r in enumerate(seed)}
    indices = np.array([lookup[(r['chirality'], r['spin'], r['l'])] for r in metadata])
    x = gauss(order)[0]
    ss = np.concatenate([(hi+lo)/2+(hi-lo)*x/2 for lo, hi in RADIAL_PANELS])
    first = _radial_factors(ss, [1])[0][:, 0]
    if np.any(first <= np.finfo(float).tiny):
        raise ArithmeticError('Cannot recover profiles from a nonpositive/subnormal first sine.')
    angles, consumed = 3*order, 0
    for old in old_bridge_chunks(order, slices_per_chunk):
        ns, remainder = divmod(len(old['P']), angles)
        if remainder or old['value'].shape[-1] != OLD_COLUMNS:
            raise ValueError('Frozen158 field ordering changed.')
        hi = consumed+ns
        denominator = first[consumed:hi, None, None, None]
        profile = old['value'][:, :, seed_columns].reshape(ns, angles, 4, 22)/denominator
        action = old['angular'][:, :, seed_columns].reshape(ns, angles, 4, 22)/denominator
        h, hp = radial_factors(ss[consumed:hi])
        selected = profile[:, :, :, indices]
        candidates = dict(value=selected*h[:, None, None, :],
            radial=-1j*SIGMA[None, None, :, None]*selected*hp[:, None, None, :],
            angular=action[:, :, :, indices]*h[:, None, None, :])
        result = dict(old)
        for key in FIELD_KEYS:
            result[key] = np.concatenate((old[key], candidates[key].reshape(ns*angles, 4, -1)), axis=2)
        consumed = hi
        yield result
    if consumed != len(ss):
        raise ArithmeticError('Incomplete frozen radial stream.')


def self_check():
    """Focused pointwise controls only; no form, exterior or source solve."""
    from dihole_radial_stream import bridge_chunks as reference_stream, radial_recovery_bounds

    start, metadata = perf_counter(), column_metadata()
    reference_labels = base_metadata(.5, 5, 26, (1, 1))['columns']
    lookup = {(r['chirality'], r['spin'], r['l'], r['radial_mode']): r['index']
              for r in reference_labels if r['family'] == 'compact'}
    selected = [lookup[(r['chirality'], r['spin'], r['l'], r['radial_mode'])]
                for r in metadata if r['group'] == 'global']
    def pair_key(r):
        return tuple(r[k] for k in ('group', 'side', 'width', 'spin', 'l', 'radial_mode'))
    paired = {pair_key(r): r['index'] for r in metadata if r['chirality'] == 'L'}
    right = [r['index'] for r in metadata if r['chirality'] == 'R']
    left = [paired[pair_key(r)] for r in metadata if r['chirality'] == 'R']
    old_exact, scalar_exact, pairing, points, largest = True, True, 0., 0, 0
    errors = {key: 0. for key in FIELD_KEYS}
    order = 28
    streams = (bridge_chunks(order, 3), old_bridge_chunks(order, 3),
               reference_stream(.5, 0., 5, 26, order, (1, 1), slices_per_chunk=3))
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
            errors[key] = max(errors[key], float(np.max(abs(candidate[:, :, 132:]-expected))/(1+np.max(abs(expected)))))
            sign = 1 if key == 'value' else -1
            pairing = max(pairing, float(np.max(abs(candidate[:, :2, right]-sign*candidate[:, 2:, left]))))
    trace, _ = radial_factors(np.array([-.1, .1]))
    support_error, cutoff_error, mouth_derivative_error = 0., 0., 0.
    samples = np.array([.0003, .0015, .003, .006, .014, .028, .043, .049])
    derivative_error, reflection_error = 0., 0.
    h = 2e-7
    for side, sign in (('left', 1.), ('right', -1.)):
        sample_s = sign*(samples-.1)
        for width in WIDTHS:
            values = [collar_factor(sample_s+j*h, side, width)[0] for j in (-2, -1, 1, 2)]
            fd = (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
            exact = collar_factor(sample_s, side, width)[1]
            derivative_error = max(derivative_error, float(np.max(abs(fd-exact))/(1+np.max(abs(exact)))))
            value, derivative = collar_factor(np.array([-sign*.1, -sign*.05, 0., sign*.1]), side, width)
            cutoff_error = max(cutoff_error, float(abs(value[1])), float(abs(derivative[1])))
            support_error = max(support_error, float(np.max(abs(value[2:]))), float(np.max(abs(derivative[2:]))))
            mouth_derivative_error = max(mouth_derivative_error, float(abs(derivative[0]-sign*2/width**1.5)))
            other = 'right' if side == 'left' else 'left'
            mirrored = collar_factor(-sample_s, other, width)
            own = collar_factor(sample_s, side, width)
            reflection_error = max(reflection_error, float(np.max(abs(own[0]-mirrored[0]))), float(np.max(abs(own[1]+mirrored[1]))))
    modes = np.arange(21, 27)
    s = np.array([-.079, -.013, .046, .093])
    values = [_radial_factors(s+j*h, modes)[0] for j in (-2, -1, 1, 2)]
    fd = (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
    exact = _radial_factors(s, modes)[1]
    sine_derivative_error = float(np.max(abs(fd-exact))/(1+np.max(abs(exact))))
    recovery, mapping = radial_recovery_bounds(64, 26), local_to_global()
    rejected = 0
    for bad in (True, 26, 0, 32.5):
        try:
            next(bridge_chunks(bad))
        except ValueError:
            rejected += 1
    tests = dict(old596_bitwise=old_exact, scalar_data_bitwise=scalar_exact,
        global_stream_agreement=max(errors.values()) < 3e-13, exact_chiral_pairing=pairing == 0.,
        exact_mouth_traces=np.all(trace == 0.), flat_inner_cutoff=cutoff_error == 0.,
        exact_collar_support=support_error == 0., mouth_derivative=mouth_derivative_error < 4e-12,
        collar_derivative=derivative_error < 2e-9, sine_derivative=sine_derivative_error < 2e-9,
        left_right_orientation=reflection_error == 0.,
        counts=len(metadata) == 264 and sum(r['group'] == 'collar' for r in metadata) == 132,
        unique_labels=len({(r['chirality'],)+pair_key(r) for r in metadata}) == 264,
        old_mapping=np.array_equal(mapping[:596], np.r_[np.arange(312), np.arange(344, 628)]),
        new_mapping=np.array_equal(mapping[596:], np.arange(628, 892)),
        all_points=points == 9*order**2,
        production_recovery=recovery['minimum_sine1'] > 0 and recovery['maximum_sine_ratio_divided_by_mode'] < 1+2e-12,
        invalid_orders=rejected == 4)
    return dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        seconds=perf_counter()-start, order=order, points=points, largest_yielded_chunk_bytes=largest,
        global_relative_errors=errors, chiral_error=pairing, collar_derivative_error=derivative_error,
        sine_derivative_error=sine_derivative_error, mouth_derivative_error=mouth_derivative_error,
        cutoff_error=cutoff_error, support_error=support_error, reflection_error=reflection_error,
        production_recovery=recovery, columns=metadata, local_to_global=mapping.tolist(),
        scope=['Frozen158 bridge, mu=1/2 and phase0; no new geometry, angular action or exterior fields.',
               'Collar fields vanish at mouths; their nonzero normal derivatives are permitted in the first-order form domain.',
               'The inner collar cutoff is flat and the value/first derivative match the zero extension.',
               'Pointwise floating checks only; no form assembly, source value, residual gain or convergence certificate.'])


if __name__ == '__main__':
    result = self_check()
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['input_sha256'] = {name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
        for name in ('dihole_residual_candidates.py', 'dihole_residual_candidates_checks.json',
                     'dihole_mouth_fields.py', 'dihole_radial_stream.py', 'dihole_enriched_bridge.py',
                     'dihole_stretched_collar.py', 'project_paths.py')}
    project_path('dihole_collar_candidates_checks.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('columns', 'local_to_global')}, indent=2))
    if not result['passed']:
        raise SystemExit('Collar candidate field checks failed.')
