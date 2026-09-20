"""Refine only independent radial modes in the dominant sewn Dirac sector.

Frozen139 angular/exterior families are retained. Streaming reuses the Nr=1
pointwise data; nested raw subspaces and both chiralities retain every mixed
term. Frequency bounds enclose the saved finite dyadic models only.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson

from dihole_bridge_enrichment import canonical, raw_polynomial, select_columns
from dihole_enriched_bridge import column_metadata
from dihole_euclidean_form import relative_error
from dihole_finite_frequency_tail import coefficient_bounds, high_tail
from dihole_radial_schur import source_change
from dihole_radial_stream import bridge_chunks
from dihole_sewn_form import exterior_data, fibre_response
from dihole_stretched_collar import gauss
from project_paths import project_path


def stream_form(phase, radial_modes=12, order=64):
    result = None
    for chunk in bridge_chunks(.5, phase, 3, radial_modes, order, (3, 4)):
        value = raw_polynomial(chunk)
        if result is None:
            result = value
        else:
            for key in result:
                result[key] += value[key]
    return result


def chiral_blocks(raw, metadata, radial_modes):
    indices = select_columns(metadata, (3, 4), 3, radial_modes)
    right, left = [], []
    for index in indices:
        row = metadata['columns'][index]
        is_right = row['component'] < 2 if row['family'] == 'global' else row['chirality'] == 'R'
        (right if is_right else left).append(index)
    for key in ('norm', 'q', 'derivative'):
        if np.any(raw[key][..., right, :][..., left] != 0):
            raise ArithmeticError('Nonzero cross-chirality block; splitting would change this matrix model.')
    blocks = [canonical(raw, np.array(which)) for which in (right, left)]
    pairing = max(relative_error(blocks[0][key][j], (-1)**j*blocks[1][key][j])
                  for key in ('q', 'derivative') for j in range(3))
    exactly_paired = all(np.array_equal(blocks[0][key][j], (-1)**j*blocks[1][key][j])
                        for key in ('q', 'derivative') for j in range(3))
    return indices, blocks, dict(pairing_error=pairing, exactly_paired=exactly_paired,
                                block_dimensions=[len(right), len(left)])


def even_response(blocks, kappa, exactly_paired=False):
    if exactly_paired:
        return sum(fibre_response(block, kappa) for block in blocks)
    # Do not silently impose an exact symmetry on separately rounded blocks.
    return .5*sum(fibre_response(block, sign*kappa) for block in blocks for sign in (-1, 1))


def frequency_integral(blocks, exactly_paired=False, target=1e-7):
    bounds = [coefficient_bounds(block) for block in blocks]
    endpoint = 10.
    for _ in range(30):
        try:
            tails = [high_tail(bound, endpoint) for bound in bounds]
        except ValueError:
            tails = [np.inf, np.inf]
        if sum(tails) <= target:
            break
        endpoint *= 2
    else:
        raise ArithmeticError('Finite-matrix frequency tail target not reached.')
    cuts = [0., .01, .1, 1., 10.]
    while cuts[-1] < endpoint:
        cuts.append(min(10*cuts[-1], endpoint))
    values = []
    for order in (16, 24, 32, 48):
        x, w = gauss(order)
        total = 0.
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            points = (hi+lo)/2+(hi-lo)*x/2
            total += (hi-lo)/2*sum(weight*even_response(blocks, k, exactly_paired)
                                  for weight, k in zip(w, points))
        values.append(dict(order=order, integral=float(total/np.pi)))
        if len(values) > 1 and abs(values[-1]['integral']-values[-2]['integral']) < 1e-10:
            break
    return dict(endpoint=endpoint, coefficient_bounds=bounds, block_tail_bounds=tails,
                certified_matrix_tail=float(np.nextafter(sum(tails), np.inf)), middle_values=values,
                middle_integral=values[-1]['integral'],
                middle_change=abs(values[-1]['integral']-values[-2]['integral']),
                zero_frequency=even_response(blocks, 0., exactly_paired),
                scope='Exact saved canonical dyadic blocks only; middle integration, assembly/whitening, and continuum errors are not enclosed.')


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), layout='constrained')
    for phase, label in ((0., 'Phase 0'), (float(np.pi), 'Phase pi')):
        rows = [r for r in result['records'] if r['phase'] == phase]
        axes[0].plot([r['radial_modes'] for r in rows], [r['integral_33'] for r in rows], marker='o', label=label)
        complete = [r for r in rows if 'full_frequency' in r]
        axes[1].plot([r['radial_modes'] for r in complete],
                     [r['full_frequency']['middle_integral'] for r in complete], marker='o', label=label)
    axes[0].set(title='Same finite window as139', xlabel='Radial sine modes',
                ylabel='Response over kappa D = .01 to 10')
    axes[1].set(title='Including the low-frequency contribution', xlabel='Radial sine modes',
                ylabel='Numerical response from zero to K')
    for ax in axes:
        ax.axhline(0, color='grey', lw=.7)
        ax.legend()
    fig.suptitle('Radial refinement at fixed angular/exterior space; supporting stress remains unknown', fontsize=12)
    fig.supxlabel('Responses in hbar c / D, MD=1. Finite-matrix frequency tails are separate; field-space errors remain.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_radial_refinement_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    baseline = json.loads(project_path('dihole_bridge_enrichment_results.json').read_text())
    exterior = raw_polynomial(exterior_data(.5, 24, (3, 4)))
    metadata = column_metadata(.5, 3, 12, (3, 4))
    records, controls, matrices = [], [], {}
    kappas = np.geomspace(.01, 10., 33)
    for pi, phase in enumerate((0., np.pi)):
        raw = stream_form(phase)
        for key in raw:
            raw[key][..., :48, :48] += exterior[key]
            matrices[f'p{pi}_{key}'] = raw[key]
        weak = raw['weak'][0]+.7*raw['weak'][1]
        controls.append(dict(name='global_sewing', phase=phase, skew_error=relative_error(weak, -weak.conj().T)))
        print(json.dumps(dict(assembled_phase=phase, dimension=len(raw['norm']), seconds=perf_counter()-start)), flush=True)
        for number in (3, 4, 6, 8, 10, 12):
            indices, blocks, pairing = chiral_blocks(raw, metadata, number)
            values = [even_response(blocks, k, pairing['exactly_paired']) for k in kappas]
            full = canonical(raw, indices)
            comparison = max(abs(even_response(blocks, k, pairing['exactly_paired'])
                                 -.5*(fibre_response(full, k)+fibre_response(full, -k))) for k in (.1, .7))
            row = dict(phase=phase, radial_modes=number, dimension=len(indices), indices=indices.tolist(),
                       pairing=pairing, split_response_error=comparison,
                       kappas=kappas.tolist(), responses=values,
                       integral_33=float(simpson(values, x=kappas)/np.pi),
                       integral_17=float(simpson(values[::2], x=kappas[::2])/np.pi))
            if number == 3:
                old = next(r for r in baseline['records'] if r['mu'] == .5 and r['phase'] == phase and r['name'] == 'E48_B3x3')
                row['baseline_relative_error'] = relative_error(np.array(values), np.array(old['responses']))
            if number in (6, 8, 10, 12):
                for bi, block in enumerate(blocks):
                    for key in block:
                        matrices[f'p{pi}_n{number}_b{bi}_{key}'] = np.array(block[key])
                row['full_frequency'] = frequency_integral(blocks, pairing['exactly_paired'])
            records.append(row)
            print(json.dumps(dict(phase=phase, radial_modes=number, window_response=row['integral_33'],
                                  seconds=perf_counter()-start)), flush=True)
        keep = select_columns(metadata, (3, 4), 3, 3)
        for k in (.1, .7):
            controls.append(dict(name='radial_schur', phase=phase, kappa=k, **source_change(raw, k, keep)))
    refined = stream_form(0., order=80)
    ext = raw_polynomial(exterior_data(.5, 32, (3, 4)))
    for key in refined:
        refined[key][..., :48, :48] += ext[key]
    original = {key: matrices['p0_'+key] for key in refined}
    controls.append(dict(name='coupled_quadrature', **{key+'_change': relative_error(refined[key], original[key])
                                                       for key in ('q', 'norm', 'derivative')}))
    indices = select_columns(metadata, (3, 4), 3, 12)
    controls[-1]['fibre_change'] = abs(fibre_response(canonical(refined, indices), .1)
                                      -fibre_response(canonical(original, indices), .1))
    tests = dict(global_sewing=max(r['skew_error'] for r in controls if r['name'] == 'global_sewing') < 2e-9,
        frozen_radial3_recovery=max(r['baseline_relative_error'] for r in records if r['radial_modes'] == 3) < 2e-8,
        chirality_split=max(r['split_response_error'] for r in records) < 2e-10,
        radial_schur=max(r['error'] for r in controls if r['name'] == 'radial_schur') < 2e-9,
        coupled_quadrature=max(controls[-1][k+'_change'] for k in ('q', 'norm', 'derivative')) < 2e-7,
        matrix_frequency_tail=all(r['full_frequency']['certified_matrix_tail'] <= 1e-7 for r in records if 'full_frequency' in r),
        middle_frequency_quadrature=all(r['full_frequency']['middle_change'] < 1e-10 for r in records if 'full_frequency' in r))
    np.savez_compressed(project_path('dihole_radial_refinement_state.npz'), **matrices)
    names = ('dihole_radial_stream.py', 'dihole_radial_schur.py', 'dihole_bridge_enrichment.py',
             'dihole_enriched_bridge.py', 'dihole_finite_frequency_tail.py', 'dihole_sewn_form.py',
             'dihole_euclidean_form.py', 'dihole_stretched_collar.py', 'dihole_interval_arithmetic.py')
    result = dict(complete=True, passed=all(tests.values()), tests=tests, controls=controls, records=records,
        metadata=metadata, seconds=perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
        baseline_sha256=hashlib.sha256(project_path('dihole_bridge_enrichment_results.json').read_bytes()).hexdigest(),
        state_sha256=hashlib.sha256(project_path('dihole_radial_refinement_state.npz').read_bytes()).hexdigest(),
        scope=['Only radial modes refined, with angular3/exterior48 in central mu=1/2; other sectors remain139.',
               'Exact cross terms, fixed Euclidean measure and physical proper-mass prescription retained.',
               'Outward tails cover saved dyadic chiral blocks only; numerical middle integration is not enclosed.',
               'Radial convergence differences and finite Schur accounting are not infinite radial remainder bounds.',
               'No continuum source, regulator removal, physical supporting material, stability or backward passage proved.'])
    project_path('dihole_radial_refinement_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, seconds=result['seconds']), indent=2))
    if not result['passed']:
        raise RuntimeError('Radial refinement checks failed; evidence retained.')
    plot(result)
    return result


if __name__ == '__main__':
    run()
