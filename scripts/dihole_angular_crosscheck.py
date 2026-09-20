"""Angular 3/4/5 crosscheck at fixed radial12 and exterior48.

The 140 angular3 frequency result is reused with authenticated input bytes.
New raw subspaces are restricted before whitening, and both chiralities are
evaluated. Finite dyadic-matrix frequency tails do not enclose assembly,
middle quadrature, continuum, omitted-sector or regulator-limit errors.
"""
import gc
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson
from scipy.linalg import eigvalsh

from dihole_bridge_enrichment import canonical, raw_polynomial, select_columns
from dihole_enriched_bridge import column_metadata
from dihole_euclidean_form import relative_error
from dihole_radial_refinement import even_response, frequency_integral
from dihole_radial_stream import bridge_chunks
from dihole_sewn_form import exterior_data, fibre_response
from project_paths import project_path


NAME = 'dihole_angular_crosscheck'
INPUTS = ('dihole_radial_refinement.py', 'dihole_radial_stream.py',
          'dihole_bridge_enrichment.py', 'dihole_enriched_bridge.py',
          'dihole_finite_frequency_tail.py', 'dihole_sewn_form.py',
          'dihole_sewn_quadrature.py', 'dihole_euclidean_form.py',
          'dihole_stretched_collar.py', 'dihole_interval_arithmetic.py',
          'dihole_azimuthal_bridge.py', 'dihole_bridge_scattering.py',
          'dihole_curved_exterior.py', 'balanced_dihole.py', 'project_paths.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stream_form(phase, order=64):
    """Only this angular argument differs from the frozen 140 assembler."""
    result = None
    for chunk in bridge_chunks(.5, phase, 5, 12, order, (3, 4)):
        value = raw_polynomial(chunk)
        if result is None:
            result = value
        else:
            for key in result:
                result[key] += value[key]
    return result


def restrict_chiral(raw, metadata, angular_modes):
    """Select labelled raw fields before either Cholesky normalization."""
    indices = select_columns(metadata, (3, 4), angular_modes, 12)
    right, left = [], []
    for index in indices:
        row = metadata['columns'][index]
        is_right = (row['component'] < 2 if row['family'] == 'global'
                    else row['chirality'] == 'R')
        (right if is_right else left).append(index)
    cross = max(float(np.max(abs(raw[key][..., right, :][..., left])))
                for key in ('norm', 'q', 'derivative'))
    if cross != 0:
        raise ArithmeticError('Nonzero chiral cross block; split would change the model.')
    blocks = [canonical(raw, np.asarray(which)) for which in (right, left)]
    pairing = max(relative_error(blocks[0][key][j], (-1)**j*blocks[1][key][j])
                  for key in ('q', 'derivative') for j in range(3))
    exact = all(np.array_equal(blocks[0][key][j], (-1)**j*blocks[1][key][j])
                for key in ('q', 'derivative') for j in range(3))
    gram = eigvalsh(raw['norm'][np.ix_(indices, indices)])
    return indices, blocks, dict(exact_cross_zero=True, cross_maximum=cross,
        pairing_error=pairing, exactly_paired=exact,
        block_dimensions=[len(right), len(left)], gram_minimum=float(gram[0]),
        gram_condition=float(gram[-1]/gram[0]))


def save_blocks(matrices, prefix, blocks):
    for bi, block in enumerate(blocks):
        for key in block:
            matrices[f'{prefix}_b{bi}_{key}'] = np.asarray(block[key])


def checkpoint(result, matrices, start, stage):
    result.update(stage=stage, seconds=perf_counter()-start)
    if matrices:
        target = project_path(NAME+'_state.npz')
        np.savez_compressed(target, **matrices)
        result['state_sha256'] = sha(target)
    project_path(NAME+'_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def angular_differences(records, baseline):
    answer = []
    for phase in (0., float(np.pi)):
        rows = {r['angular_modes']: r for r in records if r['phase'] == phase}
        old10 = next(r for r in baseline['records']
                     if r['phase'] == phase and r['radial_modes'] == 10)
        old12 = next(r for r in baseline['records']
                     if r['phase'] == phase and r['radial_modes'] == 12)
        radial = abs(old12['full_frequency']['middle_integral']
                     - old10['full_frequency']['middle_integral'])
        radial_tail = float(np.nextafter(sum(r['full_frequency']['certified_matrix_tail']
                                             for r in (old10, old12)), np.inf))
        for low, high in ((3, 4), (4, 5)):
            before, after = [rows[a]['full_frequency'] for a in (low, high)]
            change = after['middle_integral']-before['middle_integral']
            tail = float(np.nextafter(before['certified_matrix_tail']
                                      + after['certified_matrix_tail'], np.inf))
            lower = max(0., abs(change)-tail)
            upper = abs(change)+tail
            radial_lower = max(0., radial-radial_tail)
            radial_upper = radial+radial_tail
            comparison = ('angular_change_exceeds_last_radial' if lower > radial_upper
                else 'angular_change_below_last_radial' if upper < radial_lower
                else 'finite_frequency_tail_comparison_unresolved')
            answer.append(dict(phase=phase, angular_from=low, angular_to=high,
                signed_change=change, absolute_change=abs(change),
                last_radial_change=radial, ratio_to_last_radial=abs(change)/radial,
                difference_frequency_tail=tail, last_radial_frequency_tail=radial_tail,
                middle_quadrature_change_sum=before['middle_change']+after['middle_change'],
                last_radial_middle_change_sum=sum(r['full_frequency']['middle_change']
                                                   for r in (old10, old12)),
                finite_tail_only_absolute_change_interval=[lower, upper],
                finite_tail_only_radial_interval=[radial_lower, radial_upper],
                comparison=comparison,
                scope='Comparison includes finite-matrix frequency tails only; numerical middle and field-space errors are not enclosed.'))
    return answer


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout='constrained')
    for phase, label, color in ((0., 'Phase 0', '#136f8a'),
                                (float(np.pi), 'Phase pi', '#b85825')):
        rows = [r for r in result['records'] if r['phase'] == phase]
        axes[0].plot([r['angular_modes'] for r in rows],
                     [r['full_frequency']['middle_integral'] for r in rows],
                     'o-', label=label, color=color)
        differences = [r for r in result['comparisons'] if r['phase'] == phase]
        offset = -.04 if phase == 0 else .04
        axes[1].errorbar(np.asarray([1., 2.])+offset,
            [r['absolute_change'] for r in differences],
            yerr=[r['difference_frequency_tail'] for r in differences],
            marker='o', capsize=4, label=label, color=color)
    radial = result['comparisons'][0]['last_radial_change']
    axes[1].axhline(radial, color='black', ls='--', lw=1,
                    label='Previous radial 10 to 12 change')
    axes[0].set(xticks=[3, 4, 5], xlabel='Compact angular modes',
                ylabel='Numerical response in hbar c / D',
                title='Radial 12 and exterior 48 held fixed')
    axes[1].set(xticks=[1, 2], xticklabels=['Angular 3 to 4', 'Angular 4 to 5'],
                ylabel='Absolute response change', yscale='log',
                title='Tail bars concern the saved finite matrices')
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle('Angular sensitivity at MD=1; physical supporting stress remains undetermined', fontsize=12)
    fig.supxlabel('Both chiralities, mu=1/2. Middle quadrature, field truncation and regulator removal remain separate.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    source_sha = sha(__file__)
    input_hashes = {name: sha(project_path(name)) for name in INPUTS}
    baseline_path = project_path('dihole_radial_refinement_results.json')
    state_path = project_path('dihole_radial_refinement_state.npz')
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    provenance = (baseline['passed']
        and sha(project_path('dihole_radial_refinement.py')) == baseline['source_sha256']
        and sha(state_path) == baseline['state_sha256']
        and all(sha(project_path(name)) == value for name, value in baseline['input_sha256'].items()))
    if not provenance:
        raise ArithmeticError('Frozen140 input provenance mismatch; no new calculation started.')
    with np.load(state_path) as archive:
        old_raw = [{key: archive[f'p{pi}_{key}'] for key in ('norm', 'q', 'derivative', 'weak')}
                   for pi in (0, 1)]
        old_blocks = [[{key: archive[f'p{pi}_n12_b{bi}_{key}'] for key in ('q', 'derivative')}
                       for bi in (0, 1)] for pi in (0, 1)]
    metadata = column_metadata(.5, 5, 12, (3, 4))
    result = dict(complete=False, passed=False, tests={}, controls=[], records=[], comparisons=[],
        metadata=metadata, parameters=dict(mu=.5, proper_mass_times_D=1., radial_modes=12,
            exterior_columns=48, angular_modes=[3, 4, 5], phases=[0., float(np.pi)],
            bridge_order=64, exterior_order=24, control_bridge_order=80, control_exterior_order=32),
        source_sha256=source_sha, input_sha256=input_hashes,
        baseline_sha256=sha(baseline_path), baseline_state_sha256=sha(state_path),
        scope=['Only compact angular modes are refined at radial12/exterior48 in mu=1/2, with both actual chiralities.',
               'Frozen140 angular3 frequency integrals and their dyadic-matrix tail bounds are reused, after fresh raw/canonical/fibre recovery.',
               'All mixed terms and fixed physical proper-mass prescription are retained; restriction precedes normalization.',
               'Finite frequency tails cover exact saved canonical dyadic matrices only; numerical middle integration and assembly/whitening are not enclosed.',
               'Angular differences and the selected exterior family are not bounds on omitted fields, other azimuthal sectors, or the continuum.',
               'No regulator removal, renormalized supporting source, stability, new physical law or backward passage is established.'])
    matrices = {}
    checkpoint(result, matrices, start, 'source_frozen')
    try:
        exterior = raw_polynomial(exterior_data(.5, 24, (3, 4)))
        kappas = np.geomspace(.01, 10., 33)
        for pi, phase in enumerate((0., float(np.pi))):
            raw = stream_form(phase)
            for key in raw:
                raw[key][..., :48, :48] += exterior[key]
                matrices[f'p{pi}_{key}'] = raw[key]
            sewing = max(relative_error(raw['weak'][0]+k*raw['weak'][1],
                -(raw['weak'][0]+k*raw['weak'][1]).conj().T) for k in (0., .1, .7))
            result['controls'].append(dict(name='global_sewing', phase=phase, skew_error=sewing))
            print(json.dumps(dict(stage='assembled', phase=phase, dimension=len(raw['norm']),
                                  seconds=perf_counter()-start)), flush=True)
            for angular in (3, 4, 5):
                indices, blocks, pairing = restrict_chiral(raw, metadata, angular)
                values = [even_response(blocks, k, pairing['exactly_paired']) for k in kappas]
                row = dict(phase=phase, angular_modes=angular, radial_modes=12,
                    dimension=len(indices), indices=indices.tolist(), pairing=pairing,
                    kappas=kappas.tolist(), responses=values,
                    integral_33=float(simpson(values, x=kappas)/np.pi),
                    integral_17=float(simpson(values[::2], x=kappas[::2])/np.pi))
                row['sampled_form_minimum'] = min(float(eigvalsh(sum((sign*k)**j*a
                    for j, a in enumerate(block['q'])))[0])
                    for block in blocks for k in (0., .1, .7, 1., 10., 160.) for sign in (-1, 1))
                full = canonical(raw, indices)
                row['split_response_error'] = max(abs(even_response(blocks, k, pairing['exactly_paired'])
                    -.5*(fibre_response(full, k)+fibre_response(full, -k))) for k in (.1, .7))
                if angular == 3:
                    prior = next(r for r in baseline['records'] if r['phase'] == phase and r['radial_modes'] == 12)
                    row['baseline_recovery'] = dict(
                        raw_relative_errors={key: relative_error(raw[key][..., indices, :][..., indices], old_raw[pi][key])
                                             for key in raw},
                        canonical_relative_error=max(relative_error(blocks[bi][key], old_blocks[pi][bi][key])
                            for bi in (0, 1) for key in ('q', 'derivative')),
                        fibre_absolute_error=float(np.max(abs(np.asarray(values)-np.asarray(prior['responses'])))),
                        fibre_relative_error=relative_error(np.asarray(values), np.asarray(prior['responses'])))
                    save_blocks(matrices, f'p{pi}_a3_recovered', blocks)
                    save_blocks(matrices, f'p{pi}_a3', old_blocks[pi])
                    row['full_frequency'] = prior['full_frequency']
                    row['frequency_provenance'] = 'Reused exact frozen140 Nr12 canonical blocks/integral/tail; newly recovered blocks are saved separately.'
                else:
                    save_blocks(matrices, f'p{pi}_a{angular}', blocks)
                    row['full_frequency'] = frequency_integral(blocks, pairing['exactly_paired'])
                    row['frequency_provenance'] = 'Fresh integral and dyadic tail from the saved new canonical blocks.'
                result['records'].append(row)
                print(json.dumps(dict(stage='frequency_complete', phase=phase, angular_modes=angular,
                    integral=row['full_frequency']['middle_integral'],
                    tail=row['full_frequency']['certified_matrix_tail'], seconds=perf_counter()-start)), flush=True)
            checkpoint(result, matrices, start, f'phase_{pi}_complete')
            del raw, blocks, full
            gc.collect()
        result['comparisons'] = angular_differences(result['records'], baseline)
        print(json.dumps(dict(stage='quadrature_control_start', seconds=perf_counter()-start)), flush=True)
        refined = stream_form(0., order=80)
        higher_exterior = raw_polynomial(exterior_data(.5, 32, (3, 4)))
        for key in refined:
            refined[key][..., :48, :48] += higher_exterior[key]
            matrices['q80_'+key] = refined[key]
        original = {key: matrices['p0_'+key] for key in refined}
        _, blocks80, pairing80 = restrict_chiral(refined, metadata, 5)
        _, blocks64, pairing64 = restrict_chiral(original, metadata, 5)
        save_blocks(matrices, 'q80_a5', blocks80)
        sample = [.01, .1, .7, 2.]
        quad = dict(name='coupled_quadrature', phase=0., angular_modes=5,
            bridge_orders=[64, 80], exterior_orders=[24, 32],
            matrix_relative_changes={key: relative_error(refined[key], original[key]) for key in refined},
            kappas=sample,
            responses64=[even_response(blocks64, k, pairing64['exactly_paired']) for k in sample],
            responses80=[even_response(blocks80, k, pairing80['exactly_paired']) for k in sample],
            pairing=pairing80)
        quad['fibre_maximum_change'] = float(np.max(abs(np.asarray(quad['responses64'])-np.asarray(quad['responses80']))))
        result['controls'].append(quad)
        rows = result['records']
        recoveries = [r['baseline_recovery'] for r in rows if r['angular_modes'] == 3]
        result['tests'] = dict(
            frozen_provenance=provenance and sha(__file__) == source_sha
                and all(sha(project_path(n)) == h for n, h in input_hashes.items())
                and sha(baseline_path) == result['baseline_sha256'] and sha(state_path) == result['baseline_state_sha256'],
            labelled_nested_dimensions=[r['dimension'] for r in rows] == [216, 264, 312]*2,
            frozen_raw_recovery=max(v for r in recoveries for v in r['raw_relative_errors'].values()) < 2e-11,
            frozen_canonical_recovery=max(r['canonical_relative_error'] for r in recoveries) < 2e-9,
            frozen_33_fibre_recovery=max(r['fibre_absolute_error'] for r in recoveries) < 2e-10
                and max(r['fibre_relative_error'] for r in recoveries) < 2e-8,
            positive_sampled_forms=all(r['pairing']['gram_minimum'] > 0 and r['sampled_form_minimum'] > 0 for r in rows),
            exact_chiral_cross_zero=all(r['pairing']['exact_cross_zero'] for r in rows),
            beta_pairing=max(r['pairing']['pairing_error'] for r in rows) < 2e-11,
            chiral_split=max(r['split_response_error'] for r in rows) < 2e-10,
            global_sewing=max(r['skew_error'] for r in result['controls'] if r['name'] == 'global_sewing') < 2e-9,
            coupled_quadrature=max(quad['matrix_relative_changes'].values()) < 2e-7
                and quad['fibre_maximum_change'] < 2e-9,
            finite_matrix_frequency_tail=all(r['full_frequency']['certified_matrix_tail'] <= 1e-7 for r in rows),
            middle_frequency_quadrature=all(r['full_frequency']['middle_change'] < 1e-10 for r in rows))
        result.update(complete=True, passed=all(result['tests'].values()))
        checkpoint(result, matrices, start, 'complete')
        print(json.dumps(dict(passed=result['passed'], tests=result['tests'],
                              comparisons=result['comparisons'], seconds=result['seconds']), indent=2), flush=True)
        if not result['passed']:
            raise ArithmeticError('Angular crosscheck failed; completed diagnostics retained.')
        plot(result)
        return result
    except Exception as exc:
        if not result['complete']:
            result['failure'] = dict(type=type(exc).__name__, message=str(exc))
            checkpoint(result, matrices, start, 'failed_partial')
        raise


if __name__ == '__main__':
    run()
