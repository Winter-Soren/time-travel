"""One exterior-only source challenge, reusing the frozen radial-12 bridge.

Only new exterior columns and their cross terms are assembled. This tests
one finite trial extension at MD=1; it does not complete the exterior.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigvalsh

from dihole_bridge_enrichment import canonical, raw_polynomial
from dihole_euclidean_form import relative_error
from dihole_exterior_extension import extension_data, extension_metadata
from dihole_radial_refinement import even_response, frequency_integral
from dihole_radial_schur import source_change
from dihole_sewn_form import fibre_response
from project_paths import project_path


def append_exterior(old, extra):
    """The old 216 block stays exact; only old global48 overlap new fields."""
    n, added = len(old['norm']), len(extra['norm'])-48
    result = {}
    if np.any(extra['derivative'] != 0):
        raise ArithmeticError('Exterior fields overlap the metric probe.')
    for key in old:
        value = np.zeros(old[key].shape[:-2]+(n+added, n+added), dtype=old[key].dtype)
        value[..., :n, :n] = old[key]
        value[..., :48, n:] = extra[key][..., :48, 48:]
        value[..., n:, :48] = extra[key][..., 48:, :48]
        value[..., n:, n:] = extra[key][..., 48:, 48:]
        if not np.array_equal(value[..., :n, :n], old[key]):
            raise ArithmeticError('Frozen old block changed.')
        result[key] = value
    return result


def chiral_blocks(raw, old_columns, new_columns):
    right, left = [], []
    for index, row in enumerate(old_columns+new_columns):
        is_right = row['component'] < 2 if 'component' in row else row['chirality'] == 'R'
        (right if is_right else left).append(index)
    for key in ('norm', 'q', 'derivative'):
        if np.any(raw[key][..., right, :][..., left] != 0):
            raise ArithmeticError('Nonzero raw cross-chirality block.')
    blocks = [canonical(raw, np.array(indices)) for indices in (right, left)]
    paired = all(np.array_equal(blocks[0][key][j], (-1)**j*blocks[1][key][j])
                 for key in ('q', 'derivative') for j in range(3))
    return blocks, paired


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), layout='constrained')
    ax = axes[0]
    ax.add_patch(Rectangle((.005, -.85), .495, 1.7, color='#9bbdd8', alpha=.5, label='Existing exterior support'))
    for lo in (-1.2, .6):
        ax.add_patch(Rectangle((.005, lo), .495, .6, fill=False, edgecolor='#bd5429', lw=2,
                               label='Added exterior fields' if lo > 0 else None))
    for z in (-.5, .5):
        ax.add_patch(Circle((0, z), .05, color='white', ec='black'))
    ax.set(xlim=(-.04, .58), ylim=(-1.3, 1.3), xlabel='rho / D', ylabel='z / D',
           title='New fields extend beyond the old axial support')
    ax.legend(loc='upper right', fontsize=8)
    ax = axes[1]
    for i, row in enumerate(result['records']):
        ax.plot([0, 1], [row['baseline_integral'], row['full_frequency']['middle_integral']],
                marker='o', label='Phase 0' if i == 0 else 'Phase pi')
    ax.set(xticks=[0, 1], xticklabels=['Radial 12 baseline', '+ exterior fields'],
           ylabel='Numerical response in hbar c / D', title='Change comes entirely through field coupling')
    ax.legend()
    fig.suptitle('One exterior extension at MD=1; physical supporting stress remains undetermined', fontsize=12)
    fig.supxlabel('Both ends are extended. Finite frequency-tail bounds do not enclose continuum or assembly errors.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_exterior_feedback_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    baseline_path = project_path('dihole_radial_refinement_results.json')
    baseline = json.loads(baseline_path.read_text())
    state_path = project_path('dihole_radial_refinement_state.npz')
    if hashlib.sha256(state_path.read_bytes()).hexdigest() != baseline['state_sha256']:
        raise ArithmeticError('Frozen140 state hash mismatch.')
    with np.load(state_path) as saved:
        old_raw = [{key: saved[f'p{pi}_{key}'] for key in ('norm', 'q', 'derivative', 'weak')}
                   for pi in (0, 1)]
    new_columns = extension_metadata()
    extra = raw_polynomial(extension_data())
    refined = raw_polynomial(extension_data(order=32))
    controls = dict(extension_quadrature={key: relative_error(extra[key][..., 48:, :], refined[key][..., 48:, :])
                                          for key in ('norm', 'q', 'weak')},
                    probe_disjoint=bool(np.all(extra['derivative'] == 0)),
                    cross_norm=float(np.linalg.norm(extra['norm'][:48, 48:])))
    records, matrices = [], {}
    for pi, phase in enumerate((0., float(np.pi))):
        raw = append_exterior(old_raw[pi], extra)
        blocks, paired = chiral_blocks(raw, baseline['metadata']['columns'], new_columns)
        gram_eigenvalues = eigvalsh(raw['norm'])
        sample_kappas = np.r_[0., np.geomspace(.01, 160., 17)]
        minimum_form = min(float(eigvalsh(sum((sign*k)**j*a for j, a in enumerate(block['q'])))[0])
                           for block in blocks for k in sample_kappas for sign in (-1, 1))
        sewing = max(relative_error(raw['weak'][0]+k*raw['weak'][1],
                                    -(raw['weak'][0]+k*raw['weak'][1]).conj().T) for k in (0., .1, .7))
        row = dict(phase=phase, dimension=len(raw['norm']), added_columns=len(new_columns),
                   exactly_paired=paired, gram_minimum=float(gram_eigenvalues[0]),
                   gram_condition=float(gram_eigenvalues[-1]/gram_eigenvalues[0]),
                   sampled_form_minimum=minimum_form, sewing_skew=sewing,
                   frozen_old_blocks=all(np.array_equal(raw[key][..., :216, :216], old_raw[pi][key]) for key in raw))
        for key in raw:
            matrices[f'p{pi}_{key}'] = raw[key]
        for bi, block in enumerate(blocks):
            for key in block:
                matrices[f'p{pi}_b{bi}_{key}'] = np.array(block[key])
        row['full_frequency'] = frequency_integral(blocks, paired)
        old12 = next(r for r in baseline['records'] if r['phase'] == phase and r['radial_modes'] == 12)
        old10 = next(r for r in baseline['records'] if r['phase'] == phase and r['radial_modes'] == 10)
        row['baseline_integral'] = old12['full_frequency']['middle_integral']
        row['exterior_change'] = row['full_frequency']['middle_integral']-row['baseline_integral']
        row['last_radial_change'] = abs(row['baseline_integral']-old10['full_frequency']['middle_integral'])
        row['change_over_radial'] = abs(row['exterior_change'])/row['last_radial_change']
        row['difference_frequency_tail'] = float(np.nextafter(
            old12['full_frequency']['certified_matrix_tail']+row['full_frequency']['certified_matrix_tail'], np.inf))
        row['schur'] = [dict(kappa=k, **source_change(raw, k, np.arange(216))) for k in (.1, .7)]
        higher = append_exterior(old_raw[pi], refined)
        refined_blocks, refined_pairing = chiral_blocks(higher, baseline['metadata']['columns'], new_columns)
        row['quadrature_fibre_change'] = max(abs(even_response(blocks, k, paired)
                                                -even_response(refined_blocks, k, refined_pairing)) for k in (.1, .7))
        full = canonical(raw, np.arange(len(raw['norm'])))
        row['split_response_error'] = max(abs(even_response(blocks, k, paired)
                                              -.5*(fibre_response(full, k)+fibre_response(full, -k))) for k in (.1, .7))
        # This comparison only includes finite frequency tails; other errors
        # remain numerical controls, so it is not a physical certificate.
        lo = abs(row['exterior_change'])-row['difference_frequency_tail']
        hi = abs(row['exterior_change'])+row['difference_frequency_tail']
        row['decision'] = ('prioritize_exterior' if lo > row['last_radial_change'] else
                           'this_family_below_radial_change' if hi < row['last_radial_change'] else
                           'frequency_tail_comparison_unresolved')
        records.append(row)
        print(json.dumps(dict(phase=phase, exterior_change=row['exterior_change'],
                              ratio=row['change_over_radial'], seconds=perf_counter()-start)), flush=True)
    tests = dict(frozen_core_and_disjoint_probe=controls['probe_disjoint'] and all(r['frozen_old_blocks'] for r in records),
        positive_sampled_forms=all(r['gram_minimum'] > 0 and r['sampled_form_minimum'] > 0 for r in records),
        global_sewing=max(r['sewing_skew'] for r in records) < 2e-9,
        mixed_quadrature=max(controls['extension_quadrature'].values()) < 2e-8
            and max(r['quadrature_fibre_change'] for r in records) < 2e-9,
        chiral_split=max(r['split_response_error'] for r in records) < 2e-10,
        pure_feedback=all(s['direct_high'] == 0. and s['error'] < 2e-9 for r in records for s in r['schur']),
        nonzero_overlap=controls['cross_norm'] > 0,
        frequency_control=all(r['full_frequency']['certified_matrix_tail'] <= 1e-7
            and r['full_frequency']['middle_change'] < 1e-10 for r in records))
    np.savez_compressed(project_path('dihole_exterior_feedback_state.npz'), **matrices)
    names = ('dihole_exterior_extension.py', 'dihole_radial_refinement.py', 'dihole_radial_schur.py',
             'dihole_bridge_enrichment.py', 'dihole_sewn_form.py', 'dihole_finite_frequency_tail.py')
    result = dict(complete=True, passed=all(tests.values()), tests=tests, controls=controls, records=records,
        new_columns=new_columns, seconds=perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
        baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        baseline_state_sha256=baseline['state_sha256'],
        state_sha256=hashlib.sha256(project_path('dihole_exterior_feedback_state.npz').read_bytes()).hexdigest(),
        scope=['One axial exterior family only, fixed mu=.5, angular3, radial12, MD1; old216 block preserved.',
               'Positive sampled pencils, quadrature differences and middle integrals are numerical controls, not interval proofs.',
               'Frequency-tail bounds cover exact saved dyadic canonical matrices only; physical exterior/axis/angular/IR errors remain.',
               'No regulator removal, renormalized source, supporting material, stability, new law or backward passage established.'])
    project_path('dihole_exterior_feedback_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, seconds=result['seconds']), indent=2))
    if not result['passed']:
        raise ArithmeticError('Exterior challenge failed checks; evidence saved.')
    plot(result)
    return result


if __name__ == '__main__':
    run()
