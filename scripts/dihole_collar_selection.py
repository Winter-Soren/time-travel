"""Compare localized collars with more radial sines at one fixed frequency.

Reuse the frozen628-field forms. Only new/cross Gram, full-action Gram and
weak-action blocks are integrated. Saved matrices are explicitly evaluated
at kappa=.7, not frequency polynomials or a physical stress calculation.
"""
import json
from time import perf_counter

import numpy as np
from scipy.linalg import eigh, solve

from dihole_bridge_enrichment import raw_polynomial
from dihole_collar_candidates import bridge_chunks, column_metadata, local_to_global
from dihole_common_pv_residual import hermitian
from dihole_euclidean_form import relative_error
from dihole_mouth_enrichment import solve_case
from dihole_residual_selection import load_problem as load_previous, select_directions, sha
from project_paths import project_path


NAME, PARENT = 'dihole_collar_selection', 'dihole_residual_selection'
OLD, LOCAL_OLD, TOTAL, KAPPA = 628, 596, 892, .7
KEYS = ('gram', 'action_gram', 'weak_action')


def fixed_forms(raw):
    return dict(gram=raw['norm'].copy(),
        action_gram=sum(KAPPA**j*v for j, v in enumerate(raw['q'])),
        weak_action=raw['weak'][0]+KAPPA*raw['weak'][1])


def load_problem():
    old, exterior, indices, source, inputs = load_previous()
    parent = json.loads(project_path(PARENT+'_results.json').read_text(encoding='utf-8'))
    if not parent['passed'] or not parent['complete']:
        raise ArithmeticError('Frozen628 producer did not pass.')
    inputs.update(parent['input_sha256'])
    inputs.update({PARENT+'.py': parent['source_sha256'], PARENT+'_state.npz': parent['state_sha256'],
                   PARENT+'_results.json': sha(project_path(PARENT+'_results.json'))})
    for filename, digest in inputs.items():
        if sha(project_path(filename)) != digest:
            raise ArithmeticError('Frozen input mismatch: '+filename)
    with np.load(project_path(PARENT+'_state.npz')) as saved:
        raw = {k: saved['full_'+k].copy() for k in old}
        h = saved['source_coefficients'].copy()
    if not all(np.array_equal(raw[k][..., :356, :356], old[k]) for k in old):
        raise ArithmeticError('Frozen356 core changed.')
    if not np.array_equal(h[:356], source) or np.any(h[356:]):
        raise ArithmeticError('The retained source changed.')
    insertion = sum(KAPPA**j*v for j, v in enumerate(raw['derivative']))[np.ix_(indices, indices)]
    return fixed_forms(raw), fixed_forms(exterior), indices, insertion, h, inputs, parent


def partial_fixed(data, old_count, batch=256):
    n = data['value'].shape[-1]
    if not 0 < old_count < n:
        raise ValueError('Require old and new field columns.')
    result = {k: np.zeros((n, n), complex) for k in KEYS}
    for lo in range(0, len(data['P']), batch):
        value, radial, angular = [data[k][lo:lo+batch] for k in ('value', 'radial', 'angular')]
        p, w, weight = [data[k][lo:lo+batch, None, None] for k in ('P', 'W', 'weight')]
        z = (np.sqrt(weight*p)*value).reshape(-1, n)
        a = (np.sqrt(weight/p)*(radial+angular+KAPPA*1j*w*value))[:, [2, 3, 0, 1], :].reshape(-1, n)
        for key, left, right in (('gram', z, z), ('action_gram', a, a), ('weak_action', z, a)):
            ll, ln = left[:, :old_count], left[:, old_count:]
            rl, rn = right[:, :old_count], right[:, old_count:]
            result[key][:old_count, old_count:] += ll.conj().T@rn
            result[key][old_count:, :old_count] += ln.conj().T@rl
            result[key][old_count:, old_count:] += ln.conj().T@rn
    for key in ('gram', 'action_gram'):
        result[key] = hermitian(result[key])
    return result


def assemble_increment(order):
    accumulated, control = None, None
    start = perf_counter()
    for number, data in enumerate(bridge_chunks(order=order)):
        if control is None:
            ids = np.array([0, 48, 312, 595, 596, 662, 728, 859])
            sample = {k: (v[:31, :, ids] if k in ('value', 'radial', 'angular') else v[:31])
                      for k, v in data.items()}
            direct = fixed_forms(raw_polynomial(sample))
            partial = partial_fixed(sample, 4)
            for value in direct.values():
                value[:4, :4] = 0
            control = {k: relative_error(direct[k], partial[k]) for k in KEYS}
        part = partial_fixed(data, LOCAL_OLD)
        if accumulated is None:
            accumulated = part
        else:
            for k in KEYS:
                accumulated[k] += part[k]
        if (number+1) % 16 == 0:
            print(json.dumps(dict(stage='fixed_forms', order=order, chunks=number+1,
                                  seconds=perf_counter()-start)), flush=True)
    ids = local_to_global()
    result = {k: np.zeros((TOTAL, TOTAL), complex) for k in KEYS}
    for k in KEYS:
        result[k][np.ix_(ids, ids)] = accumulated[k]
    return result, control


def add_core(base, increment):
    out = {k: v.copy() for k, v in increment.items()}
    for k in KEYS:
        if np.any(out[k][:OLD, :OLD]):
            raise ArithmeticError('Increment changes frozen core.')
        out[k][:OLD, :OLD] = base[k]
    return out


def compress(forms, transform):
    out = {}
    for k, a in forms.items():
        out[k] = np.block([[a[:OLD, :OLD], a[:OLD, OLD:]@transform],
            [transform.conj().T@a[OLD:, :OLD], transform.conj().T@a[OLD:, OLD:]@transform]])
    return out


def adapter(forms, indices, insertion, kappa=KAPPA):
    """Ephemeral compatibility view for existing single-frequency solvers only.

    One coefficient holds the already evaluated action; weak coefficient1
    is zero. The only supplied insertion entries are V x V, the sole entries
    those solvers read. These arrays MUST NOT be used as frequency polynomials
    or as a full insertion; saved archives instead use the explicit fixed keys.
    """
    if kappa != KAPPA or insertion.shape != (len(indices), len(indices)):
        raise ValueError('This cache represents only kappa=.7 and its retained insertion.')
    d = np.zeros_like(forms['gram'])
    d[np.ix_(indices, indices)] = insertion
    return dict(norm=forms['gram'], q=forms['action_gram'][None], derivative=d[None],
                weak=np.stack((forms['weak_action'], np.zeros_like(d))))


def case(forms, exterior, indices, insertion, label):
    raw = adapter(forms, indices, insertion)
    ext = {k: np.zeros_like(v) for k, v in forms.items()}
    for k in KEYS:
        ext[k][:356, :356] = exterior[k]
    ext = adapter(ext, indices, np.zeros_like(insertion))
    regions = dict(exterior=ext, bridge={k: raw[k]-ext[k] for k in raw})
    return solve_case(raw, regions, indices, np.arange(len(raw['norm'])), label)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10.8, 4.3), layout='constrained')
    rows = result['records']+[result['full_pool']]
    values = [r['first_regions']['regions']['full']['residual_norm_squared'] for r in rows]
    ax[0].bar(range(len(rows)), values, color='#278b96')
    ax[0].set(xticks=range(len(rows)), xticklabels=[r['label'] for r in rows],
              ylabel='Squared auxiliary-solver error', title='Same geometry and retained source')
    ax[0].tick_params(axis='x', labelrotation=25)
    groups = result['group_gains']
    ax[1].bar(['Localized collars', 'More global sines', 'Both together'],
              [groups['collar'], groups['global'], result['selection']['total_gain']], color='#ce873b')
    ax[1].set(ylabel='Reduction available in each tested pool', title='Which new field shapes help?')
    ax[1].tick_params(axis='x', labelrotation=15)
    fig.suptitle('Testing corrections concentrated near the bridge openings')
    fig.supxlabel('Finite calculation at one frequency; no supporting stress, physical medium or backward travel established.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    base, exterior, indices, insertion, source, inputs, parent = load_problem()
    provider = json.loads(project_path('dihole_collar_candidates_checks.json').read_text(encoding='utf-8'))
    if not provider['passed'] or provider['source_sha256'] != sha(project_path('dihole_collar_candidates.py')):
        raise ArithmeticError('Candidate checks or provenance failed.')
    inputs.update(provider['input_sha256'])
    for n in ('dihole_collar_candidates.py', 'dihole_collar_candidates_checks.json'):
        inputs[n] = sha(project_path(n))
    result = dict(complete=False, passed=False, source_sha256=sha(__file__), input_sha256=inputs,
        parameters=dict(kappa=KAPPA, mass=.1, source_dimension=216, old_dimension=OLD,
            candidate_dimension=264, bridge_orders=[64, 80], selected_counts=[8, 16, 32, 64]),
        matrix_format='Fixed-frequency gram, action_gram, weak_action; retained insertion separately. Not polynomial coefficients.',
        candidates=column_metadata(), records=[], scope=[
            'Old628 fields and all exterior forms stay fixed; only264 candidate blocks are integrated.',
            'SVD optimizes the first weighted solve in this pool, not the four-filter response.',
            'The matrix cache and selected load apply only at kappa=.7; no other frequency may be inferred.',
            'Quadrature comparisons refine new blocks only; source complement, sectors and physical stress remain unresolved.'])
    saved = dict(retained_indices=indices, retained_insertion=insertion)
    try:
        increment, control = assemble_increment(64)
        forms = add_core(base, increment)
        saved.update({'full_'+k: v for k, v in forms.items()})
        g, q, w = [forms[k] for k in KEYS]
        normal = hermitian(q+.1*(w+w.conj().T)+.01*g)
        h = np.zeros((TOTAL, 216), complex); h[:OLD] = source
        rhs = (w.conj().T+.1*g)@h
        t, selection, schur, reduced = select_directions(normal, rhs, OLD)
        result['selection'] = selection
        saved.update(selection_transform=t, source_coefficients=h, candidate_schur=schur, candidate_reduced_rhs=reduced)
        groups = {name: [r['index'] for r in result['candidates'] if r['group']==name] for name in ('collar', 'global')}
        result['group_gains'] = {name: float(np.vdot(reduced[ids], solve(schur[np.ix_(ids, ids)],
            reduced[ids], assume_a='pos')).real) for name, ids in groups.items()}
        for count in (0, 8, 16, 32, 64):
            row = case(compress(forms, t[:, :count]), exterior, indices, insertion,
                       'Existing628' if count==0 else '+'+str(count)+' selected')
            result['records'].append(row)
            print(json.dumps(dict(stage=row['label'], residual_squared=row['first_regions']['regions']['full']['residual_norm_squared'],
                                  seconds=perf_counter()-start)), flush=True)
        result['full_pool'] = case(forms, exterior, indices, insertion, 'Full264 pool')
        higher, higher_control = assemble_increment(80)
        saved.update({'refined_increment_'+k: v for k, v in higher.items()})
        result['refined'] = case(compress(add_core(base, higher), t[:, :64]), exterior, indices, insertion, 'Same64 refined')
        result['quadrature_changes'] = {k: relative_error(increment[k], higher[k]) for k in KEYS}
        residuals = [r['first_regions']['regions']['full']['residual_norm_squared'] for r in result['records']]
        full = result['full_pool']['first_regions']['regions']['full']['residual_norm_squared']
        result['predicted_gain_errors'] = [abs(residuals[0]-residuals[i+1]-selection['cumulative_gains'][k-1])
                                          for i, k in enumerate((8, 16, 32, 64))]+[abs(residuals[0]-full-selection['total_gain'])]
        result['baseline_reproduction_error'] = abs(residuals[0]-parent['full_pool']['first_regions']['regions']['full']['residual_norm_squared'])
        result['refined_residual_relative_change'] = abs(result['refined']['first_regions']['regions']['full']['residual_norm_squared']-residuals[-1])/max(1., residuals[-1])
        result['partial_form_controls'] = [control, higher_control]
        rows = result['records']+[result['full_pool'], result['refined']]
        cases = [c for r in rows for c in r['common_filter']]
        result['tests'] = {k: bool(v) for k, v in dict(
            exact_provenance=sha(__file__)==result['source_sha256'] and all(sha(project_path(n))==v for n, v in inputs.items()),
            old_core_fixed=all(np.array_equal(forms[k][:OLD, :OLD], base[k]) for k in KEYS),
            baseline_reproduced=result['baseline_reproduction_error']<2e-8,
            partial_forms=max(v for c in (control, higher_control) for v in c.values())<2e-12,
            action_orthonormality=selection['action_metric_error']<2e-8,
            predicted_gains=max(result['predicted_gain_errors'])<2e-6,
            residual_nonincrease=all(b<=a+1e-8 for a,b in zip(residuals, residuals[1:]+[full])),
            positive_grams=all(r['gram_minimum']>0 for r in rows),
            global_weak_skew=max(c['normalized_weak_skew_relative_error'] for c in cases)<2e-8,
            source_overlap=max(c['overlap_isometry_error'] for c in cases)<2e-8,
            residual_identity=max(s['gram_identity_relative_error'] for c in cases for s in c['stages'])<2e-8,
            candidate_quadrature=max(result['quadrature_changes'].values())<2e-6,
            selected_quadrature=result['refined_residual_relative_change']<2e-6).items()}
        result.update(complete=True, passed=all(result['tests'].values()))
    except Exception as exc:
        result['failure'] = dict(type=type(exc).__name__, message=str(exc))
        raise
    finally:
        result['seconds'] = perf_counter()-start
        path = project_path(NAME+'_state.npz'); np.savez_compressed(path, **saved)
        result['state_sha256'] = sha(path)
        project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], groups=result['group_gains'],
        residuals=residuals, full_pool=full, seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Collar-selection checks failed; evidence saved.')
    plot(result)


if __name__ == '__main__':
    run()
