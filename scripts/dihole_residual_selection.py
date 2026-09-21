"""Choose additional bridge fields by their actual first-residual reduction.

The source216 and frozen157 core356 stay fixed. A finite pool of272 new
radial/angular fields is scored by an action-Gram Schur complement. SVD
selects the best k combinations for this one weighted least-squares step.
This is an auxiliary linear solve, not an Einstein equation or new matter.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh, solve, solve_triangular, svd

from dihole_bridge_enrichment import raw_polynomial, select_columns
from dihole_common_pv_residual import hermitian
from dihole_euclidean_form import relative_error
from dihole_mouth_enrichment import partial_polynomial, augment, solve_case
from dihole_residual_regions import first_residual_regions, load_baseline
from project_paths import project_path


NAME = 'dihole_residual_selection'
KEYS = ('norm', 'q', 'derivative', 'weak')
OLD = 356


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_problem():
    old, metadata, old_exterior, inputs = load_baseline()
    name = 'dihole_mouth_enrichment'
    record = json.loads(project_path(name+'_results.json').read_text(encoding='utf-8'))
    if not record['passed']:
        raise ArithmeticError('The opening-field baseline did not pass.')
    checks = dict(record['input_sha256'])
    checks.update({name+'.py': record['source_sha256'], name+'_state.npz': record['state_sha256']})
    for path, digest in checks.items():
        if sha(project_path(path)) != digest:
            raise ArithmeticError('Frozen157 input mismatch: '+path)
    inputs.update(checks)
    inputs[name+'_results.json'] = sha(project_path(name+'_results.json'))
    with np.load(project_path(name+'_state.npz')) as archive:
        base = {k: archive['full_'+k].copy() for k in KEYS}
        exterior = augment(old_exterior, {k: archive['exterior_increment_'+k] for k in KEYS})
    if not all(np.array_equal(base[k][..., :344, :344], old[k]) for k in KEYS):
        raise ArithmeticError('Frozen344 core changed in the356 baseline.')
    regions = dict(exterior=exterior, bridge={k: base[k]-exterior[k] for k in KEYS})
    source_indices = select_columns(metadata, (3, 4), 3, 12)
    first, columns = first_residual_regions(base, regions, source_indices, return_columns=True)
    expected = record['records'][-1]['first_regions']['regions']['full']['residual_norm_squared']
    if abs(first['regions']['full']['residual_norm_squared']-expected) > 1e-8:
        raise ArithmeticError('Frozen157 residual did not reproduce.')
    return base, exterior, source_indices, columns['source'], inputs


def assemble_increment(order):
    from dihole_residual_candidates import bridge_chunks
    accumulated, control = None, None
    for data in bridge_chunks(order=order):
        if control is None:
            ids = np.r_[[0, 1, 48, 49, 312, 313], [324, 325, 500, 501]]
            sample = {k: (v[:31, :, ids] if k in ('value', 'radial', 'angular') else v[:31])
                      for k, v in data.items()}
            direct = raw_polynomial(sample)
            for value in direct.values():
                value[..., :6, :6] = 0
            partial = partial_polynomial(sample, 6)
            control = {k: relative_error(direct[k], partial[k]) for k in KEYS}
        part = partial_polynomial(data, 324)
        if accumulated is None:
            accumulated = part
        else:
            for k in KEYS:
                accumulated[k] += part[k]
    # Exterior-only old columns312:344 vanish everywhere on this bridge.
    ids = np.r_[np.arange(312), np.arange(344, 356), np.arange(356, 628)]
    result = {}
    for key, value in accumulated.items():
        out = np.zeros(value.shape[:-2]+(628, 628), complex)
        out[..., ids[:, None], ids] = value
        result[key] = out
    return result, control


def add_to_core(base, increment):
    out = {k: a.copy() for k, a in increment.items()}
    for key in KEYS:
        if np.any(out[key][..., :OLD, :OLD]):
            raise ArithmeticError('Candidate increment modifies the old core.')
        out[key][..., :OLD, :OLD] = base[key]
    return out


def select_directions(normal, rhs, old_count):
    """Best k directions in this finite candidate pool; no small modes dropped."""
    m, cross = normal[:old_count, :old_count], normal[:old_count, old_count:]
    schur = hermitian(normal[old_count:, old_count:]-cross.conj().T@solve(m, cross, assume_a='pos'))
    reduced_rhs = rhs[old_count:]-cross.conj().T@solve(m, rhs[:old_count], assume_a='pos')
    diagonal = schur.diagonal().real
    if np.min(diagonal) <= 0:
        raise ArithmeticError('Candidate action Schur diagonal is not positive.')
    scale = 1/np.sqrt(diagonal)
    balanced = hermitian(scale[:, None]*schur*scale[None, :])
    spectrum = eigh(balanced, eigvals_only=True)
    if spectrum[0] <= 1e-11*spectrum[-1]:
        raise ArithmeticError('Candidate Schur rank is numerically unresolved; no modes silently discarded.')
    lower = np.linalg.cholesky(balanced)
    whitened = solve_triangular(lower, scale[:, None]*reduced_rhs, lower=True)
    vectors, singular, _ = svd(whitened, full_matrices=False)
    transform = scale[:, None]*solve_triangular(lower.conj().T, vectors, lower=False)
    metric_error = float(np.linalg.norm(transform.conj().T@schur@transform-np.eye(len(singular)), 2))
    return transform, dict(gains=(singular**2).tolist(), cumulative_gains=np.cumsum(singular**2).tolist(),
        total_gain=float(np.linalg.norm(whitened, 'fro')**2), balanced_schur_minimum=float(spectrum[0]),
        balanced_schur_maximum=float(spectrum[-1]), action_metric_error=metric_error,
        discarded_schur_modes=0), schur, reduced_rhs


def compress(raw, transform):
    """Append candidate combinations; old356 entries remain bitwise intact."""
    count = transform.shape[1]
    result = {}
    for key, value in raw.items():
        out = np.zeros(value.shape[:-2]+(OLD+count, OLD+count), complex)
        out[..., :OLD, :OLD] = value[..., :OLD, :OLD]
        out[..., :OLD, OLD:] = value[..., :OLD, OLD:]@transform
        out[..., OLD:, :OLD] = transform.conj().T@value[..., OLD:, :OLD]
        out[..., OLD:, OLD:] = transform.conj().T@value[..., OLD:, OLD:]@transform
        result[key] = out
    return result


def regions_for(raw, old_exterior):
    exterior = {k: np.zeros_like(v) for k, v in raw.items()}
    for k in KEYS:
        exterior[k][..., :OLD, :OLD] = old_exterior[k]
    return dict(exterior=exterior, bridge={k: raw[k]-exterior[k] for k in KEYS})


def selection_control():
    """Independent direct least-squares test of the optimal rank-k gain."""
    rng = np.random.default_rng(158)
    f = rng.normal(size=(30, 14))+1j*rng.normal(size=(30, 14))
    h = rng.normal(size=(30, 5))+1j*rng.normal(size=(30, 5))
    t, info, _, _ = select_directions(f.conj().T@f, f.conj().T@h, 6)
    before = np.linalg.norm(h-f[:, :6]@np.linalg.lstsq(f[:, :6], h, rcond=None)[0], 'fro')**2
    errors = []
    for k in (1, 2, 4, 5):
        candidate = np.column_stack((f[:, :6], f[:, 6:]@t[:, :k]))
        residual = h-candidate@np.linalg.lstsq(candidate, h, rcond=None)[0]
        errors.append(abs(before-np.linalg.norm(residual, 'fro')**2-info['cumulative_gains'][k-1]))
    return dict(maximum_direct_gain_error=float(max(errors)), action_metric_error=info['action_metric_error'])


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), layout='constrained')
    rows = result['records']
    labels = [r['label'] for r in rows]
    values = [r['first_regions']['regions']['full']['residual_norm_squared'] for r in rows]
    axes[0].bar(range(len(rows)), values, color='#288990')
    axes[0].set(xticks=range(len(rows)), xticklabels=labels, ylabel='Squared residual of the auxiliary solve',
                title='Same source and existing opening fields')
    axes[0].tick_params(axis='x', labelrotation=30)
    gains = np.array(result['selection']['cumulative_gains'])
    axes[1].plot(np.arange(1, len(gains)+1), 100*gains/result['selection']['total_gain'], color='#c67b30')
    axes[1].set(xlabel='Selected field combinations', ylabel='Available finite-pool improvement captured (%)',
                title='Optimal for this first solve only', ylim=(0, 102))
    fig.suptitle('Letting the remaining residual choose the added fields')
    fig.supxlabel('Numerical approximation of a prescribed geometry; no physical support or time passage established.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    from dihole_residual_candidates import column_metadata
    start = perf_counter()
    base, exterior, source_indices, source, inputs = load_problem()
    for name in ('dihole_residual_candidates.py', 'dihole_residual_candidates_checks.json',
                 'dihole_residual_regions.py', 'dihole_common_pv_residual.py', 'dihole_retained_resolvent.py'):
        inputs[name] = sha(project_path(name))
    checks = json.loads(project_path('dihole_residual_candidates_checks.json').read_text(encoding='utf-8'))
    if not checks['passed'] or checks['source_sha256'] != inputs['dihole_residual_candidates.py']:
        raise ArithmeticError('Candidate helper checks/provenance failed.')
    result = dict(passed=False, complete=False, source_sha256=sha(__file__), input_sha256=inputs,
        parameters=dict(kappa=.7, mass=.1, phase=0., mu=.5, source_dimension=216, old_dimension=OLD,
                        candidate_dimension=272, bridge_orders=[64, 80], selected_counts=[8, 16, 32, 64]),
        candidates=column_metadata(), records=[], scope=[
            'The old356 fields, exterior, and retained216-column source are fixed.',
            'SVD maximizes total first-step residual reduction only within this finite pool.',
            'Source eigen-sign labels are unchanged; the four-filter source is not the optimized objective.',
            'Unused singular-value gain is not a continuum or physical-source error bound.',
            'Floating geometry integration is not enclosed; no new law, material or supported passage.'])
    saved = {}
    try:
        result['independent_control'] = selection_control()
        increment, control = assemble_increment(64)
        raw = add_to_core(base, increment)
        saved.update({'full_'+k: v for k, v in raw.items()})
        kappa, mass = .7, .1
        weak = raw['weak'][0]+kappa*raw['weak'][1]
        normal = hermitian(sum(kappa**j*v for j, v in enumerate(raw['q']))
                           +mass*(weak+weak.conj().T)+mass**2*raw['norm'])
        h = np.zeros((len(normal), source.shape[1]), complex); h[:OLD] = source
        rhs = (weak.conj().T+mass*raw['norm'])@h
        transform, selection, schur, reduced_rhs = select_directions(normal, rhs, OLD)
        saved.update(selection_transform=transform, source_coefficients=h,
                     candidate_schur=schur, candidate_reduced_rhs=reduced_rhs)
        result['selection'] = selection
        groups = dict(radial=[r['index'] for r in result['candidates'] if r['group']=='radial'],
                      angular=[r['index'] for r in result['candidates'] if r['group']=='angular'],
                      nested92=[r['index'] for r in result['candidates']
                          if (r['group']=='radial' and r['radial_mode']<=14) or (r['group']=='angular' and r['l']==6)])
        result['group_gains'] = {name: float(np.vdot(reduced_rhs[ids],
            solve(schur[np.ix_(ids, ids)], reduced_rhs[ids], assume_a='pos')).real)
            for name, ids in groups.items()}
        spaces = [('baseline', np.zeros((272, 0), complex))]+[(str(k)+' selected', transform[:, :k]) for k in (8, 16, 32, 64)]
        for label, t in spaces:
            small = compress(raw, t)
            row = solve_case(small, regions_for(small, exterior), source_indices, np.arange(len(small['norm'])), label)
            result['records'].append(row)
            print(json.dumps(dict(stage=label, residual_squared=row['first_regions']['regions']['full']['residual_norm_squared'],
                seconds=perf_counter()-start)), flush=True)
        result['full_pool'] = solve_case(raw, regions_for(raw, exterior), source_indices, np.arange(628), 'full272 pool')
        higher, higher_control = assemble_increment(80)
        saved.update({'refined_increment_'+k: v for k, v in higher.items()})
        refined_raw = add_to_core(base, higher)
        selected_refined = compress(refined_raw, transform[:, :64])
        result['refined'] = solve_case(selected_refined, regions_for(selected_refined, exterior), source_indices,
                                      np.arange(420), 'same64 refined')
        result['quadrature_changes'] = {k: relative_error(increment[k], higher[k]) for k in KEYS}
        baseline = result['records'][0]['first_regions']['regions']['full']['residual_norm_squared']
        residuals = [r['first_regions']['regions']['full']['residual_norm_squared'] for r in result['records']]
        full_residual = result['full_pool']['first_regions']['regions']['full']['residual_norm_squared']
        result['predicted_gain_errors'] = [abs(baseline-residuals[i+1]-selection['cumulative_gains'][k-1])
                                          for i, k in enumerate((8, 16, 32, 64))]
        result['full_gain_error'] = abs(baseline-full_residual-selection['total_gain'])
        result['refined_residual_relative_change'] = abs(result['refined']['first_regions']['regions']['full']['residual_norm_squared']
                                                        -residuals[-1])/max(1., residuals[-1])
        result['partial_form_controls'] = [control, higher_control]
        rows = result['records']+[result['full_pool'], result['refined']]
        cases = [c for row in rows for c in row['common_filter']]
        result['tests'] = {k: bool(v) for k,v in dict(
            exact_provenance=sha(__file__)==result['source_sha256'] and all(sha(project_path(n))==v for n,v in inputs.items()),
            old_core_fixed=all(np.array_equal(raw[k][..., :OLD, :OLD], base[k]) for k in KEYS),
            partial_form_agreement=max(v for c in (control, higher_control) for v in c.values())<2e-12,
            independent_selection=result['independent_control']['maximum_direct_gain_error']<2e-10,
            action_orthonormality=selection['action_metric_error']<2e-9,
            predicted_gain=max(result['predicted_gain_errors']+[result['full_gain_error']])<2e-7,
            residual_nonincrease=all(b<=a+1e-8 for a,b in zip(residuals, residuals[1:]+[full_residual])),
            positive_grams=all(r['gram_minimum']>0 for r in rows),
            global_weak_skew=max(c['normalized_weak_skew_relative_error'] for c in cases)<2e-8,
            source_overlap=max(c['overlap_isometry_error'] for c in cases)<2e-9,
            residual_identity=max(s['gram_identity_relative_error'] for c in cases for s in c['stages'])<2e-9,
            candidate_quadrature=max(result['quadrature_changes'].values())<2e-6,
            selected_residual_quadrature=result['refined_residual_relative_change']<2e-6).items()}
        result.update(complete=True, passed=all(result['tests'].values()))
    except Exception as exc:
        result['failure'] = dict(type=type(exc).__name__, message=str(exc))
        raise
    finally:
        result['seconds'] = perf_counter()-start
        if saved:
            path = project_path(NAME+'_state.npz'); np.savez_compressed(path, **saved); result['state_sha256'] = sha(path)
        project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], groups=result['group_gains'],
                          residuals=residuals, full_pool=full_residual, seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Residual-selection checks failed; evidence saved.')
    plot(result)


if __name__ == '__main__':
    run()
