"""Targeted global mouth-trace enrichment of the saved sewn Dirac space.

Only new-new and old-new differential-action Gram blocks are assembled.
The old344 fields and retained216-column source remain fixed. First-order
residuals and common-PV intervals are diagnostics, not physical stress.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh

from dihole_bridge_enrichment import raw_polynomial, select_columns
from dihole_common_pv_residual import chain
from dihole_euclidean_form import relative_error
from dihole_retained_resolvent import normalized
from project_paths import project_path


NAME = 'dihole_mouth_enrichment'
KEYS = ('norm', 'q', 'derivative', 'weak')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def partial_polynomial(data, old_count, batch=256):
    """The existing form algebra, evaluated only where a new column participates."""
    n = data['value'].shape[-1]
    if not 0 < old_count < n:
        raise ValueError('Both old and new field columns are required.')
    result = dict(norm=np.zeros((n, n), complex), q=np.zeros((3, n, n), complex),
                  derivative=np.zeros((3, n, n), complex), weak=np.zeros((2, n, n), complex))
    sigma = np.array([1., -1., -1., 1.])[None, :, None]

    def gram(left, right):
        value = np.zeros((n, n), complex)
        value[:old_count, old_count:] = left[:, :old_count].conj().T @ right[:, old_count:]
        value[old_count:, :old_count] = left[:, old_count:].conj().T @ right[:, :old_count]
        value[old_count:, old_count:] = left[:, old_count:].conj().T @ right[:, old_count:]
        return value

    for lo in range(0, len(data['P']), batch):
        part = {k: v[lo:lo+batch] for k, v in data.items()}
        value, radial, angular = [part[k] for k in ('value', 'radial', 'angular')]
        p, w, f, fp = [part[k][:, None, None] for k in ('P', 'W', 'f', 'fp')]
        weight = np.sqrt(part['weight'])[:, None, None]
        x = weight*(radial+angular)/np.sqrt(p)
        y = weight*1j*w*value/np.sqrt(p)
        z = weight*np.sqrt(p)*value
        beta_x, beta_y = [a[:, [2, 3, 0, 1], :].reshape(-1, n) for a in (x, y)]
        x, y, z = [a.reshape(-1, n) for a in (x, y, z)]
        result['norm'] += gram(z, z)
        result['q'] += np.array([gram(x, x), gram(x, y)+gram(y, x), gram(y, y)])
        result['weak'] += np.array([gram(z, beta_x), gram(z, beta_y)])
        if np.any(f) or np.any(fp):
            dx = (weight*(f*radial-.5j*fp*sigma*value)/np.sqrt(p)).reshape(-1, n)
            dy = (-weight*1j*f*w*value/np.sqrt(p)).reshape(-1, n)
            result['derivative'] += np.array([gram(dx, x)+gram(x, dx),
                gram(dx, y)+gram(x, dy)+gram(dy, x)+gram(y, dx), gram(dy, y)+gram(y, dy)])
    for key in ('q', 'derivative'):
        result[key] = (result[key]+result[key].conj().transpose(0, 2, 1))/2
    result['norm'] = (result['norm']+result['norm'].conj().T)/2
    return result


def embed(partial, old_count, total_old=344):
    """Append new columns without touching any old-old entry."""
    added = len(partial['norm'])-old_count
    result = {}
    for key, value in partial.items():
        out = np.zeros(value.shape[:-2]+(total_old+added, total_old+added), complex)
        out[..., :old_count, total_old:] = value[..., :old_count, old_count:]
        out[..., total_old:, :old_count] = value[..., old_count:, :old_count]
        out[..., total_old:, total_old:] = value[..., old_count:, old_count:]
        result[key] = out
    return result


def augment(old, addition):
    result = {key: value.copy() for key, value in addition.items()}
    for key in result:
        if np.any(result[key][..., :344, :344]):
            raise ArithmeticError('Increment must not modify an old-old block.')
        result[key][..., :344, :344] = old[key]
    return result


def assemble_increment(bridge_order, exterior_order):
    from dihole_mouth_fields import bridge_chunks, exterior_data
    bridge, control = None, None
    for chunk in bridge_chunks(phase=0., order=bridge_order):
        if control is None:
            columns = np.r_[np.arange(6), np.arange(312, 324)]
            sample = {k: (v[:17, :, columns] if k in ('value', 'radial', 'angular') else v[:17])
                      for k, v in chunk.items()}
            direct = raw_polynomial(sample)
            for value in direct.values():
                value[..., :6, :6] = 0
            incremental = partial_polynomial(sample, 6)
            control = {key: relative_error(direct[key], incremental[key]) for key in KEYS}
        value = partial_polynomial(chunk, 312)
        if bridge is None:
            bridge = value
        else:
            for key in KEYS:
                bridge[key] += value[key]
    data, _ = exterior_data(phase=0., order=exterior_order)
    exterior = partial_polynomial(data, 48)
    if np.any(exterior['derivative']):
        raise ArithmeticError('Exterior extension overlaps the source variation.')
    return embed(bridge, 312), embed(exterior, 48), control


def solve_case(raw, regions, source_indices, trial_indices, label):
    from dihole_residual_regions import first_residual_regions
    _, my, derivative, _, overlap = normalized(raw, source_indices, trial_indices, .7)
    inverse = np.linalg.inv(np.linalg.cholesky(raw['norm'][np.ix_(trial_indices, trial_indices)]))
    weak = raw['weak'][0]+.7*raw['weak'][1]
    b = inverse @ weak[np.ix_(trial_indices, trial_indices)] @ inverse.conj().T
    cases = [chain(my, b, derivative, overlap, eta)[0] for eta in (.01, .0001)]
    regional = first_residual_regions(raw, regions, source_indices,
        kappa=.7, mass=.1, trial_indices=trial_indices)
    return dict(label=label, dimension=len(trial_indices), first_regions=regional,
        gram_minimum=float(eigh(raw['norm'][np.ix_(trial_indices,trial_indices)],eigvals_only=True)[0]),
        common_filter=cases)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), layout='constrained')
    rows = result['records']
    labels = ['Existing fields', '+4 opening fields', '+12 opening fields']
    residual = [r['common_filter'][0]['stages'][0]['total_residual_squared'] for r in rows]
    axes[0].bar(labels, residual, color=['#557786','#248c94','#d18135'])
    axes[0].set(ylabel='First field-equation residual squared', title='Same retained source in every case')
    for i, value in enumerate(residual):
        axes[0].text(i, value, f'{value:.3g}', ha='center', va='bottom')
    for i, row in enumerate(rows):
        lo, hi = row['common_filter'][0]['conditional_retained_source_interval']
        axes[1].plot([i, i], [lo, hi], color=['#557786','#248c94','#d18135'][i], lw=3)
        axes[1].plot(i, row['common_filter'][0]['trial_filter_response'], 'ko', ms=4)
    axes[1].axhline(0, color='gray', lw=1)
    axes[1].set(xticks=range(3), xticklabels=labels, ylabel='Conditional retained-source interval',
                title='Positive regulator shift eta = 0.01')
    fig.suptitle('Testing fields that cross the bridge openings', fontsize=13)
    fig.supxlabel('Finite trial fields on a prescribed geometry; no physical supporting stress or time passage is established.', fontsize=8)
    for ext in ('png','pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext),dpi=160)
    plt.close(fig)


def run():
    from dihole_mouth_fields import column_metadata
    from dihole_residual_regions import load_baseline
    start = perf_counter()
    base, metadata, exterior_base, inputs = load_baseline()
    field_checks = json.loads(project_path('dihole_mouth_fields_checks.json').read_text(encoding='utf-8'))
    if not field_checks['passed'] or field_checks['source_sha256'] != sha(project_path('dihole_mouth_fields.py')):
        raise ArithmeticError('Mouth-field controls or source provenance failed.')
    inputs = dict(inputs)
    for name in ('dihole_mouth_fields.py','dihole_mouth_fields_checks.json',
                 'dihole_residual_regions.py','dihole_residual_regions_results.json',
                 'dihole_residual_regions_state.npz','dihole_common_pv_residual.py',
                 'dihole_retained_resolvent.py','dihole_bridge_enrichment.py',
                 'dihole_exterior_feedback.py','dihole_euclidean_form.py','project_paths.py'):
        inputs[name] = sha(project_path(name))
    source_indices = select_columns(metadata, (3,4), 3, 12)
    result = dict(complete=False, passed=False, source_sha256=sha(__file__),input_sha256=inputs,
        new_columns=column_metadata(),parameters=dict(phase=0.,mu=.5,kappa=.7,source_dimension=216,
            old_dimension=344,new_fields=[4,12],bridge_orders=[64,80],exterior_orders=[24,32]),
        records=[],scope=['The retained216-field source is fixed while only the trial solution space is enlarged.',
            'New mouth fields are global conforming spinors; they are not new material or a changed geometry.',
            'Only new-new and old-new forms are assembled; the old344-field core is preserved exactly.',
            'Residual and PV intervals use floating integrals and concern only the retained insertion at positive eta.',
            'A residual reduction does not supply the omitted physical source, continuum error or stable support.'])
    saved = {}
    try:
        bridge, exterior, control = assemble_increment(64,24)
        addition = {key: bridge[key]+exterior[key] for key in KEYS}
        raw = augment(base,addition)
        regions = dict(exterior=augment(exterior_base,exterior),
            bridge=augment({key:base[key]-exterior_base[key] for key in KEYS},bridge))
        for prefix, data in (('full',raw),('bridge_increment',bridge),('exterior_increment',exterior)):
            for key,value in data.items(): saved[prefix+'_'+key]=value
        for label,indices in (('baseline',np.arange(344)),('zero_modes',np.arange(348)),('first_pair',np.arange(356))):
            row=solve_case(raw,regions,source_indices,indices,label)
            result['records'].append(row)
            print(json.dumps(dict(stage=label,residual_squared=row['common_filter'][0]['stages'][0]['total_residual_squared'],
                interval=row['common_filter'][0]['conditional_retained_source_interval'],seconds=perf_counter()-start)),flush=True)
        higher_bridge,higher_exterior,higher_control=assemble_increment(80,32)
        # Compare new-field integration only; the same fixed old core/source is retained.
        refined=augment(base,{key:higher_bridge[key]+higher_exterior[key] for key in KEYS})
        refined_regions=dict(exterior=augment(exterior_base,higher_exterior),
            bridge=augment({key:base[key]-exterior_base[key] for key in KEYS},higher_bridge))
        result['refined']=solve_case(refined,refined_regions,source_indices,np.arange(356),'higher_new_field_quadrature')
        for key,value in refined.items(): saved['refined_'+key]=value
        result['partial_form_control']=control
        result['refined_partial_form_control']=higher_control
        result['new_block_quadrature_change']={key:relative_error(
            addition[key],higher_bridge[key]+higher_exterior[key]) for key in KEYS}
        rows=result['records']
        residual=[r['common_filter'][0]['stages'][0]['total_residual_squared'] for r in rows]
        refined_residual=result['refined']['common_filter'][0]['stages'][0]['total_residual_squared']
        result['residual_reduction_fraction']=[1-x/residual[0] for x in residual]
        result['residual_quadrature_relative_change']=abs(refined_residual-residual[-1])/max(1.,residual[-1])
        all_cases=[c for r in rows+[result['refined']] for c in r['common_filter']]
        result['tests']={key:bool(value) for key,value in dict(
            input_provenance=sha(__file__)==result['source_sha256'] and all(sha(project_path(n))==h for n,h in inputs.items()),
            fixed_old_core=all(np.array_equal(raw[key][...,:344,:344],base[key]) for key in KEYS),
            partial_form_agreement=max(list(control.values())+list(higher_control.values()))<2e-12,
            positive_grams=all(r['gram_minimum']>0 for r in rows+[result['refined']]),
            measured_global_skew=max(c['normalized_weak_skew_relative_error'] for c in all_cases)<2e-8,
            first_residual_monotonic=all(residual[i+1]<=residual[i]+1e-8 for i in (0,1)),
            source_overlap=max(c['overlap_isometry_error'] for c in all_cases)<2e-9,
            residual_identity=max(s['gram_identity_relative_error'] for c in all_cases for s in c['stages'])<2e-9,
            new_quadrature=max(result['new_block_quadrature_change'].values())<2e-6,
            residual_quadrature=result['residual_quadrature_relative_change']<2e-6).items()}
        result.update(complete=True,passed=all(result['tests'].values()))
    except Exception as exc:
        result['failure']=dict(type=type(exc).__name__,message=str(exc))
        raise
    finally:
        result['seconds']=perf_counter()-start
        if saved:
            path=project_path(NAME+'_state.npz');np.savez_compressed(path,**saved);result['state_sha256']=sha(path)
        project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=result['tests'],reduction=result['residual_reduction_fraction'],
        quadrature=result['residual_quadrature_relative_change'],seconds=result['seconds']),indent=2),flush=True)
    if not result['passed']:
        raise ArithmeticError('Mouth-field integration checks failed; evidence saved.')
    plot(result)


if __name__=='__main__':
    run()
