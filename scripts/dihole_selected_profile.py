"""Where the selected64-field first residual remains on the unchanged bridge.

Only saved matrices and one pointwise bridge pass are used. The source is
the same retained216-column insertion; the selected trial is lifted back
into the provider's raw628-column chart before its full action is evaluated.
These numerical squared-error densities are not energy or supporting stress.
"""
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_key, '1')

import numpy as np
from scipy.linalg import eigh, solve

from dihole_residual_candidates import bridge_chunks, local_to_global
from dihole_residual_profile import bridge_grid, summarize
from dihole_residual_regions import first_residual_regions
from dihole_residual_selection import load_problem, compress, regions_for
from project_paths import project_path


NAME = 'dihole_selected_profile'
PARENT = 'dihole_residual_selection'
INITIAL = 'dihole_residual_profile'
KEYS = ('norm', 'q', 'derivative', 'weak')
INITIAL_HASHES = {
    INITIAL+'.py': 'cc08d16b1e25288f7bd8aa377037782e1c2ab9ef4930423c82913b7e8c7c91e3',
    INITIAL+'_results.json': '6436ab832fe57584fa3486776d5795f885bd0eb60224890d69ea5a59fcd3bbcf',
    INITIAL+'_state.npz': '2b8b375e326fa2f9aa9869b8efa7b9b1cf429c61512cd123bf574e6b4ac43fd1'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(inputs):
    bad = [name for name, digest in inputs.items() if sha(project_path(name)) != digest]
    if bad:
        raise ArithmeticError('Frozen input changed: '+', '.join(bad))


def selection_audit(raw, source, transform, parent):
    """Independent raw normal solves and a216-square Gram eigenvalue check.

    This does not call the production selection routine or its SVD. The
    spectrum of R* S^-1 R is the squared-singular-value gain spectrum.
    """
    g = raw['norm']
    weak = raw['weak'][0]+.7*raw['weak'][1]
    normal = sum(.7**j*c for j, c in enumerate(raw['q']))+.1*(weak+weak.conj().T)+.01*g
    normal = (normal+normal.conj().T)/2
    rhs = (weak.conj().T+.1*g)@source
    source_square = float(np.vdot(source, g@source).real)
    old, cross = normal[:356, :356], normal[:356, 356:]
    old_solution = solve(old, rhs[:356], assume_a='pos')
    schur = normal[356:, 356:]-cross.conj().T@solve(old, cross, assume_a='pos')
    schur = (schur+schur.conj().T)/2
    reduced = rhs[356:]-cross.conj().T@old_solution
    gram = reduced.conj().T@solve(schur, reduced, assume_a='pos')
    gains = eigh((gram+gram.conj().T)/2, eigvals_only=True)[::-1]
    expected_gains = np.asarray(parent['selection']['gains'])
    rows = []
    for count, saved in zip((0, 8, 16, 32, 64), parent['records']):
        t = transform[:, :count]
        block = np.block([[old, cross@t],
                          [t.conj().T@cross.conj().T, t.conj().T@normal[356:, 356:]@t]])
        load = np.concatenate((rhs[:356], t.conj().T@rhs[356:]))
        solution = solve(block, load, assume_a='pos')
        residual = source_square-float(np.vdot(load, solution).real)
        target = saved['first_regions']['regions']['full']['residual_norm_squared']
        rows.append(dict(selected=count, residual_squared=residual, saved_difference=residual-target,
            relative_normal_residual=float(np.linalg.norm(block@solution-load)/np.linalg.norm(load))))
    full_solution = solve(normal, rhs, assume_a='pos')
    full_residual = source_square-float(np.vdot(rhs, full_solution).real)
    return dict(raw_solve_rows=rows, full_pool_residual_squared=full_residual,
        full_pool_saved_difference=full_residual-parent['full_pool']['first_regions']['regions']['full']['residual_norm_squared'],
        gain_spectrum_maximum_absolute_error=float(np.max(abs(gains-expected_gains))),
        gain_spectrum_relative_error=float(np.max(abs(gains-expected_gains))/max(1., expected_gains[0])),
        smallest_computed_gain=float(gains[-1]), total_gain=float(np.trace(gram).real),
        total_gain_saved_difference=float(np.trace(gram).real)-parent['selection']['total_gain'])


def load_selected():
    """Authenticate the successful producer, recover selected trial columns."""
    parent = json.loads(project_path(PARENT+'_results.json').read_text(encoding='utf-8'))
    if not parent['complete'] or not parent['passed'] or not all(parent['tests'].values()):
        raise ArithmeticError('The residual-selection producer is not successful.')
    inputs = {**parent['input_sha256'], **INITIAL_HASHES,
              PARENT+'.py': parent['source_sha256'], PARENT+'_state.npz': parent['state_sha256'],
              PARENT+'_results.json': sha(project_path(PARENT+'_results.json'))}
    authenticate(inputs)
    base, exterior, indices, old_source, inherited = load_problem()
    inputs.update(inherited)
    with np.load(project_path(PARENT+'_state.npz')) as saved:
        raw = {key: saved['full_'+key].copy() for key in KEYS}
        transform = saved['selection_transform'].copy()
        parent_source = saved['source_coefficients'].copy()
    if not all(np.array_equal(raw[key][..., :356, :356], base[key]) for key in KEYS):
        raise ArithmeticError('The saved candidate space changed its old356 core.')
    if not np.array_equal(parent_source[:356], old_source) or np.any(parent_source[356:]):
        raise ArithmeticError('The selection source is not the unchanged retained insertion.')
    audit = selection_audit(raw, parent_source, transform, parent)
    small = compress(raw, transform[:, :64])
    first, compact = first_residual_regions(small, regions_for(small, exterior), indices,
                                            kappa=.7, mass=.1, return_columns=True)
    columns = {}
    for name in ('source', 'trial'):
        out = np.zeros((628, 216), complex)
        out[:356] = compact[name][:356]
        out[356:] = transform[:, :64]@compact[name][356:]
        columns[name] = out
    columns['source_eigenvalues'] = compact['source_eigenvalues']
    if not np.array_equal(columns['source'], parent_source):
        raise ArithmeticError('Lifting the selected basis changed the weighted source.')
    initial = json.loads(project_path(INITIAL+'_results.json').read_text(encoding='utf-8'))
    with np.load(project_path(INITIAL+'_state.npz')) as saved:
        initial_arrays = {key: saved[key].copy() for key in
            ('s', 'radial_weights', 'cosine', 'cosine_weights', 'radial_density')}
    target = next(row for row in parent['records'] if row['label'] == '64 selected')
    return columns, dict(input_sha256=inputs, audit=audit, first_solve=first,
        source_indices=indices, transform=transform[:, :64], initial=initial,
        initial_arrays=initial_arrays, target=target,
        parameters=dict(mu=.5, phase=0., kappa=.7, first_mass=.1, eta=.01,
                        source_dimension=216, old_dimension=356, selected_dimension=64,
                        trial_dimension=420, raw_dimension=628))


def selected_profile(columns, order=64, progress=False):
    grid = bridge_grid(order)
    ns, nc = len(grid['s']), len(grid['cosine'])
    mapping = local_to_global()
    trial = np.ascontiguousarray(columns['trial'][mapping])
    combination = np.ascontiguousarray((.1*columns['trial']-columns['source'])[mapping])
    positive, negative = columns['source_eigenvalues'] > 0, columns['source_eigenvalues'] < 0
    if not np.all(positive | negative):
        raise ArithmeticError('Unexpected zero source eigenvalue.')
    density = np.empty((ns, nc, 3))
    offset, weight_error, largest, start = 0, 0., 0, perf_counter()
    for chunk_number, data in enumerate(bridge_chunks(order=order)):
        points, spin, local = data['value'].shape
        if spin != 4 or local != len(mapping) or points % nc:
            raise ArithmeticError('Candidate stream chart changed.')
        count = points//nc
        expected = (grid['radial_weights'][offset:offset+count, None]*grid['cosine_weights'][None, :]).ravel()
        weight_error = max(weight_error, float(np.max(abs(data['weight']-expected))))
        action = data['radial']+data['angular']+.7j*data['W'][:, None, None]*data['value']
        residual = (action.reshape(-1, local)@trial).reshape(points, 4, 216)
        residual = residual[:, [2, 3, 0, 1], :]/np.sqrt(data['P'])[:, None, None]
        residual += (data['value'].reshape(-1, local)@combination).reshape(points, 4, 216) \
                    *np.sqrt(data['P'])[:, None, None]
        squared = abs(residual)**2
        pos, neg = [np.sum(squared[:, :, at], axis=(1, 2)) for at in (positive, negative)]
        density[offset:offset+count] = np.stack((pos+neg, pos, neg), axis=1).reshape(count, nc, 3)
        offset += count
        largest = max(largest, sum(v.nbytes for v in data.values())+action.nbytes+residual.nbytes+squared.nbytes)
        if progress and chunk_number == 0:
            print(json.dumps(dict(stage='first_chunk', seconds=perf_counter()-start,
                                  working_chunk_bytes=largest)), flush=True)
        if progress and offset % order == 0:
            integral = np.einsum('s,c,scj->j', grid['radial_weights'][:offset], grid['cosine_weights'], density[:offset])
            print(json.dumps(dict(stage='radial_panel', slices=offset,
                cumulative_residual_squared=integral.tolist(), seconds=perf_counter()-start)), flush=True)
    if offset != ns:
        raise ArithmeticError('The whole bridge was not evaluated.')
    radial = np.einsum('c,scj->sj', grid['cosine_weights'], density)
    angular = np.einsum('s,scj->cj', grid['radial_weights'], density)
    return {**grid, 'residual_density': density, 'radial_density': radial,
            'cosine_density': angular, 'theta_density': angular*np.sin(grid['theta'])[:, None]}, \
        dict(seconds=perf_counter()-start, working_chunk_bytes=largest,
             weight_maximum_absolute_error=weight_error, radial_slices=offset)


def draw(result, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), layout='constrained')
    colors = ('#547887', '#b86d2c')
    for key, label, color in (('initial_radial_density', 'Before selection:356 fields', colors[0]),
                              ('radial_density', 'After selection:420 fields', colors[1])):
        axes[0].plot(arrays['s'], arrays[key][:, 0], label=label, color=color, lw=1.6)
    for lo, hi in ((-.1, -.05), (.05, .1)):
        axes[0].axvspan(lo, hi, color='#e5bb84', alpha=.15)
    axes[0].set(xlabel='Bridge coordinate s', ylabel='Squared residual per unit s',
                title='Full first-order auxiliary residual')
    axes[0].legend(fontsize=8)
    before = [r['residual_squared'][0] for r in result['initial_summary']['radial_regions']]
    after = [r['residual_squared'][0] for r in result['summary']['radial_regions']]
    positions = np.arange(3)
    axes[1].bar(positions-.18, before, .36, color=colors[0], label='Before')
    axes[1].bar(positions+.18, after, .36, color=colors[1], label='After')
    axes[1].set(xticks=positions, xticklabels=['Left collar', 'Center', 'Right collar'],
                ylabel='Integrated squared residual', title='Where the remaining error lies')
    axes[1].legend(fontsize=8)
    fig.suptitle('Locating the error after residual-guided field selection', fontsize=13)
    fig.supxlabel('Same retained source and prescribed geometry. Auxiliary error, not energy density or supporting matter.', fontsize=8)
    for extension in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+extension), dpi=160)
    plt.close(fig)


def run():
    start, own_hash = perf_counter(), sha(__file__)
    columns, problem = load_selected()
    print(json.dumps(dict(stage='saved_matrix_audit', **problem['audit'])), flush=True)
    profile, timing = selected_profile(columns, progress=True)
    summary = summarize(profile, 64)
    total = np.array(summary['total_residual_squared'])
    target = problem['target']['first_regions']['regions']['bridge']
    expected = np.array([target['residual_norm_squared'], target['positive_source_residual_squared'],
                         target['negative_source_residual_squared']])
    reproduction = float(np.max(abs(total-expected))/max(1., expected[0]))
    audit = problem['audit']
    tests = dict(input_provenance=sha(__file__) == own_hash,
        independent_raw_normal_solves=max(abs(row['saved_difference']) for row in audit['raw_solve_rows']) < 2e-7
            and abs(audit['full_pool_saved_difference']) < 2e-7,
        independent_gain_eigenvalues=audit['gain_spectrum_relative_error'] < 2e-10
            and abs(audit['total_gain_saved_difference']) < 2e-7,
        unchanged_comparison_grid=all(np.array_equal(profile[key], problem['initial_arrays'][key])
            for key in ('s', 'radial_weights', 'cosine', 'cosine_weights')),
        pointwise_selected_bridge=reproduction < 3e-11,
        quadrature_coordinates=timing['weight_maximum_absolute_error'] < 2e-18,
        nonnegative_finite_density=np.isfinite(profile['residual_density']).all()
            and np.min(profile['residual_density']) >= 0,
        radial_and_angular_marginals=max(np.max(abs(profile['radial_weights']@profile['radial_density']-total)),
            np.max(abs(profile['cosine_weights']@profile['cosine_density']-total)))/max(1., total[0]) < 3e-13,
        exact_radial_partition=abs(sum(r['residual_squared'][0] for r in summary['radial_regions'])-total[0])
            /max(1., total[0]) < 3e-13)
    authenticate(problem['input_sha256'])
    arrays = {key: value for key, value in profile.items() if key != 'angular_panels'}
    arrays.update(columns)
    arrays.update(initial_radial_density=problem['initial_arrays']['radial_density'],
                  source_indices=problem['source_indices'], local_to_raw=local_to_global(),
                  selected_transform=problem['transform'])
    state_path = project_path(NAME+'_state.npz')
    np.savez_compressed(state_path, **arrays)
    result = dict(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        parameters={**problem['parameters'], 'quadrature_order':64, 'measure':'ds*dcos(theta)',
                    'categories':['total', 'positive_source_directions', 'negative_source_directions']},
        summary=summary, initial_summary=problem['initial']['summary'], target_bridge=target,
        pointwise_bridge_relative_error=reproduction, saved_matrix_audit=audit,
        first_solve=problem['first_solve'], evaluation=timing,
        equations=dict(coefficient_lift='C628[:356]=C420[:356]; C628[356:]=T64*C420[356:]',
            residual='P^(-1/2)*beta*(radial+angular+i*kappa*W*value)*trial + P^(1/2)*value*(m*trial-source)',
            beta='component permutation [2,3,0,1]; raw local map excludes exterior-only312:344'),
        scope=['Same retained216-column insertion; the trial has64 selected additions to the frozen356 fields.',
               'Selection is optimal only for the first weighted solve inside the tested272-candidate pool.',
               'Both eigen-sign groups contain nonnegative squared errors; their labels are not physical stress signs.',
               'No exterior or form assembly; exactly one bridge pointwise pass at the frozen nominal order64.',
               'Sampled densities and floating integrals do not enclose continuum error or supply supporting energy.'],
        source_sha256=own_hash, input_sha256=problem['input_sha256'], state_sha256=sha(state_path),
        seconds=perf_counter()-start)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], summary=summary,
                         reproduction=reproduction, seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Selected residual profile failed; evidence saved.')
    draw(result, arrays)
    return result


if __name__ == '__main__':
    run()
