"""Pointwise bridge localization of checkpoint157's first resolvent residual.

The nominal356 trial and its retained216-column source are fixed. This
evaluates the full first-order action, not a compressed Dirac square. The
bridge measure is ds*dcos(theta); no extra radius, sine, or 2*pi is inserted.
Saved densities are numerical diagnostics, not physical stress or bounds.
"""
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_key, '1')

import numpy as np

from dihole_bridge_enrichment import select_columns
from dihole_mouth_fields import RC, bridge_chunks
from dihole_radial_stream import RADIAL_PANELS
from dihole_residual_regions import first_residual_regions, load_baseline
from dihole_sewn_form import RHO_BOUNDS
from dihole_stretched_collar import gauss
from project_paths import project_path


NAME = 'dihole_residual_profile'
PARENT = 'dihole_mouth_enrichment'
KEYS = ('norm', 'q', 'derivative', 'weak')
EXPECTED = {
    PARENT+'.py': 'f3a4d0d6062e1793717e57bb8f642ce70295428fadb60b763c53faee7ebbfb42',
    PARENT+'_results.json': '490a62e74e47edd9627f566cf7d00e72f5ecfade5ca32e7f4c07c9ecdc0fb423',
    PARENT+'_state.npz': '9c7774e9ac746a520dcaa7fc82ba16cac4f51f815167d320245afccda14d0767'}
LOCAL_TO_RAW = np.r_[np.arange(312), np.arange(344, 356)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(inputs):
    bad = [name for name, digest in inputs.items() if sha(project_path(name)) != digest]
    if bad:
        raise ArithmeticError('Frozen input changed: '+', '.join(bad))


def load_problem():
    """Return raw first-solve columns and authenticated problem metadata.

    Only cached matrices are loaded; the one least-squares first step is
    repeated to recover its actual coefficient columns. No assembly occurs.
    """
    authenticate(EXPECTED)
    parent = json.loads(project_path(PARENT+'_results.json').read_text(encoding='utf-8'))
    if not parent['complete'] or not parent['passed']:
        raise ArithmeticError('Checkpoint157 did not complete successfully.')
    inputs = {**parent['input_sha256'], **EXPECTED}
    authenticate(inputs)
    base, metadata, exterior, baseline_inputs = load_baseline()
    inputs.update(baseline_inputs)
    source_indices = select_columns(metadata, (3, 4), 3, 12)
    with np.load(project_path(PARENT+'_state.npz')) as saved:
        raw = {key: saved['full_'+key].copy() for key in KEYS}
        regions = {}
        for label, old in (('exterior', exterior),
                           ('bridge', {key: base[key]-exterior[key] for key in KEYS})):
            region = {key: saved[label+'_increment_'+key].copy() for key in KEYS}
            for key in KEYS:
                region[key][..., :344, :344] = old[key]
            regions[label] = region
    if not all(np.array_equal(raw[key][..., :344, :344], base[key]) for key in KEYS):
        raise ArithmeticError('The nominal trial no longer contains the frozen old344 block.')
    record, columns = first_residual_regions(raw, regions, source_indices,
        kappa=.7, mass=.1, return_columns=True)
    target = next(row for row in parent['records'] if row['label'] == 'first_pair')
    difference = max(abs(record['regions'][name]['residual_norm_squared']
        -target['first_regions']['regions'][name]['residual_norm_squared'])
        for name in ('full', 'bridge', 'exterior'))
    inputs.update({name: sha(project_path(name)) for name in
        ('dihole_radial_stream.py', 'dihole_mouth_fields.py', 'dihole_enriched_bridge.py',
         'dihole_residual_regions.py', 'dihole_stretched_collar.py', 'dihole_sewn_form.py')})
    return columns, dict(input_sha256=inputs, source_indices=source_indices,
        first_solve=record, saved_first_solve_difference=difference,
        target_bridge=target['first_regions']['regions']['bridge'],
        parameters=dict(mu=.5, phase=0., kappa=.7, first_mass=.1, eta=.01,
                        source_dimension=216, trial_dimension=356))


def bridge_grid(order):
    """Exactly the panel ordering consumed by mouth_fields.bridge_chunks."""
    x, w = gauss(order)
    cmax = np.sqrt(1-(RHO_BOUNDS[0]/RC)**2)
    angular_panels = ((-1., -cmax), (-cmax, cmax), (cmax, 1.))
    def nodes(panels):
        return (np.concatenate([(hi+lo)/2+(hi-lo)*x/2 for lo, hi in panels]),
                np.concatenate([(hi-lo)*w/2 for lo, hi in panels]))
    s, ws = nodes(RADIAL_PANELS)
    cosine, wc = nodes(angular_panels)
    return dict(s=s, radial_weights=ws, cosine=cosine, cosine_weights=wc,
                theta=np.arccos(cosine), angular_panels=angular_panels)


def bridge_profile(columns, order=64, slices_per_chunk=4, progress=False):
    """Evaluate positive/negative source-column residuals on the actual bridge.

    Returned residual_density has axes (radial node, cosine node, category),
    with category [total, positive insertion eigenvalue, negative eigenvalue].
    Each category is a sum of nonnegative squared residuals, not a signed
    stress contribution. Raw coefficients retain all216 source directions.
    """
    if isinstance(order, bool) or not isinstance(order, int) or order <= 12:
        raise ValueError('An integer quadrature order greater than12 is required.')
    grid = bridge_grid(order)
    ns, nc = len(grid['s']), len(grid['cosine'])
    values = columns['source_eigenvalues']
    if columns['source'].shape != (356, 216) or columns['trial'].shape != (356, 216):
        raise ValueError('Require the saved nominal356 first-solve coefficient chart.')
    positive, negative = values > 0, values < 0
    if np.any(~(positive | negative)):
        raise ArithmeticError('A zero source eigenvalue needs an explicit direction convention.')
    trial = np.ascontiguousarray(columns['trial'][LOCAL_TO_RAW])
    combination = np.ascontiguousarray((.1*columns['trial']-columns['source'])[LOCAL_TO_RAW])
    density = np.empty((ns, nc, 3))
    offset, weight_error, elapsed = 0, 0., perf_counter()
    largest_chunk = 0
    for chunk_number, data in enumerate(bridge_chunks(phase=0., order=order,
                                                    slices_per_chunk=slices_per_chunk)):
        points, spin, local = data['value'].shape
        if spin != 4 or local != len(LOCAL_TO_RAW) or points % nc:
            raise ArithmeticError('Pointwise bridge chart or radial ordering changed.')
        count = points//nc
        expected_weight = (grid['radial_weights'][offset:offset+count, None]
                           *grid['cosine_weights'][None, :]).ravel()
        weight_error = max(weight_error, float(np.max(abs(data['weight']-expected_weight))))
        action = data['radial']+data['angular']+.7j*data['W'][:, None, None]*data['value']
        residual = (action.reshape(-1, local)@trial).reshape(points, 4, 216)
        residual = residual[:, [2, 3, 0, 1], :]/np.sqrt(data['P'])[:, None, None]
        residual += (data['value'].reshape(-1, local)@combination).reshape(points, 4, 216) \
                    *np.sqrt(data['P'])[:, None, None]
        squared = abs(residual)**2
        pos = np.sum(squared[:, :, positive], axis=(1, 2))
        neg = np.sum(squared[:, :, negative], axis=(1, 2))
        density[offset:offset+count] = np.stack((pos+neg, pos, neg), axis=1).reshape(count, nc, 3)
        largest_chunk = max(largest_chunk, sum(v.nbytes for v in data.values())
                            +action.nbytes+residual.nbytes+squared.nbytes)
        offset += count
        if progress and chunk_number == 0:
            print(json.dumps(dict(stage='first_chunk', radial_slices=count,
                seconds=perf_counter()-elapsed, working_chunk_bytes=largest_chunk)), flush=True)
        if progress and offset % order == 0:
            total = np.einsum('s,c,scj->j', grid['radial_weights'][:offset],
                              grid['cosine_weights'], density[:offset])
            print(json.dumps(dict(stage='radial_panel_complete', completed_slices=offset,
                cumulative_residual_squared=total.tolist(), seconds=perf_counter()-elapsed)), flush=True)
    if offset != ns:
        raise ArithmeticError('The full radial grid was not consumed.')
    radial = np.einsum('c,scj->sj', grid['cosine_weights'], density)
    angular = np.einsum('s,scj->cj', grid['radial_weights'], density)
    return {**grid, 'residual_density': density, 'radial_density': radial,
            'cosine_density': angular, 'theta_density': angular*np.sin(grid['theta'])[:, None]}, \
        dict(weight_maximum_absolute_error=weight_error, largest_working_chunk_bytes=largest_chunk,
             seconds=perf_counter()-elapsed, all_radial_slices=offset)


def summarize(profile, order):
    density, ws, wc = [profile[k] for k in ('residual_density', 'radial_weights', 'cosine_weights')]
    total = np.einsum('s,c,scj->j', ws, wc, density)
    radial_names = ('left_collar', 'central', 'right_collar')
    angular_names = ('cosine_negative_polar_cap', 'middle_angles', 'cosine_positive_polar_cap')
    rows = []
    for i, name in enumerate(radial_names):
        at = slice(i*order, (i+1)*order)
        amount = np.einsum('s,c,scj->j', ws[at], wc, density[at])
        rows.append(dict(region=name, s_interval=list(RADIAL_PANELS[i]),
            residual_squared=amount.tolist(), fraction_of_bridge=float(amount[0]/total[0])))
    angular = []
    for i, name in enumerate(angular_names):
        at = slice(i*order, (i+1)*order)
        amount = np.einsum('s,c,scj->j', ws, wc[at], density[:, at])
        angular.append(dict(region=name, cosine_interval=list(profile['angular_panels'][i]),
            residual_squared=amount.tolist(), fraction_of_bridge=float(amount[0]/total[0])))
    peak = np.unravel_index(np.argmax(density[:, :, 0]), density.shape[:2])
    radial_peak = int(np.argmax(profile['radial_density'][:, 0]))
    return dict(total_residual_squared=total.tolist(), radial_regions=rows, angular_regions=angular,
        sampled_density_peak=dict(s=float(profile['s'][peak[0]]), theta=float(profile['theta'][peak[1]]),
                                  density=float(density[peak[0], peak[1], 0])),
        sampled_radial_peak=dict(s=float(profile['s'][radial_peak]),
                                density=float(profile['radial_density'][radial_peak, 0])))


def run():
    start, source_hash = perf_counter(), sha(__file__)
    columns, problem = load_problem()
    profile, timing = bridge_profile(columns, progress=True)
    summary = summarize(profile, 64)
    total = np.array(summary['total_residual_squared'])
    target = problem['target_bridge']
    expected = np.array([target['residual_norm_squared'], target['positive_source_residual_squared'],
                         target['negative_source_residual_squared']])
    reproduction = float(np.max(abs(total-expected))/max(1., expected[0]))
    radial_total = profile['radial_weights']@profile['radial_density']
    angular_total = profile['cosine_weights']@profile['cosine_density']
    native = dict(
        input_provenance=sha(__file__) == source_hash,
        saved_first_solve=problem['saved_first_solve_difference'] < 1e-7,
        saved_bridge_full_action=reproduction < 3e-11,
        quadrature_coordinates=timing['weight_maximum_absolute_error'] < 2e-18,
        nonnegative_finite_density=np.isfinite(profile['residual_density']).all()
            and np.min(profile['residual_density']) >= 0,
        signed_direction_sum=np.max(abs(profile['residual_density'][:, :, 0]
            -profile['residual_density'][:, :, 1]-profile['residual_density'][:, :, 2]))
            /max(1., np.max(profile['residual_density'])) < 3e-16,
        marginal_integrals=max(np.max(abs(radial_total-total)), np.max(abs(angular_total-total)))
            /max(1., total[0]) < 3e-13,
        exact_radial_partition=abs(sum(row['residual_squared'][0] for row in summary['radial_regions'])-total[0])
            /max(1., total[0]) < 3e-13,
        first_solve_normal_equation=problem['first_solve']['least_squares_relative_error'] < 2e-12)
    authenticate(problem['input_sha256'])
    state = {key: value for key, value in profile.items() if key != 'angular_panels'}
    state.update(columns)
    state.update(source_indices=problem['source_indices'], local_to_raw=LOCAL_TO_RAW)
    state_path = project_path(NAME+'_state.npz')
    np.savez_compressed(state_path, **state)
    result = dict(complete=True, passed=bool(all(native.values())),
        tests={key: bool(value) for key, value in native.items()},
        parameters={**problem['parameters'], 'quadrature_order':64, 'radial_slices':192,
                    'polar_nodes':192, 'measure':'ds*dcos(theta)',
                    'categories':['total', 'positive_source_directions', 'negative_source_directions']},
        summary=summary, target_bridge=target, saved_bridge_relative_error=reproduction,
        saved_first_solve_difference=problem['saved_first_solve_difference'],
        source_positive_negative_ranks=problem['first_solve']['source_positive_negative_ranks'],
        first_solve=problem['first_solve'], evaluation=timing,
        equations=dict(residual='P^(-1/2)*beta*(radial+angular+i*kappa*W*value)*trial + P^(1/2)*value*(m*trial-source)',
                       density='sum over physical spin components and retained source directions of |residual|^2',
                       bridge_norm='integral ds dcos(theta) density; beta permutes components [2,3,0,1]'),
        scope=['The fixed nominal356 trial and retained216-column insertion of checkpoint157 are unchanged.',
               'Positive/negative labels classify insertion eigenvectors; both residual squared densities are nonnegative.',
               'Only bridge fields are evaluated; saved exterior matrices enter the one first-step solve.',
               'The central/collar and polar-cap partitions coincide with quadrature panel boundaries.',
               'Sampled peaks are grid diagnostics; no interval quadrature, continuum source, physical stress or support bound follows.'],
        source_sha256=source_hash, input_sha256=problem['input_sha256'],
        state_sha256=sha(state_path), seconds=perf_counter()-start)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], summary=summary,
                         reproduction=reproduction, seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Residual localization checks failed; evidence saved.')
    return result


if __name__ == '__main__':
    run()
