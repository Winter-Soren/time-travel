"""Independent cached-matrix audit of checkpoint159's fixed-frequency solve.

No geometry, field values, exterior integration, or production SVD is run.
The first residual is recovered directly in raw coordinates; an independent
source-space Hermitian eigenproblem checks all singular-value gains.
Condition and cancellation figures are diagnostics, not roundoff enclosures.
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

from project_paths import project_path


NAME, PRODUCER, PRIOR = 'dihole_collar_audit', 'dihole_collar_selection', 'dihole_residual_selection'
OLD = 628


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(inputs):
    bad = [name for name, digest in inputs.items() if sha(project_path(name)) != digest]
    if bad:
        raise ArithmeticError('Saved input mismatch: '+', '.join(bad))


def hermitian(a):
    return (a+a.conj().T)/2


def direct_residual(normal, rhs, source_square):
    solution = solve(normal, rhs, assume_a='pos')
    recovered = float(np.vdot(rhs, solution).real)
    residual = source_square-recovered
    expanded = source_square+float(np.vdot(solution, normal@solution).real)-2*recovered
    return dict(residual_squared=residual, expanded_residual_squared=expanded,
        expanded_difference=expanded-residual,
        relative_normal_residual=float(np.linalg.norm(normal@solution-rhs)/max(1., np.linalg.norm(rhs))),
        scalar_subtraction_factor=(source_square+abs(recovered))/max(abs(residual), np.finfo(float).tiny))


def compressed_normal(normal, rhs, transform):
    return (np.block([[normal[:OLD, :OLD], normal[:OLD, OLD:]@transform],
                      [transform.conj().T@normal[OLD:, :OLD],
                       transform.conj().T@normal[OLD:, OLD:]@transform]]),
            np.concatenate((rhs[:OLD], transform.conj().T@rhs[OLD:])))


def run():
    start, own_hash = perf_counter(), sha(__file__)
    parent = json.loads(project_path(PRODUCER+'_results.json').read_text(encoding='utf-8'))
    if not parent['complete'] or not parent['passed'] or not all(parent['tests'].values()):
        raise ArithmeticError('The producer must finish successfully before this cached audit.')
    inputs = {**parent['input_sha256'], PRODUCER+'.py': parent['source_sha256'],
              PRODUCER+'_state.npz': parent['state_sha256'],
              PRODUCER+'_results.json': sha(project_path(PRODUCER+'_results.json'))}
    authenticate(inputs)
    previous = json.loads(project_path(PRIOR+'_results.json').read_text(encoding='utf-8'))
    with np.load(project_path(PRODUCER+'_state.npz')) as archive:
        saved = {key: archive[key].copy() for key in archive.files}
    with np.load(project_path(PRIOR+'_state.npz')) as archive:
        prior = {key: archive['full_'+key].copy() for key in ('norm', 'q', 'weak', 'derivative')}
        prior_source = archive['source_coefficients'].copy()
    g, q, w = [saved['full_'+key] for key in ('gram', 'action_gram', 'weak_action')]
    indices, insertion = saved['retained_indices'], saved['retained_insertion']
    source, transform = saved['source_coefficients'], saved['selection_transform']
    kappa, mass = parent['parameters']['kappa'], parent['parameters']['mass']
    if kappa != .7 or mass != .1 or g.shape != (892, 892) or source.shape != (892, 216):
        raise ArithmeticError('This audit is for the stated892-field kappa=.7, mass=.1 model.')
    expected = (prior['norm'], sum(kappa**j*c for j, c in enumerate(prior['q'])),
                prior['weak'][0]+kappa*prior['weak'][1])
    old_exact = all(np.array_equal(a[:OLD, :OLD], b) for a, b in zip((g, q, w), expected))
    expected_insertion = sum(kappa**j*c for j, c in enumerate(prior['derivative']))[np.ix_(indices, indices)]
    source_exact = np.array_equal(source[:OLD], prior_source) and not np.any(source[OLD:])
    normal = hermitian(q+mass*(w+w.conj().T)+mass**2*g)
    rhs = (w.conj().T+mass*g)@source
    source_square = float(np.vdot(source, g@source).real)
    old, cross = normal[:OLD, :OLD], normal[:OLD, OLD:]
    polarization = cross.conj().T@solve(old, cross, assume_a='pos')
    schur = hermitian(normal[OLD:, OLD:]-polarization)
    old_load = solve(old, rhs[:OLD], assume_a='pos')
    load_correction = cross.conj().T@old_load
    reduced = rhs[OLD:]-load_correction
    gain_gram = hermitian(reduced.conj().T@solve(schur, reduced, assume_a='pos'))
    gains = eigh(gain_gram, eigvals_only=True)[::-1]
    expected_gains = np.asarray(parent['selection']['gains'])
    rows = []
    for count, target in zip((0, 8, 16, 32, 64), parent['records']):
        block, load = compressed_normal(normal, rhs, transform[:, :count])
        row = direct_residual(block, load, source_square)
        row.update(selected=count, saved_difference=row['residual_squared']
            -target['first_regions']['regions']['full']['residual_norm_squared'])
        rows.append(row)
    full = direct_residual(normal, rhs, source_square)
    full['saved_difference'] = full['residual_squared']-parent['full_pool']['first_regions']['regions']['full']['residual_norm_squared']
    refined = {key: saved['refined_increment_'+key].copy() for key in ('gram', 'action_gram', 'weak_action')}
    refined_old_zero = all(not np.any(a[:OLD, :OLD]) for a in refined.values())
    for key, old_value in zip(('gram', 'action_gram', 'weak_action'), expected):
        refined[key][:OLD, :OLD] = old_value
    rg, rq, rw = [refined[key] for key in ('gram', 'action_gram', 'weak_action')]
    rn = hermitian(rq+mass*(rw+rw.conj().T)+mass**2*rg)
    rrhs = (rw.conj().T+mass*rg)@source
    block, load = compressed_normal(rn, rrhs, transform[:, :64])
    refined_check = direct_residual(block, load, float(np.vdot(source, rg@source).real))
    refined_check['saved_difference'] = refined_check['residual_squared']-parent['refined']['first_regions']['regions']['full']['residual_norm_squared']
    groups = {}
    for label in ('collar', 'global'):
        at = np.array([r['index'] for r in parent['candidates'] if r['group'] == label])
        value = float(np.vdot(reduced[at], solve(schur[np.ix_(at, at)], reduced[at], assume_a='pos')).real)
        groups[label] = dict(gain=value, saved_difference=value-parent['group_gains'][label])
    diagonal = 1/np.sqrt(np.diag(schur).real)
    balanced = hermitian(diagonal[:, None]*schur*diagonal[None, :])
    se = eigh(balanced, eigvals_only=True)
    normal_diagonal = 1/np.sqrt(np.diag(normal).real)
    ne = eigh(hermitian(normal_diagonal[:, None]*normal*normal_diagonal[None, :]), eigvals_only=True)
    condition = dict(balanced_schur_minimum=float(se[0]), balanced_schur_maximum=float(se[-1]),
        balanced_schur_condition=float(se[-1]/se[0]), balanced_full_normal_condition=float(ne[-1]/ne[0]),
        schur_frobenius_subtraction_factor=float((np.linalg.norm(normal[OLD:, OLD:])+np.linalg.norm(polarization))/np.linalg.norm(schur)),
        reduced_rhs_subtraction_factor=float((np.linalg.norm(rhs[OLD:])+np.linalg.norm(load_correction))/np.linalg.norm(reduced)),
        machine_epsilon=float(np.finfo(float).eps),
        scope='Basis-dependent floating condition and cancellation diagnostics, not certified arithmetic errors.')
    baseline_chain_errors = []
    for current, old_chain in zip(parent['records'][0]['common_filter'], previous['full_pool']['common_filter']):
        for key in ('trial_filter_response', 'conditional_retained_source_interval', 'positive_negative_error_bounds'):
            a, b = np.asarray(current[key]), np.asarray(old_chain[key])
            baseline_chain_errors.append(float(np.max(abs(a-b))/max(1., float(np.max(abs(b))))))
        baseline_chain_errors.extend(abs(a['total_residual_squared']-b['total_residual_squared'])
            /max(1., abs(b['total_residual_squared'])) for a, b in zip(current['stages'], old_chain['stages']))
    gain_error = float(np.max(abs(gains-expected_gains)))
    all_solves = rows+[full, refined_check]
    tests = dict(exact_provenance=sha(__file__) == own_hash,
        fixed_frequency_archive=all('full_'+key in saved for key in ('gram', 'action_gram', 'weak_action'))
            and not any(key in saved for key in ('q', 'full_q', 'derivative', 'full_derivative')),
        exact_old628=old_exact, exact_source=source_exact,
        exact_retained_insertion=np.array_equal(insertion, expected_insertion), refined_old_block_zero=refined_old_zero,
        raw_residual_reproduction=max(abs(r['saved_difference']) for r in all_solves) < 2e-6,
        normal_equation=max(r['relative_normal_residual'] for r in all_solves) < 2e-10,
        independent_gain_spectrum=gain_error/max(1., expected_gains[0]) < 2e-9,
        gain_total_and_groups=max([abs(float(np.trace(gain_gram).real)-parent['selection']['total_gain'])]
                                  +[abs(r['saved_difference']) for r in groups.values()]) < 2e-6,
        adapter_baseline_chain=max(baseline_chain_errors) < 2e-9,
        expanded_residual_identity=max(abs(r['expanded_difference']) for r in all_solves) < 2e-7,
        positive_balanced_matrices=se[0] > 0 and ne[0] > 0)
    authenticate(inputs)
    result = dict(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        parameters=dict(kappa=kappa, first_mass=mass, old_dimension=OLD, total_dimension=892, source_dimension=216),
        raw_solve_rows=rows, full_pool=full, same64_refined=refined_check, group_gains=groups,
        gain_spectrum_maximum_absolute_error=gain_error, smallest_computed_gain=float(gains[-1]),
        independent_total_gain=float(np.trace(gain_gram).real), condition=condition,
        adapter_baseline_chain_maximum_relative_error=max(baseline_chain_errors),
        equations=dict(normal='N=M+m(B+B*)+m^2 G; rhs=(B*+m G)h',
            residual='tr(h*Gh)-Re tr(rhs*N^-1 rhs)',
            gain='eigenvalues of Q* S^-1 Q, S=N_HH-N_HL N_LL^-1 N_LH, Q=rhs_H-N_HL N_LL^-1 rhs_L'),
        scope=['Independent raw saved-matrix solves and source-space eigenvalues; no assembly or production SVD rerun.',
               'The ephemeral adapter supplies only the actual VxV insertion read by the fixed-frequency callers.',
               'First-step gains optimize absolute-eigenvalue-weighted errors, not the signed four-filter source.',
               'Floating quadrature, omitted fields, other frequencies/sectors, and physical source errors are not enclosed.'],
        source_sha256=own_hash, input_sha256=inputs, seconds=perf_counter()-start)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], rows=rows, full=full,
                         groups=groups, condition=condition, seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Independent cached collar audit failed; evidence saved.')
    return result


if __name__ == '__main__':
    run()
