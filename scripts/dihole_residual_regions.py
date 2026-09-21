"""First common-filter residual, split over the bridge and exterior.

The source V |D_V|^(1/2) is fixed while the trial space changes. Regional
weak matrices need not be skew: their mass cross terms are retained.
Saved forms use floating quadrature, not a continuum error enclosure.
Only the old exterior48 forms are assembled anew and cached for reuse.
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

from dihole_bridge_enrichment import raw_polynomial, select_columns
from dihole_exterior_feedback import append_exterior
from dihole_sewn_form import exterior_data
from project_paths import project_path


NAME = 'dihole_residual_regions'
KEYS = ('norm', 'q', 'derivative', 'weak')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _authenticate(inputs):
    failed = [name for name, digest in inputs.items() if sha(project_path(name)) != digest]
    if failed:
        raise ArithmeticError('Saved input hash mismatch: '+', '.join(failed))


def assemble_exterior_polynomial(order=24):
    """One explicit regional assembly, not the old full calculation."""
    return raw_polynomial(exterior_data(.5, order, (3, 4)))


def load_exterior_polynomial():
    """Load the frozen48-column cache after checking its complete provenance."""
    result_path = project_path(NAME+'_results.json')
    record = json.loads(result_path.read_text(encoding='utf-8'))
    path = project_path(NAME+'_state.npz')
    if not record['passed'] or record['source_sha256'] != sha(__file__) or record['state_sha256'] != sha(path):
        raise ArithmeticError('Regional cache provenance failed.')
    _authenticate(record['input_sha256'])
    with np.load(path) as archive:
        return {key: archive[key].copy() for key in KEYS}


def load_baseline(exterior=None):
    """Return (raw344, metadata142, exterior344, input_hashes).

    Default calls reuse the authenticated cache. Production supplies its one
    freshly assembled exterior48 form explicitly. The old312 block is exact;
    only the saved141 new/cross blocks are appended. No old-old block is added.
    """
    inputs, saved, records = {}, {}, {}
    for name in ('dihole_angular_crosscheck', 'dihole_exterior_feedback'):
        record_path = project_path(name+'_results.json')
        record = json.loads(record_path.read_text(encoding='utf-8'))
        if not record['passed']:
            raise ArithmeticError('A baseline did not pass: '+name)
        checks = dict(record['input_sha256'])
        checks.update({name+'.py': record['source_sha256'], name+'_state.npz': record['state_sha256']})
        _authenticate(checks)
        inputs.update(checks)
        inputs[record_path.name] = sha(record_path)
        records[name] = record
        with np.load(project_path(name+'_state.npz')) as archive:
            saved[name] = {key: archive['p0_'+key] for key in KEYS}
    if exterior is None:
        exterior = load_exterior_polynomial()
        for name in (NAME+'.py', NAME+'_results.json', NAME+'_state.npz'):
            inputs[name] = sha(project_path(name))
    old = saved['dihole_angular_crosscheck']
    extension = saved['dihole_exterior_feedback']
    indices = np.r_[np.arange(48), np.arange(216, 248)]
    extra = {key: extension[key][..., indices, :][..., indices].copy() for key in KEYS}
    for key in KEYS:
        if np.any(extension[key][..., 48:216, 216:]) or np.any(extension[key][..., 216:, 48:216]):
            raise ArithmeticError('Compact bridge fields overlap the exterior extension.')
        extra[key][..., :48, :48] = 0
    combined = append_exterior(old, extra)
    regional = {key: np.zeros_like(combined[key]) for key in KEYS}
    for key in KEYS:
        if exterior[key].shape != old[key].shape[:-2]+(48, 48):
            raise ValueError('The cache must contain the old exterior48 polynomial.')
        regional[key][..., :48, :48] = exterior[key]
        regional[key][..., :48, 312:] = extra[key][..., :48, 48:]
        regional[key][..., 312:, :48] = extra[key][..., 48:, :48]
        regional[key][..., 312:, 312:] = extra[key][..., 48:, 48:]
    return combined, records['dihole_angular_crosscheck']['metadata'], regional, inputs


def first_residual_regions(raw, regions, source_indices, kappa=.7, mass=.1,
                           trial_indices=None, return_columns=False):
    """First least-squares solve with one fixed, source-weighted load.

    All region forms use the SAME ambient raw field coordinates as raw.
    source_indices must lie within trial_indices. With return_columns=True,
    also return raw source/trial coefficients and signed source eigenvalues;
    these permit independent pointwise residual localization.
    """
    if not np.isfinite(kappa) or not np.isfinite(mass) or mass <= 0:
        raise ValueError('Finite frequency and positive first-filter mass required.')
    dimension = len(raw['norm'])
    v = np.asarray(source_indices)
    y = np.arange(dimension) if trial_indices is None else np.asarray(trial_indices)
    for ids in (v, y):
        if ids.ndim != 1 or ids.dtype.kind not in 'iu' or len(set(ids.tolist())) != len(ids) or not len(ids) or np.min(ids) < 0 or np.max(ids) >= dimension:
            raise ValueError('Distinct valid integer field indices required.')
    if not set(v.tolist()).issubset(y.tolist()) or 'full' in regions or not regions:
        raise ValueError('The trial must contain the source; supply named non-full regions.')
    g = raw['norm']
    q, derivative = [sum(kappa**j*c for j, c in enumerate(raw[key])) for key in ('q', 'derivative')]
    weak = raw['weak'][0]+kappa*raw['weak'][1]
    iv = np.linalg.inv(np.linalg.cholesky(g[np.ix_(v, v)]))
    iy = np.linalg.inv(np.linalg.cholesky(g[np.ix_(y, y)]))
    dv = iv@derivative[np.ix_(v, v)]@iv.conj().T
    values, vectors = eigh((dv+dv.conj().T)/2)
    weighted = vectors*np.sqrt(abs(values))
    positive, negative = values > 0, values < 0
    source = np.zeros((dimension, len(v)), complex)
    source[v] = iv.conj().T@weighted
    my = iy@q[np.ix_(y, y)]@iy.conj().T
    b = iy@weak[np.ix_(y, y)]@iy.conj().T
    u = iy@g[np.ix_(y, v)]@iv.conj().T@weighted
    normal = my+mass*(b+b.conj().T)+mass**2*np.eye(len(y))
    normal = (normal+normal.conj().T)/2
    rhs = (b.conj().T+mass*np.eye(len(y)))@u
    new = solve(normal, rhs, assume_a='pos')
    trial = np.zeros_like(source)
    trial[y] = iy.conj().T@new
    regional = {}
    for name, rr in dict(full=raw, **regions).items():
        rg = rr['norm']
        rq = sum(kappa**j*c for j, c in enumerate(rr['q']))
        rw = rr['weak'][0]+kappa*rr['weak'][1]
        rn = rq+mass*(rw+rw.conj().T)+mass**2*rg
        source_squares = np.sum(source.conj()*(rg@source), axis=0).real
        residual_squares = (np.sum(trial.conj()*(rn@trial), axis=0).real+source_squares
            -2*np.sum(trial.conj()*((rw.conj().T+mass*rg)@source), axis=0).real)
        regional[name] = dict(source_norm_squared=float(sum(source_squares)),
            residual_norm_squared=float(sum(residual_squares)),
            positive_source_residual_squared=float(sum(residual_squares[positive])),
            negative_source_residual_squared=float(sum(residual_squares[negative])))
    inside = (b+mass*np.eye(len(y)))@new-u
    in_y = float(np.linalg.norm(inside, 'fro')**2)
    outside_y = float(np.vdot(new, (my-b.conj().T@b)@new).real)
    full = regional['full']['residual_norm_squared']
    if full < 0:
        raise ArithmeticError('Negative total residual square.')
    recombination = max(float(np.max(abs(raw[key]-sum(rr[key] for rr in regions.values()))))
        / max(1., float(np.max(abs(raw[key])))) for key in ('norm', 'q', 'weak'))
    result = dict(kappa=float(kappa), mass=float(mass), eta=float(mass**2),
        source_dimension=len(v), trial_dimension=len(y), ambient_dimension=dimension,
        source_positive_negative_ranks=[int(sum(positive)), int(sum(negative))],
        residual_norm=float(np.sqrt(full)), regions=regional,
        in_trial_residual_squared=in_y, unrepresented_action_residual_squared=outside_y,
        regional_form_sum_relative_error=recombination,
        regional_residual_sum_error=full-sum(regional[name]['residual_norm_squared'] for name in regions),
        regional_source_sum_error=regional['full']['source_norm_squared']
            -sum(regional[name]['source_norm_squared'] for name in regions),
        full_residual_identity_error=full-in_y-outside_y,
        least_squares_relative_error=float(np.linalg.norm(normal@new-rhs, 'fro')/max(1., np.linalg.norm(rhs, 'fro'))))
    if return_columns:
        return result, dict(source=source, trial=trial, source_eigenvalues=values)
    return result


def regional_control():
    """Independent strong-action check, including non-skew regional weak forms."""
    rng = np.random.default_rng(157)
    z = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    a = (z-z.conj().T)/2
    derivative = np.diag([2., -1., .3, 0., 0., 0., 0., 0.])
    def forms(mask):
        return dict(norm=mask, q=np.array([a.conj().T@mask@a, np.zeros((8, 8)), np.zeros((8, 8))]),
                    derivative=np.array([derivative, np.zeros((8, 8)), np.zeros((8, 8))]),
                    weak=np.array([mask@a, np.zeros((8, 8))]))
    masks = dict(first=np.diag([1.]*3+[0.]*5), second=np.diag([0.]*3+[1.]*5))
    record, columns = first_residual_regions(forms(np.eye(8)), {n: forms(p) for n, p in masks.items()},
        np.arange(3), trial_indices=np.arange(5), return_columns=True)
    residual = (a+.1*np.eye(8))@columns['trial']-columns['source']
    errors = [abs(float(np.vdot(residual, p@residual).real)-record['regions'][n]['residual_norm_squared'])
              for n, p in dict(full=np.eye(8), **masks).items()]
    return dict(maximum_direct_action_error=max(errors),
                non_skew_regional_weak_norm=float(np.linalg.norm(masks['first']@a-a@masks['first'])))


def run():
    start, source_hash = perf_counter(), sha(__file__)
    exterior = assemble_exterior_polynomial()
    raw, metadata, exterior_raw, inputs = load_baseline(exterior)
    prior_path = project_path('dihole_common_pv_residual_results.json')
    prior = json.loads(prior_path.read_text(encoding='utf-8'))
    if not prior['passed']:
        raise ArithmeticError('The frozen common-filter control did not pass.')
    inputs.update(prior['input_sha256'])
    inputs.update({prior_path.name: sha(prior_path), 'dihole_common_pv_residual.py': prior['source_sha256']})
    _authenticate(inputs)
    bridge = {key: raw[key]-exterior_raw[key] for key in KEYS}
    regions = dict(bridge=bridge, exterior=exterior_raw)
    source = select_columns(metadata, (3, 4), 3, 12)
    trials = [(f'angular{angular}_exterior48', select_columns(metadata, (3, 4), angular, 12))
              for angular in (3, 4, 5)]+[('angular5_exterior80', np.arange(344))]
    rows = []
    for label, indices in trials:
        rows.append(dict(label=label, **first_residual_regions(raw, regions, source, trial_indices=indices)))
    baseline = next(r for r in prior['records'] if r['eta'] == .01)['stages'][0]
    measured = rows[2]
    recovery = max(abs(measured[a]-baseline[b])/max(1., abs(baseline[b])) for a, b in
        [('in_trial_residual_squared', 'in_space_residual_squared'),
         ('unrepresented_action_residual_squared', 'complement_residual_squared')])
    recovery = max(recovery, abs(measured['regions']['full']['residual_norm_squared']
        -baseline['total_residual_squared'])/baseline['total_residual_squared'])
    control = regional_control()
    tests = dict(exact_input_provenance=sha(__file__) == source_hash,
        disjoint_bridge_extension=all(np.all(bridge[k][..., 312:, :] == 0)
            and np.all(bridge[k][..., :, 312:] == 0) for k in KEYS),
        exterior_probe_zero=bool(np.all(exterior_raw['derivative'] == 0)),
        frozen156_recovery=recovery < 2e-12,
        regional_and_full_identities=max(max(abs(r['regional_residual_sum_error']), abs(r['regional_source_sum_error']),
            abs(r['full_residual_identity_error']))/r['regions']['full']['source_norm_squared'] for r in rows) < 2e-12
            and max(r['regional_form_sum_relative_error'] for r in rows) < 2e-14,
        least_squares=max(r['least_squares_relative_error'] for r in rows) < 2e-12,
        same_source=max(abs(r['regions']['full']['source_norm_squared']-rows[0]['regions']['full']['source_norm_squared']) for r in rows) < 1e-10,
        nested_residual_nonincrease=all(rows[i+1]['regions']['full']['residual_norm_squared'] <= rows[i]['regions']['full']['residual_norm_squared']+1e-8 for i in range(3)),
        independent_regional_action=control['maximum_direct_action_error'] < 2e-12 and control['non_skew_regional_weak_norm'] > .1)
    _authenticate(inputs)
    state_path = project_path(NAME+'_state.npz')
    np.savez_compressed(state_path, **exterior)
    result = dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        parameters=dict(phase=0., mu=.5, regulator_mass=1., kappa=.7, first_mass=.1, eta=.01,
            source_dimension=216, source_radial_modes=12, source_angular_modes=3,
            exterior_assembly_order=24, cache_columns=48),
        records=rows, independent_control=control, frozen156_relative_error=recovery,
        source_sha256=source_hash, input_sha256=inputs, state_sha256=sha(state_path),
        cache_arrays={k: dict(shape=list(exterior[k].shape), dtype=str(exterior[k].dtype),
            sha256_C_bytes=hashlib.sha256(np.ascontiguousarray(exterior[k]).tobytes()).hexdigest()) for k in KEYS},
        seconds=perf_counter()-start,
        scope=['Same finite-rank source, first common-filter step only; no full source or completed filter bound.',
               'Regional squared norms retain the non-skew regional weak cross terms.',
               'Bridge/exterior split is geometric; it does not identify radial versus mouth-trace deficiencies.',
               'Floating quadrature and arithmetic are not enclosed; no continuum, omitted-insertion or azimuthal-tail certificate.',
               'The cached48-column exterior form is newly assembled once; frozen140-142 files are unchanged.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], seconds=result['seconds'],
        source_sha256=source_hash, state_sha256=result['state_sha256'],
        residual_squared=[r['regions']['full']['residual_norm_squared'] for r in rows]), indent=2))
    if not result['passed']:
        raise ArithmeticError('Regional residual checks failed; evidence saved.')
    return result


if __name__ == '__main__':
    run()
