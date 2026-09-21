"""First-order residual bound for the common positive PV filter, at fixed eta>0.

For skew-adjoint A, P=A* A and m_j=sqrt(eta+j M^2), the commuting filter
T=sqrt(6) M^3 prod_j(A+m_j)^-1 satisfies T* T=6 M^6/prod_j(P+m_j^2).
For orthonormal conforming Y, retain BOTH MY=(AY)*AY and B=Y*AY.
Delta=MY-B*B is the unrepresented action Gram, not an optional correction.
Each least-squares step has residual Gram
  ((B+m)u_new-u_old)*((B+m)u_new-u_old)+u_new*Delta u_new.
The error recurrence e_new=(e_old+||residual||F)/m uses only first derivatives.
Source-positive/negative column groups give a sign-aware squared-norm interval.

Numbers below use saved floating integrals, not continuum enclosures. This
bounds only a finite-rank insertion conditionally on exact forms; no physical
stress, eta=0 limit, omitted source insertion or angular tail is certified.
"""
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '1')

import numpy as np
from scipy.linalg import eigh, solve

from dihole_bridge_enrichment import select_columns
from dihole_retained_resolvent import normalized
from project_paths import project_path


NAME = 'dihole_common_pv_residual'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hermitian(matrix):
    return (matrix + matrix.conj().T) / 2


def chain(my, b, derivative, overlap, eta, mass=1.):
    """Return conditional scalar bounds and final weighted trial columns."""
    if not np.isfinite(eta) or eta <= 0 or not np.isfinite(mass) or mass <= 0:
        raise ValueError('The common filter requires strictly positive eta and regulator mass.')
    values, vectors = eigh(hermitian(derivative))
    u = (overlap @ vectors) * np.sqrt(abs(values))
    masks = (values > 0, values < 0)
    delta = hermitian(my - b.conj().T @ b)
    delta_min = float(eigh(delta, eigvals_only=True, subset_by_index=(0, 0))[0])
    delta_scale = max(1., float(np.linalg.norm(delta, 'fro')))
    if delta_min < -2e-10 * delta_scale:
        raise ArithmeticError('The measured full-action complement is materially indefinite.')
    errors = np.zeros(2)
    records = []
    eye = np.eye(len(my))
    for j in range(4):
        m = float(np.sqrt(eta + j * mass**2))
        # The tiny measured skew defect is retained rather than silently zeroed.
        normal = hermitian(my + m * (b + b.conj().T) + m*m * eye)
        rhs = (b.conj().T + m * eye) @ u
        new = solve(normal, rhs, assume_a='pos')
        inside = (b + m * eye) @ new - u
        inside_gram = inside.conj().T @ inside
        complement_gram = new.conj().T @ delta @ new
        gram = hermitian(inside_gram + complement_gram)
        direct = hermitian(new.conj().T @ normal @ new
            - new.conj().T @ (b.conj().T + m * eye) @ u
            - u.conj().T @ (b + m * eye) @ new + u.conj().T @ u)
        squares = np.array([float(np.trace(gram[np.ix_(mask, mask)]).real) for mask in masks])
        tolerance = 2e-11 * max(1., float(np.linalg.norm(gram, 'fro')))
        if np.min(squares) < -tolerance:
            raise ArithmeticError('A residual squared norm is materially negative.')
        residuals = np.sqrt(np.maximum(squares, 0.))
        errors = (errors + residuals) / m
        records.append(dict(index=j, squared_mass=m*m, mass=m,
            inverse_coercivity=1/m, pv_coefficient=(1, -3, 3, -1)[j],
            source_coefficient=(-.5, 1.5, -1.5, .5)[j],
            positive_negative_residual_norms=residuals.tolist(),
            in_space_residual_squared=float(np.trace(inside_gram).real),
            complement_residual_squared=float(np.trace(complement_gram).real),
            total_residual_squared=float(np.trace(gram).real),
            gram_identity_relative_error=float(np.linalg.norm(gram-direct, 'fro'))
                / max(1., float(np.linalg.norm(gram, 'fro'))),
            least_squares_relative_error=float(np.linalg.norm(normal @ new-rhs, 'fro'))
                / max(1., float(np.linalg.norm(rhs, 'fro'))),
            negative_square_roundoff_clamps=int(np.count_nonzero(squares < 0))))
        u = new
    scale = np.sqrt(6.) * mass**3
    output = scale * u
    norms = np.array([np.linalg.norm(output[:, mask], 'fro') for mask in masks])
    errors *= scale
    lower = np.maximum(0., norms-errors)**2
    upper = (norms+errors)**2
    interval = [.5*(lower[1]-upper[0]), .5*(upper[1]-lower[0])]
    return dict(eta=eta, regulator_mass=mass, filter_scale=float(scale),
        source_positive_negative_ranks=[int(np.count_nonzero(mask)) for mask in masks],
        trial_filter_response=float(-.5*(norms[0]**2-norms[1]**2)),
        positive_negative_output_norms=norms.tolist(), positive_negative_error_bounds=errors.tolist(),
        positive_negative_squared_norm_intervals=list(zip(lower.tolist(), upper.tolist())),
        conditional_retained_source_interval=[float(x) for x in interval],
        sign_resolved=bool(interval[0] > 0 or interval[1] < 0),
        delta_minimum_eigenvalue=delta_min, delta_frobenius_norm=delta_scale,
        normalized_weak_skew_relative_error=float(np.linalg.norm(b+b.conj().T, 'fro'))
            / max(1., float(np.linalg.norm(b, 'fro'))),
        overlap_isometry_error=float(np.linalg.norm(overlap.conj().T @ overlap
            - np.eye(overlap.shape[1]), 2)), stages=records), output


def controls():
    """Complete and incomplete trial spaces for an actual finite skew operator."""
    rng = np.random.default_rng(156)
    random = rng.normal(size=(13, 13)) + 1j*rng.normal(size=(13, 13))
    vectors = np.linalg.qr(random)[0]
    a = 1j*(vectors*np.arange(-6., 7.)) @ vectors.conj().T
    derivative = hermitian(random[:5, :5])
    d, z = eigh(derivative)
    source = np.eye(13, 5) @ (z*np.sqrt(abs(d)))
    records = []
    for eta, mass in ((.01, 1.), (.03, .7)):
        exact = source.copy()
        for j in range(4):
            exact = solve(a + np.sqrt(eta+j*mass**2)*np.eye(13), exact)
        exact *= np.sqrt(6.)*mass**3
        exact_response = float(-.5*(np.linalg.norm(exact[:, d > 0], 'fro')**2
                                   - np.linalg.norm(exact[:, d < 0], 'fro')**2))
        p = a.conj().T @ a
        rational = 6*mass**6*np.eye(13, dtype=complex)
        for j in range(4):
            rational = solve(p+(eta+j*mass**2)*np.eye(13), rational, assume_a='pos')
        rational_response = float(-.5*np.trace(derivative @ rational[:5, :5]).real)
        for count in (8, 13):
            y = np.eye(13, count)
            ay = a @ y
            record, output = chain(ay.conj().T @ ay, y.conj().T @ ay,
                                    derivative, np.eye(count, 5), eta, mass)
            lo, hi = record['conditional_retained_source_interval']
            exact_error = [float(np.linalg.norm((y@output-exact)[:, mask], 'fro'))
                           for mask in (d > 0, d < 0)]
            # Independent strong residual check uses A directly, not Delta.
            old = source[:count].copy()
            residual_error = 0.
            for stage in record['stages']:
                m = stage['mass']
                action = ay+m*y
                new = np.linalg.lstsq(action, y@old, rcond=None)[0]
                squared = float(np.linalg.norm(action@new-y@old, 'fro')**2)
                residual_error = max(residual_error,
                    abs(squared-stage['total_residual_squared'])/max(1., squared))
                old = new
            scale = max(1., abs(exact_response))
            records.append(dict(eta=eta, regulator_mass=mass, flux_dimension=count,
                rational_response_error=abs(rational_response-exact_response)/scale,
                complete_response_error=(abs(record['trial_filter_response']-exact_response)/scale
                                         if count == 13 else None),
                source_interval_margin=min(exact_response-lo, hi-exact_response)/scale,
                error_bound_margin=min(np.array(record['positive_negative_error_bounds'])-exact_error),
                independent_residual_relative_error=residual_error,
                gram_identity_relative_error=max(s['gram_identity_relative_error'] for s in record['stages'])))
    rejected = False
    try:
        chain(np.eye(2), 1j*np.eye(2), np.eye(1), np.eye(2, 1), 0.)
    except ValueError:
        rejected = True
    return dict(records=records, zero_eta_rejected=rejected,
        passed=bool(rejected and all(r['rational_response_error'] < 2e-10
            and (r['complete_response_error'] is None or r['complete_response_error'] < 2e-10)
            and r['source_interval_margin'] > -2e-10 and r['error_bound_margin'] > -2e-10
            and r['independent_residual_relative_error'] < 2e-10
            and r['gram_identity_relative_error'] < 2e-10 for r in records)))


def run():
    start = perf_counter()
    baseline_path = project_path('dihole_angular_crosscheck_results.json')
    state_path = project_path('dihole_angular_crosscheck_state.npz')
    previous_path = project_path('dihole_retained_resolvent_results.json')
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    previous = json.loads(previous_path.read_text(encoding='utf-8'))
    inputs = dict(baseline['input_sha256'])
    inputs.update({'dihole_angular_crosscheck.py': baseline['source_sha256'],
        state_path.name: baseline['state_sha256'], baseline_path.name: sha(baseline_path),
        'dihole_retained_resolvent.py': previous['source_sha256'], previous_path.name: sha(previous_path)})
    if not baseline['passed'] or not previous['passed'] or not all(
            sha(project_path(name)) == digest for name, digest in inputs.items()):
        raise ArithmeticError('Strict saved-form provenance failed; no calculation started.')
    retained = select_columns(baseline['metadata'], (3, 4), 3, 12)
    flux = select_columns(baseline['metadata'], (3, 4), 5, 12)
    with np.load(state_path) as archive:
        raw = {key: archive['p0_'+key] for key in ('norm', 'q', 'derivative', 'weak')}
    mv, my, derivative, _, overlap = normalized(raw, retained, flux, .7)
    iy = np.linalg.inv(np.linalg.cholesky(raw['norm'][np.ix_(flux, flux)]))
    weak = raw['weak'][0] + .7*raw['weak'][1]
    b = iy @ weak[np.ix_(flux, flux)] @ iy.conj().T
    result = dict(source_sha256=sha(__file__), input_sha256=inputs, controls=controls(),
        parameters=dict(phase=0., mu=.5, kappa=.7, regulator_mass=1., radial_modes=12,
            retained_angular_modes=3, flux_angular_modes=5, retained_dimension=len(retained),
            flux_dimension=len(flux), mass_order='ascending j=0,1,2,3'), records=[],
        scope=['The common rational PV filter is represented by one squared first-order product.',
            'MY is the full differential-action form; the weak matrix B is not squared in its place.',
            'Bounds apply to the finite-rank source V D_V V* for positive eta, with exact continuum forms.',
            'Saved numbers use unenclosed floating quadrature and arithmetic; they are conditional numerical bounds.',
            'The filter trial centre and the old r(MV) centre are different approximations; their sign difference is not physical.',
            'No omitted insertion, other sector, frequency integral, eta=0 or regulator-removal error is bounded.',
            'The large first-step residual and its small-mass amplification leave the retained sign unresolved.'])
    for eta in (.01, .0001):
        record, _ = chain(my, b, derivative, overlap, eta)
        old = next(row for row in previous['records'] if row['eta'] == eta and row['flux_dimension'] == len(flux))
        record['previous_independent_shift_interval'] = old['conditional_retained_pv_interval']
        record['previous_finite_positive_form_response'] = old['finite_pv_response']
        record['interval_radius_over_previous_finite_signal'] = max(abs(x-record['trial_filter_response'])
            for x in record['conditional_retained_source_interval']) / abs(old['finite_pv_response'])
        result['records'].append(record)
    rows = result['records']
    result['tests'] = dict(complete_and_partial_operator_controls=result['controls']['passed'],
        exact_input_provenance=all(sha(project_path(n)) == h for n, h in inputs.items()),
        retained_dimensions=len(retained) == 216 and len(flux) == 312,
        nested_overlap=max(r['overlap_isometry_error'] for r in rows) < 2e-10,
        measured_weak_skew=max(r['normalized_weak_skew_relative_error'] for r in rows) < 2e-8,
        sampled_complement_positive=min(r['delta_minimum_eigenvalue'] for r in rows) >= -2e-9,
        residual_gram_identity=max(s['gram_identity_relative_error'] for r in rows for s in r['stages']) < 2e-10,
        least_squares_equations=max(s['least_squares_relative_error'] for r in rows for s in r['stages']) < 2e-10)
    result.update(passed=all(result['tests'].values()), seconds=perf_counter()-start,
        retained_sign_certified=False, physical_source_computed=False)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], seconds=result['seconds'],
        records=[{key:r[key] for key in ('eta', 'trial_filter_response', 'conditional_retained_source_interval',
                                        'interval_radius_over_previous_finite_signal')} for r in rows]), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Common-filter control failed; diagnostics saved.')


if __name__ == '__main__':
    run()
