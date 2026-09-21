"""Positive-shift feedback bounds from the saved, globally sewn Dirac forms.

The analytic bracket includes the exterior for a finite-rank insertion only.
Numbers use floating matrix integrals, not certified continuum enclosures.
No geometry assembly or old simulation is rerun.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh, solve

from dihole_bridge_enrichment import select_columns
from project_paths import project_path


NAME = 'dihole_retained_resolvent'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hermitian(matrix):
    return (matrix + matrix.conj().T) / 2


def inverse(matrix, shift):
    if not np.isfinite(shift) or shift <= 0:
        raise ValueError('The retained bracket requires a strictly positive shift.')
    return solve(matrix + shift * np.eye(len(matrix)), np.eye(len(matrix)), assume_a='pos')


def bracket(mv, my, cross, shift):
    lower = inverse(mv, shift)
    upper = hermitian((np.eye(len(mv)) - cross.conj().T @ inverse(my, shift) @ cross) / shift)
    return lower, upper


def source_parts(derivative):
    values, vectors = eigh(derivative)
    return tuple((vectors * np.maximum(sign * values, 0)) @ vectors.conj().T for sign in (1, -1))


def trace_error(derivative, gap):
    positive, negative = source_parts(derivative)
    return [-float(np.trace(negative @ gap).real), float(np.trace(positive @ gap).real)]


def controls():
    """Compare with a complete noncommuting skew operator, including its nullspace."""
    rng = np.random.default_rng(155)
    random = rng.normal(size=(13, 13)) + 1j * rng.normal(size=(13, 13))
    vectors = np.linalg.qr(random)[0]
    a = 1j * (vectors * np.arange(-6., 7.)) @ vectors.conj().T
    p = hermitian(a.conj().T @ a)
    derivative = hermitian(random[:5, :5])
    errors, margins = [], []
    for shift in (1e-3, 1., 3.):
        exact = inverse(p, shift)[:5, :5]
        for count in (5, 9, 13):
            lower, upper = bracket(p[:5, :5], p[:count, :count], a[:count, :5], shift)
            scale = max(1., np.linalg.norm(upper, 2))
            margins.extend(float(eigh(hermitian(x), eigvals_only=True)[0]) / scale
                           for x in (exact - lower, upper - exact))
            lo, hi = trace_error(derivative, upper - lower)
            truth = float(np.trace(derivative @ (exact - lower)).real)
            margins.extend(((truth - lo) / scale, (hi - truth) / scale))
            if count == 13:
                errors.append(float(np.linalg.norm(upper - exact, 2)) / scale)
    rejected_zero = False
    try:
        inverse(p, 0.)
    except ValueError:
        rejected_zero = True
    return dict(minimum_scaled_margin=min(margins), complete_flux_relative_error=max(errors),
                zero_shift_rejected=rejected_zero,
                passed=bool(min(margins) > -2e-10 and max(errors) < 2e-10 and rejected_zero))


def normalized(raw, retained, flux, kappa):
    iv = np.linalg.inv(np.linalg.cholesky(raw['norm'][np.ix_(retained, retained)]))
    iy = np.linalg.inv(np.linalg.cholesky(raw['norm'][np.ix_(flux, flux)]))
    q, derivative = [sum(kappa**j * c for j, c in enumerate(raw[key])) for key in ('q', 'derivative')]
    weak = raw['weak'][0] + kappa * raw['weak'][1]
    congruence = lambda x, i, j, l, r: l @ x[np.ix_(i, j)] @ r.conj().T
    mv = hermitian(congruence(q, retained, retained, iv, iv))
    my = hermitian(congruence(q, flux, flux, iy, iy))
    dv = hermitian(congruence(derivative, retained, retained, iv, iv))
    cross = congruence(weak, flux, retained, iy, iv)
    overlap = congruence(raw['norm'], flux, retained, iy, iv)
    return mv, my, dv, cross, overlap


def evaluate(raw, retained, flux, eta, kappa=.7):
    mv, my, derivative, cross, overlap = normalized(raw, retained, flux, kappa)
    values, vectors = eigh(mv)
    inserted = np.einsum('ij,ij->j', vectors.conj(), derivative @ vectors).real
    # Keep the common rational form for the central PV value; do not subtract large traces.
    central = -3 * float(inserted @ (1 / np.prod(values[:, None] + eta + np.arange(4), axis=1)))
    low, high, direct, rows = 0., 0., 0., []
    for j, coefficient in enumerate((-.5, 1.5, -1.5, .5)):
        shift = eta + j
        lower, upper = bracket(mv, my, cross, shift)
        gap = hermitian(upper - lower)
        err = trace_error(derivative, gap)
        low += min(coefficient * err[0], coefficient * err[1])
        high += max(coefficient * err[0], coefficient * err[1])
        direct += coefficient * float(np.trace(derivative @ lower).real)
        # A larger Galerkin solve supplies a separate lower bracket for the same loads V.
        enlarged_lower = hermitian(overlap.conj().T @ inverse(my, shift) @ overlap)
        scale = max(1., np.linalg.norm(upper, 2))
        rows.append(dict(shift=shift, trace_error=err,
            scaled_order_margin=min(float(eigh(hermitian(x), eigvals_only=True)[0]) / scale
                for x in (gap, enlarged_lower - lower, upper - enlarged_lower))))
    interval = [central + low, central + high]
    return dict(kappa=kappa, eta=eta, proper_mass=1., retained_dimension=len(retained),
        flux_dimension=len(flux), finite_pv_response=central,
        conditional_retained_pv_interval=interval,
        error_radius_over_signal=max(abs(low), abs(high)) / abs(central),
        sign_resolved=bool(interval[0] > 0 or interval[1] < 0),
        rational_identity_error=abs(direct - central),
        overlap_isometry_error=float(np.linalg.norm(overlap.conj().T @ overlap - np.eye(len(retained)), 2)),
        shifts=rows)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3), layout='constrained')
    for eta, color in ((.01, '#156487'), (.0001, '#bc6627')):
        rows = [r for r in result['records'] if r['eta'] == eta]
        axes[0].semilogy([r['flux_dimension'] for r in rows],
                         [r['error_radius_over_signal'] for r in rows], 'o-', color=color,
                         label=f'Positive infrared shift = {eta:g}')
    axes[0].axhline(1., color='gray', ls='--', label='Error equals the finite signal')
    axes[0].set(xlabel='Fields used to estimate feedback', ylabel='Error radius / |finite response|',
                title='This bound is too wide to decide the sign', xticks=[216, 264, 312])
    axes[0].legend(fontsize=8)
    axes[1].axis('off')
    text = ('Saved bridge + exterior fields\n216 retained patterns, one angular sector\n\n'
            'What the calculation adds\nA bound on feedback from the full operator\ninto this finite-rank source insertion\n\n'
            'What is still missing\nThe omitted source insertion and sectors;\ncontrolled integration and infrared errors\n\n'
            'Outcome\nNo supporting stress or time passage established')
    axes[1].text(.03, .96, text, va='top', fontsize=10, linespacing=1.4)
    fig.suptitle('A useful error test rejects a shortcut to the physical source', fontsize=12)
    fig.supxlabel('Conditional floating-matrix bounds at phase 0, kappa D = 0.7, MD = 1; not a continuum certificate.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME + '_progress.' + ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    baseline_path = project_path('dihole_angular_crosscheck_results.json')
    state_path = project_path('dihole_angular_crosscheck_state.npz')
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    inputs = dict(baseline['input_sha256'])
    inputs.update({'dihole_angular_crosscheck.py': baseline['source_sha256'],
                   state_path.name: baseline['state_sha256'], baseline_path.name: sha(baseline_path)})
    if not baseline['passed'] or not all(sha(project_path(n)) == h for n, h in inputs.items()):
        raise ArithmeticError('Saved-form provenance mismatch; no calculation started.')
    result = dict(source_sha256=sha(__file__), input_sha256=inputs, controls=controls(),
        records=[], parameters=dict(phase=0., mu=.5, radial_modes=12, retained_angular_modes=3,
                                    flux_angular_modes=[3, 4, 5]),
        scope=['The analytic theorem brackets retained resolvent matrix elements for sigma>0 on the full sewn operator.',
               'Displayed numbers are conditional floating-matrix evaluations; geometry quadrature is not enclosed.',
               'PV intervals bound only the finite-rank insertion V D_V V*, not the full physical metric variation.',
               'Individual-shift bounds lose PV cancellations and deteriorate as eta decreases.',
               'One phase, one sector and one frequency do not determine an integrated or renormalized source.'])
    retained = select_columns(baseline['metadata'], (3, 4), 3, 12)
    with np.load(state_path) as archive:
        raw = {k: archive['p0_' + k] for k in ('norm', 'q', 'derivative', 'weak')}
        refined = {k: archive['q80_' + k] for k in raw}
        for eta in (.01, .0001):
            for angular in (3, 4, 5):
                flux = select_columns(baseline['metadata'], (3, 4), angular, 12)
                result['records'].append(evaluate(raw, retained, flux, eta))
            control = evaluate(refined, retained, flux, eta)
            original = result['records'][-1]
            control['interval_relative_change'] = float(np.max(abs(np.array(control['conditional_retained_pv_interval'])
                - original['conditional_retained_pv_interval'])) / max(abs(x) for x in original['conditional_retained_pv_interval']))
            result.setdefault('quadrature_controls', []).append(control)
    rows = result['records'] + result['quadrature_controls']
    result['tests'] = dict(
        complete_operator_control=result['controls']['passed'],
        exact_input_provenance=all(sha(project_path(n)) == h for n, h in inputs.items()),
        nested_overlap=max(r['overlap_isometry_error'] for r in rows) < 2e-10,
        rational_pv_identity=max(r['rational_identity_error'] for r in rows) < 2e-10,
        sampled_matrix_order=min(s['scaled_order_margin'] for r in rows for s in r['shifts']) > -2e-9,
        flux_enlargement_tightens=all(group[i+1]['conditional_retained_pv_interval'][0] >= group[i]['conditional_retained_pv_interval'][0]
            and group[i+1]['conditional_retained_pv_interval'][1] <= group[i]['conditional_retained_pv_interval'][1]
            for eta in (.01, .0001) for group in [[r for r in result['records'] if r['eta'] == eta]] for i in (0, 1)),
        quadrature_sensitivity=max(r['interval_relative_change'] for r in result['quadrature_controls']) < 2e-7)
    result.update(passed=all(result['tests'].values()), seconds=perf_counter() - start,
                  retained_sign_certified=False, physical_source_computed=False)
    project_path(NAME + '_results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], records=[{
        k: r[k] for k in ('eta', 'flux_dimension', 'finite_pv_response', 'conditional_retained_pv_interval', 'error_radius_over_signal')}
        for r in result['records']], seconds=result['seconds']), indent=2), flush=True)
    if not result['passed']:
        raise ArithmeticError('Retained-feedback check failed; diagnostics were saved.')
    plot(result)


if __name__ == '__main__':
    run()
