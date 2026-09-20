"""Independent bridge fields coupled to the globally sewn physical form.

One maximal form yields nested variational spaces by principal submatrices.
Matrices are assembled from pointwise differential actions in small batches;
the unchanged exterior is added only to its global block. Frequency-tail
bounds concern the saved finite matrices, not the physical continuum source.
"""
import gc
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson

from dihole_enriched_bridge import bridge_data, column_metadata
from dihole_euclidean_form import assemble, relative_error
from dihole_finite_frequency_tail import coefficient_bounds, high_tail
from dihole_sewn_form import exterior_data, fibre_response
from dihole_stretched_collar import gauss
from project_paths import project_path


def raw_polynomial(data, batch=256):
    """Exact kappa-polynomial algebra, numerical pointwise Gram integration."""
    n = data['value'].shape[-1]
    result = dict(norm=np.zeros((n, n), complex), q=np.zeros((3, n, n), complex),
                  derivative=np.zeros((3, n, n), complex), weak=np.zeros((2, n, n), complex))
    sigma = np.array([1., -1., -1., 1.])[None, :, None]
    for begin in range(0, len(data['P']), batch):
        part = {k: v[begin:begin+batch] for k, v in data.items()}
        v, radial, angular = [part[k] for k in ('value', 'radial', 'angular')]
        p, w, f, fp = [part[k][:, None, None] for k in ('P', 'W', 'f', 'fp')]
        weight = np.sqrt(part['weight'])[:, None, None]
        x = weight*(radial+angular)/np.sqrt(p)
        y = weight*1j*w*v/np.sqrt(p)
        z = weight*np.sqrt(p)*v
        beta_x, beta_y = [a[:, [2, 3, 0, 1], :].reshape(-1, n) for a in (x, y)]
        x, y, z = [a.reshape(-1, n) for a in (x, y, z)]
        xy = x.conj().T@y
        result['q'] += np.array([x.conj().T@x, xy+xy.conj().T, y.conj().T@y])
        result['norm'] += z.conj().T@z
        result['weak'] += np.array([z.conj().T@beta_x, z.conj().T@beta_y])
        if np.any(f) or np.any(fp):
            dx = (weight*(f*radial-.5j*fp*sigma*v)/np.sqrt(p)).reshape(-1, n)
            dy = (-weight*1j*f*w*v/np.sqrt(p)).reshape(-1, n)
            for index, cross in enumerate((dx.conj().T@x, dx.conj().T@y+dy.conj().T@x,
                                            dy.conj().T@y)):
                result['derivative'][index] += cross+cross.conj().T
    for key in ('q', 'derivative'):
        result[key] = (result[key]+result[key].conj().transpose(0, 2, 1))/2
    result['norm'] = (result['norm']+result['norm'].conj().T)/2
    return result


def assembled_batches(data, kappa, epsilon=0., mass=0., batch=256):
    """Reuse137's independent single-frequency implementation as a control."""
    result = None
    for begin in range(0, len(data['P']), batch):
        row = assemble({k: v[begin:begin+batch] for k, v in data.items()}, kappa, epsilon, mass)
        if result is None:
            result = row
        else:
            for key in result:
                result[key] += row[key]
    return result


def select_columns(metadata, degrees, angular_modes=0, radial_modes=0):
    first = max(1, abs(metadata['a']))
    return np.array([r['index'] for r in metadata['columns']
        if ((r['family'] == 'global' and r['rho_degree'] < degrees[0] and r['z_degree'] < degrees[1])
            or (r['family'] == 'compact' and angular_modes > 0 and r['radial_mode'] <= radial_modes
                and (r['l'] == 0 or first <= r['l'] < first+angular_modes)))], int)


def canonical(raw, indices):
    norm = raw['norm'][np.ix_(indices, indices)]
    inverse = np.linalg.inv(np.linalg.cholesky(norm))
    result = {}
    for key in ('q', 'derivative'):
        result[key] = []
        for coefficient in raw[key]:
            matrix = inverse@coefficient[np.ix_(indices, indices)]@inverse.conj().T
            result[key].append((matrix+matrix.conj().T)/2)
    return result


def full_frequency(polynomial, tail_target=1e-7):
    """Numerical middle integral plus outward tails of the stored matrix model."""
    bounds = coefficient_bounds(polynomial)
    endpoint = 10.
    for _ in range(30):
        try:
            tail = high_tail(bounds, endpoint)
        except ValueError:
            tail = np.inf
        if tail <= tail_target:
            break
        endpoint *= 2
    else:
        raise ArithmeticError('Finite-matrix tail target not reached.')
    cuts = [0., .01, .1, 1., 10.]
    while cuts[-1] < endpoint:
        cuts.append(min(10*cuts[-1], endpoint))
    values = []
    for order in (16, 24):
        x, w = gauss(order)
        total = 0.
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            points = (hi+lo)/2+(hi-lo)*x/2
            total += (hi-lo)/2*sum(weight*fibre_response(polynomial, k) for weight, k in zip(w, points))
        values.append(float(total/np.pi))
    return dict(endpoint=endpoint, certified_matrix_tail=tail, coefficient_bounds=bounds,
                integral_order16=values[0], integral_order24=values[1],
                middle_quadrature_change=abs(values[1]-values[0]),
                zero_frequency_value=fibre_response(polynomial, 0.),
                scope='Outward tail for exact dyadic canonical matrices only; middle integration and physical projection/assembly error are not certified.')


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout='constrained')
    labels = ['24 exterior', '+6 bridge', '+20 bridge', '+42 bridge', '48 ext +42 bridge']
    wanted = ['E24', 'E24_B1x1', 'E24_B2x2', 'E24_B3x3', 'E48_B3x3']
    for phase, label in ((0., 'Phase 0'), (float(np.pi), 'Phase pi')):
        records = [r for r in result['records'] if r['mu'] == .5 and r['phase'] == phase]
        values = [next(r['integral_33'] for r in records if r['name'] == key) for key in wanted]
        axes[0].plot(range(len(labels)), values, marker='o', label=label)
    axes[0].set_xticks(range(len(labels)), labels, rotation=25, ha='right')
    axes[0].set(ylabel='Finite-window response (hbar c / D)', title='Independent interior fields change the model')
    axes[0].legend()
    axes[0].axhline(0, color='grey', lw=.7)
    central = [r for r in result['records'] if r['mu'] == .5 and r['name'] == 'E48_B3x3']
    for r in central:
        axes[1].plot(r['kappas'], r['responses'], label='Phase '+('0' if r['phase'] == 0 else 'pi'))
    axes[1].set(xscale='log', xlabel='Frequency kappa D', ylabel='PV fibre response', title='Largest tested central-sector space')
    axes[1].axhline(0, color='grey', lw=.7)
    axes[1].legend()
    fig.suptitle('Bridge fields can now vary independently; physical supporting stress is still unknown', fontsize=12)
    fig.supxlabel('Finite common physical form, MD=1. Field-space and continuum errors remain unbounded.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_bridge_enrichment_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    baseline = json.loads(project_path('dihole_sewn_form_results.json').read_text())
    records, controls, matrices = [], [], {}
    kappa = np.geomspace(.01, 10., 33)
    for sector, mu in enumerate((-.5, .5, 1.5)):
        exterior = exterior_data(mu, order=24, degrees=(3, 4))
        exterior_raw = raw_polynomial(exterior)
        del exterior
        metadata = column_metadata(mu, 3, 3, (3, 4))
        count = metadata['global_columns']
        for phase_index, phase in enumerate((0., np.pi)):
            data = bridge_data(mu, phase, angular_modes=3, radial_modes=3, order=48, degrees=(3, 4))
            raw = raw_polynomial(data)
            if mu == .5 and phase == 0:
                check = assembled_batches(data, .7, mass=1.3)
                plus, minus = [assembled_batches(data, .7, epsilon=e) for e in (1e-5, -1e-5)]
                controls.append(dict(name='pointwise_polynomial',
                    q_error=relative_error(sum(.7**j*raw['q'][j] for j in range(3)), check['q']),
                    derivative_error=relative_error(sum(.7**j*raw['derivative'][j] for j in range(3)), check['derivative']),
                    derivative_fd_error=relative_error(check['derivative'], (plus['q']-minus['q'])/2e-5)))
            del data
            gc.collect()
            for key in raw:
                raw[key][..., :count, :count] += exterior_raw[key]
            weak = raw['weak'][0]+.7*raw['weak'][1]
            controls.append(dict(name='global_sewing', mu=mu, phase=phase,
                                 skew_error=relative_error(weak, -weak.conj().T)))
            for key in raw:
                matrices[f's{sector}_p{phase_index}_{key}'] = raw[key]
            cases = [('E24', (2, 3), 0, 0), ('E48', (3, 4), 0, 0),
                     ('E24_B1x1', (2, 3), 1, 1), ('E24_B2x2', (2, 3), 2, 2),
                     ('E24_B2x3', (2, 3), 2, 3), ('E24_B3x2', (2, 3), 3, 2),
                     ('E24_B3x3', (2, 3), 3, 3), ('E48_B3x3', (3, 4), 3, 3)]
            for name, degrees, angular, radial in cases:
                indices = select_columns(metadata, degrees, angular, radial)
                polynomial = canonical(raw, indices)
                response = [fibre_response(polynomial, k) for k in kappa]
                row = dict(mu=mu, phase=phase, name=name, dimension=len(indices), indices=indices.tolist(),
                           kappas=kappa.tolist(), responses=response,
                           integral_33=float(simpson(response, x=kappa)/np.pi),
                           integral_17=float(simpson(response[::2], x=kappa[::2])/np.pi))
                if name == 'E24':
                    old = next(r for r in baseline['records'] if r['mu'] == mu and r['phase'] == phase)
                    row['baseline_relative_error'] = relative_error(np.array(response), np.array(old['responses']))
                if name == 'E48_B3x3':
                    # Preserve the exact canonical dyadic inputs enclosed by
                    # the frequency-tail helper, as well as the raw forms.
                    for key in polynomial:
                        matrices[f's{sector}_p{phase_index}_canonical_{key}'] = np.array(polynomial[key])
                    row['full_frequency'] = full_frequency(polynomial)
                records.append(row)
            print(json.dumps(dict(assembled_mu=mu, phase=phase, dimension=metadata['total_columns'],
                                  seconds=perf_counter()-start)), flush=True)
    # Stronger quadrature for the entire coupled central matrix, not merely
    # its separate diagonal blocks. No physical-domain remainder follows.
    data = bridge_data(.5, 0., 3, 3, order=64, degrees=(3, 4))
    refined = raw_polynomial(data)
    del data
    exterior = exterior_data(.5, order=32, degrees=(3, 4))
    exterior_raw = raw_polynomial(exterior)
    del exterior
    for key in refined:
        refined[key][..., :48, :48] += exterior_raw[key]
    controls.append(dict(name='coupled_quadrature', **{key+'_change': relative_error(
        refined[key], matrices['s1_p0_'+key]) for key in ('q', 'norm', 'derivative')}))
    indices = np.arange(len(refined['norm']))
    refined_response = fibre_response(canonical(refined, indices), .7)
    original = {key: matrices['s1_p0_'+key] for key in refined}
    controls[-1]['fibre_change'] = abs(refined_response-fibre_response(canonical(original, indices), .7))
    polynomial_control = next(r for r in controls if r['name'] == 'pointwise_polynomial')
    refinement = controls[-1]
    tests = dict(actual_polynomial=polynomial_control['q_error'] < 1e-12 and polynomial_control['derivative_error'] < 1e-12,
        metric_variation=polynomial_control['derivative_fd_error'] < 2e-8,
        global_sewing=max(r['skew_error'] for r in controls if r['name'] == 'global_sewing') < 2e-9,
        frozen_global_recovery=max(r['baseline_relative_error'] for r in records if r['name'] == 'E24') < 2e-8,
        coupled_quadrature=max(refinement[k+'_change'] for k in ('q', 'norm', 'derivative')) < 2e-7,
        finite_matrix_tail=all(r['full_frequency']['certified_matrix_tail'] <= 1e-7 for r in records if 'full_frequency' in r),
        finite_responses=all(np.all(np.isfinite(r['responses'])) for r in records))
    np.savez_compressed(project_path('dihole_bridge_enrichment_state.npz'), **matrices)
    names = ('dihole_enriched_bridge.py', 'dihole_finite_frequency_tail.py', 'dihole_euclidean_form.py',
             'dihole_sewn_form.py', 'dihole_sewn_quadrature.py', 'dihole_azimuthal_bridge.py',
             'dihole_interval_arithmetic.py', 'dihole_stretched_collar.py', 'balanced_dihole.py')
    result = dict(complete=True, passed=all(tests.values()), tests=tests, controls=controls, records=records,
        seconds=perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
        baseline_sha256=hashlib.sha256(project_path('dihole_sewn_form_results.json').read_bytes()).hexdigest(),
        state_sha256=hashlib.sha256(project_path('dihole_bridge_enrichment_state.npz').read_bytes()).hexdigest(),
        scope=['Independent zero-trace bridge modes and exact global trace lifts assembled together, including cross terms.',
               'Fixed exterior support and axis exclusion, three azimuthal sectors, finite angular/radial trial dimensions.',
               'Outward frequency tails cover exact saved canonical finite-matrix coefficients, not assembly or continuum errors.',
               'No physical angular/exterior/IR remainder, covariant subtraction or regulator-removal limit is supplied.',
               'No physical source sign, obtainable supporting matter, stability, new law or backward passage established.'])
    project_path('dihole_bridge_enrichment_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, seconds=result['seconds']), indent=2))
    if not result['passed']:
        raise RuntimeError('Bridge enrichment checks failed; evidence retained.')
    plot(result)
    return result


if __name__ == '__main__':
    run()
