"""Proper heat response and a combined saved angular/exterior trial space.

Uses the squared differential action on the sewn geometry, both chiralities,
and the exact frequency polynomial. Finite-space comparisons are sensitivity
tests, not bounds on the physical continuum or a renormalized source.
"""
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh

import dihole_interval_arithmetic as ia
from dihole_exterior_feedback import append_exterior, chiral_blocks
from dihole_finite_frequency_tail import coefficient_bounds, high_tail, PI
from dihole_stretched_collar import gauss
from project_paths import project_path


NAME = 'dihole_proper_heat'
KEYS = ('norm', 'q', 'derivative', 'weak')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def heat_tail(bounds, cutoff, heat_time):
    """Outward Gaussian tail of the exact supplied dyadic finite matrices.

    e^-x <= 2^-floor(x) avoids an unenclosed floating exponential. The
    elementary Gaussian-integral majorants include positive kappa and 1/pi.
    """
    if cutoff <= 0 or heat_time <= 0 or bounds['a'] <= 0:
        raise ValueError('Positive cutoff, heat time and quadratic coefficient required.')
    p = ia.point
    k, s = p(float(cutoff)), p(float(heat_time))
    chi = ia.sub(ia.sub(p(bounds['a']), ia.div(p(bounds['b']), k)),
                 ia.div(p(bounds['c']), ia.square(k)))
    if chi[0] <= 0:
        raise ValueError('No positive high-frequency quadratic bound at this cutoff.')
    alpha = ia.mul(s, chi)
    exponent = float(ia.mul(alpha, ia.square(k))[0])
    exponential = p(math.ldexp(1., -min(1074, math.floor(exponent))))
    two_alpha = ia.mul(p(2.), alpha)
    terms = ia.add(ia.add(ia.div(p(bounds['d0']), k), p(bounds['d1'])),
        ia.mul(p(bounds['d2']), ia.add(k, ia.div(p(1.), ia.mul(two_alpha, k)))))
    value = ia.div(ia.mul(ia.mul(ia.mul(s, p(float(bounds['n']))), exponential), terms),
                   ia.mul(PI, two_alpha))
    return float(value[1])


def uv_coefficient(blocks):
    total = 0.
    for block in blocks:
        values, vectors = eigh(block['q'][2])
        if values[0] <= 0:
            raise ArithmeticError('The temporal coefficient must be positive.')
        diagonal = np.einsum('ij,ij->j', vectors.conj(), block['derivative'][2] @ vectors).real
        total -= float(diagonal @ values**-1.5) / (4 * np.sqrt(np.pi))
    return total


def spectrum(blocks, nodes):
    eigenvalues, inserted = [], []
    for k in nodes:
        es, ds = [], []
        for block in blocks:
            q, dq = [sum(k**j * c for j, c in enumerate(block[key])) for key in ('q', 'derivative')]
            values, vectors = eigh(q, check_finite=False)
            if values[0] <= 0:
                raise ArithmeticError('Nonpositive sampled finite form; do not assume a gap.')
            es.extend(values)
            ds.extend(np.einsum('ij,ij->j', vectors.conj(), dq @ vectors).real)
        eigenvalues.append(es)
        inserted.append(ds)
    return np.asarray(eigenvalues), np.asarray(inserted)


def integrate(blocks, heat_times, cutoff, order):
    cuts = [0., .01, .1, 1.]
    while cuts[-1] < cutoff:
        cuts.append(min(cutoff, 10 * cuts[-1]))
    x, w = gauss(order)
    nodes = np.concatenate([(lo + hi) / 2 + (hi - lo) * x / 2 for lo, hi in zip(cuts[:-1], cuts[1:])])
    weights = np.concatenate([(hi - lo) * w / 2 for lo, hi in zip(cuts[:-1], cuts[1:])])
    eigenvalues, inserted = spectrum(blocks, nodes)
    samples = dict(nodes=nodes, weights=weights, eigenvalues=eigenvalues, inserted=inserted)
    heat, pv = sampled_integrals(samples, heat_times)
    return heat, pv, samples


def sampled_integrals(samples, heat_times):
    weights, eigenvalues, inserted = [samples[k] for k in ('weights', 'eigenvalues', 'inserted')]
    heat = np.array([-s * np.einsum('kn,kn,k->', np.exp(-s * eigenvalues), inserted, weights) / np.pi
                     for s in heat_times])
    rational = 6 / np.prod(eigenvalues[..., None] + np.arange(4), axis=-1)
    pv = -.5 * float(np.einsum('kn,kn,k->', rational, inserted, weights)) / np.pi
    return heat, pv


def benchmark():
    """Exact noncommuting Gaussian pencil and an independent trace derivative."""
    sx = np.array([[0., 1.], [1., 0.]])
    sy = np.array([[0., -1j], [1j, 0.]])
    sz = np.diag([1., -1.])
    eye = np.eye(2)
    a = np.array([2.6, 1.4])
    block = dict(q=np.array([2*eye-.6*sz, 1.2*sx, 2*eye+.6*sz]),
                 derivative=np.array([eye+sy, np.zeros((2, 2)), eye]))
    times = np.array([.001, .03, 1.])
    values, _, _ = integrate([block], times, 256., 40)
    exact = np.array([-np.sum(np.exp(-s*a) * (s/np.sqrt(a) + .5/a**1.5))
                      / (2*np.sqrt(np.pi*s)) for s in times])
    bounds = coefficient_bounds(block)
    from scipy.integrate import quad
    s, cutoff = .03, 16.
    actual = quad(lambda k: s * (1+k*k) * np.sum(np.exp(-s*a*(1+k*k))) / np.pi,
                  cutoff, np.inf, epsabs=1e-16)[0]
    q = np.array([[2., .3+.2j], [.3-.2j, 4.]])
    dq = np.array([[.7, -.4j], [.4j, -.2]])
    ev, vec = eigh(q)
    trace_derivative = -.3 * np.sum(np.exp(-.3*ev)*np.einsum('ij,ij->j',vec.conj(),dq@vec).real)
    eps = 1e-5
    fd = (np.exp(-.3*eigh(q+eps*dq, eigvals_only=True)).sum()
          - np.exp(-.3*eigh(q-eps*dq, eigvals_only=True)).sum())/(2*eps)
    return dict(gaussian_maximum_error=float(np.max(abs(values-exact))),
                noncommuting_derivative_error=float(abs(fd-trace_derivative)),
                sampled_tail=actual, outward_tail=heat_tail(bounds, cutoff, s))


def load_models():
    baselines, states, inputs = {}, {}, {}
    for name in ('dihole_radial_refinement', 'dihole_angular_crosscheck', 'dihole_exterior_feedback'):
        path = project_path(name+'_results.json')
        baseline = json.loads(path.read_text(encoding='utf-8'))
        if not baseline['passed']:
            raise ArithmeticError('A saved calculation did not pass.')
        checks = dict(baseline['input_sha256'])
        checks.update({name+'.py': baseline['source_sha256'], name+'_state.npz': baseline['state_sha256']})
        if not all(sha(project_path(n)) == h for n, h in checks.items()):
            raise ArithmeticError('Saved-input authentication failed: '+name)
        inputs.update(checks)
        inputs[path.name] = sha(path)
        baselines[name] = baseline
        with np.load(project_path(name+'_state.npz')) as archive:
            states[name] = {key: archive[key] for key in archive.files}
    radial, angular, exterior = [states[n] for n in ('dihole_radial_refinement', 'dihole_angular_crosscheck', 'dihole_exterior_feedback')]
    blocks = lambda state, prefix: [{key: state[f'{prefix}_b{b}_{key}'] for key in ('q', 'derivative')} for b in (0, 1)]
    models = {f'r{n}_a3_e48': blocks(radial, f'p0_n{n}') for n in (8, 10, 12)}
    models.update({f'r12_a{n}_e48': blocks(angular, f'p0_a{n}') for n in (4, 5)})
    models['r12_a3_e80'] = blocks(exterior, 'p0')
    ids = np.r_[np.arange(48), np.arange(216, 248)]
    extra = {key: exterior['p0_'+key][..., ids, :][..., ids].copy() for key in KEYS}
    for key in KEYS:
        extra[key][..., :48, :48] = 0
    if any(np.any(exterior['p0_'+key][..., 48:216, 216:]) for key in KEYS):
        raise ArithmeticError('The supposedly disjoint compact/exterior blocks are nonzero.')
    columns = baselines['dihole_angular_crosscheck']['metadata']['columns']
    new = baselines['dihole_exterior_feedback']['new_columns']
    for prefix, label in (('p0_', 'r12_a5_e80'), ('q80_', 'q80_r12_a5_e80')):
        raw = {key: angular[prefix+key] for key in KEYS}
        combined = append_exterior(raw, extra)
        if not all(np.array_equal(combined[key][..., :312, :312], raw[key]) for key in KEYS):
            raise ArithmeticError('The old angular space was modified during combination.')
        models[label], paired = chiral_blocks(combined, columns, new)
        if not paired:
            raise ArithmeticError('Unexpected loss of chiral pairing in the combined space.')
    return models, baselines, inputs


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    times = np.array(result['heat_times'])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    records = {r['model']: r for r in result['records']}
    for name, label in (('r8_a3_e48', 'Radial 8, angular 3'), ('r12_a3_e48', 'Radial 12, angular 3'),
                        ('r12_a5_e48', 'Radial 12, angular 5'), ('r12_a5_e80', 'Angular 5 + larger exterior')):
        axes[0].semilogx(times, np.sqrt(times)*np.array(records[name]['heat']), label=label)
    axes[0].set(xlabel='Heat parameter s / D squared', ylabel='sqrt(s) times finite heat variation',
                title='Short scales still depend on omitted fields')
    axes[0].legend(fontsize=8)
    for row in result['comparisons']:
        if row['name'] == 'earlier_radial':
            continue
        axes[1].loglog(times, np.maximum(row['relative_change'], 1e-16), label=row['name'].replace('_', ' '))
    axes[1].axhline(.01, color='gray', ls='--', label='1% sensitivity threshold')
    axes[1].set(xlabel='Heat parameter s / D squared', ylabel='Change / |larger-space response|',
                title='Refinement sensitivity is not a continuum bound')
    axes[1].legend(fontsize=8)
    fig.suptitle('Quantum source: combined angular and exterior fields, one azimuthal sector', fontsize=12)
    fig.supxlabel('Proper heat response on the prescribed sewn geometry. No renormalized support or time passage established.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run(recover=False):
    start = perf_counter()
    models, baselines, inputs = load_models()
    for name in ('dihole_proper_heat.py', 'dihole_exterior_feedback.py', 'dihole_finite_frequency_tail.py',
                 'dihole_interval_arithmetic.py', 'dihole_stretched_collar.py', 'project_paths.py'):
        inputs[name] = sha(project_path(name))
    cache = None
    if recover:
        manifest_path = project_path(NAME+'_serialization_failure.json')
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        source_path, state_path = [project_path(manifest[k]) for k in ('source_archive', 'state_file')]
        if sha(source_path) != manifest['source_sha256'] or sha(state_path) != manifest['state_sha256']:
            raise ArithmeticError('Failed-run source or saved spectra do not match their provenance.')
        for path in (manifest_path, source_path, state_path):
            inputs[path.relative_to(project_path('README.md').parent).as_posix()] = sha(path)
        with np.load(state_path) as archive:
            cache = {key: archive[key] for key in archive.files}
    times = np.geomspace(1e-6, 10., 57)
    result = dict(complete=False, passed=False, source_sha256=sha(__file__), input_sha256=inputs,
        heat_times=times.tolist(), benchmark=benchmark(), records=[], comparisons=[],
        scope=['One azimuthal sector mu=1/2, phase0, both chiralities; physical proper Euclidean form.',
               'The new344-column model combines authenticated angular and disjoint exterior extensions without reassembly.',
               'Gaussian and rational frequency-tail bounds concern exact saved dyadic canonical matrices only.',
               'Middle frequency quadrature, assembly and infinite field/sector complements are not certified.',
               'q80 combined control holds the appended exterior blocks fixed; their higher-order matrices were not saved.',
               'Fixed finite-space ultraviolet behavior cannot supply the continuum heat coefficients or renormalized source.'])
    saved = {}
    try:
        for label, blocks in models.items():
            # Exact +/- coefficient pairing gives identical norm bounds, but both spectra are evaluated.
            if not all(np.array_equal(blocks[1][key][j], (-1)**j*blocks[0][key][j])
                       for key in ('q', 'derivative') for j in range(3)):
                raise ArithmeticError('Full even-frequency trace was not established.')
            bounds = coefficient_bounds(blocks[0])
            cutoff = 1.
            for _ in range(30):
                try:
                    tails = [2*heat_tail(bounds, cutoff, s) for s in times]
                    pv_tail = 2*high_tail(bounds, cutoff)
                    if max(tails) < 1e-10 and pv_tail < 1e-10:
                        break
                except ValueError:
                    pass
                cutoff *= 2
            else:
                raise ArithmeticError('Finite-matrix frequency tails did not reach their target.')
            low, pv_low, _ = integrate(blocks, times, cutoff, 24)
            if cache is None:
                high, pv_high, spectra = integrate(blocks, times, cutoff, 40)
            else:
                for bi, block in enumerate(blocks):
                    for key, value in block.items():
                        if not np.array_equal(value, cache[f'{label}_b{bi}_{key}']):
                            raise ArithmeticError('Recovered spectra do not use these exact canonical matrices.')
                spectra = {key: cache[label+'_'+key] for key in ('nodes','weights','eigenvalues','inserted')}
                if not 0 < spectra['nodes'][-1] < cutoff or len(spectra['nodes']) % 40:
                    raise ArithmeticError('Unexpected recovered quadrature grid.')
                high, pv_high = sampled_integrals(spectra, times)
            for key, value in spectra.items():
                saved[label+'_'+key] = value
            for bi, block in enumerate(blocks):
                for key, value in block.items():
                    saved[f'{label}_b{bi}_{key}'] = np.asarray(value)
            uv = uv_coefficient(blocks)
            row = dict(model=label, dimension=sum(len(b['q'][0]) for b in blocks), cutoff=cutoff,
                heat=high.tolist(), heat_quadrature_change=abs(high-low).tolist(),
                heat_frequency_tail=tails, coefficient_bounds=bounds,
                finite_pv_integral=pv_high, pv_quadrature_change=abs(pv_high-pv_low),
                pv_frequency_tail=pv_tail, finite_uv_coefficient=uv,
                uv_at_smallest_time=float(np.sqrt(times[0])*high[0]),
                uv_relative_discrepancy=abs(np.sqrt(times[0])*high[0]/uv-1))
            result['records'].append(row)
            print(json.dumps(dict(model=label, pv=pv_high, heat_quadrature=max(abs(high-low)/(1+abs(high))),
                                  seconds=perf_counter()-start)), flush=True)
        records = {r['model']: r for r in result['records']}
        for name, lo, hi in (('earlier_radial','r8_a3_e48','r10_a3_e48'),
            ('radial','r10_a3_e48','r12_a3_e48'), ('angular','r12_a4_e48','r12_a5_e48'),
            ('exterior','r12_a5_e48','r12_a5_e80'), ('quadrature','r12_a5_e80','q80_r12_a5_e80')):
            before, after = (np.array(records[key]['heat']) for key in (lo, hi))
            relative = abs(after-before)/np.maximum(abs(after), 1e-12)
            result['comparisons'].append(dict(name=name, smaller=lo, larger=hi, relative_change=relative.tolist()))
        # This is only a sensitivity window in the tested families, never a continuum error certificate.
        changes = [r['relative_change'] for r in result['comparisons'] if r['name'] in ('radial','angular','exterior')]
        window = np.max(changes, axis=0) < .01
        window &= abs(np.array(records['r12_a5_e80']['heat'])) > 1e-10
        result['one_percent_sensitivity_times'] = times[window].tolist()
        reference = next(r for r in baselines['dihole_angular_crosscheck']['records'] if r['phase']==0 and r['angular_modes']==5)['full_frequency']
        recovery = abs(records['r12_a5_e48']['finite_pv_integral']-reference['middle_integral'])
        result['frozen_pv_recovery'] = dict(error=recovery, old_tail=reference['certified_matrix_tail'])
        rows = result['records']
        result['tests'] = dict(
            analytic_gaussian=result['benchmark']['gaussian_maximum_error'] < 2e-11,
            noncommuting_variation=result['benchmark']['noncommuting_derivative_error'] < 2e-10,
            gaussian_tail_control=result['benchmark']['sampled_tail'] <= result['benchmark']['outward_tail'],
            exact_input_hashes=all(sha(project_path(n))==h for n,h in inputs.items()),
            frequency_tail=max(max(r['heat_frequency_tail']) for r in rows) < 1e-10,
            heat_quadrature=max(max(np.array(r['heat_quadrature_change'])/(1+abs(np.array(r['heat'])))) for r in rows) < 1e-8,
            pv_quadrature=max(r['pv_quadrature_change'] for r in rows) < 1e-9,
            finite_uv_limit=max(r['uv_relative_discrepancy'] for r in rows) < .01,
            old_pv_recovered=recovery < reference['certified_matrix_tail']+records['r12_a5_e48']['pv_frequency_tail']+1e-10,
            combined_quadrature=max(next(r for r in result['comparisons'] if r['name']=='quadrature')['relative_change']) < 2e-6)
        result['tests'] = {key: bool(value) for key, value in result['tests'].items()}
        result.update(complete=True, passed=all(result['tests'].values()))
    except Exception as exc:
        result['failure'] = dict(type=type(exc).__name__, message=str(exc))
        raise
    finally:
        result['seconds'] = perf_counter()-start
        if saved:
            state_path=project_path(NAME+'_state.npz')
            np.savez_compressed(state_path, **saved)
            result['state_sha256']=sha(state_path)
        project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=result['tests'],
                          sensitivity_times=result['one_percent_sensitivity_times'],seconds=result['seconds']),indent=2),flush=True)
    if not result['passed']:
        raise ArithmeticError('Proper heat checks failed; all computed evidence was retained.')
    plot(result)


if __name__ == '__main__':
    import sys
    if sys.argv[1:] not in ([], ['--recover-serialization-failure']):
        raise SystemExit('Usage: dihole_proper_heat.py [--recover-serialization-failure]')
    run(recover=bool(sys.argv[1:]))
