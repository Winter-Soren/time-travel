"""Complete fixed-l<=4 bridge basis across all nine azimuthal sectors.

Use the physical-mass exterior and PV prescription candidate from135.
An analytic scattering variation replaces three perturbed bridge solves.
This is still a finite absorbing DG model, finite M and frequency window;
angular completeness at one cutoff is not a continuum or source proof.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import block_diag

from balanced_dihole import symbolic_fields
from dihole_azimuthal_bridge import AzimuthalBridgeMatrices
from dihole_bridge_scattering import collar_mask, embed
from dihole_curved_exterior import curved_model, graded_mesh
from dihole_direct_source import direct_response
from dihole_pv_source import PV_WEIGHTS, input_hashes, massive_return, prepare_massive
from dihole_return_pilot import cut_basis as legacy_cut_basis
from dihole_source_window import return_map
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


def differentiated_bridge(coefficients, n, kappa, rtol=2e-11, atol=2e-13):
    """Scattering, its paired-metric derivative, and direct term in one ODE."""
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError('Use positive finite imaginary frequency.')
    state = np.zeros(8*n*n+3, complex)
    blocks = state[:8*n*n].reshape(8, n, n)
    blocks[1] = blocks[2] = np.eye(n)

    def rhs(s, y):
        r, forward, reverse, left, rz, fz, vz, lz = y[:8*n*n].reshape(8, n, n)
        wp, wm, x = coefficients(s); mask = collar_mask(s)
        a, b, c, d = -kappa*wp, -1j*x, 1j*x.conj().T, kappa*wm
        az, bz, cz, dz = -2*mask*a, -mask*b, -mask*c, -2*mask*d
        ar, cd = a-r@c, c@r+d
        arz, cdz = az-rz@c-r@cz, cz@r+c@rz+dz
        values = (ar@r+b-r@d, ar@forward, -reverse@cd, -reverse@c@forward,
                  arz@r+ar@rz+bz-rz@d-r@dz,
                  arz@forward+ar@fz, -vz@cd-reverse@cdz,
                  -vz@c@forward-reverse@cz@forward-reverse@c@fz)
        trace_w = np.trace(wp)+np.trace(wm)
        direct = np.trace(c@r)+kappa*trace_w/2
        sub = np.trace(c@(rz-mask*r))
        return np.r_[np.array(values).ravel(), direct, sub, mask*trace_w]

    nfev = 0
    for lo, hi in ((-.1, -.05), (-.05, .05), (.05, .1)):
        solution = solve_ivp(rhs, (lo, hi), state, method='DOP853', rtol=rtol, atol=atol)
        if not solution.success: raise RuntimeError(solution.message)
        state = solution.y[:, -1]; nfev += solution.nfev
    if not np.all(np.isfinite(state)): raise ArithmeticError('Nonfinite scattering variation.')
    r, forward, reverse, left, rz, fz, vz, lz = state[:8*n*n].reshape(8, n, n)
    scattering = np.block([[left, reverse], [forward, r]])
    derivative = np.block([[lz, vz], [fz, rz]])
    direct, sub, volume = state[-3:]
    return scattering, derivative, dict(direct_logdet=float(direct.real),
        direct_dlog=float((sub-kappa*volume).real), power_subtracted_dlog=float(sub.real),
        volume_factor=float(volume.real), scalar_imaginary_error=float(max(abs(state[-3:].imag))),
        nfev=int(nfev))


def cut_basis(model, kind, points):
    """Evaluate declared continuum frames directly, without parity guesses."""
    center = .5 if kind == 'upper' else -.5
    theta = np.arctan2(points[:, 0], points[:, 1]-center)
    c, s = np.cos(theta/2), np.sin(theta/2)
    radial_plus, radial_minus = np.stack((c, s)), np.stack((-s, c))
    f, g = model.profiles(theta if kind == 'upper' else np.pi-theta)
    if kind == 'upper':
        return 1j*radial_minus[:, :, None]*g/.05, radial_plus[:, :, None]*f/.05
    return -radial_minus[:, :, None]*f/.05, 1j*radial_plus[:, :, None]*g/.05


def prepare_sector(bridge, triangulation):
    model, mass, geometry = curved_model(kappa=.1, mu=bridge.mu, triangulation=triangulation)
    ports = [model.modal_port(kind, lambda points, kind=kind: cut_basis(bridge, kind, points))
             for kind in ('upper', 'lower')]
    factors = []
    for port in ports:
        pair = []
        for index, projector in ((0, 'pp'), (1, 'pm')):
            raw = np.concatenate([np.einsum('ab,pbc->pac', edge[projector],
                np.moveaxis(cut_basis(bridge, port['kind'], edge['points'])[index], 0, 1))
                                  for edge in port['edges']])
            gram = np.einsum('p,pai,paj->ij', port['weights'], raw.conj(), raw)
            pair.append(np.linalg.cholesky((gram+gram.conj().T)/2).conj().T)
        factors.append(pair)
    return dict(model=model, mass=mass, geometry=geometry, ports=ports,
                zero_operator=(model.matrix_full-.1*mass).tocsc(),
                lift=np.column_stack([p['G'] for p in ports]),
                receiver=np.vstack([p['C'] for p in ports]),
                fi=block_diag(factors[0][1], factors[1][1]),
                fo=block_diag(factors[0][0], factors[1][0]))


def sew(scattering, phase, p, m, fi, fo):
    n = p+m
    value = scattering.copy()
    value[n:, :n] *= np.exp(1j*phase)
    value[:n, n:] *= np.exp(-1j*phase)
    rows = np.r_[n+np.arange(p), np.arange(m), n+p+np.arange(m), m+np.arange(p)]
    cols = np.r_[n+np.arange(m), np.arange(p), n+m+np.arange(p), p+np.arange(m)]
    raw = value[np.ix_(rows, cols)]
    full_i, full_o = block_diag(fi, fo), block_diag(fo, fi)
    return np.linalg.solve(full_o.T, (full_i@raw).T).T


def response(bridge, prepared, returned, kappa):
    scattering, derivative, direct = differentiated_bridge(bridge.coefficients, bridge.n, kappa)
    rows = []
    for phase in (0., np.pi):
        b, bz = [sew(s, phase, bridge.p, bridge.m, prepared['fi'], prepared['fo'])
                 for s in (scattering, derivative)]
        denominator = np.eye(len(b))-b@returned
        dlog = -float(np.trace(np.linalg.solve(denominator, bz@returned)).real)
        rows.append(dict(phase=phase, loop_dlog=dlog,
                         loop_logdet=float(np.linalg.slogdet(denominator)[1]),
                         denominator_condition=float(np.linalg.cond(denominator))))
    return dict(direct=direct, phases=rows)


def controls(function, parameters):
    rows = []
    for mu in (.5, -.5, 1.5):
        bridge = AzimuthalBridgeMatrices(function, parameters, 3, 48, mu=mu, mass=1.)
        n, kappa, eta = bridge.n, .1, 1e-4
        t = perf_counter()
        s, z, direct = differentiated_bridge(bridge.coefficients, n, kappa)
        analytic_seconds = perf_counter()-t
        t = perf_counter()
        base, _ = embed(bridge.coefficients, n, n, 1j*kappa)
        varied = [embed(lambda x, e=e: bridge.coefficients(x, e), n, n, 1j*kappa)[0]
                  for e in (-eta, eta)]
        old_direct, _ = direct_response(bridge.coefficients, n, n, kappa)
        reference_seconds = perf_counter()-t
        refined = differentiated_bridge(bridge.coefficients, n, kappa, 2e-12, 2e-14)
        rows.append(dict(mu=mu, scattering_error=float(np.max(abs(s-base))),
            derivative_fd_error=float(np.max(abs(z-(varied[1]-varied[0])/(2*eta)))),
            direct_error=max(abs(direct[k]-old_direct[k]) for k in
                             ('direct_logdet', 'direct_dlog', 'power_subtracted_dlog')),
            tighter_derivative_change=float(np.max(abs(z-refined[1]))),
            analytic_seconds=analytic_seconds, reference_seconds=reference_seconds))
    bridge = AzimuthalBridgeMatrices(function, parameters, 4, 48)
    theta = np.linspace(.001, np.pi-.001, 83)
    legacy_errors = []
    for kind, center in (('upper', .5), ('lower', -.5)):
        points = np.column_stack((.05*np.sin(theta), center+.05*np.cos(theta)))
        current, old = cut_basis(bridge, kind, points), legacy_cut_basis(kind, points, 4)
        legacy_errors.append(float(max(np.max(abs(a-b)) for a, b in zip(current, old))))
    return dict(derivatives=rows, legacy_cut_errors=legacy_errors)


def hashes():
    values = input_hashes()
    for name in ('dihole_azimuthal_bridge.py', 'complete_angular_modes.py', 'dihole_pv_source.py',
                 'dihole_pv_source_results.json'):
        p = project_path(name); values[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return values


def integrals(records, stride=1):
    selected = sorted(records, key=lambda r: r['grid_index'])[::stride]
    k = np.array([r['kappa'] for r in selected])
    apply = lambda values: (-simpson(k[:, None]*np.asarray(values), x=np.log(k), axis=0)/np.pi).tolist()
    return dict(points=len(selected), full=apply([r['pv_full_dlog'] for r in selected]),
                direct=apply([[r['pv_direct_dlog']] for r in selected])[0],
                loop=apply([r['pv_loop_dlog'] for r in selected]))


def run():
    started = perf_counter()
    output = project_path('dihole_pv_sectors_results.json')
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest(); inputs = hashes()
    result = dict(complete=False, source_sha256=source_hash, input_sha256=inputs,
                  lmax=4, pv_scale=1., window=[.01, 10.], points=17, sectors={})
    if output.exists():
        saved = json.loads(output.read_text())
        if saved['source_sha256'] != source_hash or saved['input_sha256'] != inputs:
            raise RuntimeError('Changed source/inputs: preserve earlier evidence before rerunning.')
        if saved['complete']:
            print('Exact sector run already complete; no recomputation.'); return saved
        result = saved
    function, _ = symbolic_fields(); parameters = exact_target_parameters()
    if 'controls' not in result:
        result['controls'] = controls(function, parameters)
    checks = result['controls']
    tests = dict(legacy_frames=max(checks['legacy_cut_errors']) < 2e-10,
                 analytic_scattering=max(r['scattering_error'] for r in checks['derivatives']) < 2e-9,
                 analytic_derivative=max(r['derivative_fd_error'] for r in checks['derivatives']) < 3e-7,
                 direct_term=max(r['direct_error'] for r in checks['derivatives']) < 2e-8,
                 derivative_tolerance=max(r['tighter_derivative_change'] for r in checks['derivatives']) < 2e-9)
    result['tests'] = tests
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(controls=checks, tests=tests)), flush=True)
    if not all(tests.values()): raise RuntimeError('New derivative/frame controls failed; evidence retained.')
    previous = json.loads(project_path('dihole_pv_source_results.json').read_text())
    result['sectors']['0.5'] = dict(mu=.5, a=0, size=4, reused_from135=True,
        plus_labels=list(range(5)), minus_labels=list(range(1, 5)), records=previous['records'])
    triangulation = graded_mesh(Rtarget=16.)
    grid = np.linspace(np.log(.01), np.log(10.), 17)
    # Complete each magnetic angular multiplet, not equal numbers of modes
    # at each mu. For q=1, a=mu-.5 ranges -l..l at angular label l.
    for a in (-1, 1, -2, 2, -3, 3, -4, 4):
        mu, size = a+.5, 5-abs(a); key = str(mu)
        if key not in result['sectors']:
            result['sectors'][key] = dict(mu=mu, a=a, size=size, records=[], reused_from135=False)
        sector = result['sectors'][key]
        done = {r['grid_index'] for r in sector['records']}
        if done == set(range(17)): continue
        models = [AzimuthalBridgeMatrices(function, parameters, size, 48, mu=mu, mass=np.sqrt(j)) for j in range(4)]
        sector.update(plus_labels=models[0].plus_labels.tolist(), minus_labels=models[0].minus_labels.tolist())
        prepared = prepare_sector(models[0], triangulation); system = prepare_massive(prepared)
        sector['geometry'] = prepared['geometry']
        for index in list(range(0, 17, 2))+list(range(1, 17, 2)):
            if index in done: continue
            t = perf_counter(); kappa = float(np.exp(grid[index]))
            row = dict(grid_index=index, kappa=kappa, masses=[])
            for j, bridge in enumerate(models):
                if j == 0:
                    right, diag = return_map(prepared, kappa)
                    returned = block_diag(right, right.conj().T)
                else:
                    returned, diag = massive_return(system, kappa, bridge.mass,
                                                    ledger=index in (0, 16) and j == 3)
                row['masses'].append(dict(mass=bridge.mass, diagnostics=diag,
                                         **response(bridge, prepared, returned, kappa)))
            direct = np.array([m['direct']['power_subtracted_dlog'] for m in row['masses']])
            loop = np.array([[p['loop_dlog'] for p in m['phases']] for m in row['masses']])
            row.update(pv_direct_dlog=float(PV_WEIGHTS@direct), pv_loop_dlog=(PV_WEIGHTS@loop).tolist(),
                       pv_full_dlog=(PV_WEIGHTS@(direct[:, None]+loop)).tolist(),
                       runtime_seconds=perf_counter()-t)
            sector['records'].append(row)
            output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            print(json.dumps(dict(mu=mu, done=len(sector['records']), kappa=kappa,
                                  pv_full_dlog=row['pv_full_dlog'], seconds=row['runtime_seconds'])), flush=True)
        sector['integrals'] = {label: integrals(sector['records'], stride) for label, stride in (('fine', 1), ('coarse', 2))}
    for sector in result['sectors'].values():
        sector['integrals'] = {label: integrals(sector['records'], stride) for label, stride in (('fine', 1), ('coarse', 2))}
    # Check the actual retained labels, including the lone upper l=0 mode.
    degeneracies = {str(l): dict(upper=sum(l in s['plus_labels'] for s in result['sectors'].values()),
                                lower=sum(l in s['minus_labels'] for s in result['sectors'].values())) for l in range(5)}
    result['angular_degeneracies'] = degeneracies
    result['totals'] = {name: {quantity: np.sum([s['integrals'][name][quantity] for s in result['sectors'].values()], axis=0).tolist()
                               for quantity in ('full', 'direct', 'loop')} for name in ('fine', 'coarse')}
    new = [m for s in result['sectors'].values() if not s['reused_from135'] for r in s['records'] for m in r['masses']]
    ledgers = [m['diagnostics'] for m in new if 'ledger' in m['diagnostics']]
    tests.update(complete_l_multiplets=all(d['upper'] == 2*int(l)+1 and d['lower'] == (0 if l == '0' else 2*int(l)+1)
                                          for l, d in degeneracies.items()),
                 all_sector_grids=all({r['grid_index'] for r in s['records']} == set(range(17)) for s in result['sectors'].values()),
                 passive_exteriors=max(m['diagnostics']['maximum_singular'] for m in new) < 1+1e-8,
                 linear_solves=max(m['diagnostics']['relative_linear_residual'] for m in new) < 2e-10,
                 nonnegative_unassigned_loss=min(m['diagnostics'].get('minimum_unassigned_loss', 0.) for m in new) > -2e-9,
                 independent_current_accounting=max(d[k]['identity_error'] for d in ledgers for k in ('ledger', 'off_shell_ledger')) < 2e-9,
                 frozen_inputs=hashes() == inputs)
    result.update(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                  runtime_seconds=perf_counter()-started,
                  scope=['Complete reference l<=4 basis, all nine mu sectors, 49 one-Weyl or98 full-Dirac angular components.',
                         'Actual angular coefficients mix l within each mu; a finite spherical basis is not a spectral cutoff of the curved operator.',
                         'Proper MD=1, finite frequency window and P1 absorbing exterior; neither angular nor spatial nor frequency tails are bounded.',
                         'Fine/coarse grids show sensitivity, not certified errors. Baseline mu=.5 is reused with its finite-difference derivative.',
                         'No local large-M subtraction or physical continuum/removal limit is performed. Do not add the scale-only convention shift to this pilot.',
                         'No supporting medium, stability, formation, new physical law or backward-time channel established.'])
    render(result)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'], totals=result['totals'], seconds=result['runtime_seconds']), indent=2))
    if not result['passed']: raise RuntimeError('Sector pilot checks failed; evidence retained.')
    return result


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sectors = sorted(result['sectors'].values(), key=lambda s: s['mu'])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), layout='constrained')
    for j, label in enumerate(('phase 0', 'phase pi')):
        axes[0].plot([s['mu'] for s in sectors], [s['integrals']['fine']['full'][j] for s in sectors], 'o-', label=label)
        axes[1].plot([1, 9], [result['sectors']['0.5']['integrals']['fine']['full'][j], result['totals']['fine']['full'][j]], 'o-', label=label)
    axes[0].set(xlabel='Azimuthal sector mu', ylabel='Finite-M response [1/D]')
    axes[0].set_yscale('symlog', linthresh=1e-7)
    axes[1].set(xlabel='Included sectors', ylabel='Summed finite-window response [1/D]')
    axes[1].set_xticks([1, 9], ['Previous one sector', 'All sectors at l <= 4'])
    for ax in axes: ax.axhline(0, color='gray', lw=.7); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Complete finite angular multiplets; physical supporting stress remains unknown', fontsize=11)
    for ext in ('png', 'pdf'): fig.savefig(project_path('dihole_pv_sectors_progress.'+ext), dpi=180)
    plt.close(fig)


if __name__ == '__main__': run()
