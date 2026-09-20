"""Physical-mass Pauli--Villars pilot on the existing connected finite model.

Both Weyls are explicit. The exterior weak operator has diagonal blocks
A0+kappa*W and A0.adjoint()+kappa*W, and off-diagonal i*m*P. This preserves
the discrete total-current identity; m is a fixed proper mass, not a
frequency shift. The absorbing outer boundary and finite angular sector
remain numerical approximations, not a covariant four-dimensional limit.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import quad, simpson
from scipy.linalg import block_diag
from scipy.sparse import block_diag as sparse_block_diag, bmat
from scipy.sparse.linalg import splu

from balanced_dihole import symbolic_fields
from dihole_bridge_scattering import complex_record, embed
from dihole_direct_source import direct_response
from dihole_massive_bridge import MassiveBridgeMatrices
from dihole_return_pilot import sew
from dihole_source_window import prepare_exterior, return_map
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


PV_WEIGHTS = np.array([1., -3., 3., -1.])


def prepare_massive(prepared):
    model = prepared['model']
    active = model.active
    local = [np.kron(np.einsum('p,pa,pb->ab', cell['weights']*cell['lapse'],
                              cell['basis'], cell['basis']), np.eye(2))
             for cell in model.cells]
    proper = sparse_block_diag(local, format='csc')[active][:, active]
    weak = prepared['mass'][active][:, active]
    zero = prepared['zero_operator'][active][:, active]
    lift = prepared['lift'][active]
    receiver = prepared['receiver'][:, active]
    # Dual ports reverse the physical characteristic directions of L.
    full_lift = block_diag(lift, receiver.conj().T)
    full_receiver = block_diag(receiver, lift.conj().T)
    return dict(prepared=prepared, zero=zero, proper=proper, weak=weak,
                lift=full_lift, receiver=full_receiver,
                weight=sparse_block_diag((weak, weak), format='csc'),
                minimum_local_proper_eigenvalue=float(min(np.linalg.eigvalsh(v)[0] for v in local)))


def operator(system, kappa, mass):
    if not np.isfinite(kappa) or kappa <= 0 or not np.isfinite(mass) or mass < 0:
        raise ValueError('Positive imaginary frequency and nonnegative proper mass required.')
    a, w, p = system['zero'], system['weak'], system['proper']
    return bmat([[a+kappa*w, 1j*mass*p],
                 [1j*mass*p, a.conj().T+kappa*w]], format='csc')


def physical_ledger(system, matrix, fields, coefficients, kappa):
    """Independent boundary/jump/volume quadrature for both physical Weyls."""
    prepared = system['prepared']; model = prepared['model']
    active = model.active; count = len(active)
    vector = fields@coefficients
    solution = np.zeros((2, model.ndof), complex)
    solution[:, active] = vector.reshape(2, count)
    ports = prepared['ports']; modes = prepared['lift'].shape[1]
    incoming = [{}, {}]
    for chirality, key in enumerate(('minus', 'plus')):
        start = chirality*modes
        for port in ports:
            size = port[key].shape[-1]
            values = port[key]@coefficients[start:start+size]
            start += size
            offset = 0
            for edge in port['edges']:
                end = offset+len(edge['weights'])
                incoming[chirality][id(edge)] = values[offset:end]
                offset = end
    outgoing = dict.fromkeys(('upper', 'lower', 'outer', 'axis'), 0.)
    input_norm = mismatch = jump = bulk = 0.
    for chirality, u in enumerate(solution):
        for edge in model.boundary:
            trace = model.trace(u, edge)
            out = trace@edge['pp' if chirality == 0 else 'pm'].T
            inc = trace@edge['pm' if chirality == 0 else 'pp'].T
            prescribed = incoming[chirality].get(id(edge), np.zeros_like(trace))
            norm = lambda v: float(np.einsum('p,ps,ps->', edge['weights'], v.conj(), v).real)
            outgoing[edge['kind']] += norm(out)
            input_norm += norm(prescribed)
            mismatch += norm(inc-prescribed)
        for edge in model.interior:
            left = model.trace(u, edge)
            j = edge['other']
            right = edge['other_basis']@u[6*j:6*j+6].reshape(3, 2)
            delta = left-right
            jump += float(np.einsum('p,ps,ps->', edge['weights'], delta.conj(), delta).real)
        for cell in model.cells:
            j = cell['triangle']
            value = cell['basis']@u[6*j:6*j+6].reshape(3, 2)
            bulk += float(np.einsum('p,ps,ps->', cell['weights'], value.conj(), value).real)
    rhs = system['lift']@coefficients
    residual = matrix@vector-rhs
    balance = sum(outgoing.values())+mismatch+jump+2*kappa*bulk-input_norm
    algebraic = float(2*np.vdot(vector, residual).real)
    return dict(input_norm=input_norm, outgoing=outgoing, incoming_mismatch=mismatch,
                interior_jump=jump, attenuation=2*kappa*bulk, balance=balance,
                algebraic_residual=algebraic, identity_error=abs(balance-algebraic),
                volume_quadrature_error=abs(bulk-np.vdot(vector, system['weight']@vector).real))


def massive_return(system, kappa, mass, ledger=False):
    matrix = operator(system, kappa, mass)
    fields = splu(matrix).solve(system['lift'])
    returned = system['receiver']@fields
    residual = matrix@fields-system['lift']
    defect = np.eye(len(returned))-returned.conj().T@returned-2*kappa*(fields.conj().T@(system['weight']@fields))
    diagnostics = dict(maximum_singular=float(np.linalg.svd(returned, compute_uv=False)[0]),
        relative_linear_residual=float(np.linalg.norm(residual)/np.linalg.norm(system['lift'])),
        minimum_unassigned_loss=float(np.linalg.eigvalsh((defect+defect.conj().T)/2)[0]))
    if ledger:
        rng = np.random.default_rng(135)
        coef = rng.normal(size=len(returned))+1j*rng.normal(size=len(returned))
        coef /= np.linalg.norm(coef)
        diagnostics['ledger'] = physical_ledger(system, matrix, fields, coef, kappa)
        # Non-solutions catch an erroneous integration-by-parts identity that
        # a small solved residual could conceal. One column is sufficient.
        trial = rng.normal(size=(matrix.shape[0], 1))+1j*rng.normal(size=(matrix.shape[0], 1))
        trial /= np.linalg.norm(trial)
        padded = np.zeros_like(fields); padded[:, 0] = trial[:, 0]
        one = np.eye(len(returned))[0]
        diagnostics['off_shell_ledger'] = physical_ledger(system, matrix, padded, one, kappa)
    return returned, diagnostics


def sew_full(scattering, phase, size, fi, fo):
    """Same charged loop phase for both Weyls, in exterior R/L ordering."""
    p, m, n = size+1, size, 2*size+1
    s = scattering.copy()
    s[n:, :n] *= np.exp(1j*phase)
    s[:n, n:] *= np.exp(-1j*phase)
    rows = np.r_[n+np.arange(p), np.arange(m), n+p+np.arange(m), m+np.arange(p)]
    cols = np.r_[n+np.arange(m), np.arange(p), n+m+np.arange(p), p+np.arange(m)]
    raw = s[np.ix_(rows, cols)]
    full_i, full_o = block_diag(fi, fo), block_diag(fo, fi)
    return np.linalg.solve(full_o.T, (full_i@raw).T).T


def bridge_response(model, prepared, returned, kappa, phases=(0., np.pi)):
    eta = 1e-4
    n = model.plus
    scattering = [embed(lambda s, e=e: model.coefficients(s, e), n, n, 1j*kappa)[0]
                  for e in (0., -eta, eta)]
    direct, _ = direct_response(model.coefficients, n, n, kappa)
    rows = []
    for phase in phases:
        b, bm, bp = [sew_full(s, phase, model.size, prepared['fi'], prepared['fo']) for s in scattering]
        denominator = np.eye(len(b))-b@returned
        dlog = -float(np.trace(np.linalg.solve(denominator, ((bp-bm)/(2*eta))@returned)).real)
        fd = (np.linalg.slogdet(np.eye(len(b))-bp@returned)[1]
              -np.linalg.slogdet(np.eye(len(b))-bm@returned)[1])/(2*eta)
        rows.append(dict(phase=float(phase), loop_dlog=dlog,
                         loop_logdet=float(np.linalg.slogdet(denominator)[1]),
                         determinant_derivative_error=float(abs(dlog-fd))))
    return dict(direct=direct, phases=rows)


def flat_control():
    """Analytic PV moments and an independent stable bulk quadrature."""
    scale, transverse, meridian, weight = 1.7, .8, 2.3, 3.1
    j = np.arange(4)
    t2 = transverse**2+meridian**2*j*scale**2
    # Rationalized form of 2*x-(2*x*x+1)/sqrt(1+x*x).
    integral, error = quad(lambda x: -1/(np.hypot(x, 1)*(np.hypot(x, 1)+x)**2), 0, np.inf,
                           epsabs=2e-12, epsrel=2e-12)
    residues = t2/weight*integral
    return dict(mass_scale=scale, moment_sums=[float(PV_WEIGHTS@(j**r)) for r in range(3)],
                dimensionless_bulk_integral=integral, quadrature_error=error,
                individual_residues=residues.tolist(), combined_residue=float(PV_WEIGHTS@residues))


def control_run(prepared, system, zero_model):
    kappa = .1
    right, _ = return_map(prepared, kappa)
    expected = block_diag(right, right.conj().T)
    actual, diagnostics = massive_return(system, kappa, 0., ledger=True)
    full = bridge_response(zero_model, prepared, actual, kappa, (0., .7, np.pi))
    p, m = zero_model.size+1, zero_model.size
    single = [embed(lambda s, e=e: zero_model.weyl.coefficients(s, e), p, m, 1j*kappa)[0]
              for e in (0., -1e-4, 1e-4)]
    doubling = []
    for phase, row in zip((0., .7, np.pi), full['phases']):
        b, bm, bp = [sew(s, phase, zero_model.size, prepared['fi'], prepared['fo']) for s in single]
        denominator = np.eye(len(b))-b@right
        derivative = -float(np.trace(np.linalg.solve(denominator, ((bp-bm)/2e-4)@right)).real)
        doubling.append(dict(phase=phase, log_error=abs(row['loop_logdet']-2*np.linalg.slogdet(denominator)[1]),
                             derivative_error=abs(row['loop_dlog']-2*derivative)))
    t = perf_counter()
    massive, info = massive_return(system, kappa, 1., ledger=True)
    return dict(flat=flat_control(), minimum_proper_eigenvalue=system['minimum_local_proper_eigenvalue'],
                massless_return_error=float(np.max(abs(actual-expected))),
                massless_diagnostics=diagnostics, massless_determinant_doubling=doubling,
                massive_diagnostics=info, massive_return=complex_record(massive),
                massive_solve_seconds=perf_counter()-t,
                maximum_ledger_identity_error=max(d[key]['identity_error'] for d in (diagnostics, info)
                                                  for key in ('ledger', 'off_shell_ledger')))


def input_hashes():
    names = ('dihole_massive_bridge.py', 'dihole_source_window.py', 'dihole_direct_source.py',
             'dihole_curved_exterior.py', 'dihole_return_pilot.py', 'dihole_bridge_scattering.py',
             'dihole_cylinder_transfer.py', 'dihole_stretched_collar.py', 'balanced_dihole.py',
             'exterior_boundary_map.py', 'two_mouth_mesh.py', 'spherical_weyl_leads.py', 'project_paths.py')
    return {str(project_path(name)): hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names}


def run(probe=False):
    started = perf_counter()
    output = project_path('dihole_pv_probe_results.json' if probe else 'dihole_pv_source_results.json')
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    hashes = input_hashes()
    report = dict(complete=False, source_sha256=source_hash, input_sha256=hashes, records=[],
                  parameters=dict(radius=16., angular_size=4, q=1, mu=.5, pv_scale=1.,
                                  squared_mass_multipliers=[0, 1, 2, 3], weights=PV_WEIGHTS.tolist(),
                                  window=[.01, 10.], points=17))
    if output.exists():
        saved = json.loads(output.read_text())
        if saved['source_sha256'] != source_hash or saved['input_sha256'] != hashes:
            raise RuntimeError('Changed source/inputs: preserve the old checkpoint before a changed run.')
        if saved['complete']:
            print('Exact run already complete; not recomputing.'); return saved
        report = saved
    function, _ = symbolic_fields()
    models = [MassiveBridgeMatrices(function, exact_target_parameters(), 4, 48, mass=np.sqrt(j)) for j in range(4)]
    prepared = prepare_exterior(16.)
    system = prepare_massive(prepared)
    report['geometry'] = prepared['geometry']
    if 'controls' not in report:
        report['controls'] = control_run(prepared, system, models[0])
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        print(json.dumps(dict(controls=report['controls'], elapsed=perf_counter()-started)), flush=True)
    controls = report['controls']
    tests = dict(flat_moments=controls['flat']['moment_sums'] == [0., 0., 0.],
                 flat_bulk_cancellation=abs(controls['flat']['combined_residue']) < 2e-12,
                 flat_bulk_exact=abs(controls['flat']['dimensionless_bulk_integral']+.5) < 2e-12,
                 positive_mass_weak_form=controls['minimum_proper_eigenvalue'] > 0,
                 massless_adjoint_return=controls['massless_return_error'] < 2e-10,
                 massless_determinant_doubling=max(max(r['log_error'], r['derivative_error'])
                                                    for r in controls['massless_determinant_doubling']) < 3e-6,
                 independent_current_identity=controls['maximum_ledger_identity_error'] < 2e-9)
    if not all(tests.values()):
        report.update(passed=False, tests=tests)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        raise RuntimeError('Physical-mass integration controls failed; evidence retained.')
    if not probe:
        grid = np.linspace(np.log(.01), np.log(10.), 17)
        done = {r['grid_index'] for r in report['records']}
        for index in list(range(0, 17, 2))+list(range(1, 17, 2)):
            if index in done: continue
            t = perf_counter(); kappa = float(np.exp(grid[index]))
            row = dict(grid_index=index, kappa=kappa, masses=[])
            for j, model in enumerate(models):
                if j == 0:
                    right, diag = return_map(prepared, kappa)
                    returned = block_diag(right, right.conj().T)
                else:
                    returned, diag = massive_return(system, kappa, model.mass, ledger=index in (0, 16) and j == 3)
                response = bridge_response(model, prepared, returned, kappa)
                row['masses'].append(dict(mass=model.mass, diagnostics=diag, **response))
            # The identical leading-volume derivative cancels analytically
            # since sum(c)=0; use the remaining values to avoid cancellation.
            sub = np.array([r['direct']['power_subtracted_dlog'] for r in row['masses']])
            raw = np.array([r['direct']['direct_dlog'] for r in row['masses']])
            loop = np.array([[p['loop_dlog'] for p in r['phases']] for r in row['masses']])
            row.update(pv_direct_dlog=float(PV_WEIGHTS@sub), pv_loop_dlog=(PV_WEIGHTS@loop).tolist(),
                       pv_full_dlog=(PV_WEIGHTS@(sub[:, None]+loop)).tolist(),
                       raw_combination_difference=float(abs(PV_WEIGHTS@(raw-sub))),
                       runtime_seconds=perf_counter()-t)
            report['records'].append(row)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            print(json.dumps(dict(completed=len(report['records']), kappa=kappa,
                                  pv_full_dlog=row['pv_full_dlog'], seconds=row['runtime_seconds'])), flush=True)
        rows = sorted(report['records'], key=lambda r: r['grid_index'])
        report['integrals'] = {}
        for stride, name in ((1, 'fine'), (2, 'coarse')):
            selected = rows[::stride]
            k = np.array([r['kappa'] for r in selected])
            integrate = lambda vals: (-1/np.pi*simpson(k[:, None]*np.asarray(vals), x=np.log(k), axis=0)).tolist()
            report['integrals'][name] = dict(points=len(selected),
                full=integrate([r['pv_full_dlog'] for r in selected]),
                direct=integrate([[r['pv_direct_dlog']] for r in selected])[0],
                loop=integrate([r['pv_loop_dlog'] for r in selected]))
        masses = [m for r in rows for m in r['masses']]
        tests.update(complete_grid={r['grid_index'] for r in rows} == set(range(17)),
                     passive_maps=max(m['diagnostics']['maximum_singular'] for m in masses) < 1+1e-8,
                     linear_solves=max(m['diagnostics']['relative_linear_residual'] for m in masses) < 2e-10,
                     nonnegative_remaining_loss=min(m['diagnostics'].get('minimum_unassigned_loss', 0.) for m in masses) > -2e-9,
                     determinant_derivatives=max(p['determinant_derivative_error'] for m in masses for p in m['phases']) < 3e-6,
                     common_volume_cancellation=max(r['raw_combination_difference'] for r in rows) < 2e-7,
                     extreme_frequency_ledger=max(d[key]['identity_error'] for m in masses
                        if 'ledger' in (d := m['diagnostics']) for key in ('ledger', 'off_shell_ledger')) < 2e-9)
    tests['frozen_inputs'] = input_hashes() == hashes
    report.update(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                  runtime_seconds=perf_counter()-started,
                  normalization='Full Dirac: -1/pi integral PV-weighted direct+loop metric derivative d(kappa). hbar=c=D=1.',
                  scope=['Finite M=1/D, one mu sector, four angular pairs and finite frequency window only.',
                         'Negative PV coefficients are subtraction weights, never a physical material or state.',
                         'Outer incoming-zero boundary and P1 jump loss remain; no physical exterior or spatial/angular limit.',
                         'No omitted frequency tail bound, all-mu sum, large-M counterterms, Hadamard matching or Ward validation.',
                         'No renormalized supporting stress, stability, formation or backward travel is established.'])
    if not probe: render(report)
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=report['passed'], tests=report['tests'], integrals=report.get('integrals'),
                          seconds=report['runtime_seconds']), indent=2))
    if not report['passed']: raise RuntimeError('Connected PV pilot controls failed; results retained.')
    return report


def render(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = sorted(report['records'], key=lambda r: r['kappa'])
    fig, axes = plt.subplots(1, 2, figsize=(9.3, 4.1), layout='constrained')
    k = [r['kappa'] for r in rows]
    for j, label in enumerate(('phase 0', 'phase pi')):
        axes[0].semilogx(k, [-r['pv_full_dlog'][j]/np.pi for r in rows], label=label)
    axes[0].set(xlabel='Imaginary frequency kappa D', ylabel='Finite-M PV response integrand')
    values = report['integrals']['fine']
    axes[1].bar(['direct', 'loop', 'combined'], [values['direct'], values['loop'][0], values['full'][0]],
                color=['tab:orange', 'tab:blue', 'tab:purple'])
    axes[1].set(ylabel='Finite-window response [1/D]', title='M D = 1; phase 0')
    for ax in axes: ax.axhline(0, color='gray', lw=.7); ax.grid(axis='y', alpha=.2)
    axes[0].legend()
    fig.suptitle('Physical-mass regulator pilot; supporting stress remains undetermined', fontsize=11)
    for ext in ('png', 'pdf'): fig.savefig(project_path('dihole_pv_source_progress.'+ext), dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', action='store_true', help='Run integration controls and a massive exterior timing probe only.')
    run(parser.parse_args().probe)
