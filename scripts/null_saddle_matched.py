"""Match the local Dirac node to finite, varying bilinear-arm continuations.

All powers refer to zero channels at the remote caps. This is still a finite
Hermite arm model: overlapping full-line arm corners and continuum tails are
not certified. Historical frozen-node source/results remain unchanged.
"""
import gc
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag

from null_saddle_scattering import SaddleDG, modal_port, complex_json, FACES
from null_saddle_leads import frozen_lead
from project_paths import project_path


def frozen_port(face, radius, eta, pairs):
    lead = frozen_lead(face, radius, eta, pairs)
    source, output = np.zeros(len(lead.minus_indices), complex), np.zeros(len(lead.plus_indices), complex)
    if face.startswith('v'):
        source[0] = 1
    else:
        output[0] = 1
    return dict(lead=lead, Z=lead.Z, source=source, output=output,
                metadata=dict(frozen=True, face=face))


def node_maps(model, pairs):
    """One interior solve shared by all cap/ODE comparisons at this basis."""
    leads = [frozen_lead(face, model.radius, model.eta, pairs) for face in FACES]
    ports = [modal_port(model, lead) for lead in leads]
    g = np.column_stack([p['G'] for p in ports])
    c = np.vstack([p['C'] for p in ports])
    lift = model.solve(g)
    return dict(ports=ports, G=g, C=c, lift=lift, response=c@lift, pairs=pairs)


def match(model, maps, tails, name):
    ports, g, c, lift = (maps[key] for key in ('ports', 'G', 'C', 'lift'))
    z = block_diag(*[p['Rminus']@tail['Z']@np.linalg.inv(p['Rplus']) for p, tail in zip(ports, tails)])
    drive = np.zeros((g.shape[1], 2), complex)
    offset = 0
    for port, tail in zip(ports, tails):
        count = len(port['Rminus'])
        if port['kind'].startswith('v'):
            drive[offset:offset+count, ('v+', 'v-').index(port['kind'])] = port['Rminus']@tail['source']
        offset += count
    feedback = np.eye(len(c))-maps['response']@z
    outgoing = np.linalg.solve(feedback, maps['response']@drive)
    incoming = drive+z@outgoing
    solution = lift@incoming
    total_rhs = g@incoming
    extracted = c@solution
    scattering = np.zeros((2, 2), complex)
    gram_term, raw_flux = np.zeros(2), np.zeros(2)
    incoming_by_port = []
    ap = bm = 0
    for port, tail in zip(ports, tails):
        np_, nm_ = len(port['Rplus']), len(port['Rminus'])
        a = np.linalg.solve(port['Rplus'], extracted[ap:ap+np_])
        b = np.linalg.solve(port['Rminus'], incoming[bm:bm+nm_])
        for sign, coef, gram in ((1, a, port['Gramplus']), (-1, b, port['Gramminus'])):
            gram_term += sign*np.einsum('ij,ik,kj->j', coef.conj(), gram-np.eye(len(gram)), coef).real
            raw_flux += sign*np.sum(abs(coef)**2, axis=0)
        if port['kind'].startswith('u'):
            scattering[('u+', 'u-').index(port['kind'])] = tail['output']@a
        incoming_by_port.append(incoming[bm:bm+nm_])
        ap += np_; bm += nm_
    powers = abs(scattering)**2
    tail_error = raw_flux-(np.sum(powers, axis=0)-1)
    ledger, omitted, identity = [], [], []
    for column in range(2):
        data = model.modal_incoming(ports, [x[:, column] for x in incoming_by_port])
        row = model.ledger(solution[:, column], data, total_rhs[:, column])
        missing = sum(row['outgoing'].values())-float(np.vdot(extracted[:, column], extracted[:, column]).real)
        error = 1-np.sum(powers[:, column])-row['numerical_loss']-missing-gram_term[column]-tail_error[column]
        ledger.append(row); omitted.append(missing); identity.append(float(abs(error)))
    return dict(name=name, radius=model.radius, eta=model.eta, spacing=model.spacing, pairs=maps['pairs'],
                energy=0., dofs=model.ndof, scattering=complex_json(scattering), outgoing_zero_powers=powers.tolist(),
                zero_power_sum=np.sum(powers, axis=0).tolist(), ledger=ledger, omitted_outgoing_norm=omitted,
                signed_finite_face_Gram_term=gram_term.tolist(), signed_tail_current_error=tail_error.tolist(),
                zero_power_ledger_error=identity, feedback_condition=float(np.linalg.cond(feedback)),
                maximum_face_Gram_deviation=max(p['gram_deviation'] for p in ports),
                maximum_face_Gram_condition=max(p['gram_condition'] for p in ports),
                boundary_response_error=float(np.max(abs(extracted-outgoing))),
                tails=[tail['metadata'] for tail in tails])


def run_group(radius, spacing, pairs, specifications, eta=1.):
    from null_saddle_tail import tail_port
    started = time.perf_counter()
    model = SaddleDG(radius, spacing, eta)
    maps = node_maps(model, pairs)
    results = {}
    for name, far, rtol in specifications:
        tails = ([frozen_port(face, radius, eta, pairs) for face in FACES] if far is None else
                 [tail_port(face, radius, eta, pairs=pairs, Rfar=far, rtol=rtol) for face in FACES])
        row = match(model, maps, tails, name)
        results[name] = row
        print(json.dumps(dict(case=name, powers=row['outgoing_zero_powers'], tail_error=row['signed_tail_current_error'])), flush=True)
    print(json.dumps(dict(group_seconds=time.perf_counter()-started)), flush=True)
    return results


def checks():
    from null_saddle_tail import self_check
    helper = self_check()
    cases = {}
    groups = [(2.5, .125, 4, [('coarse', 12., 1e-10)], 1.),
              (2.5, .0625, 4, [('frozen_fine', None, 1e-10), ('cap6', 6., 1e-10),
                               ('fine', 12., 1e-10), ('cap24', 24., 1e-10), ('tighter', 12., 1e-12)], 1.),
              (3., .0625, 4, [('larger_box', 12., 1e-10)], 1.),
              (2.5, .0625, 6, [('extra_pairs', 12., 1e-10)], 1.),
              (2.5, .125, 4, [('reflected', 12., 1e-10)], -1.)]
    # ponytail: a bounded set of distinct error checks; no parameter sweep.
    for args in groups:
        cases.update(run_group(*args))
        gc.collect()
    matrix = lambda row: np.asarray(row['scattering']['real'])+1j*np.asarray(row['scattering']['imag'])
    power = lambda name: np.asarray(cases[name]['outgoing_zero_powers'])
    old = json.loads(project_path('null_saddle_scattering_results.json').read_text())['cases']['extra_pair']
    frozen_error = float(np.max(abs(matrix(cases['frozen_fine'])-matrix(old))))
    symmetry = float(np.max(abs(matrix(cases['reflected'])-np.array([[0, 1], [1, 0]])@matrix(cases['coarse'])@np.diag([-1., 1.]))))
    changes = {name: float(np.max(abs(power(name)-power('fine')))) for name in
               ('frozen_fine', 'coarse', 'cap6', 'cap24', 'tighter', 'larger_box', 'extra_pairs')}
    rows = list(cases.values())
    tests = dict(frozen_regression=frozen_error < 2e-12, reflected_amplitudes=symmetry < 2e-9,
                 current_identity=max(x['identity_error'] for row in rows for x in row['ledger']) < 2e-10,
                 power_ledger=max(max(row['zero_power_ledger_error']) for row in rows) < 2e-10,
                 tail_current=max(max(abs(np.asarray(row['signed_tail_current_error']))) for row in rows) < 2e-8,
                 interior_solve=max(x['linear_residual_relative'] for row in rows for x in row['ledger']) < 2e-10,
                 boundary_feedback=max(row['boundary_response_error'] for row in rows) < 2e-10,
                 decreasing_mesh_loss=cases['fine']['ledger'][0]['numerical_loss'] < cases['coarse']['ledger'][0]['numerical_loss'],
                 ODE_control=changes['tighter'] < 2e-8,
                 outgoing_projection=all(min(row['omitted_outgoing_norm']) > -2e-10 for row in rows))
    tests = {name: bool(value) for name, value in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests, helper_checks=helper, cases=cases,
                  comparisons=dict(frozen_regression_amplitude_error=frozen_error,
                                   reflection_amplitude_error=symmetry, maximum_power_changes_from_fine=changes),
                  scope='Finite Hermite projections of varying bilinear arms with frozen remote caps attached to a finite DG node.',
                  limits=['Remote-cap, mesh, mode and box changes are sensitivities, not certified error bars.',
                          'Four full-line arms overlap in corners; the independent attachment does not certify the complete plane exterior.',
                          'No global curved-dihole matching, supporting vacuum stress, backreaction or chronology change.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in
                                ('null_saddle_tail.py', 'null_saddle_leads.py', 'null_saddle_scattering.py',
                                 'null_saddle_scattering_results.json', 'exterior_boundary_map.py', 'project_paths.py')})
    if not result['passed']:
        raise AssertionError(json.dumps(dict(tests=tests, comparisons=result['comparisons']), indent=2))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.9), layout='constrained')
    ax = axes[0]
    r = np.linspace(2.5, 12., 200)
    center, width = 1/(2*r), 1/(2*np.sqrt(r))
    ax.fill_between(r, center-width, center+width, color='#267f81', alpha=.25, label='Varying Gaussian: ±1 std. dev.')
    ax.plot(r, center, color='#267f81', lw=2)
    ax.fill_between(r, center[0]-width[0], center[0]+width[0], facecolor='none', edgecolor='#ae7b2c',
                    hatch='///', linewidth=.6, label='Previous frozen Gaussian')
    ax.axvline(2.5, color='#555', ls='--')
    ax.set(xlabel='Outward coordinate r', ylabel='Transverse coordinate', title='Outgoing-arm basis follows the field')
    ax.legend(fontsize=8)
    ax.text(.01, -.24, 'Basis envelope only, not the full quantum field.\nSource and output maps evolve as well.', transform=ax.transAxes, fontsize=9)
    ax = axes[1]
    labels = ['Frozen', 'Varying', 'Larger node', 'More modes']
    names = ['frozen_fine', 'fine', 'larger_box', 'extra_pairs']
    values = [result['cases'][name]['outgoing_zero_powers'][1][0]*100 for name in names]
    ax.bar(labels, values, color=['#ae7b2c', '#267f81', '#5470b0', '#7a608c'])
    ax.set(ylabel='Weaker outgoing channel power (%)', title=r'Actual matched routing at $\eta=1$')
    ax.tick_params(axis='x', labelrotation=20)
    for i, value in enumerate(values):
        ax.text(i, value, f'{value:.5f}', ha='center', va='bottom', fontsize=8)
    ax.set_ylim(0, max(values)*1.19)
    ax.text(.01, -.24, 'Raw powers at the remote caps.\nNumerical losses remain in the saved ledger.', transform=ax.transAxes, fontsize=9)
    ax = axes[2]
    differences = result['comparisons']['maximum_power_changes_from_fine']
    keys = ['cap6', 'cap24', 'tighter', 'coarse', 'larger_box', 'extra_pairs']
    ax.barh(['Cap 6 vs 12', 'Cap 24 vs 12', 'Tighter ODE', 'Finer mesh', 'Larger node', 'More modes'],
            [max(differences[key], 1e-17) for key in keys], color='#5470b0')
    ax.set_xscale('log'); ax.invert_yaxis()
    ax.set(xlabel='Maximum output-power change', title='Different errors checked separately')
    ax.text(.01, -.24, 'Measured sensitivities, not error bounds.\nCorners and curved exterior still unresolved.', transform=ax.transAxes, fontsize=9)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Replacing frozen outer channels | Local Dirac node with varying magnetic arms', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('null_saddle_matched_progress.'+ext), dpi=175)
    plt.close(fig)


if __name__ == '__main__':
    result = checks()
    project_path('null_saddle_matched_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result)
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'comparisons')}, indent=2))
