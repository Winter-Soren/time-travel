"""Frequency-dependent finite-node response with anchored bilinear-arm maps.

Real subthreshold energies have flux amplitudes. Upper-half-plane responses
are analytic diagnostics, not probabilities or a renormalized vacuum source.
"""
import gc
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from null_saddle_scattering import SaddleDG, AU, FACES
from null_saddle_matched import node_maps, match
from project_paths import project_path


class FrequencyDG(SaddleDG):
    """Existing complex H-E assembly with the missing absorption ledger term."""
    def __init__(self, radius, spacing, eta=1., energy=0., constant_mass=None):
        energy = complex(energy)
        if not np.isfinite([energy.real, energy.imag]).all() or energy.imag < 0:
            raise ValueError('Require finite energy in the closed upper half plane.')
        super().__init__(radius, spacing, eta, energy, constant_mass)

    def ledger(self, solution, incoming_data, rhs):
        row = super().ledger(solution, incoming_data, rhs)
        norm = 0.
        for cell in self.cells:
            it = cell['triangle']
            field = cell['basis']@solution[6*it:6*it+6].reshape(3, 2)
            norm += float(np.einsum('p,ps,ps->', cell['weights'], field.conj(), field).real)
        absorption = 2*self.energy.imag*norm
        row['bulk_norm'] = norm
        row['bulk_absorption'] = absorption
        row['balance_residual'] += absorption
        residual = self.matrix@solution-rhs
        row['identity_error'] = abs(row['balance_residual']-float(2*np.vdot(solution, residual).real))
        return row


def frequency_group(energy, spacing, specifications, radius=2.5, pairs=4, eta=1.):
    from null_saddle_frequency_tail import frequency_port
    model = FrequencyDG(radius, spacing, eta, energy)
    maps = node_maps(model, pairs)
    results = {}
    for name, cap, rtol in specifications:
        tails = [frequency_port(face, radius, eta, energy, pairs=pairs, Rfar=cap, rtol=rtol) for face in FACES]
        row = match(model, maps, tails, name)
        row['energy'] = dict(real=float(model.energy.real), imag=float(model.energy.imag))
        row['squared_amplitudes'] = row.pop('outgoing_zero_powers')
        row['squared_amplitude_sum'] = row.pop('zero_power_sum')
        row['anchored_tail_current_difference'] = row.pop('signed_tail_current_error')
        row.pop('zero_power_ledger_error')
        row['amplitude_ledger_error'] = [float(abs(1-row['squared_amplitude_sum'][j]
            -row['ledger'][j]['numerical_loss']-row['ledger'][j]['bulk_absorption']
            -row['omitted_outgoing_norm'][j]-row['signed_finite_face_Gram_term'][j]
            -row['anchored_tail_current_difference'][j])) for j in range(2)]
        row['interpretation'] = ('Real-energy flux amplitudes referred to the near arm planes.' if model.energy.imag == 0 else
            'Anchored analytic amplitudes at complex energy; squared moduli are NOT probabilities or physical far flux.')
        results[name] = row
        print(json.dumps(dict(case=name, energy=row['energy'], squared_amplitudes=row['squared_amplitudes'],
                              absorption=[x['bulk_absorption'] for x in row['ledger']])), flush=True)
    return results


def absorption_control(spacing):
    energy = .7+.4j
    _, spinors = np.linalg.eigh(AU)
    spinor = spinors[:, -1]
    exact = lambda points: np.exp(1j*energy*points[:, 0])[:, None]*spinor
    model = FrequencyDG(1., spacing, energy=energy, constant_mass=0.)
    rhs, data, _ = model.input(lambda kind, points, normal: exact(points))
    solution = model.solve(rhs)
    return dict(spacing=spacing, l2_error=model.l2_error(solution, exact),
                exact_bulk_norm=2*math.sinh(2*energy.imag)/energy.imag,
                ledger=model.ledger(solution, data, rhs))


def matrix(row):
    return np.asarray(row['scattering']['real'])+1j*np.asarray(row['scattering']['imag'])


def checks():
    from null_saddle_frequency_tail import self_check
    helper = self_check()
    plane = [absorption_control(h) for h in (.25, .125)]
    cases = {}
    for name, energy in [('zero', 0.), ('real_half', .5), ('real_minus_half', -.5), ('real_one_coarse', 1.),
                         ('imag_small', .2j), ('imag_half_coarse', .5j), ('imag_one', 1j)]:
        cases.update(frequency_group(energy, .125, [(name, 12., 1e-10)]))
        gc.collect()
    for energy, prefix in [(1., 'real_one'), (.5j, 'imag_half')]:
        specs = [(prefix, 12., 1e-10), (prefix+'_cap6', 6., 1e-10), (prefix+'_cap24', 24., 1e-10)]
        cases.update(frequency_group(energy, .0625, specs))
        gc.collect()
    cases.update(frequency_group(.5j, .0625, [('larger_box', 12., 1e-10)], radius=3.))
    gc.collect()
    cases.update(frequency_group(.5j, .0625, [('extra_pairs', 12., 1e-10)], pairs=6))
    old = json.loads(project_path('null_saddle_matched_results.json').read_text())['cases']['coarse']
    zero_error = float(np.max(abs(matrix(cases['zero'])-matrix(old))))
    conjugation = float(np.max(abs(matrix(cases['real_minus_half'])+matrix(cases['real_half']).conj())))
    imaginary_real_part = max(float(np.max(abs(matrix(row).real))) for row in cases.values() if row['energy']['imag'] > 0)
    change = lambda a, b: float(np.max(abs(matrix(cases[a])-matrix(cases[b]))))
    comparisons = dict(zero_energy_regression=zero_error, particle_hole_symmetry=conjugation,
                       imaginary_axis_real_part=imaginary_real_part,
                       real_mesh=change('real_one', 'real_one_coarse'), imag_mesh=change('imag_half', 'imag_half_coarse'),
                       real_cap6=change('real_one', 'real_one_cap6'), real_cap24=change('real_one', 'real_one_cap24'),
                       imag_cap6=change('imag_half', 'imag_half_cap6'), imag_cap24=change('imag_half', 'imag_half_cap24'),
                       imag_pairs=change('imag_half', 'extra_pairs'),
                       imag_radius_anchored=float(np.max(abs(np.exp(-2j*.5j*(3.-2.5))*matrix(cases['larger_box'])-matrix(cases['imag_half'])))))
    rows = list(cases.values())
    tests = dict(zero_energy_regression=zero_error < 2e-10, particle_hole_symmetry=conjugation < 2e-9,
                 imaginary_axis_reality=imaginary_real_part < 2e-9,
                 manufactured_refinement=plane[1]['l2_error'] < plane[0]['l2_error']/3,
                 manufactured_absorption=abs(plane[1]['ledger']['bulk_norm']-plane[1]['exact_bulk_norm']) < 2e-4,
                 current_absorption_identity=max(row['identity_error'] for item in rows for row in item['ledger']) < 2e-10,
                 amplitude_ledger=max(max(row['amplitude_ledger_error']) for row in rows) < 2e-10,
                 linear_solve=max(item['linear_residual_relative'] for row in rows for item in row['ledger']) < 2e-10,
                 nonnegative_absorption=min(item['bulk_absorption'] for row in rows for item in row['ledger']) >= 0,
                 real_tail_current=max(max(abs(np.asarray(row['anchored_tail_current_difference']))) for row in rows if row['energy']['imag'] == 0) < 2e-8)
    tests = {name: bool(value) for name, value in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests, helper_checks=helper, absorption_controls=plane,
                  cases=cases, comparisons=comparisons,
                  convention='Each arm zero amplitude is phase/attenuation referenced to the near plane r=R; trivial remote propagation removed.',
                  limits=['Imaginary-frequency squared amplitudes are not probabilities, physical far flux, or a stress tensor.',
                          'Finite DG/Hermite/cap model; corner overlap and infinite-channel errors remain uncontrolled.',
                          'Changing near-plane R changes the reference: compare exp[-2iE(R2-R1)] S_R2 to S_R1.',
                          'No global curved exterior, vacuum renormalization, self-consistent throat or chronology change.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in
                                ('null_saddle_frequency_tail.py', 'null_saddle_tail.py', 'null_saddle_matched.py',
                                 'null_saddle_matched_results.json', 'null_saddle_scattering.py', 'null_saddle_leads.py',
                                 'null_saddle_corner.py', 'null_saddle_corner_results.json', 'project_paths.py')})
    if not result['passed']:
        raise AssertionError(json.dumps(dict(tests=tests, comparisons=comparisons), indent=2))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from null_saddle_corner import corner_gap_squared
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.9), layout='constrained')
    ax = axes[0]
    names = ['zero', 'real_half', 'real_one_coarse']
    for index, label, color in [(0, 'Right output', '#267f81'), (1, 'Left output', '#b07b28')]:
        values = [result['cases'][name]['squared_amplitudes'][index][0] for name in names]
        ax.plot([0., .5, 1.], values, 'o--', color=color, label=label)
    ax.set(xlabel='Real energy E (local units)', ylabel='Raw outgoing flux fraction', ylim=(0, 1.04), title='Real-frequency quantum routing')
    ax.legend(fontsize=9)
    ax.text(.01, -.25, 'Three finite-model points; finer E=1 control saved.\nInput and output phases have fixed reference planes.', transform=ax.transAxes, fontsize=9)
    ax = axes[1]
    names = ['zero', 'imag_small', 'imag_half_coarse', 'imag_one']
    ax.semilogy([0., .2, .5, 1.], [np.linalg.norm(matrix(result['cases'][name]), 2) for name in names], 'o--', color='#536eaf')
    ax.set(xlabel=r'Imaginary energy $E=i\kappa$', ylabel=r'Anchored response norm $\|S(i\kappa)\|$', title='Response used toward a vacuum calculation')
    ax.text(.01, -.25, 'Analytic response, not a transmission probability.\nThis alone is not a quantum stress or energy.', transform=ax.transAxes, fontsize=9)
    ax = axes[2]
    radius = np.linspace(1.2, 3.2, 100)
    gap = np.sqrt([corner_gap_squared(float(r), 1.) for r in radius])
    ax.plot(radius, gap, color='#7a608c', lw=2)
    ax.plot([2.5], [math.sqrt(corner_gap_squared(2.5, 1.))], 'o', color='#7a608c')
    ax.set(xlabel='Corner begins at |u|, |v| ≥ R', ylabel='Local lower gap bound', title=r'Corner coercivity for $|\eta|\leq1$')
    ax.text(.01, -.25, 'For compactly supported fields in the corner.\nNot a boundary-map or scattering-error bound.', transform=ax.transAxes, fontsize=9)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Frequency response and the remaining exterior problem | Local Dirac model', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('null_saddle_spectrum_progress.'+ext), dpi=175)
    plt.close(fig)


if __name__ == '__main__':
    result = checks()
    project_path('null_saddle_spectrum_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result)
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'comparisons')}, indent=2))
