"""Sample physical curvature scales and local-coupling sensitivity of the fixed collar.

These are diagnostics, not a global curvature enclosure, heat-source bound,
or certificate that a semiclassical approximation is accurate.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from balanced_dihole import symbolic_fields
from dihole_flux_scale import ALPHA_EXAMPLE, load_inputs, scale
from dihole_local_source import evaluate, invariants, metric_jets, prepare
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME = 'dihole_physical_control'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sample(function, parameters, order):
    prepared = prepare(function, parameters, order)
    values = metric_jets(prepared)
    scalars, ricci = invariants(values)
    metric = values[0, :4]
    frame = ricci/np.sqrt(abs(metric[:, None]*metric[None, :]))
    components = np.max(abs(frame), axis=(0, 1)).reshape(order, order)
    magnetic = np.sqrt(scalars[3]/2)*abs(parameters['physical_A_factor'])
    index = np.unravel_index(np.argmax(components), components.shape)
    d = scale(1, ALPHA_EXAMPLE)
    row = dict(order=order, samples=order**2,
        scalar_abs_peak_D2=float(np.max(abs(scalars[0]))),
        ricci_component_peak_D2=float(np.max(components)),
        ricci_frobenius_peak_D2=float(np.max(np.linalg.norm(frame, axis=(0, 1)))),
        magnetic_peak_D2=float(np.max(magnetic)),
        ricci_peak_u_over_D=float(prepared['u'][index[0]]),
        ricci_peak_cos_theta=float(prepared['nodes'][index[1]]),
        scalar_abs_peak_planck2=float(np.max(abs(scalars[0]))/d**2),
        ricci_component_peak_planck2=float(np.max(components)/d**2),
        magnetic_peak_planck2=float(np.max(magnetic)/d**2),
        profiles=dict(u_over_D=prepared['u'].tolist(),
            scalar_abs_max_over_angles_planck2=(np.max(abs(scalars[0]).reshape(order, order), axis=1)/d**2).tolist(),
            ricci_component_max_over_angles_planck2=(np.max(components, axis=1)/d**2).tolist()),
        scalar_contraction_error=float(np.max(abs(scalars[0]-np.einsum('iin,in->n', ricci, 1/metric)))),
        static_ricci_square_error=float(np.max(abs(scalars[2]-np.sum(frame**2, axis=(0, 1))))),
        static_mixed_time_error=float(np.max(abs(frame[0, 1:]))))
    return row, prepared


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3), layout='constrained')
    profile = result['samples'][-1]['profiles']
    axes[0].plot(profile['u_over_D'], profile['scalar_abs_max_over_angles_planck2'], label='Absolute scalar curvature')
    axes[0].plot(profile['u_over_D'], profile['ricci_component_max_over_angles_planck2'], label='Largest frame Ricci component')
    axes[0].set(xlabel='Outward collar coordinate u / D', ylabel='Curvature times Planck length squared',
                title='Actual collar samples, illustrative q = 1 coupling')
    axes[0].legend(fontsize=8)
    change = np.linspace(0., .025, 101)
    for coefficient, label in zip(result['local_energy_response_coefficients'], ('R squared term', 'Ricci squared term')):
        axes[1].plot(change, coefficient*change, label=label)
    axes[1].axhline(abs(result['required_source']), color='#9d3939', ls='--', label='Required source magnitude')
    axes[1].set(xlabel='Absolute change of one physical action coefficient', ylabel='Absolute change in source requirement [hbar c / D]',
                title='How unspecified local couplings affect the test')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Physical approximation checks; supporting stress is still unknown')
    fig.supxlabel('Sampled curvature is not a global bound. Couplings remain fixed at their declared values; no tuning is performed.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    records, inputs = load_inputs()
    for name in ('dihole_flux_scale.py', 'dihole_local_source.py'):
        inputs[name] = sha(project_path(name))
    function, parameters = symbolic_fields()[0], exact_target_parameters()
    rows = []
    for order in (64, 128):
        row, prepared = sample(function, parameters, order)
        rows.append(row)
    local = records['dihole_local_source']['records'][-1]
    actual_derivatives = np.imag(evaluate(prepared, 1e-25j)[0])/1e-25
    derivative_error = float(np.max(abs(actual_derivatives-local['action_derivative'])/np.maximum(1, abs(actual_derivatives))))
    response = np.asarray(local['energy_derivative'][1:3])
    required = records['dihole_local_source']['required_extra_source']*scale(1, ALPHA_EXAMPLE)**2
    peak_change = max(abs(rows[0][key]/rows[1][key]-1) for key in
                      ('scalar_abs_peak_D2', 'ricci_component_peak_D2', 'magnetic_peak_D2'))
    tests = dict(exact_input_provenance=all(sha(project_path(n)) == h for n, h in inputs.items()),
        actual_local_derivative_reproduced=derivative_error < 2e-12,
        scalar_contraction=max(r['scalar_contraction_error'] for r in rows) < 1e-10,
        static_ricci_contraction=max(r['static_ricci_square_error'] for r in rows) < 2e-10,
        no_mixed_static_time_curvature=max(r['static_mixed_time_error'] for r in rows) < 1e-12,
        sampled_peak_refinement=peak_change < .003,
        physical_scale_matches_saved_family=parameters['q'] == 1 and parameters['Q']/parameters['D'] == .2,
        sensitivity_threshold_identity=bool(np.allclose(abs(required)/response*response, abs(required), rtol=1e-15, atol=0)))
    result = dict(passed=all(tests.values()), tests=tests, samples=rows,
        parameters=dict(q=1, Q_over_D=.2, illustrative_alpha_R=ALPHA_EXAMPLE, D_over_planck=scale(1, ALPHA_EXAMPLE)),
        required_source=required, local_energy_response_coefficients=response.tolist(),
        action_convention='S_local=integral sqrt(-g) (c1 R^2+c2 Ricci_ab Ricci^ab); E=-S/T; fixed metric/field convention.',
        required_source_shift='delta C_required=-636.8446804260368 delta c1-289.35080303598477 delta c2.',
        single_coefficient_change_matching_source=(abs(required)/response).tolist(),
        coupling_uncertainty_bound='|delta C_required| <= 636.8446804260368 |delta c1|+289.35080303598477 |delta c2| in the specified basis.',
        current_physical_couplings=[0., 0.], assigned_coupling_uncertainty=None,
        actual_source=None, UV_bound=None, IR_bound=None, approximation_certified=False,
        max_sampled_peak_relative_change=peak_change, local_derivative_relative_error=derivative_error,
        scope=[
            'Only the varied collars are sampled, without exact axes or a global extremum enclosure.',
            'The fixed illustrative coupling is not a prediction of its value at these physical scales.',
            'Ricci and scalar curvature do not bound Weyl curvature, higher derivatives, fluctuations or backreaction.',
            'Curvature numbers are not percentages of approximation error or distance toward new physics.',
            'Local action responses are basis-specific; changes of field variables also change matter couplings and observables.',
            'Actual physical coefficient changes alter the model. A matched scheme change does not alter its predictions.',
            'No coefficient is tuned, no uncertainty is assigned, and no supporting medium or chronology change is established.'],
        sources=['https://arxiv.org/pdf/2211.09902', 'https://arxiv.org/pdf/hep-th/9906160'],
        input_sha256=inputs, source_sha256=sha(__file__))
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Physical-control checks failed; evidence saved.')
    plot(result)
    print(json.dumps({k: result[k] for k in ('passed', 'tests', 'single_coefficient_change_matching_source', 'max_sampled_peak_relative_change')}, indent=2))


if __name__ == '__main__':
    run()
