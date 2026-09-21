"""Leading constant-field heat response on the actual compact paired collars.

This is a local approximation, with no inhomogeneous/global remainder bound.
The fixed physical subtraction length is distinguished from the changing D.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from mpmath import mp
from mpmath.ctx_iv import MPIntervalContext

from balanced_dihole import symbolic_fields
from dihole_flux_scale import load_inputs
from dihole_local_source import prepare
from dihole_pv_counterterms import MU_H_D0
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME = 'dihole_strong_field_response'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def response_functions(x):
    """Return g(x)=(x*coth x)' and its tower part, without large-x overflow."""
    x = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError('Finite nonnegative proper-time times field required.')
    g, tower = np.ones_like(x), np.zeros_like(x)
    small, middle = x < .1, (x >= .1) & (x < 300.)
    z = x[small]
    g[small] = z*(2/3+z*z*(-4/45+z*z*(4/315+z*z*(-8/4725+z*z*4/18711))))
    tower[small] = g[small]-1
    z = x[middle]
    e = np.exp(-2*z)
    tower[middle] = 2*e/(1-e)-4*z*e/(1-e)**2
    g[middle] = 1+tower[middle]
    # ponytail: x>=300 omits a tower smaller than binary64 relative accuracy;
    # an exponentially resolved tower needs a logarithmic representation.
    return g, tower


def fields(prepared, parameters):
    g = np.exp(2*prepared['logs'][0])
    potential = prepared['potential']*parameters['physical_A_factor']
    radial = potential[1]**2/(g[1]*g[3])
    angular = potential[2]**2/(g[2]*g[3])
    magnetic = np.sqrt(radial+angular)
    if np.any(magnetic <= 0):
        raise ArithmeticError('Sampled magnetic zero needs separate treatment.')
    delta = prepared['mask'][0, 0]*radial/magnetic
    return magnetic, delta, radial, angular


def heat_response(y, magnetic, delta, measure):
    """deltaK(s;q)/q^2 at y=s*q/D^2, with q>0 and B stored for q=1."""
    if not math.isfinite(y) or y <= 0:
        raise ValueError('Positive finite scaled proper time required.')
    g, tower = response_functions(y*magnetic)
    weight = measure*delta/(4*math.pi**2*y)
    return np.array([weight@g, np.sum(weight), weight@tower])


def local_response_per_q2(q, mu_D, magnetic, delta, measure, log_glaisher):
    if not all(math.isfinite(v) and v > 0 for v in (q, mu_D)):
        raise ValueError('Positive finite flux magnitude and mu*D required.')
    bracket = 12*log_glaisher-.5-math.log(2)-math.log(q)-np.log(magnetic)+2*math.log(mu_D)
    return float(np.dot(measure*magnetic*delta, bracket)/(12*math.pi**2))


def constant_control():
    """Independent proper-time quadrature versus the zeta constant, 65 digits."""
    with mp.workdps(65):
        tiny = mp.mpf('.001')

        def g(x):
            return mp.coth(x)-x/mp.sinh(x)**2

        def sub(x):
            if abs(x) < tiny:
                return -4*x/45+4*x**3/315-8*x**5/4725+4*x**7/18711
            return g(x)/x**2-2/(3*x)

        log_a = mp.mpf(1)/12-mp.diff(mp.zeta, -1)
        constant = mp.quad(sub, [0, tiny, 1])+mp.quad(lambda x: g(x)/x**2, [1, 5, mp.inf])
        exact = 8*log_a-mp.mpf(1)/3-mp.mpf(2)/3*(mp.log(2)+mp.euler)
        cuts = []
        for cut in map(mp.mpf, ('.1', '1', '10')):
            integral = mp.quad(sub, [0, tiny, cut])+mp.quad(lambda x: g(x)/x**2, [cut, max(cut, 5), mp.inf])
            cuts.append(float(abs(integral+mp.mpf(2)/3*mp.log(cut)-exact)))
        errors = []
        for x in map(mp.mpf, ('1e-8', '.001', '.0099', '.01', '.0999', '.1', '1', '10', '299', '300')):
            actual = response_functions(float(x))[0].item()
            errors.append(float(abs(mp.mpf(actual)-g(x))/max(abs(g(x)), mp.mpf('1e-100'))))
        return dict(log_glaisher=float(log_a), integral_constant=mp.nstr(constant, 55),
            zeta_constant=mp.nstr(exact, 55), integral_error=float(abs(constant-exact)),
            cutoff_identity_errors=cuts, response_function_max_relative_error=max(errors))


def ideal_loop_bound(q, geometry, required):
    """Conditional fixed-holonomy circle bound, not an exact global LLL bound."""
    enclosure = geometry['interval_zero_exclusion']['records'][-1]
    if not enclosure['passed']:
        raise ArithmeticError('Closed collar enclosure did not pass.')
    iv = MPIntervalContext(); iv.dps = 65
    sigma = iv.mpf(enclosure['Sigma_lower'])
    rmax, dmax, width = iv.mpf(1)/20, iv.mpf(21)/20, iv.mpf(1)/20
    length = 2*width*sigma**2/(rmax*dmax)**2
    lower = float(np.nextafter(float(length.a), -np.inf))
    per_channel = iv.pi/(3*iv.mpf(lower))
    upper = float(np.nextafter(float(per_channel.b), np.inf))
    errors = []
    for f in (0., .3, 1.):
        for radial_fraction in (0., .2, 1.):
            eps = 1e-25j
            varied = np.sqrt(radial_fraction*np.exp(-2*eps*f)+1-radial_fraction)*np.exp(-eps*f)
            errors.append(abs(varied.imag/eps.imag+f*(1+radial_fraction)))
    return dict(two_full_collar_optical_length_lower_over_D=lower,
        negative_Casimir_response_upper_magnitude_per_channel=upper,
        nominal_q_channels_upper_magnitude_per_q2=upper/q,
        ratio_to_required_magnitude=upper/q/abs(required),
        optical_variation_max_error=float(max(errors)),
        formula='deltaL=-J; J=integral f*(1+radial_tangent_fraction) dl/N <=2L; Q_Casimir>=-pi*N_closed/(3*Lmin).',
        assumptions=['Every counted loop crosses both full prescribed collars.',
            'At most |q| independent mouth channels, each a massless complex 1+1 Dirac circle with fixed holonomy.',
            'Only the topology-dependent flat-circle Casimir term; no curved conformal anomaly or projector/mixing contribution included.'],
        scope='The optical length lower bound is geometric. Applying the Casimir bound to the exact global lowest Landau sector remains unproved.')


def run():
    records, inputs = load_inputs()
    atom_path = project_path('atomic_aperture_budget_results.json')
    atom = json.loads(atom_path.read_text(encoding='utf-8'))
    if not atom['passed'] or sha(project_path('atomic_aperture_budget.py')) != atom['source_sha256']:
        raise ArithmeticError('Atomic reference provenance failed.')
    for filename, digest in atom['input_sha256'].items():
        if sha(project_path(filename)) != digest:
            raise ArithmeticError('Atomic input changed: '+filename)
        inputs[filename] = digest
    inputs[atom_path.name] = sha(atom_path)
    for name in ('atomic_aperture_budget.py', 'dihole_flux_scale.py', 'dihole_pv_counterterms.py'):
        inputs[name] = sha(project_path(name))
    control = constant_control()
    log_a = control['log_glaisher']
    q_atom = atom['geometry']['continuous_flux_requirement_at_illustrative_alpha']
    function, parameters = symbolic_fields()[0], exact_target_parameters()
    grid = np.geomspace(1e-6, 100., 161)
    rows = []
    for order in (64, 128):
        prepared = prepare(function, parameters, order)
        magnetic, delta, radial, angular = fields(prepared, parameters)
        measure = prepared['measure']
        curves = np.array([heat_response(float(y), magnetic, delta, measure) for y in grid])
        log_coefficient = float(measure@(magnetic*delta)/(12*math.pi**2))
        base = local_response_per_q2(1., MU_H_D0, magnetic, delta, measure, log_a)
        atom_fixed = local_response_per_q2(q_atom, MU_H_D0*q_atom, magnetic, delta, measure, log_a)
        atom_scaled = local_response_per_q2(q_atom, MU_H_D0, magnetic, delta, measure, log_a)
        rows.append(dict(order=order, sampled_B_min_D2=float(magnetic.min()), sampled_B_max_D2=float(magnetic.max()),
            integral_deltaB_D=float(measure@delta), integral_B_deltaB_over_D=float(measure@(magnetic*delta)),
            magnetic_a4_derivative=log_coefficient*32*math.pi**2,
            log_q_coefficient=log_coefficient, q1_local_response=base,
            atomic_fixed_physical_reference_response_per_q2=atom_fixed,
            atomic_fixed_muD_loop_component_per_q2=atom_scaled,
            atomic_compensating_Maxwell_change_per_q2=2*log_coefficient*math.log(q_atom),
            scaled_proper_time_y=grid.tolist(), heat_response_per_q2=curves.tolist()))
    reference = rows[-1]
    curve_error = float(np.max(abs(np.array(rows[0]['heat_response_per_q2'])[:, 0]/curves[:, 0]-1)))
    eps = 1e-25j
    varied_b = np.sqrt(radial*np.exp(2*eps*prepared['mask'][0, 0])+angular)
    metric_errors = []
    for y in (.03, .3, 3.):
        x = y*varied_b
        varied_trace = measure@(x/np.tanh(x))/(4*math.pi**2*y*y)
        expected = heat_response(y, magnetic, delta, measure)[0]
        metric_errors.append(float(abs(varied_trace.imag/eps.imag/expected-1)))
    saved_a4 = -records['dihole_local_source']['records'][-1]['charge_included_F2_energy_derivative']*2/3
    physical = json.loads(project_path('dihole_physical_control_results.json').read_text())
    if not physical['passed'] or sha(project_path('dihole_physical_control.py')) != physical['source_sha256']:
        raise ArithmeticError('Physical source requirement provenance failed.')
    for filename, digest in physical['input_sha256'].items():
        if sha(project_path(filename)) != digest:
            raise ArithmeticError('Physical reference input changed: '+filename)
        inputs[filename] = digest
    source_per_q2 = physical['required_source']
    inputs['dihole_physical_control_results.json'] = sha(project_path('dihole_physical_control_results.json'))
    inputs['dihole_physical_control.py'] = sha(project_path('dihole_physical_control.py'))
    geometry_path = project_path('dihole_strong_field_geometry_results.json')
    geometry = json.loads(geometry_path.read_text(encoding='utf-8'))
    if not geometry['passed'] or sha(project_path('dihole_strong_field_geometry.py')) != geometry['source_sha256']:
        raise ArithmeticError('Strong-field geometry provenance failed.')
    for filename, digest in geometry['input_sha256'].items():
        if sha(project_path(filename)) != digest:
            raise ArithmeticError('Strong-field geometry input changed: '+filename)
        inputs[filename] = digest
    inputs[geometry_path.name] = sha(geometry_path)
    inputs['dihole_strong_field_geometry.py'] = sha(project_path('dihole_strong_field_geometry.py'))
    ideal = ideal_loop_bound(q_atom, geometry, source_per_q2)
    scale_error = abs(reference['atomic_fixed_muD_loop_component_per_q2']
        +reference['atomic_compensating_Maxwell_change_per_q2']
        -reference['atomic_fixed_physical_reference_response_per_q2'])
    functions = response_functions(np.geomspace(1e-8, 1e4, 300))
    tests = dict(exact_current_inputs=all(sha(project_path(n)) == h for n, h in inputs.items()),
        independent_proper_time_constant=control['integral_error'] < 1e-30,
        proper_time_cutoff_independence=max(control['cutoff_identity_errors']) < 1e-30,
        stable_response_function=control['response_function_max_relative_error'] < 2e-11,
        actual_paired_metric_variation=max(metric_errors) < 3e-13,
        magnetic_heat_coefficient=abs(reference['magnetic_a4_derivative']/saved_a4-1) < 3e-13,
        collar_quadrature_sensitivity=curve_error < 1e-8,
        heat_response_signs=bool(np.all(functions[0] > 0) and np.all(functions[0] <= 1) and np.all(functions[1] <= 0)),
        matched_Maxwell_scale_change=scale_error < 1e-20,
        fixed_path_optical_variation=ideal['optical_variation_max_error'] < 1e-15,
        conditional_loop_length_positive=ideal['two_full_collar_optical_length_lower_over_D'] > 0,
        fixed_physical_reference_scaling=abs(reference['atomic_fixed_physical_reference_response_per_q2']
            -reference['q1_local_response']-reference['log_q_coefficient']*math.log(q_atom)) < 1e-20)
    result = dict(passed=all(tests.values()), tests=tests, samples=rows, constants=control,
        atom_flux_continuous_estimate=q_atom, q1_reference_muD=MU_H_D0,
        atomic_fixed_reference_muD=MU_H_D0*q_atom,
        required_source_per_q2=source_per_q2,
        local_to_required_magnitude_ratio=abs(reference['atomic_fixed_physical_reference_response_per_q2']/source_per_q2),
        actual_metric_variation_relative_errors=metric_errors, quadrature_max_relative_change=curve_error,
        scale_cancellation_absolute_error=scale_error,
        conditional_ideal_returning_loops=ideal,
        column_order=['full_local', 'lowest_Landau', 'gapped_tower'],
        units='D is current Weyl scale; heat response per Killing time in D^-1, energy response in hbar*c/D.',
        fixed_reference='D_ref is the q=1 physical scale. Hold physical mu=sqrt(2)*exp(-gamma)/D_ref and the declared coupling fixed; then mu*D=q*mu*D_ref.',
        local_formula='Q_local/q^2=integral dV B1 deltaB1 [12 log(A_Glaisher)-1/2-log(2*q*B1/(mu*D)^2)]/(12*pi^2).',
        scale_change='The fixed-muD loop number is a component in a different subtraction convention. The matched Maxwell term cancels its change; do not tune mu to obtain support.',
        approximation_remainder_bound=None, actual_full_source=None, atom_transit=None,
        limitations=['Leading constant-field density applied locally to the varying collar; no quantitative inhomogeneous or global remainder.',
            'The local lowest Landau term is not the exact global lowest-sector trace of the curved connected geometry.',
            'The atomic flux is a continuous size estimate, not a resolved integer-sector source calculation.',
            'Declared illustrative coupling at the fixed reference is not a measured high-energy or hydrogen interaction.',
            'Sampled magnetic positivity excludes no unsampled or exterior zero. Curvature, gradients and longitudinal physics remain omitted.',
            'Small local response cannot exclude full support without bounding the omitted contribution.',
            'No supporting material, stable opening, atomic survival or time difference is established.'],
        failed_pilot=dict(reason='Cancellation at x=.01 exceeded the response-function accuracy target; the x^9 series now covers x<.1.',
            source='provenance/dihole_strong_field_response_small_x_trial.py',
            source_sha256=sha(Path(__file__).resolve().parents[1]/'provenance/dihole_strong_field_response_small_x_trial.py'),
            result='results/dihole_strong_field_response_small_x_trial.json',
            result_sha256=sha(project_path('dihole_strong_field_response_small_x_trial.json'))),
        sources=['https://arxiv.org/pdf/0906.2430', 'https://arxiv.org/pdf/0810.4889', 'https://arxiv.org/pdf/1807.04726'],
        source_sha256=sha(__file__), input_sha256=inputs)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Strong-field controls failed; result preserved.')
    plot(result)
    print(json.dumps({k: result[k] for k in ('passed', 'tests', 'local_to_required_magnitude_ratio', 'quadrature_max_relative_change')}, indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.5), layout='constrained')
    row = result['samples'][-1]
    x, curves = np.array(row['scaled_proper_time_y']), np.array(row['heat_response_per_q2'])
    axes[0].loglog(x, curves[:, 0], label='Full leading local response', lw=2)
    axes[0].loglog(x, curves[:, 1], '--', label='Local lowest Landau term')
    axes[0].loglog(x, -curves[:, 2], ':', label='Magnitude of negative tower term')
    axes[0].set(xlabel='Scaled proper time y = s q / D squared', ylabel='Local heat variation / q squared', ylim=(1e-12, 1e-2),
                title='Both parts matter in the short-time cancellation')
    axes[0].legend(fontsize=8)
    values = [abs(row['atomic_fixed_physical_reference_response_per_q2']), abs(result['required_source_per_q2'])]
    axes[1].bar(['Computed local\napproximation (+)', 'Required supporting\nresponse (-)'], values, color=['#268a96', '#ab4f48'])
    axes[1].set(yscale='log', ylim=(1e-8, 100), ylabel='Magnitude per q squared [hbar c / D]',
                title='Atomic-scale comparison at one fixed physical reference')
    for i, v in enumerate(values):
        axes[1].text(i, v*1.4, f'{v:.3g}', ha='center')
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
    fig.suptitle('Strong-field local contribution calculated; full supporting source remains unknown')
    fig.supxlabel('Unbounded omitted curvature, gradient and global contributions prevent a support verdict. No atom transit simulated.', fontsize=8)
    for extension in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+extension), dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()
