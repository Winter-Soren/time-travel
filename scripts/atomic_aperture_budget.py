"""One-hydrogen-atom target: size, flux and an instantaneous tidal diagnostic.

No transit or state survival is simulated. The atom is momentarily at rest
relative to the static frame; a moving trajectory needs its own Fermi frame.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from balanced_dihole import symbolic_fields
from dihole_flux_scale import ALPHA_EXAMPLE, PLANCK_LENGTH_M, implied_alpha, load_inputs, scale
from dihole_local_source import invariants, metric_jets, prepare
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME = 'atomic_aperture_budget'
BOHR_M = 5.29177210544e-11
ELECTRON_KG = 9.1093837139e-31
PROTON_KG = 1.67262192595e-27
HARTREE_J = 4.3597447222060e-18
C_M_S = 299792458.
MASS_REDUCED = ELECTRON_KG*PROTON_KG/(ELECTRON_KG+PROTON_KG)
HYDROGEN_RADIUS_M = BOHR_M*ELECTRON_KG/MASS_REDUCED
BINDING_J = HARTREE_J*MASS_REDUCED/ELECTRON_KG/2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tidal_frame(prepared):
    """R_0i0j = Hess_h(N)_ij / N in the orthonormal static frame."""
    values = metric_jets(prepared)
    metric = values[0, 1:4]
    dg = np.zeros((3, 3, metric.shape[1]))
    dg[0], dg[1] = values[1, 1:4], values[2, 1:4]
    logs = prepared['logs'][:, 0]
    dn = np.stack((logs[1], logs[2], np.zeros_like(logs[1])))
    hessian = np.einsum('in,jn->ijn', dn, dn)
    hessian[0, 0] += logs[3]
    hessian[0, 1] += logs[4]
    hessian[1, 0] += logs[4]
    hessian[1, 1] += logs[5]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                connection = ((dg[i, k] if k == j else 0)+(dg[j, k] if k == i else 0)
                              -(dg[k, i] if i == j else 0))/(2*metric[k])
                hessian[i, j] -= connection*dn[k]
    tidal = hessian/np.sqrt(metric[:, None]*metric[None, :])
    ricci = invariants(values)[1]
    trace_error = np.max(abs(np.trace(tidal, axis1=0, axis2=1)-ricci[0, 0]/(-values[0, 0])))
    return tidal, float(trace_error)


def outside_probability(radius):
    x = radius/HYDROGEN_RADIUS_M
    return math.exp(-2*x)*(1+2*x+2*x*x)


def run():
    records, inputs = load_inputs()
    for name in ('dihole_flux_scale.py', 'dihole_local_source.py'):
        inputs[name] = sha(project_path(name))
    geometry = records['balanced_dihole']
    areal_ratio = next(r['area_radius'] for r in geometry['cut_summaries']
                       if r['epsilon'] == .025 and r['angular_order'] == 96)
    diameter = 1e-9
    d_m = diameter/(2*areal_ratio)
    d_q1 = scale(1, ALPHA_EXAMPLE)*PLANCK_LENGTH_M
    function, parameters = symbolic_fields()[0], exact_target_parameters()
    rows = []
    for order in (64, 128):
        prepared = prepare(function, parameters, order)
        tidal, error = tidal_frame(prepared)
        eigenvalues = np.linalg.eigvalsh(tidal.transpose(2, 0, 1))
        trace = eigenvalues.sum(axis=1)
        absolute_sum = abs(eigenvalues).sum(axis=1)
        factor = MASS_REDUCED*C_M_S**2*HYDROGEN_RADIUS_M**2/(2*d_m**2*BINDING_J)
        rows.append(dict(order=order, trace_identity_error=error,
            eigenvalue_abs_peak_D2=float(np.max(abs(eigenvalues))),
            trace_abs_peak_D2=float(np.max(abs(trace))),
            absolute_eigenvalue_sum_peak_D2=float(np.max(absolute_sum)),
            hydrogen_1s_first_order_shift_over_binding_peak=float(np.max(abs(trace))*factor),
            hydrogen_1s_absolute_expectation_bound_over_binding_peak=float(np.max(absolute_sum)*factor),
            profile_u_over_D=prepared['u'].tolist(),
            profile_shift_ratio=(np.max(abs(trace).reshape(order, order), axis=1)*factor).tolist()))
    ratio = rows[-1]['hydrogen_1s_first_order_shift_over_binding_peak']
    refinement = abs(rows[0]['trace_abs_peak_D2']/rows[1]['trace_abs_peak_D2']-1)
    tests = dict(exact_input_hashes=all(sha(project_path(n)) == h for n, h in inputs.items()),
        tidal_trace_matches_independent_ricci=max(r['trace_identity_error'] for r in rows) < 1e-10,
        sampled_tidal_refinement=refinement < .003,
        hydrogen_probability_normalization=outside_probability(0.) == 1.,
        diameter_uses_areal_radius=math.isclose(2*areal_ratio*d_m, diameter, rel_tol=1e-15),
        reduced_mass_less_than_electron_mass=0 < MASS_REDUCED < ELECTRON_KG,
        atomic_target_exceeds_q1_geometry=diameter > 2*areal_ratio*d_q1,
        energy_indicator_is_not_survival_probability=ratio > 1.)
    result = dict(passed=all(tests.values()), tests=tests,
        target=dict(payload='One neutral hydrogen atom; isolated nonrelativistic 1s state as input reference.',
            nominal_throat_areal_diameter_m=diameter, human_transport='Deferred by user.',
            acceptance=['Self-consistent source and current', 'Stability during transit',
                'Propagated atomic center and internal state with controlled errors',
                'Outside-clock comparison establishing the specified time shift']),
        geometry=dict(central_areal_radius_over_D=areal_ratio, required_D_m=d_m,
            continuous_flux_requirement_at_illustrative_alpha=d_m/d_q1,
            coupling_if_q_stays_one=implied_alpha(1, d_m/PLANCK_LENGTH_M),
            q1_illustrative_areal_diameter_m=2*areal_ratio*d_q1,
            flux_scope='Continuous estimate only. Allowed flux is integer; the full source must be recomputed in that sector.'),
        atom=dict(reduced_mass_kg=MASS_REDUCED, hydrogen_bohr_radius_m=HYDROGEN_RADIUS_M,
            binding_energy_J=BINDING_J, flat_1s_probability_outside_nominal_radius=outside_probability(diameter/2)),
        samples=rows, sampled_refinement=refinement,
        diameter_for_unit_sampled_shift_ratio_m=diameter*math.sqrt(ratio),
        diameter_for_one_percent_sampled_shift_ratio_m=diameter*math.sqrt(ratio/.01),
        tidal_model='H_tidal=(mu c^2/2) R_0i0j x_i x_j; <x_i x_j>_1s=a_H^2 delta_ij.',
        interpretation='The 1nm candidate does not permit a small-tidal-perturbation treatment of an initially ordinary 1s atom at all sampled collar points. No ionization or survival probability follows.',
        caveats=['Nominal areal diameter is a geometric target, not a hard aperture or proof of wavepacket fitting.',
            'Flat 1s probability omits center-of-mass localization and all supporting-field interactions.',
            'The static-frame tidal tensor applies instantaneously to a comoving geodesic; it is not the history of a transit.',
            'The leading local expansion and first-order atomic perturbation are diagnostics only when the predicted shift is large.',
            'Larger diameters reduce this indicator at fixed shape, but change flux and require a new supporting-field calculation.',
            'The electromagnetic or hidden nature of the supporting U(1) and its coupling to hydrogen remain unspecified.',
            'No safe minimum diameter, usable medium, atom survival or time travel is established.'],
        support=None, transit=None, atomic_survival=None, backward_passage=False,
        backward_scope='The current prescribed static geometry has a global time function; no backward causal passage in this metric.',
        constants=dict(source='https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
            retrieved_on='2026-09-22', treatment='CODATA2022 central values; no metrology enclosure.'),
        sources=['https://arxiv.org/pdf/gr-qc/0605052', 'https://arxiv.org/pdf/hep-th/9906160'],
        source_sha256=sha(__file__), input_sha256=inputs)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Atomic aperture checks failed; evidence retained.')
    plot(result)
    print(json.dumps({k: result[k] for k in ('passed', 'tests', 'geometry', 'atom', 'diameter_for_one_percent_sampled_shift_ratio_m')}, indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4), layout='constrained')
    profile = result['samples'][-1]
    axes[0].plot(profile['profile_u_over_D'], profile['profile_shift_ratio'], color='#a44840')
    axes[0].axhline(1., ls='--', color='#555', label='Binding-energy scale')
    axes[0].set(xlabel='Outward collar coordinate u / D', ylabel='Leading tidal shift / hydrogen binding energy',
                title='One atom, nominal 1nm throat diameter')
    axes[0].legend(fontsize=8)
    diameter_nm = np.geomspace(1., 1000., 200)
    ratio = profile['hydrogen_1s_first_order_shift_over_binding_peak']/diameter_nm**2
    axes[1].loglog(diameter_nm, ratio, color='#167b93')
    axes[1].axhline(.01, ls='--', color='#555', label='1% diagnostic target; not a safety bound')
    axes[1].set(xlabel='Nominal areal diameter [nm]', ylabel='Peak sampled tidal shift / binding energy',
                title='Smaller payload does not guarantee a gentler passage')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Atomic transport target: geometry and integrity must both be tested')
    fig.supxlabel('No atom transit is simulated. Curves use a local instantaneous tidal diagnostic; the supporting field remains unknown.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    run()
