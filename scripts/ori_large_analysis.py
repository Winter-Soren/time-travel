"""Compare completed finite Ori evolutions and plot their measured uncertainty.

Run after the large-grid cases documented in RESEARCH.md. No onset extrapolation.
"""
from project_paths import project_path
import json
import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(project_path('.cache')/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ori_evolution import coordinate_grid, diagnostics, self_check


def main(smooth=False, strong=False, strong_mid=False, strong_late=False):
    strong = strong or strong_mid or strong_late
    period = 'late' if strong_late else ('mid' if strong_mid else 'short')
    stem = f'ori_smooth_strong_{period}' if strong else ('ori_smooth' if smooth else 'ori_large')
    label = stem if strong else ('smooth_weak' if smooth else 'polynomial_weak')
    fine_stem = f'{stem}_129' if strong else ('ori_smooth_weak_129' if smooth else 'ori_large_polynomial_129')
    coefficients = ('ori_smooth_accepted_coefficients.npz' if strong else
                    ('ori_smooth_weak_coefficients.npz' if smooth else 'ori_large_constraint_coefficients.npz'))
    end = -.1 if strong_late else (-.5 if strong_mid else (-.8 if strong else -.1))
    result_path = project_path(f'{stem}_evolution_results.json')
    runs = json.loads(result_path.read_text()) if result_path.exists() else {}
    if strong:
        controls = json.loads((project_path(f'{stem}_controls.json')).read_text())
        runs.update({f'{label}_{name}': record for name, record in controls.items()})
    else:
        runs.update(json.loads((project_path(f'{stem}_guarded_results.json')).read_text()))
    runs[f'{label}_129'] = json.loads((project_path(f'{fine_stem}.json')).read_text())
    runs[f'{label}_linear_129'] = json.loads((project_path(f'{stem}_linear_129.json')).read_text())
    files = {f'{label}_{n}': f'{label}_{n}.npz' for n in (33, 65)}
    files.update({f'{label}_129': f'{fine_stem}.npz',
                  f'{label}_linear_65': f'{label}_linear_65.npz',
                  f'{label}_linear_129': f'{stem}_linear_129.npz',
                  f'{label}_half_dt_65': f'{label}_half_dt_65.npz'})
    states, masks = {}, {}
    for name, filename in files.items():
        assert Path(runs[name]['initial_data_file']).name == coefficients, (name, 'wrong initial-data family')
        assert (runs[name]['half_width'], runs[name]['stretch'], runs[name]['interior_radius']) == (2.1, 2.5, .2)
        saved = np.load(project_path(filename))
        assert float(saved['t']) == end and abs(runs[name]['samples'][-1]['t']-end) < 1e-12
        assert runs[name]['strong'] == strong
        states[name] = saved['state']
        n = runs[name]['n']
        grid, h = coordinate_grid(n, 2.1, 2.5)
        x, y = np.meshgrid(grid, grid, indexing='ij')
        masks[n] = x*x+y*y <= .2**2
        updated = diagnostics(states[name], h, masks[n])
        for key, value in updated.items():
            old = runs[name]['samples'][-1].get(key)
            if old is not None and value is not None:
                assert abs(old-value) < 1e-10*(1+abs(old)), (name, key, old, value)
        runs[name]['samples'][-1].update(updated)
        runs[name]['endpoint_diagnostics_recomputed'] = True
    differences = []
    for coarse, fine in ((33, 65), (65, 129)):
        differences.append(float(np.max(abs(states[f'{label}_{coarse}'][masks[coarse]]-
            states[f'{label}_{fine}'][::2, ::2][masks[coarse]]))))
    boundary = {str(n): float(np.max(abs(states[f'{label}_{n}'][masks[n]]-
        states[f'{label}_linear_{n}'][masks[n]]))) for n in (65, 129)}
    time_error = float(np.max(abs(states[f'{label}_65'][masks[65]]-
        states[f'{label}_half_dt_65'][masks[65]])))
    finals = [runs[f'{label}_{n}']['samples'][-1] for n in (33, 65, 129)]
    constraint_ratios = {key: [finals[i][key]/finals[i+1][key] for i in (0, 1)]
                         for key in ('max_H', 'max_M', 'max_gauge')}
    comparisons = dict(nested_state_differences=differences,
        state_refinement_ratio=differences[0]/differences[1],
        observed_state_order=float(np.log2(differences[0]/differences[1])),
        boundary_state_differences=boundary, time_step_difference=time_error,
        boundary_endpoint_diagnostic_changes={str(n): {key: abs(
            runs[f'{label}_{n}']['samples'][-1][key]-
            runs[f'{label}_linear_{n}']['samples'][-1][key])
            for key in ('max_H', 'max_M', 'max_gauge', 'max_einstein_component_residual',
                        'center_R', 'center_Ricci_squared', 'center_regular_tidal_norm',
                        'center_freefall_energy', 'center_freefall_tidal_norm')}
            for n in (65, 129)},
        constraint_refinement_ratios=constraint_ratios,
        ricci_squared_source_disagreements=[abs(f['center_Ricci_squared']-
            f['center_source_Ricci_squared']) for f in finals],
        sampled_boundary_margins={str(n): runs[f'{label}_{n}']['nearest_wall_to_core']-
            runs[f'{label}_{n}']['sampled_cartesian_travel_envelope'] for n in (33, 65, 129)})
    checks = dict(state_refines=differences[1] < differences[0],
        constraints_refine=all(v > 1 for pair in constraint_ratios.values() for v in pair),
        fine_boundary_error_subdominant=boundary['129'] < differences[1]/10,
        time_error_subdominant=time_error < differences[1]/10,
        configured_stage_CFL_observed=all(runs[name]['stage_cfl_enforced'] and
            runs[name]['max_stage_courant'] <= runs[name]['cfl']*(1+1e-12) for name in files),
        positive_sampled_spatial_metric=all(s['min_spatial_eigenvalue'] > 0
            for name in files for s in runs[name]['samples']))
    runs['algebraic_checks'] = self_check()
    runs['input_coefficients_sha256'] = hashlib.sha256((project_path(coefficients)).read_bytes()).hexdigest()
    runs['initial_data_family'] = label
    samples = runs[f'{label}_129']['samples']
    # Older histories recorded the regular tidal norm only at the recomputed endpoint.
    scale_samples = [s for s in samples if s.get('center_regular_tidal_norm') is not None
                     and s.get('center_regular_energy') is not None]
    c, G, separation, budget = 299792458., 6.67430e-11, 2., 1.
    tidal = max(s['center_regular_tidal_norm'] for s in scale_samples)
    density = max(s['center_regular_energy'] for s in scale_samples)
    sufficient_scale = float(np.sqrt(c*c*separation*tidal/budget))
    energy_density = c**4*density/(8*np.pi*G*sufficient_scale**2)
    runs['local_scale_illustration'] = dict(separation_m=separation,
        relative_acceleration_budget_m_s2=budget, max_sampled_center_tidal_norm=tidal,
        max_sampled_center_regular_energy=density,
        included_sample_times=[s['t'] for s in scale_samples],
        all_diagnostic_times_included=len(scale_samples) == len(samples),
        necessary_length_scale_m=sufficient_scale/3**.25,
        sufficient_length_scale_m=sufficient_scale,
        energy_density_at_sufficient_scale_J_m3=float(energy_density),
        energy_equivalent_density_kg_m3=float(energy_density/c**2),
        scope='Sampled center, instantaneous regular frame, all separation directions, linear tidal term only; '
              'length scale is not an apparatus radius; no trajectory or human-safety bound.')
    runs['comparisons'], runs['checks'] = comparisons, checks
    runs['scope'] = ('Finite r<=0.2 interior evolution of one radius-three initial-data family '
        f'from t=-1 to {end}; sampled numerical checks, no chronology extension or stability theorem.')
    result_path.write_text(json.dumps(runs, indent=2)+'\n')
    print(json.dumps(dict(checks=checks, comparisons=comparisons, finest=finals[-1]), indent=2))

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7), layout='constrained')
    times = [s['t'] for s in samples]
    energy_traces = [('center_regular_energy', 'Prescribed observer', '#167b86', '-'),
                     ('center_energy', 'Slice-normal observer', '#7753a6', '--')]
    if all('center_freefall_energy' in s for s in samples):
        energy_traces.append(('center_freefall_energy', 'Free fall, p = -0.5', '#bd6a23', '-.'))
    for key, legend_label, color, style in energy_traces:
        axes[0].plot(times, [s[key] for s in samples], style, color=color, marker='o', ms=3, label=legend_label)
    axes[0].set(title='Central energy along the computed interval', xlabel='Model time t', ylabel='Energy density (model units)')
    axes[0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    for key, legend_label, color in [('max_H', 'Hamiltonian', '#167b86'), ('max_M', 'Momentum', '#bd6a23'),
                             ('max_gauge', 'Gauge', '#7753a6')]:
        axes[1].loglog([32, 64, 128], [f[key] for f in finals], 'o-', color=color, label=legend_label)
    axes[1].set(title='Equation errors under refinement', xlabel='Grid intervals per direction', ylabel=f'Maximum core residual at t={end}')
    axes[1].set_xticks([32, 64, 128], ['32', '64', '128'])
    axes[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[1].legend(frameon=False, fontsize=9)
    names = ['Grid 33 → 65', 'Grid 65 → 129', 'Wall choice, 65', 'Wall choice, 129', 'Half time step, 65']
    values = differences+[boundary['65'], boundary['129'], time_error]
    axes[2].scatter(values, np.arange(len(names)), s=65,
        color=['#899ca8', '#566b7d', '#69a7af', '#167b86', '#7753a6'])
    axes[2].set_yticks(np.arange(len(names)), names)
    axes[2].set_xscale('log')
    axes[2].invert_yaxis()
    axes[2].set(title='Sensitivity of the final state', xlabel='Maximum absolute core state difference')
    for ax in axes:
        ax.grid(alpha=.18)
        ax.set_axisbelow(True)
    family_title = 'smooth strong pulse' if strong else label.replace('_', ' ')
    fig.suptitle(f'Coupled Ori–scalar model: {family_title} finite evolution checks\n'
                 'Same initial data in every comparison · stops before t=0 · does not establish a stable time-travel passage', fontsize=12)
    for extension in ('png', 'pdf'):
        fig.savefig(project_path(f'{stem}_evolution.{extension}'), dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smooth', action='store_true', help='Compare the separate infinitely smooth pulse family')
    parser.add_argument('--strong', action='store_true', help='Compare the short strong smooth pulse study')
    parser.add_argument('--strong-mid', action='store_true', help='Compare the strong smooth pulse extension to t=-0.5')
    parser.add_argument('--strong-late', action='store_true', help='Compare the strong smooth pulse extension to t=-0.1')
    main(**vars(parser.parse_args()))
