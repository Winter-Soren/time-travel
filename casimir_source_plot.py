"""Plot the saved source checks; does not rerun or modify the calculations."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parent
    result = json.loads((root / 'casimir_cavity_results.json').read_text())
    spatial, noise = result['spatial_source'], result['temporal_density_noise']
    floor = noise['general_normalized_sampling_bound']['relative_standard_deviation_lower_bound']
    assert spatial['center']['maximum_componentwise_error_over_A'] < 1e-12
    assert all(row['relative_standard_deviation'] >= floor - 1e-12 for row in noise['samples'])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
    fig.subplots_adjust(left=.085, right=.975, top=.76, bottom=.24, wspace=.33)
    fig.suptitle('Testing the proposed cavity source', fontsize=21, fontweight='bold', y=.96)
    fig.text(.5, .87, 'Ideal closed box • fixed quantum state • center mean matches the chosen target',
             ha='center', color='#475569')
    for name, color in zip(('x', 'y', 'z'), ('#1565c0', '#c56c00', '#008477')):
        rows = spatial['positive_axis_samples'][name]
        axes[0].loglog([100 * max(row['offset_L']) for row in rows],
                      [100 * row['maximum_componentwise_error_over_A'] for row in rows],
                      'o-', label=f'{name} axis', color=color, lw=1.7, markersize=5)
    axes[0].axhline(1, color='#9f3345', ls='--', lw=1.2, label='1% comparison level')
    axes[0].set(title='Spatial mismatch', xlabel='Distance from center (% of L)',
                ylabel='Largest component error / A (%)', xlim=(.085, 6), ylim=(.003, 30))
    axes[0].set_xticks([.1, .5, 1, 2, 5], ['0.1', '0.5', '1', '2', '5'])
    axes[0].legend(fontsize=9, loc='upper left', frameon=False)
    rows = noise['samples']
    axes[1].semilogx([row['sigma'] for row in rows],
                    [row['relative_standard_deviation'] for row in rows],
                    'o-', color='#7951b0', lw=1.8, label='Gaussian time averaging')
    axes[1].axhline(floor, color='#9f3345', ls='--', lw=1.3, label=f'Lower bound: {floor:.2f}×')
    axes[1].set(title='Time averaging leaves a noise floor', xlabel='Averaging width σ (in L/c)',
                ylabel='Density standard deviation / |mean|', xlim=(.2, 40), ylim=(10, 48))
    axes[1].set_xticks([.25, 1, 4, 16, 32], ['0.25', '1', '4', '16', '32'])
    axes[1].legend(fontsize=9, frameon=False, loc='upper right')
    for axis in axes:
        axis.grid(alpha=.18, which='major')
    fig.text(.085, .13, 'Left: sampled axes only; constant flat target; lines guide the eye, not a 3D error guarantee.',
             fontsize=9, color='#475569')
    fig.text(.085, .09, 'Right: this state and center observable; the fluctuation ratio alone does not determine gravitational stability.',
             fontsize=9, color='#475569')
    fig.text(.085, .04, '13 September 2027 → 12 September 2027: no demonstrated route.',
             fontsize=10, fontweight='bold', color='#9f3345')
    for suffix in ('png', 'pdf'):
        fig.savefig(root / f'casimir_source_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root / 'casimir_source_progress.png')


if __name__ == '__main__':
    main()
