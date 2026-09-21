"""Test the geometric truncation error without inventing a nonlinear solution.

Exponentiating/rescaling the saved linear fields is only a reconstruction:
it does not solve the nonlinear lapse, shift or evolution equations.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from project_paths import project_path


def run():
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    result_path = project_path('localized_expanding_packet_results.json')
    prior = json.loads(result_path.read_text(encoding='utf-8'))
    state_path = project_path('localized_expanding_packet_state.npz')
    if (not prior['passed'] or sha(state_path) != prior['state_sha256']
            or sha(project_path('localized_expanding_packet.py')) != prior['source_sha256']):
        raise ArithmeticError('Packet source/state provenance failed.')
    with np.load(state_path) as state:
        h = prior['parameters']['H_over_cutoff']
        radius = prior['summary']['initial_width_times_cutoff']
        ar = np.exp(state['tau'])[:, None]*radius*state['x'][None, :]
        rows = []
        for factor in (1., .1, .01, .001):
            z, n, v, beta_over_r, u = [factor*state[key] for key in
                ('zeta','lapse','velocity','shift_over_r','radial_gradient_times_r')]
            F = (v-beta_over_r*(1+u)-h*n)/(1+n)
            # Algebra avoids subtracting large background compactness terms.
            reconstructed = (ar*np.exp(z))**2*(2*h*F+F*F)-2*u-u*u
            linear = factor*state['compactness_excess']
            difference = reconstructed-linear
            peak = np.unravel_index(np.argmax(abs(difference)), difference.shape)
            reference = float(np.max(abs(linear)))
            remainder = float(np.max(abs(difference)))
            rows.append(dict(amplitude=prior['parameters']['epsilon']*factor,
                maximum_linear_compactness=reference, maximum_geometric_remainder=remainder,
                remainder_to_maximum_linear=remainder/reference,
                initial_remainder_to_maximum_initial_linear=float(np.max(abs(difference[0]))/np.max(abs(linear[0]))),
                peak_time_Ht=float(state['tau'][peak[0]]), peak_comoving_radius_over_width=float(state['x'][peak[1]]),
                minimum_areal_radius_derivative_factor=float(np.min(1+u)),
                scope='Exact geometric expression evaluated on reconstructed linear fields; not a nonlinear on-shell result.'))
    scaling = [rows[i]['maximum_geometric_remainder']/rows[i+1]['maximum_geometric_remainder'] for i in range(3)]
    checks = dict(quadratic_remainder_scaling=all(99. < x < 101. for x in scaling),
                  initial_amplitude_rejected=rows[0]['remainder_to_maximum_linear'] > 1,
                  weak_amplitude_geometric_control=rows[-1]['remainder_to_maximum_linear'] < .002,
                  no_reconstructed_spatial_throat=all(x['minimum_areal_radius_derivative_factor'] > 0 for x in rows))
    result = dict(complete=True, passed=all(checks.values()), checks=checks, amplitudes=rows,
        remainder_reduction_factors=scaling,
        decision='Reject epsilon=.001 as a quantitatively controlled linear compactness prediction. Lower amplitude suppresses this geometric remainder but is not a nonlinear solution or stronger source.',
        source_sha256=sha(__file__), baseline_sha256=sha(result_path), state_sha256=sha(state_path),
        scope='No bound on missing second-order matter/gravity fields is supplied; geometric convergence alone cannot validate their omission.')
    project_path('localized_packet_geometry_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Geometric truncation checks failed; evidence retained.')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    run()
