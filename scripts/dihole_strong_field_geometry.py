"""Local strong-field applicability diagnostics on the prescribed paired collars.

Sampled curvature/field-gradient ratios are numerical, not remainder bounds.
A separate interval calculation excludes magnetic zeros on this closed collar
rectangle only. The full exterior and massless lowest sector are not controlled.
"""
import hashlib
import json
import math
import os
from pathlib import Path
from time import perf_counter

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_key, '1')

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from dihole_uv_bound import (curvature, diagonal_metric_jets, prepare_third,
                            symbolic_fields, exact_target_parameters, gauss)
from dihole_stretched_collar import radius_profile, log_speed_derivatives
from project_paths import project_path


NAME = 'dihole_strong_field_geometry'
FAILURE_SOURCE = 'provenance/dihole_strong_field_geometry_relative_minimum_failure.py'
FAILURE_SHA256 = '7ad41d3ab43785e1c761f617f6d2610fc731856da6a4509ed00c31298d7a272b'
FROZEN = {
    'dihole_uv_bound.py': '19802f8cc95d67e9bd127254ec98d2fcfb8ce863f1acdc17312f14cb5f4dffd0',
    'dihole_uv_bound_results.json': 'bf36064e5b95401faf1abc2cf0eb5d39a9946bdd822ed7e5765c797d86cfba6c',
    'atomic_aperture_budget.py': '1b38240e8051fdbd6bf264093627a01ad3933f021ef3d50aa240bc2e714ca7ab',
    'atomic_aperture_budget_results.json': 'c1ea3dca1ef9e1a0f896434f7f0d56f426347d4521f658d2080d1ccda53c6f33',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs():
    hashes, records = {}, {}
    if sha(Path(__file__).resolve().parents[1]/FAILURE_SOURCE) != FAILURE_SHA256:
        raise ArithmeticError('Reconstructed failure provenance mismatch.')
    for filename, expected in FROZEN.items():
        path = project_path(filename)
        actual = sha(path)
        if actual != expected:
            raise ArithmeticError('Frozen input mismatch: '+filename)
        hashes[filename] = actual
        if filename.endswith('.json'):
            record = json.loads(path.read_text(encoding='utf-8'))
            if record.get('passed') is not True:
                raise ArithmeticError('Input calculation has not passed: '+filename)
            for dependency, digest in record.get('input_sha256', {}).items():
                name = Path(dependency).name
                current = sha(project_path(name))
                if current != digest:
                    raise ArithmeticError('Current dependency mismatch: '+name)
                hashes[name] = current
            records[filename] = record
    return records, hashes


def field_derivatives(function, p, radius, mu):
    """Pole-safe A_R and A_mu from the actual global potential derivatives."""
    radius, mu = np.broadcast_arrays(radius, mu)
    local = log_speed_derivatives(radius, mu, p)
    values = np.asarray(function(local['global_r'], local['global_cosine'],
                                 p['M'], p['a'], p['B']))
    distance = np.sqrt(radius**2+4*p['k']**2+4*p['k']*radius*mu)
    dr = (radius+2*p['k']*mu)/distance
    ar = p['physical_A_factor']*(values[5]*(1+dr)/2+values[8]*(dr-1)/(2*p['k']))
    am = p['physical_A_factor']*radius/distance*(p['k']*values[5]+values[8])
    return dict(A=values[2]*p['physical_A_factor'], radial=ar, mu=am,
                lapse=values[0], meridian=values[1], distance=distance)


def field_squared_from_potential(fields, radius, radial_rate, mu, epsilon=0., mask=0.):
    denominator = 1-np.asarray(mu)**2
    radial = np.divide(fields['radial']**2, denominator,
                       out=np.zeros_like(fields['radial']), where=denominator > 0)
    # The exact axis values A=0 or 1 are independent of R; hence A_R^2/(1-mu^2)
    # tends to zero at either pole. This limit is used only at mu=+-1.
    bracket = fields['mu']**2+np.asarray(radius)**2*np.asarray(radial_rate)**2*np.exp(2*epsilon*mask)*radial
    return (fields['lapse']/fields['meridian']/np.asarray(radius)**2)**2*bracket


def interval_null_exclusion(radial_panels=16, angular_panels=32, precision=70):
    """Exact-calibration interval cover, independent of the floating samples.

    MPIntervalContext performs directed interval arithmetic. All tile endpoints
    and calibration parameters are rational values or interval square roots.
    This proves a sign statement, not the sampled maxima elsewhere in the file.
    """
    iv = MPIntervalContext(); iv.dps = precision
    k, charge = iv.mpf(1)/2, iv.mpf(1)/5
    mass = (iv.sqrt(k*k+4*charge*k)-k)/2
    ag = iv.sqrt(k*k-mass*mass)
    field = mass/((k+ag)*(mass+k))
    factor = -1/(2*charge)
    am_lower, am_upper, lower_b, lower_v, lower_sigma = math.inf, -math.inf, math.inf, math.inf, math.inf
    failed = 0
    for i in range(radial_panels):
        lo = iv.mpf(1)/40+iv.mpf(i)/(40*radial_panels)
        hi = iv.mpf(1)/40+iv.mpf(i+1)/(40*radial_panels)
        R = iv.mpf([lo.a, hi.b])
        for j in range(angular_panels):
            lo = -1+iv.mpf(2*j)/angular_panels
            hi = -1+iv.mpf(2*(j+1))/angular_panels
            mu = iv.mpf([lo.a, hi.b])
            distance = iv.sqrt(R*R+4*k*k+4*k*R*mu)
            r, c = mass+(R+distance)/2, (distance-R)/(2*k)
            s = 1-c*c
            delta = r*r-2*mass*r-ag*ag
            w = (r*r-ag*ag)**2+delta*ag*ag*s
            wr = 4*r*(r*r-ag*ag)+(2*r-2*mass)*ag*ag*s
            wc = -2*delta*ag*ag*c
            h = 2*mass*r*ag+field*w/2
            hr, hc = 2*mass*ag+field*wr/2, field*wc/2
            z = ag*ag+2*field*mass*r*ag+field*field*w/4
            v = delta+s*z
            vr = 2*r-2*mass+s*(2*field*mass*ag+field*field*wr/4)
            vc = -2*c*z+s*field*field*wc/4
            ar = -s*(hr*v-h*vr)/(v*v)
            ac = -((-2*c*h+s*hc)*v-s*h*vc)/(v*v)
            am = factor*R/distance*(k*ar+ac)
            sigma = r*r-ag*ag*c*c
            good = bool(am.b < 0 and v.a > 0 and sigma.a > 0)
            failed += int(not good)
            am_lower = min(am_lower, np.nextafter(float(am.a), -np.inf))
            am_upper = max(am_upper, np.nextafter(float(am.b), np.inf))
            lower_v = min(lower_v, np.nextafter(float(v.a), -np.inf))
            lower_sigma = min(lower_sigma, np.nextafter(float(sigma.a), -np.inf))
            if good:
                # h_geometry=R*distance and N/P=h_geometry^2/Sigma^2 exactly.
                angular_b = distance**2*(-am)/sigma**2
                lower_b = min(lower_b, np.nextafter(float(angular_b.a), -np.inf))
    return dict(radial_panels=radial_panels, angular_panels=angular_panels,
        decimal_precision=precision, failed_tiles=failed,
        A_mu_enclosure=[float(am_lower), float(am_upper)],
        potential_denominator_lower=float(lower_v), Sigma_lower=float(lower_sigma),
        q1_B_D0_squared_lower=float(lower_b) if failed == 0 else None,
        passed=bool(failed == 0 and lower_b > 0))


def location(index, u, mu):
    return dict(u_over_D0=float(u[index]), cos_theta=float(mu[index]),
                theta=float(np.arccos(mu[index])))


def summary(values, u, mu, weights):
    normal = weights/np.sum(weights)
    return dict(minimum=float(np.min(values)), maximum=float(np.max(values)),
        minimum_location=location(np.argmin(values), u, mu),
        maximum_location=location(np.argmax(values), u, mu),
        mask_volume_mean=float(normal@values), mask_volume_RMS=float(np.sqrt(normal@(values**2))))


def geometry_sample(function, parameters, prepared, atomic_q, epsilon=0., chunk=256):
    n = prepared['order']; nodes, _ = gauss(n)
    u, mu = np.broadcast_arrays(.025*(nodes[:, None]+1), nodes[None, :])
    u, mu = u.ravel(), mu.ravel()
    number = len(u)
    B, riem_norm = np.empty(number), np.empty(number)
    grad_by_direction = np.empty((4, number))
    signs = np.array([1., -1., 0., 0.])
    for lo in range(0, number, chunk):
        at = slice(lo, lo+chunk)
        logs = prepared['logs'][:, :, at]+epsilon*prepared['mask'][:, None, at]*signs[None, :, None]
        third = prepared['third'][:, :, at].copy()
        third[0] += epsilon*signs[:, None]*prepared['mask3'][None, at]
        tensors = curvature(*diagonal_metric_jets(logs, third), prepared['potential'][:, at])
        B[at] = np.sqrt(np.sum(tensors['magnetic']**2, axis=(0, 1))/2)
        riem_norm[at] = np.sqrt(np.sum(tensors['riemann']**2, axis=(0, 1, 2, 3)))
        grad_by_direction[:, at] = np.sum(tensors['grad_magnetic']**2, axis=(1, 2))
    radius, rate, _ = radius_profile(u, order=128)
    field = field_derivatives(function, parameters, radius, mu)
    independent = np.sqrt(field_squared_from_potential(field, radius, rate, mu,
                                                     epsilon, prepared['mask'][0]))
    curvature_ratio = riem_norm/B
    gradient_ratio = np.sum(grad_by_direction, axis=0)/B**3
    weights = prepared['measure']*prepared['mask'][0]
    directions = grad_by_direction@weights
    directions /= np.sum(directions)
    q1 = dict(B_times_D0_squared=summary(B, u, mu, weights),
              Riemann_norm_over_B=summary(curvature_ratio, u, mu, weights),
              covariant_gradient_F_squared_over_B_cubed=summary(gradient_ratio, u, mu, weights))
    atom = dict(B_times_D0_squared=summary(atomic_q*B, u, mu, weights),
                Riemann_norm_over_B=summary(curvature_ratio/atomic_q, u, mu, weights),
                covariant_gradient_F_squared_over_B_cubed=summary(gradient_ratio/atomic_q, u, mu, weights))
    return dict(order=n, paired_metric_epsilon=epsilon, q1=q1,
        nominal_atomic_continuous_flux=atom, mask_volume=float(np.sum(weights)),
        weighted_covariant_gradient_direction_fractions=directions.tolist(),
        independent_B_relative_error=float(np.max(abs(B-independent)/B)))


def self_checks(function, p):
    logs = np.zeros((6, 4, 1)); third = np.zeros((4, 4, 1))
    rate = 1/3; logs[1, 0, 0] = rate
    potential = np.zeros((6, 1)); potential[1, 0] = 1
    t = curvature(*diagonal_metric_jets(logs, third), potential)
    # Euclidean ds2=exp(2x/3)dt2+dx2+dy2+dz2, F_xz=1 at x=0.
    time_gradient = float(t['grad_magnetic'][0, 0, 3, 0])
    riem = float(np.sqrt(np.sum(t['riemann']**2)))
    R = np.array([.025, .031, .042, .05])[:, None]
    mu = np.array([-1., -.4, 0., .4, 1.])[None, :]
    fields = field_derivatives(function, p, R, mu)
    tiny = 1e-25
    def potential_value(radius, cosine):
        distance = np.sqrt(radius**2+4*p['k']**2+4*p['k']*radius*cosine)
        rb, c = p['M']+(radius+distance)/2, (distance-radius)/(2*p['k'])
        return np.asarray(function(rb, c, p['M'], p['a'], p['B'])[2])*p['physical_A_factor']
    derivative = np.imag(potential_value(R, mu+tiny*1j))/tiny
    error = float(np.max(abs(derivative-fields['mu'])))
    axis_error = float(max(np.max(abs(fields['A'][:, 0]-1)), np.max(abs(fields['A'][:, -1]))))
    return dict(static_covariant_time_gradient=time_gradient, expected_time_gradient=rate,
        static_curvature_norm=riem, expected_curvature_norm=2*rate**2,
        pole_safe_derivative_complex_step_error=error, exact_axis_potential_error=axis_error,
        passed=bool(abs(time_gradient-rate) < 2e-14 and abs(riem-2*rate**2) < 2e-14
                    and error < 2e-11 and axis_error < 2e-11))


def closure_field_grid(function, p, radial_count=65, angular_count=257):
    u, mu = np.broadcast_arrays(np.linspace(0, .05, radial_count)[:, None],
                               np.linspace(-1, 1, angular_count)[None, :])
    radius, rate, _ = radius_profile(u, order=128)
    fields = field_derivatives(function, p, radius, mu)
    B = np.sqrt(field_squared_from_potential(fields, radius, rate, mu))
    at = np.argmin(B)
    return dict(radial_count=radial_count, angular_count=angular_count,
        q1_B_D0_squared_sampled_minimum=float(B.ravel()[at]),
        minimum_location=location(at, u.ravel(), mu.ravel()),
        sampled_A_mu_range=[float(np.min(fields['mu'])), float(np.max(fields['mu']))],
        note='Includes u=0,.05 and both polar limits. The minimum remains sampled, not enclosed.')


def run():
    start = perf_counter(); records, hashes = inputs()
    atom = records['atomic_aperture_budget_results.json']['geometry']
    atomic_q = atom['continuous_flux_requirement_at_illustrative_alpha']
    if not math.isfinite(atomic_q) or atomic_q <= 0:
        raise ValueError('Atomic flux estimate must be finite and positive.')
    function, _ = symbolic_fields(); p = exact_target_parameters()
    checks = self_checks(function, p)
    samples = []
    for order in (32, 64, 128):
        prepared = prepare_third(function, p, order)
        samples.append(geometry_sample(function, p, prepared, atomic_q))
        if order == 64:
            varied = [geometry_sample(function, p, prepared, atomic_q, e) for e in (-.001, .001)]
    metrics = ('B_times_D0_squared', 'Riemann_norm_over_B',
               'covariant_gradient_F_squared_over_B_cubed')
    sensitivity = {}
    for metric in metrics:
        coarse, fine = samples[-2]['q1'][metric], samples[-1]['q1'][metric]
        sensitivity[metric] = {key: abs(coarse[key]-fine[key])/abs(fine[key])
                               for key in ('minimum', 'maximum', 'mask_volume_mean', 'mask_volume_RMS')}
    intervals = [interval_null_exclusion(8, 16), interval_null_exclusion(16, 32)]
    closure = closure_field_grid(function, p)
    all_extrema_relative = max(r[key] for r in sensitivity.values() for key in ('minimum', 'maximum'))
    applicable_extrema_relative = max([r['maximum'] for r in sensitivity.values()]
                                      +[sensitivity['B_times_D0_squared']['minimum']])
    tests = dict(exact_current_input_hashes=True, independent_static_and_field_controls=checks['passed'],
        field_tensor_matches_potential=bool(max(s['independent_B_relative_error'] for s in samples+varied) < 2e-11),
        all_covariant_direction_fractions_normalized=bool(max(abs(sum(s['weighted_covariant_gradient_direction_fractions'])-1) for s in samples+varied) < 2e-14),
        sampled_maxima_and_positive_B_minimum_sensitivity=bool(applicable_extrema_relative < .01),
        mask_integral_sensitivity=bool(max(r[key] for r in sensitivity.values() for key in ('mask_volume_mean', 'mask_volume_RMS')) < 2e-5),
        closed_collar_magnetic_zero_exclusion=bool(all(r['passed'] for r in intervals)),
        closure_grid_respects_separate_lower_bound=bool(closure['q1_B_D0_squared_sampled_minimum'] >= intervals[-1]['q1_B_D0_squared_lower']))
    result = dict(passed=bool(all(tests.values())), tests=tests, source_sha256=sha(__file__),
        input_sha256=hashes, parameters=dict(geometry=p, q_reference=1,
            nominal_atomic_continuous_flux=atomic_q, nominal_atomic_D_m=atom['required_D_m'],
            nominal_atomic_areal_diameter_m=records['atomic_aperture_budget_results.json']['target']['nominal_throat_areal_diameter_m']),
        equations=dict(B='sqrt(F_charge_ab F_charge_ab/2), Euclidean orthonormal components',
            curvature='sqrt(R_abcd R_abcd)/B', gradient='(nabla_a F_bc)(nabla_a F_bc)/B^3; a,b,c=0..3',
            charge_scaling='B(q)=|q| B(1); both ratios scale as 1/|q| on the fixed geometric shape.',
            paired_variation='g_tt(epsilon)=exp(2 epsilon f)g_tt; g_uu(epsilon)=exp(-2 epsilon f)g_uu; angular metric and charge one-form fixed.',
            weight='Normalized actual four-volume-per-Killing-time measure times f>=0, including both collars and full azimuth.',
            field_identity='B^2=(N0/P/R^2)^2[A_mu^2+R^2 Rprime^2 exp(2 epsilon f) A_R^2/(1-mu^2)].',
            zero_exclusion='A_mu=factor*(R/distance)*(k A_r+A_c)<0 on R in[1/40,1/20],mu in[-1,1]; B>=distance^2 |A_mu|/Sigma^2.',
            pole_limits='A(mu=-1)=1,A(mu=+1)=0 independent of R; A_R^2/(1-mu^2)->0 at both poles.'),
        controls=checks, samples=samples, paired_variation_samples=varied,
        sensitivity_64_to_128=sensitivity, closed_collar_field_samples=closure,
        initial_relative_minimum_guard=dict(passed=bool(all_extrema_relative < .01),
            maximum_relative_difference=all_extrema_relative, threshold=.01,
            archived_source=FAILURE_SOURCE, archived_source_sha256=FAILURE_SHA256,
            archive_origin='Reconstructed by reversing the recorded sole correction patch; it was not copied before editing.',
            failed_JSON_was_written=False,
            measured_gradient_minimum_64=samples[-2]['q1']['covariant_gradient_F_squared_over_B_cubed']['minimum'],
            measured_gradient_minimum_128=samples[-1]['q1']['covariant_gradient_F_squared_over_B_cubed']['minimum'],
            explanation='The initial all-extrema relative guard fails at the gradient-ratio minimum tending to zero near the poles. It is retained here; no positive gradient minimum or relative convergence of that minimum is claimed. The acceptance gate checks the three upper indicators and positive B minimum.'),
        interval_zero_exclusion=dict(records=intervals,
            scope='Closed collar rectangle only, exact linked balanced constants and q!=0. This sign enclosure is independent of the floating diagnostic extrema.',
            calibration='k=1/2,Q=1/5,M=(sqrt(k^2+4Qk)-k)/2,a=sqrt(k^2-M^2),Bext=M/((k+a)(M+k)).',
            regularity='N=V/Sigma>0, P=N Sigma^2/(R distance)^2>0; axes are regular bundle patches. Paired variations leave the nonzero angular field component intact.'),
        unavailable=dict(covariant_second_field_derivative_bound=None,
                         inhomogeneous_heat_remainder_bound=None, full_source_bound=None),
        limitations=[
            'Curvature/gradient extrema and weighted diagnostics are sampled numerical values with quadrature sensitivity, not certified maxima or remainder bounds.',
            'The near-zero sampled gradient minimum fails a relative-error guard; it is not used as a nonzero lower bound or convergence claim.',
            'The interval result excludes field zeros only on the closed prescribed collars; it supplies no blanket exterior control.',
            'Small atomic-flux ratios are local applicability indicators, not a proof of a controlled strong-field approximation, source sign, support or atomic survival.',
            'The saved atomic flux is a continuous estimate, not a selected exact integer sector. Changing flux changes the Dirac operator.',
            'The q=1 reference is not atomic size. The fixed s0=.01 weak-field heat truncation is not justified for the nominal atomic flux.',
            'No covariant second derivative of F, higher inhomogeneous remainder or global massless lowest-sector contribution is evaluated.'],
        elapsed_seconds=perf_counter()-start)
    if not result['passed']:
        raise ArithmeticError(json.dumps(tests))
    return result


if __name__ == '__main__':
    result = run(); destination = project_path(NAME+'_results.json')
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], checks=len(result['tests']),
        fine_reference=result['samples'][-1]['q1'], interval=result['interval_zero_exclusion']['records'][-1],
        source_sha256=result['source_sha256'], results_sha256=sha(destination),
        elapsed_seconds=result['elapsed_seconds']), allow_nan=False))
