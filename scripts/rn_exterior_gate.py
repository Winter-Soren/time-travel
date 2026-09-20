"""An isolated-RN compatibility gate and a saved-profile stress diagnostic.

No collar, Maxwell, Dirac, or other field ODE is solved. The classical
magnetic-radius relation b=r_e is an explicit assumption of this gate.
"""
from project_paths import project_path, source_matches
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp
import mpmath as mp


def symbolic_checks():
    x = sp.Symbol('x', real=True)
    b = sp.Symbol('b', positive=True)
    r, n = sp.Function('r')(x), sp.Function('N')(x)
    rp, rpp = sp.diff(r, x), sp.diff(r, x, 2)
    u = sp.diff(n, x)/n
    rho = (1-rp**2-2*r*rpp)/r**2
    pr = (-1+rp**2+2*r*rp*u)/r**2
    pt = sp.diff(n, x, 2)/n+rp*u/r+rpp/r
    ricci = -2*sp.diff(n, x, 2)/n-4*rp*u/r-4*rpp/r+2*(1-rp**2)/r**2
    monopole = b*b/r**4
    other = (rho-monopole, pr+monopole, pt-monopole)
    ward = sp.simplify(sp.diff(other[1], x)+u*(other[0]+other[1])+2*rp/r*(other[1]-other[2]))
    trace = sp.simplify(-other[0]+other[1]+2*other[2]+ricci)
    null = sp.simplify(other[0]+other[1]-2*(rp*u-rpp)/r)
    R, mu = sp.symbols('R mu', positive=True)
    F = 1-2*mu/R+b*b/R**2
    FR, FRR = sp.diff(F, R), sp.diff(F, R, 2)
    rn = ((1-F-R*FR)/R**2-b*b/R**4,
          (-1+F+R*FR)/R**2+b*b/R**4,
          FRR/2+FR/R-b*b/R**4)
    boundary_identity = sp.simplify(
        (rp*u-rpp)/(n*r)+sp.diff(rp/(n*r), x)+rp**2/(n*r**2))
    return dict(
        radial_Ward_identity_zero=ward == 0,
        trace_is_minus_Ricci_scalar=trace == 0,
        radial_null_identity_zero=null == 0,
        RN_Einstein_minus_monopole_zero=all(sp.simplify(v) == 0 for v in rn),
        integrated_null_identity_zero=boundary_identity == 0,
        formulas=dict(
            scaled_energy_density='8*pi*G*rho=(1-r_prime^2-2*r*r_second)/r^2',
            scaled_radial_pressure='8*pi*G*p_r=(-1+r_prime^2+2*r*r_prime*N_prime/N)/r^2',
            scaled_tangential_pressure='8*pi*G*p_t=N_second/N+r_prime*N_prime/(r*N)+r_second/r',
            monopole='8*pi*G*(rho_EM,p_r_EM,p_t_EM)=(b^2/r^4,-b^2/r^4,b^2/r^4)',
            Ward='p_r_prime+(N_prime/N)*(rho+p_r)+2*(r_prime/r)*(p_r-p_t)=0',
            trace='-rho+p_r+2*p_t=-Ricci_scalar/(8*pi*G)',
            null='4*pi*G*Tgeom_kk/E^2=(r_prime*N_prime/N-r_second)/(N^2*r)',
            integrated_null='Q=4*pi*G*I/E=-[r_prime/(N*r)]_in^out-integral(r_prime^2/(N*r^2),dx)'))


def rn_table(b, outer_radius, delta_radius, width):
    maximum_mass = (outer_radius**2+b**2)/(2*outer_radius)
    margin = maximum_mass-b
    rows = []
    for fraction in (0., .25, .5, .75):
        mu = b+fraction*margin
        F = 1-2*mu/outer_radius+b*b/outer_radius**2
        v = math.sqrt(F)
        horizon_offset = math.sqrt(max(0., mu*mu-b*b))
        rplus, rminus = mu+horizon_offset, mu-horizon_offset
        FR = 2*mu/outer_radius**2-2*b*b/outer_radius**3
        rows.append(dict(mass_excess_fraction=fraction, mass_parameter_mu=mu,
                         rplus=rplus, rminus=rminus, outer_horizon_clearance=outer_radius-rplus,
                         F_outer=F, lapse_C1=v, radius_proper_derivative=v,
                         log_lapse_proper_derivative=FR/(2*v),
                         angular_extrinsic_curvature=v/outer_radius,
                         strict_minimum_nonovershooting_width=delta_radius/v,
                         width_ratio_to_saved=delta_radius/(v*width),
                         saved_width_compatible=width > delta_radius/v,
                         C_required_to_match_saved_outer_lapse=1/v,
                         hypothetical_collar_Q_boundary_term_C1=-1/outer_radius,
                         rn_lapse_derivative_identity_error=abs(F-(1-b/outer_radius)**2*(1-fraction))))
    return dict(classical_magnetic_radius=b, outer_radius=outer_radius,
                minimum_mass=b, strict_upper_mass=maximum_mass, available_mass_margin=margin,
                allowed_interval='b <= mu < (Rout^2+b^2)/(2*Rout), assuming Rout>b and mu>=b.',
                explanation='For mu>=b, rminus<=b<Rout. Thus F(Rout)>0 selects the exterior region beyond rplus, not the region inside rminus. Equality at the upper mass puts the horizon on the cut and is excluded.',
                extremal_maximum_outer_speed=1-b/outer_radius,
                universal_strict_minimum_nonovershooting_width=delta_radius/(1-b/outer_radius),
                rows=rows)


def tensor_from_saved(profile, b):
    r, rp, rpp = [np.asarray(profile[k]) for k in ('radius', 'radius_prime', 'radius_second')]
    n, np_, npp = [np.asarray(profile[k]) for k in ('lapse', 'lapse_prime', 'lapse_second')]
    u = np_/n
    up = npp/n-u*u
    rho = (1-rp*rp-2*r*rpp)/(r*r)
    pr = (-1+rp*rp+2*r*rp*u)/(r*r)
    pt = npp/n+rp*u/r+rpp/r
    em = b*b/r**4
    other_rho, other_pr, other_pt = rho-em, pr+em, pt-em
    # Differentiate the explicit radial-pressure formula; no profile differencing.
    other_pr_prime = (2*rp*rpp/r**2-2*rp*(-1+rp*rp)/r**3
                      +2*(rpp*u+rp*up)/r-2*rp*rp*u/r**2-4*b*b*rp/r**5)
    ward_terms = [other_pr_prime, u*(other_rho+other_pr), 2*rp/r*(other_pr-other_pt)]
    ward = sum(ward_terms)
    ward_scale = 1+sum(abs(a) for a in ward_terms)
    ricci = -2*npp/n-4*rp*u/r-4*rpp/r+2*(1-rp*rp)/r**2
    trace = -other_rho+other_pr+2*other_pt
    trace_error = trace+ricci
    trace_scale = 1+abs(other_rho)+abs(other_pr)+2*abs(other_pt)+abs(ricci)
    radial_null = other_rho+other_pr
    direct_null = 2*(rp*u-rpp)/r
    return dict(scaled_Einstein_rho=rho.tolist(), scaled_Einstein_pr=pr.tolist(),
                scaled_Einstein_pt=pt.tolist(), scaled_monopole_rho=em.tolist(),
                scaled_other_rho=other_rho.tolist(), scaled_other_pr=other_pr.tolist(),
                scaled_other_pt=other_pt.tolist(), scaled_other_radial_null=radial_null.tolist(),
                scaled_other_trace=trace.tolist(), Ricci_scalar=ricci.tolist(),
                scaled_other_pr_prime=other_pr_prime.tolist(),
                pointwise_scaled_Ward_residual=ward.tolist(),
                pointwise_normalized_Ward_residual=(ward/ward_scale).tolist(),
                pointwise_normalized_trace_residual=(trace_error/trace_scale).tolist(),
                maximum_Ward_absolute_error=float(np.max(abs(ward))),
                maximum_Ward_normalized_error=float(np.max(abs(ward)/ward_scale)),
                maximum_trace_absolute_error=float(np.max(abs(trace_error))),
                maximum_trace_normalized_error=float(np.max(abs(trace_error)/trace_scale)),
                maximum_radial_null_absolute_error=float(np.max(abs(radial_null-direct_null))))


def high_precision_tensor_controls(profile, b):
    """Independent expanded arithmetic on the same exact binary input values."""
    maximum_ward = maximum_trace = mp.mpf(0)
    with mp.workdps(75):
        charge = mp.mpf(b)**2
        arrays = [profile[k] for k in ('radius', 'radius_prime', 'radius_second',
                                       'lapse', 'lapse_prime', 'lapse_second')]
        for values in zip(*arrays):
            r, rp, rpp, n, np_, npp = map(mp.mpf, values)
            rho = 1/r**2-rp**2/r**2-2*rpp/r-charge/r**4
            pr = -1/r**2+rp**2/r**2+2*rp*np_/(r*n)+charge/r**4
            pt = npp/n+rp*np_/(r*n)+rpp/r-charge/r**4
            prp = (2*rp/r**3+2*rp*rpp/r**2-2*rp**3/r**3
                   +2*rpp*np_/(r*n)+2*rp*npp/(r*n)-2*rp*np_**2/(r*n*n)
                   -2*rp*rp*np_/(r*r*n)-4*charge*rp/r**5)
            terms = (prp, np_*(rho+pr)/n, 2*rp*(pr-pt)/r)
            ward = abs(sum(terms))/(1+sum(abs(z) for z in terms))
            ricci = -2*npp/n-4*rp*np_/(r*n)-4*rpp/r+2/r**2-2*rp**2/r**2
            trace = abs(-rho+pr+2*pt+ricci)/(1+abs(rho)+abs(pr)+2*abs(pt)+abs(ricci))
            maximum_ward, maximum_trace = max(maximum_ward, ward), max(maximum_trace, trace)
        return dict(decimal_digits=75, evaluated_points=len(arrays[0]),
                    maximum_normalized_Ward_residual=float(maximum_ward),
                    maximum_normalized_trace_residual=float(maximum_trace),
                    purpose='Distinguishes cancellation roundoff in float checks from a formula defect; not an enclosure of the saved physical profiles.')


def run():
    inputs = ('smooth_magnetic_collar.py', 'smooth_magnetic_collar_results.json')
    hashes = {name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in inputs}
    saved = json.loads(project_path(inputs[1]).read_text(encoding='utf-8'))
    expected = dict(q=1., d=1., core_radius=.2, collar_proper_width=.1, exterior_cut_radius=.25, rho_c=32.)
    if not saved['passed'] or any(saved['parameters'].get(k) != v for k, v in expected.items()):
        raise ValueError('This gate requires the successful fixed-parameter increment108 profiles.')
    if not source_matches(inputs[0], saved['source_sha256'], hashes[inputs[0]]):
        raise ValueError('Saved108 source identity does not match the current source.')
    b, outer_radius, width = .2, .25, .1
    symbolic = symbolic_checks()
    branch = rn_table(b, outer_radius, outer_radius-b, width)
    x = np.asarray(saved['profiles']['x'])
    profiles = saved['profiles']['geometry']
    tensor = tensor_from_saved(profiles, b)
    high_precision = high_precision_tensor_controls(profiles, b)
    samples = []
    for fraction in (0., .25, .5, .75, 1.):
        index = int(np.argmin(abs(x-fraction*width)))
        samples.append(dict(x=float(x[index]), width_fraction=fraction,
            radius=profiles['radius'][index], lapse=profiles['lapse'][index],
            scaled_other_rho=tensor['scaled_other_rho'][index],
            scaled_other_pr=tensor['scaled_other_pr'][index],
            scaled_other_pt=tensor['scaled_other_pt'][index],
            scaled_other_radial_null=tensor['scaled_other_radial_null'][index],
            scaled_other_trace=tensor['scaled_other_trace'][index]))
    rn_errors = []
    for row in branch['rows']:
        R, mu, v = outer_radius, row['mass_parameter_mu'], row['radius_proper_derivative']
        FR, FRR = 2*mu/R**2-2*b*b/R**3, -4*mu/R**3+6*b*b/R**4
        exact_rn = dict(radius=[R], radius_prime=[v], radius_second=[FR/2],
                        lapse=[v], lapse_prime=[FR/2], lapse_second=[FRR*v/2])
        check = tensor_from_saved(exact_rn, b)
        rn_errors.append(max(abs(check[key][0]) for key in ('scaled_other_rho', 'scaled_other_pr', 'scaled_other_pt')))
    target_outer = b*b/outer_radius**4*np.array([-1., 1., -1.])
    actual_outer = np.array([tensor[k][-1] for k in ('scaled_other_rho', 'scaled_other_pr', 'scaled_other_pt')])
    actual_inner = np.array([tensor[k][0] for k in ('scaled_other_rho', 'scaled_other_pr', 'scaled_other_pt')])
    # Identity is checked directly against the saved108 independently differentiated density.
    radial_density = np.asarray(tensor['scaled_other_radial_null'])/(2*np.asarray(profiles['lapse']))
    saved_density = np.asarray(saved['profiles']['radial_null_densities']['Qgeom_density'])
    tests = dict(
        generic_symbolic_identities=all(v for k, v in symbolic.items() if k != 'formulas'),
        horizon_interval=all(row['mass_parameter_mu'] >= b and row['outer_horizon_clearance'] > 0 for row in branch['rows']) and abs(branch['strict_upper_mass']-.205) < 2e-15,
        endpoint_mass_dependence=max(row['rn_lapse_derivative_identity_error'] for row in branch['rows']) < 3e-15,
        width_gate_rejects_saved_profile=all(row['strict_minimum_nonovershooting_width'] > width for row in branch['rows']) and abs(branch['universal_strict_minimum_nonovershooting_width']-.25) < 2e-15,
        exact_RN_tensor_control=max(rn_errors) < 3e-13,
        # Subtraction of almost identical monopole/core terms costs relative
        # accuracy near the inner endpoint; retain the errors and check the
        # same input independently at75 digits instead of using a 3e-15 guard.
        pointwise_Ward=tensor['maximum_Ward_normalized_error'] < 256*np.finfo(float).eps and high_precision['maximum_normalized_Ward_residual'] < 1e-60,
        pointwise_trace=tensor['maximum_trace_normalized_error'] < 256*np.finfo(float).eps and high_precision['maximum_normalized_trace_residual'] < 1e-60,
        inner_product_geometry_matches_monopole=float(np.max(abs(actual_inner))) < 3e-13,
        flat_outer_requires_nonzero_other_stress=float(np.max(abs(actual_outer-target_outer))) < 3e-13,
        saved_radial_null_density_consistency=float(np.max(abs(radial_density-saved_density))) < 3e-13,
        input_preservation=all(hashlib.sha256(project_path(name).read_bytes()).hexdigest() == digest for name, digest in hashes.items()))
    tests = {k: bool(v) for k, v in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests,
        parameters=expected, assumption='This gate assumes b=r_e and r_e^2=pi*G*q^2/g_R^2. mu=G*Mphysical is a length; it is not the magnetic-null budget M of increment108.',
        RN_branch=branch, symbolic_controls=symbolic,
        endpoint_derivation=dict(F='1-2*mu/R+b^2/R^2', lapse='C*sqrt(F)',
            radius_proper_derivative='sqrt(F)', log_lapse_proper_derivative='F_R/(2*sqrt(F))',
            angular_extrinsic_curvature='sqrt(F)/R', asymptotic_clock='C=1 in the compact table. C is a positive constant for general derivations.',
            clock_does_not_repair_geometry='Changing C changes lapse normalization but leaves r_prime, N_prime/N, and angular extrinsic curvature unchanged. Saved108 outer r_prime=1 cannot match any allowed branch r_prime<=.2.'),
        width_gate=dict(saved_width=width, radius_increase=outer_radius-b,
            hypothesis='A nonovershooting opening: r_prime(0)=0, 0<=r_prime(x)<=vout, r_prime(w)=vout, with continuous derivative.',
            strict_inequality='DeltaR=integral r_prime dx < w*vout, hence w>DeltaR/vout.',
            all_allowed_black_hole_RN_branches_reject_saved_width=True,
            limitation='This is not a no-go for all collars. Within the black-hole branch, a radius-monotone collar of width .1 could evade this bound only by overshooting r_prime>.5 somewhere before ending at r_prime<=.2. No such alternative geometry is constructed here.'),
        below_extremal_mass_escape=dict(
            assumption='Formal RN exterior with mu<b outside an excised/replaced core; no black-hole horizon branch is assumed here.',
            required_strict_outer_speed=(outer_radius-b)/width,
            strict_upper_mass_for_nonovershooting_saved_width=(outer_radius**2+b*b-outer_radius**2*((outer_radius-b)/width)**2)/(2*outer_radius),
            strict_minimum_mass_deficit_from_b=b-(outer_radius**2+b*b-outer_radius**2*((outer_radius-b)/width)**2)/(2*outer_radius),
            strict_minimum_fractional_mass_deficit=1-(outer_radius**2+b*b-outer_radius**2*((outer_radius-b)/width)**2)/(2*outer_radius*b),
            interpretation='For w=.1 and DeltaR=.05, the nonovershooting width condition requires vout>.5, equivalently mu<.17375: a deficit exceeding13.125% of b. This is an algebraic necessary condition, not a constructed collar or a calculation of achievable quantum mass reduction.',
            global_scope='An RN formula with mu<b need not represent a physical naked singularity when its central region is excised and replaced. Global regularity, field equations and the actual two-mouth exterior must be solved separately.'),
        null_integral=dict(
            formula='Q=4*pi*G*I/E=-[r_prime/(N*r)]_in^out-integral(r_prime^2/(N*r^2),dx)',
            RN_matching_boundary_term='With r_prime(inner)=0 and Nout=C*vout, the boundary contribution is -1/(C*Rout); the unknown collar bulk integral is nonnegative.',
            strict_sign='Any positive-lapse collar with these endpoints and a nonzero radius increase has Q<-1/(C*Rout). This does not specify the remaining tensor or its source.',
            clock_scaling='For t_new=alpha*t and the SAME affinely parameterized null geodesic: N_new=N/alpha, E_new=E/alpha, I_new=I, Q_new=alpha*Q. E*Q=4*pi*G*I and the sign are invariant; Q alone is not clock invariant.',
            affine_scaling='If lambda_new=beta*lambda, k_new=k/beta, E_new=E/beta and I_new=I/beta, so Q is invariant under affine rescaling.',
            saved108_Q=saved['geometry_budget']['Qgeom_integration_by_parts'],
            warning='The saved108 number belongs to its original flat-matching collar; it is not the integral for an unconstructed RN-matching collar.'),
        tensor_diagnostic=dict(
            basis='Static orthonormal frame, metric signature(-,+,+,+): covariant tensor diag(rho,p_r,p_t,p_t), with all off-diagonal components zero for this spherical monopole diagnostic.',
            normalization='Every energy density, pressure, null sum, and trace is multiplied by 8*pi*G and has dimension inverse length squared (d=1). Ward residual has dimension inverse length cubed.',
            numerical_normalization='Ward error divided by1+sum of absolute values of its three terms; trace error divided by1+abs(rho)+abs(p_r)+2abs(p_t)+abs(Ricci), using the scaled tensor.',
            field_subtracted='Only the radial monopole Maxwell tensor is subtracted. The angular multipole Maxwell stress from108 and any higher-curvature effective terms are not subtracted; this is not a unique quantum-source target.',
            angular_scope='Spherical monopole diagnostic only. The actual corrected magnetic field is nonspherical and also requires nonspherical Einstein equations or compensating stress.',
            x=x.tolist(), geometry=profiles, **tensor, samples=samples,
            expected_flat_outer_scaled_other=target_outer.tolist(),
            exact_RN_numeric_maximum_residual=max(rn_errors),
            independent_high_precision_arithmetic=high_precision),
        conclusions=dict(saved108_is_an_Einstein_monopole_solution=False,
                         saved108_outer_endpoint_matches_allowed_isolated_RN=False,
                         saved_width_admits_nonovershooting_black_hole_RN_collar=False,
                         any_other_collar_ruled_out=False,
                         exact_two_mouth_exterior_solved=False,
                         quantum_source_or_time_travel_established=False),
        scope=[
            'A rejection and parameter gate for an isolated black-hole RN exterior under an explicit classical magnetic-radius interpretation. It is not an exact two-mouth or global dihole solution.',
            'The compact table includes extremal mu=b and strictly subextremal black holes b<mu<.205. A separate analytic mu<b escape is recorded; the black-hole branch gate does not exclude every possible quantum exterior.',
            'The RN clock may be rescaled, but extrinsic geometry cannot be repaired by a clock choice. A collar must be designed against physically compatible boundary data.',
            'Saved108 geometry remains a legitimate prescribed smooth test background. The nonzero tensor residual identifies stress that would be required to interpret it gravitationally.',
            'Neither a new collar nor a supporting source is constructed. No original smooth, Maxwell, quantum, Dirac, or other field ODE is rerun.'],
        input_sha256=hashes, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    project_path('rn_exterior_gate_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, branch=branch,
                         tensor_samples=samples, Ward_normalized=tensor['maximum_Ward_normalized_error'],
                         trace_normalized=tensor['maximum_trace_normalized_error'],
                         RN_tensor_residual=max(rn_errors)), indent=2))
    if not result['passed']:
        raise AssertionError(tests)
    return result


if __name__ == '__main__':
    run()
