"""A smooth prescribed collar with source-free, angular test-Maxwell matching.

Only the collar ODEs are integrated. The original core's unit-center bases
are read from increment107, preserving their chi normalization. None of the
metric, Maxwell, or quantum fields is claimed to solve Einstein's equations.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from magnetic_mouth_matching import local_controls, solve_matching


def smooth_step(s):
    """C-infinity step, its first two s derivatives, and accurate complement."""
    s = np.asarray(s, dtype=float)
    interior = (s > 0) & (s < 1)
    step = np.where(s >= 1, 1., 0.)
    complement = np.array(1. - step)
    first = np.zeros_like(s)
    second = np.zeros_like(s)
    z = s[interior]
    h = -1/z + 1/(1-z)
    h1 = 1/z**2 + 1/(1-z)**2
    h2 = -2/z**3 + 2/(1-z)**3
    p, q = expit(h), expit(-h)
    step[interior], complement[interior] = p, q
    first[interior] = p*q*h1
    second[interior] = p*q*((q-p)*h1*h1+h2)
    return step, first, second, complement


class CollarGeometry:
    def __init__(self, b=.2, width=.1, rho_c=32., tighter=False):
        self.b, self.width, self.rho_c = b, width, rho_c
        self.y0 = math.asinh(rho_c)
        self.outer_radius = b + width/2
        self.ell = b*math.hypot(1., rho_c)
        self.tighter = tighter
        self.radius_integral = solve_ivp(
            lambda s, state: [float(smooth_step(s)[0])], (0., 1.), [0.],
            method='DOP853', rtol=2e-13 if tighter else 2e-12,
            atol=2e-15 if tighter else 2e-14,
            max_step=1/(512 if tighter else 256), dense_output=True)
        if not self.radius_integral.success:
            raise RuntimeError(self.radius_integral.message)

    def evaluate(self, x):
        x = np.asarray(x, dtype=float)
        s = x/self.width
        if np.any(s < -1e-14) or np.any(s > 1+1e-14):
            raise ValueError('Collar coordinates must lie in [0,width].')
        s = np.clip(s, 0., 1.)
        step, ss, sss, complement = smooth_step(s)
        sx, sxx = ss/self.width, sss/self.width**2
        y = self.y0+x/self.b
        f = np.logaddexp(y, -y)-np.logaddexp(self.y0, -self.y0)
        fp = np.tanh(y)/self.b
        fpp = (1/np.cosh(y)**2)/self.b**2
        logn = complement*f
        lognp = complement*fp-sx*f
        lognpp = complement*fpp-2*sx*fp-sxx*f
        lapse = np.exp(logn)
        radius = self.b+self.width*self.radius_integral.sol(s)[0]
        return dict(radius=radius, radius_prime=step, radius_second=sx,
                    lapse=lapse, lapse_prime=lapse*lognp,
                    lapse_second=lapse*(lognpp+lognp**2),
                    log_lapse_prime=lognp, step=step,
                    step_first_s=ss, step_second_s=sss)


def propagate_modes(saved_modes, geometry):
    """Propagate normalized (a,p=N a_x) without rescaling the old center data."""
    rows, solutions = [], []
    for old in saved_modes:
        degree = old['degree']
        le = old['e_prime_end']/(geometry.ell*old['e_end'])
        lo = old['o_prime_end']/(geometry.ell*old['o_end'])

        def rhs(x, state):
            g = geometry.evaluate(x)
            n, r = float(g['lapse']), float(g['radius'])
            e, pe, o, po = state
            potential = n*degree*(degree+1)/r**2
            return [pe/n, potential*e, po/n, potential*o]

        solution = solve_ivp(rhs, (0., geometry.width), [1., le, 1., lo],
                             method='DOP853', dense_output=True,
                             rtol=2e-13 if geometry.tighter else 2e-12,
                             atol=2e-14 if geometry.tighter else 2e-13,
                             max_step=geometry.width/(256 if geometry.tighter else 128))
        if not solution.success:
            raise RuntimeError(solution.message)
        e, pe, o, po = map(float, solution.y[:, -1])
        actual_e, actual_o = old['e_end']*e, old['o_end']*o
        diagonal = (pe/e+po/o)/2  # N_outer=1.
        off = -1/(2*geometry.ell*actual_e*actual_o)
        rows.append(dict(degree=degree, e_end=actual_e, o_end=actual_o,
                         physical_diagonal=diagonal, physical_offdiagonal=off,
                         even_eigenvalue=diagonal+off,
                         odd_eigenvalue=diagonal-off,
                         lower_reflection_parity=(-1)**(degree+1),
                         inner_even_value=old['e_end'], inner_odd_value=old['o_end'],
                         inner_even_log_derivative=le, inner_odd_log_derivative=lo,
                         normalized_even_outer=e, normalized_even_flux_outer=pe,
                         normalized_odd_outer=o, normalized_odd_flux_outer=po,
                         original_chi_center_wronskian=1.,
                         conserved_N_wronskian=1/geometry.ell,
                         ode_evaluations=solution.nfev))
        solutions.append(solution)
    return rows, solutions


def quadrature(width, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return width*(nodes+1)/2, width*weights/2


def geometric_budget(geometry, order):
    x, weights = quadrature(geometry.width, order)
    g = geometry.evaluate(x)
    r, rp, rpp = g['radius'], g['radius_prime'], g['radius_second']
    n, np_ = g['lapse'], g['lapse_prime']
    positive_integral = float(weights @ (rp*rp/(n*r*r)))
    integrated = -1/geometry.outer_radius-positive_integral
    direct = float(weights @ ((rp*np_/n-rpp)/(n*r)))
    return dict(positive_radius_gradient_integral=positive_integral,
                Qgeom_integration_by_parts=integrated, Qgeom_direct_derivatives=direct,
                direct_identity_absolute_error=abs(direct-integrated),
                optical_length=float(weights @ (1/n)),
                proper_length=geometry.width,
                formula='Qgeom=(4*pi*G/E)*Igeom=-1/rout-integral(r_prime^2/(N*r^2),dx), one collar.',
                direct_formula='integral[(r_prime*N_prime/N-r_second)/(N*r),dx]',
                quadrature_order=order)


def matched_mode_profile(mode, solution, geometry, x):
    """Both outward collars in the common core angular orientation."""
    state = solution.sol(x)
    end = solution.y[:, -1]
    lapse = geometry.evaluate(x)['lapse']
    left, right = mode['core_left_value'], mode['core_right_value']
    ce, co = (left+right)/2, (right-left)/2
    even, odd = state[0]/end[0], state[2]/end[2]
    evenp, oddp = state[1]/(lapse*end[0]), state[3]/(lapse*end[2])
    return dict(upper_value=ce*even+co*odd, upper_prime=ce*evenp+co*oddp,
                lower_core_value=ce*even-co*odd,
                lower_core_prime=ce*evenp-co*oddp,
                normalized_even=state[0], normalized_even_prime=state[1]/lapse,
                normalized_odd=state[2], normalized_odd_prime=state[3]/lapse)


def magnetic_budget(record, solutions, geometry, order):
    x, weights = quadrature(geometry.width, order)
    g = geometry.evaluate(x)
    rows = []
    for mode, solution in zip(record['core_modes'], solutions):
        p = matched_mode_profile(mode, solution, geometry, x)
        l = mode['degree']
        gram = 2*l*(l+1)/(2*l+1)
        upper = .5*gram*float(weights @ (p['upper_prime']**2/(g['lapse']*g['radius']**2)))
        lower = .5*gram*float(weights @ (p['lower_core_prime']**2/(g['lapse']*g['radius']**2)))
        # Static Maxwell energy on one collar, distinct from the radial-null integral.
        collar_energy = math.pi*gram*float(weights @ (g['lapse']*(
            p['upper_prime']**2+l*(l+1)*p['upper_value']**2/g['radius']**2)))
        boundary = matched_mode_profile(mode, solution, geometry, np.array([0., geometry.width]))
        endpoint_g = geometry.evaluate(np.array([0., geometry.width]))
        endpoint_energy = math.pi*gram*float(np.diff(
            endpoint_g['lapse']*boundary['upper_value']*boundary['upper_prime'])[0])
        rows.append(dict(degree=l, upper_collar_M=upper, lower_collar_M=lower,
                         upper_collar_energy_g2=collar_energy,
                         upper_collar_energy_boundary_g2=endpoint_energy,
                         energy_identity_absolute_error=abs(collar_energy-endpoint_energy)))
    cutoff = record['induced_cutoff']
    return dict(retained_upper_M=sum(m['upper_collar_M'] for m in rows[:cutoff]),
                monitored_upper_M=sum(m['upper_collar_M'] for m in rows),
                monitored_lower_M=sum(m['lower_collar_M'] for m in rows),
                retained_upper_collar_energy_g2=sum(m['upper_collar_energy_g2'] for m in rows[:cutoff]),
                monitored_upper_collar_energy_g2=sum(m['upper_collar_energy_g2'] for m in rows),
                formula='M=(g^2/E)*angle_average(I_EM)=0.5*sum_l G_l*integral(a_l_prime^2/(N*r^2),dx), one collar.',
                monopole_radial_null_contribution=0., modes=rows, quadrature_order=order)


def solve_record(nmax, order, modes, solutions, geometry, budget_order):
    record = solve_matching(nmax, order, modes, b=geometry.outer_radius, d=1., q=1.)
    # These reused energy names meant the isolated core in increment107.
    record['retained_throat_plus_two_collars_energy_g2'] = record.pop('retained_core_energy_g2')
    record['monitored_throat_plus_two_collars_energy_g2'] = record.pop('monitored_core_energy_g2')
    for mode in record['core_modes']:
        mode['gauge_coupling_squared_throat_plus_two_collars_energy'] = mode.pop('gauge_coupling_squared_core_energy')
    record['magnetic_radial_null_budget'] = magnetic_budget(record, solutions, geometry, budget_order)
    return record


def differential_controls(geometry, solutions):
    x = geometry.width*np.array([.17, .37, .61, .83])
    h = geometry.width*1e-4
    exact = geometry.evaluate(x)
    samples = [geometry.evaluate(x+j*h) for j in (-2, -1, 1, 2)]
    errors = {}
    for value, derivative in [('radius', 'radius_prime'), ('radius_prime', 'radius_second'),
                              ('lapse', 'lapse_prime'), ('lapse_prime', 'lapse_second')]:
        fd = (samples[0][value]-8*samples[1][value]+8*samples[2][value]-samples[3][value])/(12*h)
        errors[value+'_derivative_absolute_error'] = float(np.max(abs(fd-exact[derivative])))
    ode_errors = []
    for l, solution in enumerate(solutions, 1):
        y = solution.sol(x)
        ys = [solution.sol(x+j*h) for j in (-2, -1, 1, 2)]
        fd = (ys[0]-8*ys[1]+8*ys[2]-ys[3])/(12*h)
        potential = exact['lapse']*l*(l+1)/exact['radius']**2
        rhs = np.array([y[1]/exact['lapse'], potential*y[0],
                        y[3]/exact['lapse'], potential*y[2]])
        ode_errors.append(float(np.max(abs(fd-rhs)/(1+abs(rhs)))))
    errors['modal_first_order_ODE_five_point_relative_error'] = max(ode_errors)
    # Only l=1 is numerically suitable for subtraction of its nearly parallel bases.
    first = solutions[0].sol(np.linspace(0., geometry.width, 51))
    w = first[0]*first[3]-first[1]*first[2]
    errors['degree1_normalized_N_wronskian_relative_drift'] = float(np.max(abs(w/w[0]-1)))
    return errors


def run():
    inputs = ('magnetic_mouth_matching.py', 'magnetic_mouth_matching_results.json')
    hashes = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in inputs}
    saved = json.loads(Path(inputs[1]).read_text(encoding='utf-8'))
    if not saved['passed'] or len(saved['core_modes']) < 15:
        raise ValueError('The successful increment107 core endpoint data are required.')
    expected = dict(q=1., b=.2, d=1., rho_c=32.)
    if any(saved['parameters'].get(key) != value for key, value in expected.items()):
        raise ValueError('Saved107 parameters do not match this collar/core normalization.')
    if saved['source_sha256'] != hashes['magnetic_mouth_matching.py']:
        raise ValueError('Saved107 results do not identify the current matching source.')
    old = saved['core_modes'][:15]
    geometry = CollarGeometry()
    modes, solutions = propagate_modes(old, geometry)
    records = [solve_record(n, order, modes, solutions, geometry, 128)
               for n in (3, 6, 9) for order in (32, 64)]
    fine = [r for r in records if r['quadrature_order'] == 64]
    selected = fine[-1]
    budget = geometric_budget(geometry, 128)
    local = local_controls(b=geometry.outer_radius)
    # The reused helper's sharp-core clock used its b argument as core radius.
    # Here that argument is the OUTER sphere radius, so keep only exterior data.
    for key in ('core_cut_lapse', 'cut_proper_to_chi_derivative', 'physical_dtn_clock_factor'):
        local.pop(key)
    derivatives = differential_controls(geometry, solutions)
    tighter_geometry = CollarGeometry(tighter=True)
    tighter_modes, tighter_solutions = propagate_modes(old, tighter_geometry)
    tighter_record = solve_record(9, 64, tighter_modes, tighter_solutions, tighter_geometry, 256)
    tighter_budget = geometric_budget(tighter_geometry, 256)
    refinement = dict(
        maximum_physical_DtN_diagonal_change=max(abs(a['physical_diagonal']-b['physical_diagonal']) for a, b in zip(modes, tighter_modes)),
        maximum_beta_change=float(np.max(abs(np.asarray(selected['beta_upper'])-tighter_record['beta_upper']))),
        monitored_mismatch_norm_change=abs(selected['monitored_total_mismatch_weighted_norm']-tighter_record['monitored_total_mismatch_weighted_norm']),
        Qgeom_change=abs(budget['Qgeom_integration_by_parts']-tighter_budget['Qgeom_integration_by_parts']),
        optical_length_change=abs(budget['optical_length']-tighter_budget['optical_length']),
        monitored_upper_M_change=abs(selected['magnetic_radial_null_budget']['monitored_upper_M']-tighter_record['magnetic_radial_null_budget']['monitored_upper_M']),
        description='One tighter collar integration: half maximum step, tenfold tighter rtol/atol, and radial quadrature128 to256; no original core integration.')
    end = geometry.evaluate(np.array([0., geometry.width]))
    endpoints = {k: np.asarray(v).tolist() for k, v in end.items()}
    endpoint_errors = dict(
        radius_outer=abs(end['radius'][1]-.25),
        lapse= float(np.max(abs(end['lapse']-1))),
        radius_slopes=float(np.max(abs(end['radius_prime']-[0., 1.]))),
        radius_second=float(np.max(abs(end['radius_second']))),
        inner_lapse_prime=abs(end['lapse_prime'][0]-32/(.2*math.hypot(1, 32))),
        inner_lapse_second=abs(end['lapse_second'][0]-1/.2**2),
        outer_lapse_prime=abs(end['lapse_prime'][1]),
        outer_lapse_second=abs(end['lapse_second'][1]))
    quadrature_refinement = []
    for n in (3, 6, 9):
        a, b = [r for r in records if r['induced_cutoff'] == n]
        quadrature_refinement.append(dict(induced_cutoff=n,
            maximum_beta_change=float(np.max(abs(np.asarray(a['beta_upper'])-b['beta_upper']))),
            monitored_current_norm_change=abs(a['monitored_total_mismatch_weighted_norm']-b['monitored_total_mismatch_weighted_norm'])))
    cutoff_refinement = []
    for a, b in zip(fine[:-1], fine[1:]):
        n = a['induced_cutoff']
        cutoff_refinement.append(dict(lower_cutoff=n, upper_cutoff=b['induced_cutoff'],
            maximum_shared_beta_change=float(np.max(abs(np.asarray(a['beta_upper'])-b['beta_upper'][:n]))),
            monitored_current_norm_reduction=a['monitored_total_mismatch_weighted_norm']/b['monitored_total_mismatch_weighted_norm']))
    x = np.linspace(0., geometry.width, 513)
    gg = geometry.evaluate(x)
    profiles = dict(x=x.tolist(), geometry={k: np.asarray(v).tolist() for k, v in gg.items()},
        modal=[dict(degree=m['degree'], **{k: np.asarray(v).tolist() for k, v in matched_mode_profile(m, s, geometry, x).items()})
               for m, s in zip(selected['core_modes'], solutions)],
        selected_record=dict(induced_cutoff=9, monitor_degree=15, quadrature_order=64))
    magnetic_density = {}
    for side in ('upper', 'lower_core'):
        magnetic_density[side] = sum(
            m['degree']*(m['degree']+1)/(2*m['degree']+1)*np.asarray(m[side+'_prime'])**2
            for m in profiles['modal'])/(gg['lapse']*gg['radius']**2)
    profiles['radial_null_densities'] = dict(
        Qgeom_density=((gg['radius_prime']*gg['lapse_prime']/gg['lapse']-gg['radius_second'])/(gg['lapse']*gg['radius'])).tolist(),
        Maxwell_M_density_upper=magnetic_density['upper'].tolist(),
        Maxwell_M_density_lower=magnetic_density['lower_core'].tolist(),
        normalization='Integrate each density over proper x to obtain its separately normalized one-collar Qgeom or monitored M; these columns have different dimensions.')
    classical_radius_diagnostic = dict(
        assumed_in_general_symbolic_results=False,
        assumption='ONLY IF core radius b is the classical magnetic radius r_e: r_e^2=pi*G*q^2/g_R^2, with q=1 and d=1 here.',
        conditional_dimensionless_4piG_over_gR2_d2=4*geometry.b**2,
        conditional_r_e_squared_over_outer_radius_squared=(geometry.b/geometry.outer_radius)**2,
        conditional_dimensionless_remaining_non_Maxwell_Q=budget['Qgeom_integration_by_parts']-4*geometry.b**2*selected['magnetic_radial_null_budget']['monitored_upper_M'],
        isolated_extremal_RN_comparator=dict(
            applicability='Isolated extremal RN illustration, not the two-mouth exterior solution.',
            F=(1-geometry.b/geometry.outer_radius)**2,
            lapse=1-geometry.b/geometry.outer_radius,
            radius_proper_derivative=1-geometry.b/geometry.outer_radius,
            log_lapse_proper_derivative=geometry.b/geometry.outer_radius**2,
            angular_extrinsic_curvature=(1-geometry.b/geometry.outer_radius)/geometry.outer_radius),
        actual_flat_endpoint=dict(lapse=1., radius_proper_derivative=1.,
                                  log_lapse_proper_derivative=0., angular_extrinsic_curvature=1/geometry.outer_radius),
        interpretation='The conditional charge-curvature ratio is .64, not a small electrovac expansion parameter. A time-coordinate rescaling does not remove the extrinsic-curvature mismatch. Gravitational exterior matching remains unresolved.')
    tests = dict(
        smooth_endpoint_matching=max(endpoint_errors.values()) < 5e-12,
        positive_radius_and_lapse=float(np.min(gg['radius'])) >= .2-1e-14 and float(np.min(gg['lapse'])) > 0,
        geometry_derivative_control=max(v for k, v in derivatives.items() if 'derivative_absolute' in k) < 2e-7,
        propagated_modal_ODE_control=derivatives['modal_first_order_ODE_five_point_relative_error'] < 2e-8,
        conserved_normalization_control=derivatives['degree1_normalized_N_wronskian_relative_drift'] < 2e-7,
        positive_two_end_DtN=min(min(m['even_eigenvalue'], m['odd_eigenvalue']) for m in modes) > 0,
        exact_incident_and_self_projections=max(max(r['incident_value_projection_error'], r['incident_derivative_projection_error'], r['self_value_projection_error'], r['self_derivative_projection_error']) for r in records) < 3e-11,
        retained_Maxwell_matching=max(r['retained_mismatch_maximum'] for r in records) < 3e-14,
        reflected_mouth_symmetry=max(r['reflection_symmetry_error'] for r in records) < 3e-14,
        charge_unchanged=local['extra_magnetic_flux_at_either_mouth'] < 2e-14,
        exterior_derivative_control=max(local['induced_normal_derivative_five_point_error'], local['incident_normal_derivative_five_point_error']) < 3e-7,
        geometric_null_identity=budget['direct_identity_absolute_error'] < 2e-11 and budget['Qgeom_integration_by_parts'] < -4,
        nonnegative_Maxwell_null_budget=min(m['upper_collar_M'] for r in records for m in r['magnetic_radial_null_budget']['modes']) >= 0,
        positive_energy_and_identity=min(m['boundary_canonical_quadratic_form'] for r in records for m in r['core_modes']) >= 0 and max(m['energy_identity_absolute_error'] for r in records for m in r['magnetic_radial_null_budget']['modes']) < 2e-11,
        angular_quadrature_stability=max(r['maximum_beta_change'] for r in quadrature_refinement) < 2e-12,
        monitored_omitted_modes_decrease=all(r['monitored_current_norm_reduction'] > 20 for r in cutoff_refinement),
        tighter_integration_stability=refinement['maximum_physical_DtN_diagonal_change'] < 2e-9 and refinement['maximum_beta_change'] < 2e-12 and refinement['Qgeom_change'] < 2e-11 and refinement['monitored_upper_M_change'] < 2e-11,
        inputs_preserved=all(hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest for name, digest in hashes.items()))
    tests = {k: bool(v) for k, v in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests,
        parameters=dict(q=1., d=1., core_radius=.2, collar_proper_width=.1, exterior_cut_radius=.25,
                        rho_c=32., ell=geometry.ell, induced_cutoffs=[3, 6, 9], monitor_extra_degrees=6,
                        angular_quadrature_orders=[32, 64], radial_quadrature_order=128),
        conventions=dict(
            geometry='ds^2=-N(x)^2 dt^2+dx^2+r(x)^2 dOmega^2 on each outward collar; the two collars are reflected copies.',
            smooth_construction='S(s)=exp(-1/s)/(exp(-1/s)+exp(-1/(1-s))) with flat endpoint extensions; r=b+w*integral_0^(x/w)S and N=exp((1-S)*log(cosh(asinh(rho_c)+x/b)/cosh(asinh(rho_c)))). Symmetry S(s)+S(1-s)=1 gives rout=b+w/2.',
            endpoint_jets='At x=0 the jets equal r=b, N=cosh(asinh(rho_c)+x/b)/cosh(asinh(rho_c)); at x=w they equal r=rout+x-w, N=1. Flatness of S at both endpoints proves C-infinity matching, independently of finite endpoint checks.',
            Maxwell_ODE="(N*a_prime)_prime=N*l*(l+1)*a/r^2; prime denotes proper outward x derivative.",
            DtN='Unit original chi-center bases obey N*(e*o_prime-o*e_prime)=1/ell. At Nouter=1, Ddiag=(e_prime/e+o_prime/o)/2 and Doff=-1/(2*ell*e_outer*o_outer).',
            reused_center_normalization='solve_matching receives actual propagated original-center e_end/o_end. Its center_chi_derivative remains da/dchi at the old core center, not an x derivative; center_value is unchanged in meaning.',
            reused_boundary_names='record core_left_value/core_right_value now mean the boundary values at the OUTER ends of the two collars in the common core angular orientation. Boundary energy is for the old core PLUS BOTH collars.',
            lower_reflection='Lower core cosine=-global cosine. Both value and outward derivative receive (-1)^(l+1).',
            current_matching='The exterior starts at radius .25, where N=1 and dr/dx=1. Equal proper normal derivatives exclude a sheet current in retained modes.',
            geometric_budget='Qgeom=(4*pi*G/E)*integral(Tgeom_kk dlambda), one collar, radial null with conserved Killing energy E.',
            magnetic_budget='M=(g^2/E)*angle_average integral(TEM_kk dlambda), one collar; G_l=2*l*(l+1)/(2*l+1). Monopole Ftheta_phi gives zero radial-null contribution.',
            units='d=1 length unit. Qgeom has dimension inverse length, M inverse length cubed. G and g^2 remain symbolic; these numerical columns are not added.',
            relation='The remaining non-Maxwell effective stress in Einstein form must supply Qgeom-(4*pi*G/g^2)*M in Qgeom units. Fixed higher-curvature terms must be carried consistently with a specified gravitational renormalization scheme; this is not a unique Dirac-source target.'),
        endpoint_profiles=endpoints, endpoint_errors=endpoint_errors,
        geometry_budget=budget, core_plus_collar_modes=modes, records=records,
        profiles=profiles, local_controls=local, differential_controls=derivatives,
        tighter_integration_control=refinement, angular_quadrature_sensitivity=quadrature_refinement,
        angular_cutoff_sensitivity=cutoff_refinement,
        conditional_classical_magnetic_radius_diagnostic=classical_radius_diagnostic,
        conclusions=dict(smooth_prescribed_geometry_constructed=True,
                         retained_source_free_Maxwell_matching_solved=True,
                         monitored_Maxwell_radial_null_is_nonnegative=True,
                         metric_requires_negative_integrated_radial_null_source=True,
                         Einstein_equations_solved=False, all_angular_matching_certified=False),
        scope=[
            'One explicit smooth, prescribed, spherically symmetric collar per mouth. Smoothness removes the sharp metric interface; the stress required by this geometry has not been physically supplied.',
            'Source-free static Maxwell ODEs are solved on this geometry and matched in a finite exterior angular projection. The full core-plus-collar field means the full Dirichlet extension of the actual exterior trace; only modes through L+6 are evaluated.',
            'Higher Maxwell interface modes and numerical errors remain unenclosed. Sensitivity controls and dense floating profiles are not proof of an infinite-mode solution.',
            'The Maxwell field is angularly varying while the prescribed geometry is spherical. Its angle-averaged null budget is a necessary diagnostic only; the nonspherical Einstein equations remain unsolved.',
            'The classical b=r_e coupling relation is supplied only as an explicit conditional diagnostic. It would make the near-mouth flat exterior an uncontrolled electrovac approximation; it is not assumed by the symbolic test-field calculation.',
            'The Maxwell contribution is nonnegative and cannot supply the negative radial-null integral required here. G/g^2 and any physical length scale remain unspecified.',
            'No quantum source, Einstein balance, material production, formation, stability, causality, or time-travel result follows. Previous Dirac-source numbers do not transfer to this new geometry/field.',
            'Only new collar ODEs and small exterior multipole systems were solved. Saved107 core endpoints and historical files were preserved.'],
        input_sha256=hashes, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=saved.get('sources', []))
    destination = Path(__file__).with_name('smooth_magnetic_collar_results.json')
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, geometry_budget=budget,
        differential_controls=derivatives, tighter_integration_control=refinement,
        fine_records=[dict(induced_cutoff=r['induced_cutoff'],
                           bare_mismatch_weighted_norm=r['bare_mismatch_weighted_norm'],
                           retained_mismatch_maximum=r['retained_mismatch_maximum'],
                           higher_mismatch_weighted_norm=r['higher_mismatch_weighted_norm'],
                           leading_beta=r['beta_upper'][0],
                           magnetic_budget_M=r['magnetic_radial_null_budget']['monitored_upper_M'],
                           total_energy_g2=r['monitored_throat_plus_two_collars_energy_g2']) for r in fine]), indent=2))
    if not result['passed']:
        raise AssertionError(tests)
    return result


if __name__ == '__main__':
    run()
