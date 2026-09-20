"""Analytic quadrature-remainder component for the frozen117 compact trial.

This does NOT enclose the floating quadrature sum, nor the full physical
residual.  It bounds |I-Q| for EXACT tensor Gauss--Legendre quadrature on nine
width-two panels, using the exact dyadic parameters and Hermite coefficients
saved by117.  That is a different quadrature rule from117's five panels.
No new PDE solve or integral evaluation is performed here.

Proof ingredients (our deliberately conservative combination):

* Mehler, https://dlmf.nist.gov/18.18.E28, for normalized Hermite functions:
  sum t**j |phi_j(x+iy)|**2 = [pi(1-t*t)]**(-1/2)
  * exp[-(1-t)/(1+t)*x*x+(1+t)/(1-t)*y*y].  Thus the first K
  terms are bounded by t**(-(K-1)) times this expression.
* Trefethen, SIAM Review50(2008), Theorems4.2--4.3, equations4.7/4.11:
  https://people.maths.ox.ac.uk/~trefethen/publication/PDF/2008_127.pdf .
  Chebyshev coefficients obey |a_j|<=2*M*rho**(-j).  Truncating after
  degree2*n-1 has error<=2*M*rho**(1-2*n)/(rho-1).  Positivity and
  polynomial exactness of GLn imply |I-Q|<=h*delta*M on an interval of
  length h, delta=4*rho**(1-2*n)/(rho-1).  Tensorization needs only
  one complex coordinate at a time: error<=area*delta*(M_u+M_v).
* Balanced-dihole algebra is the frozen117 differential expression, based
  on Emparan hep-th/9906160 equations8--16.  Its mass can be evaluated as
  F*(mu*F-q*g*Pnum)/(rho*E*Sigma**2), without a spurious lapse pole.

Each cutoff is continued using its panel's polynomial, NOT complex abs or
clipping.  All real cutoff knots are panel boundaries.  F_res=sum r_j**2
is the holomorphic continuation of the real squared residual, rather than
the nonholomorphic sum |r_j(z)|**2.  Its modulus is bounded by ||r(z)||**2.

The certified component assumes mpmath interval arithmetic behaves as
documented.  Binary64 trial coefficients/parameters are exact inputs here;
the exact-balanced physical-parameter discrepancy is a separate error.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from project_paths import project_path, source_matches


def dyadic(iv, value):
    """Exact binary64 input, with no decimal reinterpretation."""
    numerator, denominator = float(value).as_integer_ratio()
    return iv.mpf(numerator)/denominator


def outward(value):
    """Serializable binary64 enclosure of an interval-arithmetic constant."""
    lower = math.nextafter(float(value.a), -math.inf)
    upper = math.nextafter(float(value.b), math.inf)
    return dict(lower=lower, upper=upper,
                lower_hex=lower.hex(), upper_hex=upper.hex())


def mehler_bound(iv, count, t, imaginary_extent):
    """Uniform in the real coordinate; discard a nonpositive exponent."""
    return (t**(-(count-1))/iv.sqrt(iv.pi*(1-t*t))
            *iv.exp((1+t)/(1-t)*imaginary_extent**2))


def analytic_geometry_bounds(iv):
    """Triangle/reverse-triangle inequalities on an explicit complex box.

    rho.real in [.515,.820], |z.real|<=.152, and both imaginary parts
    at most .0063.  These encompass either one-coordinate Bernstein strip.
    All quantities below are scalar formula enclosures, not sampled maxima.
    """
    p = iv.mpf
    rho_lo, rho_hi, z_re, y = map(p, ('.515', '.820', '.152', '.0063'))
    mass_lo, mass_hi, ag, field, k = map(p, ('.153', '.154', '.476', '.241', '.5'))
    rho_mod = iv.sqrt(rho_hi**2+y**2)
    z_mod = iv.sqrt(z_re**2+y**2)
    displaced_mod = iv.sqrt((k+z_re)**2+y**2)
    # Re[rho^2+(z+-k)^2]>0 fixes the analytic principal square roots.
    radicand_real_lo = rho_lo**2+(k-z_re)**2-2*y**2
    root_lo = iv.sqrt(radicand_real_lo)
    root_hi = iv.sqrt(rho_mod**2+displaced_mod**2)
    r_lo, r_hi = mass_lo+root_lo, mass_hi+root_hi
    c_hi = z_mod/root_lo  # c=2*z/(Rplus+Rminus).
    sigma_lo = r_lo**2-ag**2*c_hi**2
    sigma_hi = r_hi**2+ag**2*c_hi**2
    h_lo = root_lo**2-k**2*c_hi**2
    h_hi = root_hi**2+k**2*c_hi**2
    assert radicand_real_lo > 0 and sigma_lo > 0 and h_lo > 0

    r_rho, r_z = rho_mod/root_lo, displaced_mod/root_lo
    c_rho = c_hi*r_rho/root_lo
    c_z = (1+c_hi*r_z)/root_lo
    speed = h_hi**2/sigma_lo**2
    def speed_derivative(r_j, c_j):
        h_j = 2*root_hi*r_j+2*k**2*c_hi*c_j
        sigma_j = 2*r_hi*r_j+2*ag**2*c_hi*c_j
        return 2*h_hi*h_j/sigma_lo**2+2*h_hi**2*sigma_j/sigma_lo**3
    speed_rho, speed_z = speed_derivative(r_rho, c_rho), speed_derivative(r_z, c_z)

    # F=N*Sigma; Aphysical=g*Pnum/F, g=2.5 exactly in frozen117.
    delta_hi = r_hi**2+2*mass_hi*r_hi+ag**2
    w_hi = (r_hi**2+ag**2)**2+delta_hi*ag**2*(1+c_hi**2)
    f_hi = (sigma_hi+2*mass_hi*r_hi
            +2*field*mass_hi*r_hi*ag*(1+c_hi**2)
            +field**2*w_hi*(1+c_hi**2)/4)
    pnum_hi = (2*mass_hi*r_hi*ag+field*w_hi/2)*(1+c_hi**2)
    normalized_mass = f_hi*(400000*f_hi+512000*p('2.5')*pnum_hi)/(rho_lo*p('43.3')*sigma_lo**2)
    normalized_speed = speed/p('.509')
    connection = p('.011762')*(speed_rho+speed_z)/(2*iv.sqrt(2)*p('.509'))
    return dict(radicand_real_lower=radicand_real_lo, root_real_lower=root_lo,
        root_modulus_upper=root_hi, r_modulus_upper=r_hi, c_modulus_upper=c_hi,
        sigma_modulus_lower=sigma_lo, h_modulus_lower=h_lo,
        speed_modulus_upper=speed, speed_rho_modulus_upper=speed_rho,
        speed_z_modulus_upper=speed_z, lapse_numerator_modulus_upper=f_hi,
        potential_numerator_modulus_upper=pnum_hi,
        normalized_speed_upper=normalized_speed,
        normalized_connection_upper=connection,
        normalized_mass_upper=normalized_mass)


def self_checks(iv):
    """Independent exact-small quadrature and finite-Mehler controls."""
    p = iv.mpf
    # GL2: nodes +-1/sqrt3, weights1.  The tensor test has a known integral.
    node = 1/iv.sqrt(3)
    q1 = 1/(3-node)+1/(3+node)
    exact1 = iv.log(2)
    exact_error = exact1**2-q1**2
    # On either ellipse2 strip: |3-z|>=7/4, other coordinate denominator>=2.
    strip_bound = p(2)/7
    tensor_bound = 4*(4*p(2)**(-3))*(2*strip_bound)
    assert exact_error > 0 and exact_error < tensor_bound
    assert abs(q1-p(9)/13) < p('1e-45')

    # Independently sum normalized Hermites at a complex point.  This is
    # an implementation check of Mehler normalization, not the proof of it.
    z = iv.mpc(p('.3'), p('.2')); t = p('.25')
    phi0 = iv.exp(-z*z/2)/iv.pi**p('.25')
    phis = [phi0, iv.sqrt(2)*z*phi0]
    for j in range(1, 63):
        phis.append(iv.sqrt(p(2)/(j+1))*z*phis[-1]-iv.sqrt(p(j)/(j+1))*phis[-2])
    partial = sum((t**j*abs(value)**2 for j, value in enumerate(phis)), p(0))
    full = (iv.exp(-(1-t)/(1+t)*p('.3')**2+(1+t)/(1-t)*p('.2')**2)
            /iv.sqrt(iv.pi*(1-t*t)))
    difference = full-partial
    assert abs(difference) < p('1e-30')
    # The finite upper bound includes every term, without relying on that
    # tiny numerical difference having a definite sign.
    assert partial < mehler_bound(iv, len(phis), t, p('.2'))
    return dict(tensor_GL2_exact_error=outward(exact_error),
                tensor_GL2_remainder_bound=outward(tensor_bound),
                finite_Mehler_difference=outward(difference))


def run_checks():
    started = time.perf_counter()
    iv = MPIntervalContext(); iv.dps = 60
    p = iv.mpf
    calibration_path = project_path('dihole_null_scaling_results.json')
    previous_path = project_path('dihole_probe_residual_results.json')
    field_path = project_path('dihole_probe_residual_fields.npz')
    calibration = json.loads(calibration_path.read_text())
    previous = json.loads(previous_path.read_text())
    assert calibration['passed'] and previous['passed']
    assert source_matches('dihole_probe_residual.py', previous['source_sha256'])
    assert source_matches('dihole_null_scaling.py', calibration['source_sha256'])
    assert hashlib.sha256(calibration_path.read_bytes()).hexdigest() == previous['input_sha256']['dihole_null_scaling_results.json']
    assert hashlib.sha256(field_path.read_bytes()).hexdigest() == previous['field_sha256']
    parameters, coefficients = calibration['parameters'], calibration['coefficients']
    saved = previous['refined_corrected']['wide_q512000']
    assert saved['q'] == 512000 and saved['mu'] == 398908.5
    assert saved['kappa'] == 1 and saved['inner'] == 7 and saved['outer'] == 9
    assert parameters['k'] == .5 and parameters['physical_A_factor'] == -2.5
    for name, lower, upper in (('M', '.153', '.154'), ('a', '0', '.476'), ('B', '0', '.241')):
        assert p(lower) < dyadic(iv, parameters[name]) < p(upper)
    length = dyadic(iv, saved['length_scale'])
    energy = dyadic(iv, saved['energy_scale'])
    center = dyadic(iv, coefficients['rho_star'])
    vstar = dyadic(iv, coefficients['v_star'])
    assert length < p('.011762') and energy > p('43.3') and vstar > p('.509')

    # Either coordinate on a width2 ellipse2: |Re|<=9.25, |Im|<=.75;
    # the other is real in [-9,9].  Both physical coordinates share these
    # enclosure bounds, although their imaginary parts are correlated.
    real_extent = length*p('18.25')/iv.sqrt(2)
    imaginary_extent = length*p('.75')/iv.sqrt(2)
    assert center-real_extent > p('.515')
    assert center+real_extent < p('.820')
    assert real_extent < p('.152') and imaginary_extent < p('.0063')

    norms = {}
    with np.load(field_path) as fields:
        for key in ('wide_q512000_a', 'wide_q512000_b'):
            array = fields[key]
            assert array.shape == (128,128,2) and np.isfinite(array).all()
            # Every saved coefficient is an exact dyadic, not a rounded
            # decimal string.  Sum both spin components in the same norm.
            squared = sum((dyadic(iv,x)**2 for x in array.ravel()), p(0))
            norms[key] = iv.sqrt(squared)
            assert norms[key] < p('.65')

    geometry = analytic_geometry_bounds(iv)
    assert geometry['normalized_speed_upper'] < 8
    assert geometry['normalized_connection_upper'] < p('1.1')
    assert geometry['normalized_mass_upper'] < 500000
    kc = mehler_bound(iv, 129, p('.9'), p('.75'))
    kr = mehler_bound(iv, 129, p('.99'), p(0))
    assert kc < p('4.2e10') and kr < 15
    # Finite derivative-coefficient norm <=sqrt(2N)*||a|| =16||a||.
    field_bound = p('.65')*iv.sqrt(p('4.2e10')*15)
    assert field_bound < 520000
    # For t=(z+1)/2 on ellipse2: |t|,|1-t|<=9/8.  Continue each
    # panel polynomial (1-t)^3(1+3t+6t^2); no complex clipping.
    tmax = p(9)/8
    chi_bound = tmax**3*(1+3*tmax+6*tmax**2)
    derivative_bound = 15*tmax**4
    assert chi_bound < 18 and derivative_bound < 25
    cut_field = p(18)*520000
    cut_derivative = (18*16+25)*p(520000)
    source_bound = 18*iv.exp(p('.75')**2/2)/iv.sqrt(iv.pi)
    assert source_bound < 14
    # AU,AV,J each have Euclidean operator norm1, including complex inputs.
    residual_bound = 14+8*(2*cut_derivative)+p('2.2')*cut_field+(500000+1)*cut_field
    integrand_bound = residual_bound**2
    assert integrand_bound < p('3e25')

    remainders = []
    for order in (48, 64, 72, 80):
        delta = 4*p(2)**(1-2*order)
        # Sum of all panel areas is18^2; there are TWO coordinate errors.
        bound = 2*(18**2)*delta*p('3e25')
        remainders.append(dict(order_per_width2_panel=order,
            total_nodes_per_axis=9*order,
            quantity='Enclosure of the positive upper-bound constant, not of the signed actual quadrature error.',
            squared_residual_integral_error=outward(bound)))
    controls = self_checks(iv)
    tests = dict(input_provenance=True, dyadic_target_and_parameter_boxes=True,
        Bernstein_strips_inside_geometry_box=True, square_root_branches_and_denominators=True,
        full_coefficient_norms_enclosed=True, Mehler_and_derivative_bounds=True,
        panel_cutoff_polynomial_bounds=True, complex_residual_bound_derived=True,
        independent_tensor_quadrature_control=True, finite_Mehler_normalization_control=True)
    inputs = [calibration_path, previous_path, field_path,
        project_path('dihole_probe_residual.py'), project_path('dihole_null_scaling.py'),
        project_path('balanced_dihole.py'), project_path('dihole_return_geometry.py'),
        project_path('null_saddle_leads.py'), project_path('project_paths.py')]
    return dict(passed=all(tests.values()), tests=tests,
        target=dict(q=512000, mu=398908.5, kappa=1, field_keys=list(norms),
            parameter_interpretation='Exact saved binary64/dyadic integrand; not exact-balanced radical parameters.',
            source_cutoff=[5,7], trial_cutoff=[7,9], coefficient_shape=[128,128,2]),
        quadrature_rule=dict(panel_edges=list(range(-9,10,2)), ellipse_parameter=2,
            complex_coordinate_real_extent=9.25, complex_coordinate_imaginary_extent=.75,
            other_coordinate_interval=[-9,9], area=324,
            evaluated_here=False, same_as_frozen117_rule=False),
        enclosure=dict(interval_digits=iv.dps,
            outward_representation='binary64 endpoints rounded outward with nextafter; hex endpoints are exact.',
            physical_real_extent=outward(real_extent), physical_imaginary_extent=outward(imaginary_extent),
            coefficient_norms={k:outward(v) for k,v in norms.items()},
            geometry={k:outward(v) for k,v in geometry.items()},
            Mehler_complex=outward(kc), Mehler_real=outward(kr),
            cutoff_modulus=outward(chi_bound), cutoff_derivative_modulus=outward(derivative_bound),
            complex_source_modulus=outward(source_bound),
            complex_residual_norm=outward(residual_bound),
            complex_squared_integrand_modulus=outward(integrand_bound),
            rounded_integrand_budget='3e25'),
        remainder_components=remainders, controls=controls,
        scope=[
            'Certified analytic truncation component for exact GL nodes/weights and exact sums on the specified width2 rule.',
            'No curved residual integral or quadrature sum was evaluated or enclosed by this module.',
            'Saved117 ordinary quadrature used different panels; its numerical value does not inherit this bound.',
            'Node/weight evaluation, Hermite contractions, coefficient cancellation and sum rounding remain to be enclosed.',
            'The difference between dyadic calibration and exact-balanced physical parameters remains separate.',
            'No certified curved Green disk, quantum stress, backreaction or backward-time passage is established.'],
        formulas=dict(Mehler='sum_{j<K}|phi_j(x+iy)|^2 <= t^{-(K-1)} exp((1+t)y^2/(1-t))/sqrt(pi*(1-t^2))',
            one_coordinate_GL='|I-Q_n| <= length*4*M*rho^{1-2n}/(rho-1)',
            tensor_GL='|I2-Q2| <= area*4*rho^{1-2n}*(M_u+M_v)/(rho-1)',
            mass='F*(mu*F-q*g*Pnum)/(rho*E*Sigma^2); g=2.5'),
        sources=[dict(title='NIST DLMF18.18.28: Mehler formula',url='https://dlmf.nist.gov/18.18.E28'),
            dict(title='Trefethen2008: Is Gauss quadrature better than Clenshaw-Curtis?',
                 url='https://people.maths.ox.ac.uk/~trefethen/publication/PDF/2008_127.pdf'),
            dict(title='Emparan1999: Black diholes',url='https://arxiv.org/abs/hep-th/9906160')],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={str(path.relative_to(project_path('dihole_probe_residual.py').parent.parent)).replace('\\','/'):
                      hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
        runtime_seconds=time.perf_counter()-started)


if __name__ == '__main__':
    result = run_checks()
    path = project_path('dihole_quadrature_remainder_results.json')
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=len(result['tests']),
        remainder_components=result['remainder_components'], runtime_seconds=result['runtime_seconds']),indent=2))
