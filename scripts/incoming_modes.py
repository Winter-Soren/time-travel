"""Unit-flux monopole collar modes and an exactly positive-energy pulse.

These are reference lead modes, not the two-mouth scattering solution.
Run with the existing .venv; no grid evolution or eigensolver is needed.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, betaincinv, erfc


def collar_mode(q, k, rho, z, energy, charge=1):
    """Reduced f_mu about one mouth, north patch, unit signed radial flux.

    z is measured from the mouth; mu=charge*(k+1/2), Aphi=charge*q*t.
    The Hamiltonian is the project's +sigma.Pi Weyl block.
    """
    if not (isinstance(q, int) and q > 0 and isinstance(k, int) and 0 <= k < q):
        raise ValueError('Require integer q>0 and 0<=k<q')
    if charge not in (-1, 1) or not math.isfinite(energy) or energy <= 0:
        raise ValueError('Require charge=+/-1 and finite energy>0')
    rho, z = np.broadcast_arrays(np.asarray(rho, float), np.asarray(z, float))
    if np.any(~np.isfinite(rho)) or np.any(~np.isfinite(z)) or np.any(rho <= 0):
        raise ValueError('Use finite rho>0, z; axes require patchwise limits')
    radius = np.hypot(rho, z)
    theta = np.arctan2(rho, z)
    sh, ch = np.sin(theta / 2), np.cos(theta / 2)
    amplitude = np.exp(.5 * (np.log(rho / radius**2) - math.log(2)
                       - betaln(k + 1, q - k))
                       + k * np.log(sh) + (q - 1 - k) * np.log(ch))
    spin = np.stack((ch, sh) if charge == 1 else (sh, -ch))
    return spin * amplitude * np.exp(1j * charge * energy * radius)


def positive_pulse(u, width=1.):
    """Fourier transform of a(E)=2*width**1.5*E*exp(-width*E), E>0."""
    if not math.isfinite(width) or width <= 0:
        raise ValueError('Require finite width>0')
    u = np.asarray(u, float)
    if np.any(~np.isfinite(u)):
        raise ValueError('Require finite pulse coordinate')
    return 2 * width**1.5 / math.sqrt(2 * math.pi) / (width + 1j * u)**2


def checks():
    # Independent high precision derivatives of the analytic spinor; no stencil.
    residuals, gauge_errors, implementation_errors = [], [], []
    with mp.workdps(45):
        sx = mp.matrix([[0, 1], [1, 0]])
        sy = mp.matrix([[0, -1j], [1j, 0]])
        sz = mp.matrix([[1, 0], [0, -1]])
        for q, k in ((1, 0), (7, 0), (7, 3), (7, 6), (100000, 0), (100000, 16)):
            for sign in (-1, 1):
                E = mp.mpf('1.3')
                def f(rho, z):
                    r = mp.sqrt(rho*rho + z*z)
                    theta = mp.atan2(rho, z)
                    sh, ch = mp.sin(theta/2), mp.cos(theta/2)
                    scalar = mp.sqrt(rho/(2*r*r*mp.beta(k+1, q-k)))
                    scalar *= sh**k * ch**(q-1-k) * mp.exp(1j*sign*E*r)
                    return scalar * mp.matrix([ch, sh] if sign == 1 else [sh, -ch])
                for t in ((mp.mpf(k)+mp.mpf('.5'))/q,
                          (mp.mpf(k)+mp.mpf('.8'))/(q+1)):
                    r = mp.mpf('1.7')
                    rho, z = 2*r*mp.sqrt(t*(1-t)), r*(1-2*t)
                    v = f(rho, z)
                    reference = np.array([complex(component) for component in v])
                    numerical = collar_mode(q, k, float(rho), float(z), float(E), sign)
                    implementation_errors.append(float(np.linalg.norm(numerical-reference)/np.linalg.norm(reference)))
                    dr = mp.matrix([mp.diff(lambda x: f(x, z)[i], rho) for i in range(2)])
                    dz = mp.matrix([mp.diff(lambda x: f(rho, x)[i], z) for i in range(2)])
                    mu, A = sign*(mp.mpf(k)+mp.mpf('.5')), sign*q*t
                    Hf = -1j*sx*dr - 1j*sz*dz + sy*(mu-A)/rho*v
                    residuals.append(float(mp.norm(Hf-E*v)/(E*mp.norm(v))))
                    n = -sign*q
                    gauge_errors.append(float(abs((mu+n)-(A+n)-(mu-A))))

    flux_rows = []
    for q, k in ((1, 0), (7, 0), (7, 3), (7, 6), (100000, 0), (100000, 16), (100000, 99999)):
        for sign in (-1, 1):
            radius = 2.3
            def integrand(theta):
                v = collar_mode(q, k, radius*np.sin(theta), radius*np.cos(theta), 1.3, sign)
                current = (np.sin(theta)*2*np.real(v[0].conjugate()*v[1])
                           + np.cos(theta)*(abs(v[0])**2-abs(v[1])**2))
                return float(radius*current)
            ts = betaincinv(k+1, q-k, np.array([1e-8, .001, .1, .5, .9, .999, 1-1e-8]))
            cuts = np.r_[0., 2*np.arcsin(np.sqrt(ts)), math.pi]
            flux = sum(quad(integrand, a, b, epsabs=1e-11, epsrel=1e-11)[0]
                       for a, b in zip(cuts[:-1], cuts[1:]))
            flux_rows.append(dict(q=q, k=k, charge=sign, radial_flux=flux))

    width = .7
    spectral = lambda E: 2*width**1.5*E*math.exp(-width*E)
    moments = [quad(lambda E: E**n*spectral(E)**2, 0, np.inf)[0] for n in range(3)]
    time_norm = quad(lambda u: abs(positive_pulse(u, width))**2, -np.inf, np.inf)[0]
    transform_errors = []
    for u in (-2., 0., .7, 3.):
        real = quad(lambda E: spectral(E)*math.cos(E*u), 0, np.inf, epsabs=1e-11)[0]
        imag = quad(lambda E: -spectral(E)*math.sin(E*u), 0, np.inf, epsabs=1e-11)[0]
        transform_errors.append(abs((real+1j*imag)/math.sqrt(2*math.pi)-positive_pulse(u, width)))

    # Actual partner-field defect on an upper matching sphere, not an escape rate.
    defects = []
    q, k, separation = 31, 4, 1.
    def beta_density(t):
        return math.exp(k*math.log(t)+(q-1-k)*math.log1p(-t)-betaln(k+1, q-k))
    for b in (.02, .01, .005):
        def squared_defect(t):
            rho, z = 2*b*math.sqrt(t*(1-t)), separation+b*(1-2*t)
            other_r = math.hypot(rho, z)
            # Rationalized 1-z/other_r avoids cancellation near the axis.
            delta_W = q*rho/(2*other_r*(other_r+z))
            return beta_density(t)*delta_W**2
        defect = math.sqrt(quad(squared_defect, 0, 1, epsabs=1e-13)[0])
        leading = q*b/(2*separation**2)*math.sqrt((k+1)*(q-k)/((q+1)*(q+2)))
        defects.append(dict(q=q, k=k, radius_over_d=b/separation,
                            rms_operator_defect=defect, leading_small_radius=leading,
                            defect_over_local_gap=defect*b/math.sqrt(q+1)))

    tests = {
        'monopole_stationary_dirac_equation': max(residuals) < 1e-35,
        'implemented_mode_matches_independent_formula': max(implementation_errors) < 2e-8,
        'gauge_covariance': max(gauge_errors) < 1e-35,
        'unit_radial_probability_flux': max(abs(r['radial_flux']-r['charge']) for r in flux_rows) < 2e-8,
        'positive_spectrum_norm_and_moments': max(abs(moments[i]-v) for i, v in enumerate((1., 3/(2*width), 3/width**2))) < 1e-10,
        'pulse_time_integrated_flux': abs(time_norm-1) < 1e-10,
        'direct_fourier_transform': max(transform_errors) < 1e-10,
        'partner_defect_asymptotic': abs(defects[-1]['rms_operator_defect']/defects[-1]['leading_small_radius']-1) < .02,
    }
    tests = {name: bool(passed) for name, passed in tests.items()}
    assert all(tests.values()), tests
    return dict(passed=True, tests=tests, max_dirac_relative_residual=max(residuals),
                max_implementation_relative_error=max(implementation_errors),
                max_gauge_error=max(gauge_errors), flux_checks=flux_rows,
                pulse=dict(width=width, spectral_norm=moments[0], mean_energy=moments[1],
                           second_energy_moment=moments[2], integrated_flux=time_norm,
                           max_fourier_error=float(max(transform_errors))),
                frozen_guide_gaussian_negative_energy_weights=[
                    dict(sigma=.12, carrier_energy=E, weight=float(.5*erfc(math.sqrt(2)*.12*E)))
                    for E in (0., 5., 15.)],
                partner_field_defects=defects,
                scope='Exact monopole reference modes and positive-energy pulse; no full dipole scattering solution.',
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def saved_collar_checks():
    """One flat point-pair magnetic term at the saved finest R=100 collar."""
    source = project_path('dirac_return_closure_results.json')
    saved = json.loads(source.read_text())
    if not saved['passed']:
        raise ValueError('Saved collar calculation did not pass its controls')
    record = max(saved['records'], key=lambda row: (row['order'], row['refinement']))
    root = next(row for row in record['roots'] if row['cutoff'] == 100.)
    background = root['background']
    cut = background['cuts']['100.0']
    q = saved['parameters']['q']
    d = saved['parameters']['distance_over_re']
    b, R = cut['isotropic_radius_from_collar'], root['cutoff']
    clock, lapse = background['clock'], cut['lapse']/background['clock']
    if not (isinstance(q, int) and q > 1 and 0 < b < d and R > 0 and clock > 0 and lapse > 0):
        raise ValueError('Invalid saved collar scales')
    eps = b/d
    scale = q*b/(2*d*d)
    flat_gap = math.sqrt(q+1)/b
    tetrad_factor = lapse*b/R
    energies = math.pi/(clock*np.asarray(root['periods']))
    weights = np.asarray(root['weights'])
    energy_scales = dict(minimum_sampled=float(energies.min()),
                         maximum_sampled=float(energies.max()),
                         weighted_rms=float(np.sqrt(np.dot(weights, energies**2))))
    rows = []
    for k in (0, 16, q//2, q-1):
        # t(1-t) tilts Beta(k+1,q-k) into Beta(k+2,q-k+1).
        a, beta_b = k+2, q-k+1
        total = a+beta_b
        mode = (a-1)/(total-2)
        std = math.sqrt(a*beta_b/(total*total*(total+1)))
        with mp.workdps(50):
            m = mp.mpf(a-1)/(total-2)
            log_peak = float((a-1)*mp.log(m)+(beta_b-1)*mp.log1p(-m)-mp.log(mp.beta(a,beta_b)))
        def density(y):
            offset = std*y
            return std*math.exp(log_peak+(a-1)*math.log1p(offset/mode)
                               +(beta_b-1)*math.log1p(-offset/(1-mode)))
        def correction_scaled(y):
            c = 1-2*(mode+std*y)
            D = math.sqrt(1+2*eps*c+eps*eps)
            dm1 = (2*eps*c+eps*eps)/(D+1)
            # Stable K^2-1 even when its expectation is only ~1e-12.
            return math.expm1(-2*math.log1p(dm1)-2*math.log1p((dm1+eps*c)/2))/eps
        quantiles = betaincinv(a, beta_b, np.array([1e-9, 1e-5, .01, .1, .5, .9, .99, 1-1e-5, 1-1e-9]))
        cuts = (np.r_[0., quantiles, 1.]-mode)/std
        # Endpoint modes have an enormous far-tail interval; resolve its near end.
        tails = []
        for inner, outer in ((-cuts[1], -cuts[0]), (cuts[-2], cuts[-1])):
            count = max(2, math.ceil(math.log2(outer/inner))+1)
            tails.append(np.geomspace(inner, outer, count))
        cuts = np.unique(np.r_[cuts, -tails[0], tails[1]])
        def integral(fn):
            values = [quad(fn, left, right, epsabs=1e-12, epsrel=1e-12, limit=150)
                      for left, right in zip(cuts[:-1], cuts[1:])]
            return sum(value[0] for value in values), sum(value[1] for value in values)
        norm, norm_error = integral(density)
        correction, quadrature_error = integral(lambda y: density(y)*correction_scaled(y))
        correction *= eps
        # Independent known Beta moments through eps^2; omitted order is eps^3.
        mean_c = (beta_b-a)/total
        mean_c2 = mean_c**2+4*a*beta_b/(total*total*(total+1))
        moment_correction = -4*eps*mean_c+eps*eps*(11.5*mean_c2-1.5)
        beta_moment = (k+1)*(q-k)/((q+1)*(q+2))
        leading = scale*math.sqrt(beta_moment)
        relative_shift = math.expm1(.5*math.log1p(correction))
        defect = leading*(1+relative_shift)
        rows.append(dict(k=k, flat_rms_defect=defect, leading_beta_moment=leading,
                         relative_correction_to_leading=relative_shift,
                         analytic_lower_bound=leading/(1+eps)**2,
                         analytic_upper_bound=leading/(1-eps)**2,
                         flat_defect_over_flat_gap=defect/flat_gap,
                         conditional_tetrad_local_defect=defect*b/R,
                         conditional_tetrad_infinity_defect=tetrad_factor*defect,
                         flat_defect_over_energy={key:defect/value for key,value in energy_scales.items()},
                         conditional_tetrad_defect_over_energy={key:tetrad_factor*defect/value for key,value in energy_scales.items()},
                         quadrature_tilted_beta_norm=norm, quadrature_norm_error_estimate=norm_error,
                         squared_correction=correction, beta_moment_squared_correction=moment_correction,
                         squared_correction_difference=abs(correction-moment_correction),
                         squared_correction_quadrature_error_estimate=eps*quadrature_error))
    uniform_bound = scale/(2*(1-eps)**2)
    tests = {
        'normalized_tilted_beta_quadrature': max(abs(row['quadrature_tilted_beta_norm']-1) for row in rows) < 2e-10,
        'independent_beta_moment_correction': max(row['squared_correction_difference'] for row in rows) < 2e-18,
        'exact_channel_and_uniform_bounds': all(row['analytic_lower_bound'] <= row['flat_rms_defect'] <= row['analytic_upper_bound'] <= uniform_bound for row in rows),
        'saved_energy_rms_matches_state': abs(energy_scales['weighted_rms']/(math.sqrt(root['omega_squared'])/(2*clock))-1) < 1e-10,
        'tetrad_gap_ratio_is_consistent': all(abs(row['conditional_tetrad_infinity_defect']/(lapse*math.sqrt(q+1)/R)/row['flat_defect_over_flat_gap']-1) < 1e-14 for row in rows),
        'all_reported_scales_are_finite': all(math.isfinite(value) and value > 0 for value in (uniform_bound,flat_gap,*energy_scales.values())),
    }
    tests = {key:bool(value) for key,value in tests.items()}
    assert all(tests.values()), tests
    return dict(passed=True, tests=tests, input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                units='re=1, hbar=c=1; energies in inverse re, time at exterior infinity unless explicitly local',
                saved_record=dict(order=record['order'], refinement=record['refinement'], cutoff=R),
                parameters=dict(q=q, isotropic_b=b, areal_R=R, separation_d=d, b_over_d=eps,
                                common_clock=clock, lapse_at_infinity_normalization=lapse,
                                central_optical_x=cut['optical_x']),
                flat_angular_gap=flat_gap, conditional_tetrad_local_gap=math.sqrt(q+1)/R,
                conditional_tetrad_infinity_gap=lapse*math.sqrt(q+1)/R,
                conditional_tetrad_infinity_factor=tetrad_factor,
                uniform_flat_operator_bound=uniform_bound,
                uniform_conditional_tetrad_infinity_bound=tetrad_factor*uniform_bound,
                energy_scales=energy_scales, channels=rows,
                scope='Only the flat point-pair magnetic residual; optional spherical tetrad conversion is conditional. Sampled loop energies do not bound the continuum or assign energy to a sphere channel.')


if __name__ == '__main__':
    result = checks()
    result['saved_collar'] = saved_collar_checks()
    project_path('incoming_modes_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2))
