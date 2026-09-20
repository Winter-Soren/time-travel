"""A prescribed smooth bridge retaining the exact balanced-dihole cut germs.

This is an off-shell geometry and a required-source calculation, not a matter
construction or a semiclassical solution.  Put R=r(s), chi=R', L=log(N/P),
where N(R,theta), P(R,theta) are the actual upper-mouth exterior fields.  For

  g=-N^2 dt^2+P^2[ds^2+R^2 dtheta^2]+R^2 sin(theta)^2/N^2 dphi^2,

the extra radial null Einstein source, in units Qkk=8*pi*G*T_extra,kk, is

  Qkk = [(1-chi^2)*Delta_S2 L/R^2 - R''*(2/R-L_R)]/P^2.

Indeed K_orig=(2/R-L_R)/P, Ric_h(nn)=-nK-|K|^2-Delta_gamma(P)/P,
and static R_tt/N^2+R_nn=Ric_h(nn)+[Delta_h N-Hess_nn N]/N.
The angular remainder simplifies to
Delta_gamma(N)/N-Delta_gamma(P)/P=Delta_S2 log(N/P)/(P^2 R^2).
The pulled-back magnetic field has F_sphi=chi*F_Rphi, so its radial null
source scales by chi^2.  Its Maxwell equation need not remain source-free.

The lapse and optical longitudinal coefficient remain angle dependent even
on the constant-R central cylinder.  Its coordinate length is not a single
optical length and does not imply the old scalar AdS propagation factor.
"""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import expit

from balanced_dihole import cut_profile, symbolic_fields
from project_paths import project_path


def exact_target_parameters():
    """Floating evaluations of explicitly linked exact balanced constants."""
    ctx = mp.mp.clone()
    ctx.dps = 70
    k, charge = ctx.mpf(1)/2, ctx.mpf(1)/5
    mass = (ctx.sqrt(k*k+4*charge*k)-k)/2
    a = ctx.sqrt(k*k-mass*mass)
    rp = mass+k
    field = mass/((k+a)*rp)
    return dict(Q=float(charge), D=1., k=float(k), M=float(mass),
                a=float(a), rplus=float(rp), B=float(field), q=1.,
                physical_A_factor=float(-1/(2*charge)))


def smooth_step(x):
    """Exact endpoint values and first derivative of the flat C-infinity step."""
    x = np.asarray(x, float)
    inside = (x > 0) & (x < 1)
    value = np.where(x >= 1, 1., 0.)
    first = np.zeros_like(x)
    t = x[inside]
    h = -1/t+1/(1-t)
    value[inside] = expit(h)
    first[inside] = expit(h)*expit(-h)*(1/t**2+1/(1-t)**2)
    return value, first


@lru_cache(maxsize=8)
def gauss(order):
    return np.polynomial.legendre.leggauss(order)


def radius_profile(u, r0=.025, rc=.05, order=64):
    """One outward half-collar; all cut jets equal the exterior coordinate R."""
    if not 0 < r0 < rc < .5:
        raise ValueError('Require 0<r0<rc<k=1/2.')
    delta = 2*(rc-r0)
    u = np.asarray(u, float)
    if np.any((u < 0) | (u > delta)):
        raise ValueError('Half-collar coordinate must lie in [0,delta].')
    t = u/delta
    nodes, weights = gauss(order)
    sample = t[..., None]*(nodes+1)/2
    integral = t/2*np.sum(smooth_step(sample)[0]*weights, axis=-1)
    # S(t)+S(1-t)=1 proves the full integral is exactly 1/2.
    integral = np.where(t == 1, .5, integral)
    chi, first = smooth_step(t)
    return r0+delta*integral, chi, first/delta


def log_speed_derivatives(radius, mu, parameters):
    """Actual L=log(N/P)=2log(R D)-2log(Sigma), differentiated locally."""
    radius, mu = np.broadcast_arrays(np.asarray(radius, float), np.asarray(mu, float))
    k, mass, ag = (parameters[key] for key in ('k', 'M', 'a'))
    if np.any((radius <= 0) | (radius >= k) | (abs(mu) > 1)):
        raise ValueError('Local coordinates must have 0<R<k and -1<=mu<=1.')
    distance = np.sqrt(radius**2+4*k*k+4*k*radius*mu)
    dr = (radius+2*k*mu)/distance
    dm = 2*k*radius/distance
    dmm = -(2*k*radius)**2/distance**3
    rb = mass+(radius+distance)/2
    c = (distance-radius)/(2*k)
    rr, cr = (1+dr)/2, (dr-1)/(2*k)
    rm, cm = dm/2, dm/(2*k)
    rmm, cmm = dmm/2, dmm/(2*k)
    sigma = rb*rb-ag*ag*c*c
    sr = 2*rb*rr-2*ag*ag*c*cr
    sm = 2*rb*rm-2*ag*ag*c*cm
    smm = 2*(rm*rm+rb*rmm)-2*ag*ag*(cm*cm+c*cmm)
    lr = 2/radius+2*dr/distance-2*sr/sigma
    lm = 2*dm/distance-2*sm/sigma
    lmm = 2*(dmm/distance-(dm/distance)**2-smm/sigma+(sm/sigma)**2)
    return dict(log_speed=2*np.log(radius*distance/sigma), radial=lr,
                mu_first=lm, mu_second=lmm,
                sphere_laplacian=(1-mu*mu)*lmm-2*mu*lm,
                global_r=rb, global_cosine=c)


def field_values(function, parameters, radius, mu):
    local = log_speed_derivatives(radius, mu, parameters)
    raw = function(local['global_r'], local['global_cosine'],
                   parameters['M'], parameters['a'], parameters['B'])
    local['lapse'], local['meridian_factor'] = np.asarray(raw[0]), np.asarray(raw[1])
    local['unit_A_phi'] = np.asarray(raw[2])*parameters['physical_A_factor']
    return local


def null_parts(radius, chi, second, meridian_factor, radial_log_speed, angular_log_speed):
    """Geometric identity, usable also for independent flat-metric controls."""
    radius = np.asarray(radius)
    angular = (1-np.asarray(chi)**2)*angular_log_speed/radius**2/meridian_factor**2
    opening = -np.asarray(second)*(2/radius-radial_log_speed)/meridian_factor**2
    return angular, opening, angular+opening


def null_residual(function, parameters, radius, chi, second, mu):
    fields = field_values(function, parameters, radius, mu)
    angular, opening, total = null_parts(radius, chi, second,
                                        fields['meridian_factor'], fields['radial'],
                                        fields['sphere_laplacian'])
    return fields, angular, opening, total


def integrated_budget(parameters, order_u=64, order_mu=32, radius_order=64):
    """One half-collar, using the exact angular divergence cancellation."""
    delta = .05
    nodes, weights = gauss(order_u)
    u, wu = delta*(nodes+1)/2, delta*weights/2
    mu, wm = gauss(order_mu)
    radius, chi, second = radius_profile(u, order=radius_order)
    derivatives = log_speed_derivatives(radius[:, None], mu[None, :], parameters)
    lr_average = derivatives['radial'] @ wm/2
    weighted_slice = -4*math.pi*second*(2*radius-radius**2*lr_average)
    return dict(order_u=order_u, order_mu=order_mu, radius_order=radius_order,
                half_collar_integral=float(wu @ weighted_slice),
                both_collars_integral=float(2*wu @ weighted_slice),
                largest_slice=float(np.max(weighted_slice)),
                most_negative_slice=float(np.min(weighted_slice)))


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    profile = result['bridge_profile']
    x = np.asarray(profile['s'])
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=True,
                             layout='constrained')
    axes[0].plot(x, profile['radius'], color='#2563a6', lw=2)
    axes[0].set_ylabel('Local mouth radius R / D')
    axes[0].set_title('A smooth bridge with the actual dihole mouth data', fontsize=12)
    axes[1].plot(x, profile['weighted_null_slice'], color='#a73e39', lw=2)
    axes[1].axhline(0, color='0.45', lw=.8)
    axes[1].set_ylabel('Required angular-integrated\nnull source (dimensionless)')
    axes[1].set_xlabel('Bridge coordinate s / D (not proper or optical distance)')
    for axis in axes:
        axis.axvspan(-.05, .05, color='#dfe8ea', alpha=.7)
        axis.axvline(-.05, color='0.55', lw=.8, ls='--')
        axis.axvline(.05, color='0.55', lw=.8, ls='--')
        axis.grid(alpha=.2)
    axes[0].text(0, .037, 'variable central segment\nangle-dependent lapse',
                 ha='center', va='center', fontsize=9)
    fig.suptitle('Prescribed geometry; supporting source unconstructed',
                 color='#8c332e', fontsize=11)
    png = project_path('dihole_stretched_collar_progress.png')
    pdf = project_path('dihole_stretched_collar_progress.pdf')
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)
    return [str(png), str(pdf)]


def run():
    dependency = project_path('balanced_dihole.py')
    before_hash = hashlib.sha256(dependency.read_bytes()).hexdigest()
    parameters = exact_target_parameters()
    function, _ = symbolic_fields()
    # Two small quadrature orders are a sensitivity check, not an enclosure.
    coarse = integrated_budget(parameters)
    fine = integrated_budget(parameters, 128, 64, 128)
    mu, weights = gauss(64)
    fractions = np.array([0., .125, .25, .5, .75, .875, 1.])
    radius, chi, second = radius_profile(.05*fractions)
    rows = []
    weighted_error = []
    for fraction, r, c, rr in zip(fractions, radius, chi, second):
        fields, angular, opening, total = null_residual(function, parameters, r, c, rr, mu)
        average_lr = float(weights @ fields['radial']/2)
        exact_slice = -4*math.pi*rr*(2*r-r*r*average_lr)
        direct_slice = 2*math.pi*float(weights @ (fields['meridian_factor']**2*r*r*total))
        weighted_error.append(abs(direct_slice-exact_slice))
        sample_mu = np.array([-.75, 0., .75])
        fs, a, o, q = null_residual(function, parameters, r, c, rr, sample_mu)
        rows.append(dict(u_over_delta=float(fraction), radius=float(r),
                         radius_prime=float(c), radius_second=float(rr),
                         mu=sample_mu.tolist(), lapse=fs['lapse'].tolist(),
                         meridian_factor=fs['meridian_factor'].tolist(),
                         null_angular=a.tolist(), null_opening=o.tolist(),
                         null_extra=q.tolist(),
                         weighted_null_slice=exact_slice,
                         directly_integrated_null_slice=direct_slice))

    # Verify the derivative formulas against differences of the actual N/P.
    derivative_errors = []
    for r, m in ((.027, -.6), (.033, .1), (.047, .7)):
        f = field_values(function, parameters, r, m)
        step = 2e-5
        def log_ratio(rr, mm):
            g = field_values(function, parameters, rr, mm)
            return float(np.log(g['lapse']/g['meridian_factor']))
        vr = [log_ratio(r+j*step, m) for j in (-2, -1, 1, 2)]
        vm = [log_ratio(r, m+j*step) for j in (-2, -1, 1, 2)]
        first_r = (vr[0]-8*vr[1]+8*vr[2]-vr[3])/(12*step)
        first_m = (vm[0]-8*vm[1]+8*vm[2]-vm[3])/(12*step)
        second_m = (-vm[3]+16*vm[2]-30*log_ratio(r, m)+16*vm[1]-vm[0])/(12*step**2)
        derivative_errors.append(dict(radius=r, mu=m,
                                      radial=abs(first_r-float(f['radial'])),
                                      mu_first=abs(first_m-float(f['mu_first'])),
                                      mu_second=abs(second_m-float(f['mu_second']))))

    poles = field_values(function, parameters, radius[:, None], np.array([-1., 1.])[None, :])
    axis_error = float(np.max(abs(poles['lapse']*poles['meridian_factor']-1)))
    cut = cut_profile(function, parameters, .05, np.linspace(-1, 1, 17))
    newcut = field_values(function, parameters, float(radius[-1]), cut['mu'])
    cut_error = max(float(np.max(abs(newcut[key]-cut[old])))
                    for key, old in (('lapse', 'lapse'), ('meridian_factor', 'meridian_factor'),
                                     ('unit_A_phi', 'Aphysical')))
    coordinate = np.linspace(-.1, .1, 321)
    u = np.clip(abs(coordinate)-.05, 0., .05)
    rp, cp, rpp = radius_profile(u)
    lp = log_speed_derivatives(rp[:, None], mu[None, :], parameters)['radial']
    profile_slice = -4*math.pi*rpp*(2*rp-rp*rp*(lp @ weights/2))
    finer_radius = radius_profile(.05*fractions, order=128)[0]
    radius_error = float(np.max(abs(radius-finer_radius)))
    budget_change = abs(coarse['both_collars_integral']-fine['both_collars_integral'])
    after_hash = hashlib.sha256(dependency.read_bytes()).hexdigest()
    flat_r = np.array([.3, .6])
    flat_second = np.array([.7, -.1])
    flat_result = null_parts(flat_r, np.array([.2, .8]), flat_second, 1., 0., 0.)[-1]
    flat_error = float(np.max(abs(flat_result+2*flat_second/flat_r)))
    tests = dict(
        flat_spherical_identity=flat_error < 2e-15,
        unchanged_germ_identity=all(max(abs(x) for x in rows[-1][key]) == 0.
                                    for key in ('null_angular', 'null_opening', 'null_extra')),
        smooth_endpoint_data=(radius[0] == .025 and radius[-1] == .05
                              and chi[0] == 0 and chi[-1] == 1
                              and second[0] == 0 and second[-1] == 0),
        axis_balance=axis_error < 2e-12,
        actual_cut_values=cut_error < 2e-12,
        actual_log_speed_derivatives=max(d['radial'] for d in derivative_errors) < 2e-7
            and max(d['mu_first'] for d in derivative_errors) < 2e-8
            and max(d['mu_second'] for d in derivative_errors) < 2e-5,
        angular_divergence_cancellation=max(weighted_error) < 2e-10,
        radius_quadrature_sensitivity=radius_error < 2e-12,
        null_budget_quadrature_sensitivity=budget_change < 2e-8,
        original_dependency_preserved=(before_hash == after_hash))
    tests = {name: bool(passed) for name, passed in tests.items()}
    result = dict(
        passed=all(tests.values()), tests=tests, parameters=parameters,
        exact_target='D=1,Q=1/5,k=1/2,M=(sqrt(k^2+4Qk)-k)/2,a=sqrt(k^2-M^2),B=M/((k+a)(M+k))',
        geometry=dict(cylinder_radius=.025, cut_radius=.05, half_collar_coordinate_width=.05,
                      displayed_central_coordinate_length=.1,
                      family='Insert any nonnegative central coordinate length with the same two fixed collars.',
                      smoothness='S(t)=exp(-1/t)/(exp(-1/t)+exp(-1/(1-t))) with flat extensions; R=r0+delta integral_0^(u/delta) S.',
                      matching='At u=delta, R=rc, R_prime=1 and higher derivatives vanish; at u=0, R=r0 and every derivative vanishes.',
                      reflection='Lower mouth is the z-reflection of the upper germ; the same core theta pulls back the lower global polar angle pi-theta.'),
        half_collar_samples=rows, quadrature=dict(coarse=coarse, fine=fine,
                        both_collars_change=budget_change, radius_change=radius_error,
                        angular_identity_max_error=max(weighted_error)),
        controls=dict(axis_balance_max_error=axis_error, cut_values_max_error=cut_error,
                      flat_spherical_identity_error=flat_error,
                      finite_difference_derivative_errors=derivative_errors),
        bridge_profile=dict(s=coordinate.tolist(), radius=rp.tolist(),
                            weighted_null_slice=profile_slice.tolist()),
        formulae=dict(Qkk='[(1-chi^2)*Delta_S2 log(N/P)/R^2 - R_second*(2/R-partial_R log(N/P))]/P^2',
                      weighted_slice='integral N sqrt(h) Qkk dtheta dphi = -4pi R_second [2R-R^2 average_S2 partial_R log(N/P)]',
                      source_normalization='Qkk=8pi G T_extra,kk for orthonormal k=e_t+e_s; the trace term in Einstein cancels.',
                      units='D is the length unit: Qkk has D^-2; weighted slice is dimensionless; collar integral has units D.',
                      Maxwell_divergence='For a=Aphi, div F^phi=[(chi^2-1)*partial_R(N^2 a_R)+R_second*N^2 a_R]/[P^2 R^2 (1-mu^2)].'),
        domain_and_scope=[
            'Prescribed static smooth off-shell bridge; no supporting matter, Einstein solution, stability or formation mechanism is supplied.',
            'Axis regularity uses the exact balanced relation NP=1 at both poles. Floating evaluations are approximations to linked exact constants.',
            'A=q a_unit dphi uses integer q, with q(a_unit-1)dphi regular at the other polar patch; the opposite outward mouth orientations give opposite fluxes.',
            'The pulled-back field is smooth but generally has the displayed nonzero Maxwell current. A physical current and its stress remain unconstructed.',
            'A spin structure, fixed Wilson holonomy and static state prescription must be selected before a quantum calculation; no state or determinant is computed here.',
            'The cylinder keeps angle-dependent N and P. Its variable coordinate length is neither proper nor a single optical length; no scalar AdS phase or MMP energy term is inherited.',
            'The integrated null target concerns this geometry and magnetic field. Its angular average vanishes in the central cylinder, but its pointwise null residual need not vanish there.',
            'Saved numbers and coarse/fine changes are floating numerical diagnostics, not interval-certified bounds or proof of a sign at unsampled points.',
            'No chronological shift is present: this is a static geometric support target, not a backward-time passage.'
        ],
        sources=['https://arxiv.org/abs/hep-th/9906160'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={str(dependency):before_hash})
    result['figures'] = render(result)
    output = project_path('dihole_stretched_collar_results.json')
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests,
                         both_collars_integral=fine['both_collars_integral'],
                         quadrature_change=budget_change, result=str(output)), indent=2))
    if not result['passed']:
        raise RuntimeError('One or more focused bridge controls failed.')
    return result


if __name__ == '__main__':
    run()
