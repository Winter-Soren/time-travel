"""Finite flux contrasts of the parity-even minimal Dirac determinant on H2 x S2.

This homogeneous product benchmark is not a wormhole state or a source tensor.
Run with the project's .venv Python; only the sibling result file is written.
"""
from project_paths import project_path

from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_genlaguerre, roots_laguerre, roots_legendre


@lru_cache(maxsize=16384)
def _angular_mp(t):
    """Independent high-precision spectral sums, including their t derivatives.

    Both q=2 and q=3 share the q=0,1 subtraction. Gaussian factors are updated
    multiplicatively, with a cutoff exp(-100), in 45-digit arithmetic.
    """
    with mp.workdps(45):
        zt = mp.mpf(float(t))
        step = mp.exp(-zt)
        step2 = step * step
        gaussian = step
        gaussian_step = step**3
        linear = step
        values = [mp.mpf(q) for q in range(4)]
        derivatives = [mp.mpf(0) for _ in range(4)]
        for n in range(1, int(math.ceil(math.sqrt(100 / t))) + 3):
            factor = gaussian
            for q in range(4):
                term = 2 * (q + 2 * n) * factor
                values[q] += term
                derivatives[q] -= n * (n + q) * term
                factor *= linear
            gaussian *= gaussian_step
            gaussian_step *= step2
            linear *= step
        return tuple(
            (float(values[q] - q*q*values[1] + (q*q-1)*values[0]),
             float(derivatives[q] - q*q*derivatives[1] + (q*q-1)*derivatives[0]))
            for q in (2, 3)
        )


def angular_contrast(t, q):
    """Return D_q A and its derivative without double-precision UV cancellation."""
    if q not in (2, 3) or t <= 0:
        raise ValueError('This bounded benchmark supports q=2,3 and t>0')
    if t < .25:
        return _angular_mp(float(t))[q-2]
    n = np.arange(1, int(math.ceil(math.sqrt(100/t))) + 3, dtype=float)
    values = []
    derivatives = []
    for charge in (q, 1, 0):
        eigenvalue = n * (n + charge)
        terms = 2 * (charge + 2*n) * np.exp(-t * eigenvalue)
        values.append(charge + float(np.sum(terms)))
        derivatives.append(-float(np.dot(eigenvalue, terms)))
    return (values[0] - q*q*values[1] + (q*q-1)*values[2],
            derivatives[0] - q*q*derivatives[1] + (q*q-1)*derivatives[2])


@lru_cache(maxsize=16)
def _heat_nodes(order):
    return roots_laguerre(order), roots_genlaguerre(order, -.5)


def hyperbolic_heat(t, order=80):
    """H and H' using thermal moments or a rescaled Gaussian integral."""
    if t <= 0:
        raise ValueError('Require t>0')
    thermal, gaussian = _heat_nodes(order)
    if t <= 1:
        s, w = thermal
        factor = np.exp(-t*s*s/(4*math.pi**2)) / (-np.expm1(-s))
        h = 1/(2*t) + float(np.dot(w, s*factor))/(2*math.pi**2)
        hp = -1/(2*t*t) - float(np.dot(w, s**3*factor))/(8*math.pi**4)
    else:
        v, w = gaussian
        z = math.pi*np.sqrt(v/t)
        regular = z / np.tanh(z)
        h = float(np.dot(w, regular))/(2*math.pi*math.sqrt(t))
        hp = -float(np.dot(w, v*regular))/(2*math.pi*t**1.5)
    return h, hp


def _uv_coefficients(q):
    common = q*q*(q*q-1)
    return {3: -common/360, 4: -common/504,
            5: common*(q*q-20)/15120}


def _uv_integrals(q, cutoff, eta):
    # This truncation, not an asserted remainder bound, is tested by shrinking
    # the cutoff. H coefficients through t^3 suffice at these small cutoffs.
    hc = {-1: .5, 0: 1/12, 1: -1/120, 2: 1/504, 3: -1/1440}
    w = direct = ibp = 0.
    for m, hm in hc.items():
        for n, dn in _uv_coefficients(q).items():
            term = hm*dn*cutoff**(m+n)/(m+n)
            w += term*math.exp(2*eta*(n-m))/(8*math.pi**2)
            direct += (n-m)*term/(4*math.pi**2)
            ibp -= m*term/(2*math.pi**2)
    return np.array([w, direct, ibp])


def _log_nodes(lower, upper, panel_width, order):
    boundaries = np.linspace(math.log(lower), math.log(upper),
                             int(math.ceil(math.log(upper/lower)/panel_width))+1)
    nodes, weights = roots_legendre(order)
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        for x, w in zip((lo+hi)/2+(hi-lo)*nodes/2, (hi-lo)*weights/2):
            yield math.exp(x), w


def flux_integrals(q, config, eta=0.):
    """Action contrast plus two eta derivatives at eta=0.

    The entire IR is included by x=1/sqrt(t); there is no finite upper cutoff.
    Derivative columns always refer to eta=0, even for an action at nonzero eta.
    """
    uv = _uv_integrals(q, config['uv_cutoff'], eta)
    core = np.zeros(3)
    ir = np.zeros(3)
    ep = math.exp(2*eta)
    em = 1/ep
    def terms(t):
        h, hp = hyperbolic_heat(t, config['heat_order'])
        d, dp = angular_contrast(t, q)
        if eta:
            action_h = hyperbolic_heat(t*em, config['heat_order'])[0]
            action_d = angular_contrast(t*ep, q)[0]
        else:
            action_h, action_d = h, d
        return np.array([action_h*action_d/(8*math.pi**2),
                         t*(-hp*d+h*dp)/(4*math.pi**2),
                         -t*hp*d/(2*math.pi**2)])
    for t, weight in _log_nodes(config['uv_cutoff'], config['ir_start'],
                               config['log_panel_width'], config['time_order']):
        core += weight*terms(t)
    nodes, weights = roots_legendre(config['ir_order'])
    xmax = 1/math.sqrt(config['ir_start'])
    for x, weight in zip(xmax*(nodes+1)/2, xmax*weights/2):
        ir += weight*2*terms(1/(x*x))/x
    total = uv+core+ir
    return {'action': float(total[0]), 'derivative_direct': float(total[1]),
            'derivative_integration_by_parts': float(total[2]),
            'contributions_uv_core_ir': [uv.tolist(), core.tolist(), ir.tolist()]}


def _independent_heat(t):
    # Separate adaptive integration of the defining spectral density.
    def base(u):
        if u == 0:
            return math.sqrt(t)/math.pi
        return u/math.tanh(math.pi*u/math.sqrt(t))*math.exp(-u*u)
    h, eh = quad(base, 0, 12, epsabs=2e-13, epsrel=2e-13)
    hp, ep = quad(lambda u: u*u*base(u), 0, 12, epsabs=2e-13, epsrel=2e-13)
    return h/t, -hp/t**2, eh/t, ep/t**2


def _independent_angular(t, q):
    # Direct mp.exp sums deliberately do not use the production recurrence.
    with mp.workdps(55):
        z = mp.mpf(str(t))
        count = int(mp.ceil(mp.sqrt(140/z)))+2
        values = []
        derivatives = []
        for charge in (q, 1, 0):
            terms = [2*(charge+2*n)*mp.exp(-z*n*(n+charge)) for n in range(1,count+1)]
            values.append(charge+mp.fsum(terms))
            derivatives.append(-mp.fsum(n*(n+charge)*term for n,term in enumerate(terms,1)))
        return (float(values[0]-q*q*values[1]+(q*q-1)*values[2]),
                float(derivatives[0]-q*q*derivatives[1]+(q*q-1)*derivatives[2]))


def dirac_flux_response_checks():
    configs = {
        'coarse': dict(uv_cutoff=.003, ir_start=24., log_panel_width=1.,
                       time_order=8, heat_order=48, ir_order=24),
        'main': dict(uv_cutoff=.001, ir_start=32., log_panel_width=.75,
                     time_order=12, heat_order=80, ir_order=36),
        'fine': dict(uv_cutoff=.0005, ir_start=48., log_panel_width=.5,
                     time_order=16, heat_order=128, ir_order=48),
    }
    heat_controls = []
    for t in (.001, .03, .25, 1., 1.000001, 10., 100.):
        reference = _independent_heat(t)
        value = hyperbolic_heat(t, configs['fine']['heat_order'])
        heat_controls.append({'t': t, 'H': value[0], 'H_prime': value[1],
                              'absolute_errors': [abs(value[i]-reference[i]) for i in (0,1)],
                              'relative_errors': [abs(value[i]/reference[i]-1) for i in (0,1)],
                              'adaptive_quadrature_error_estimates': list(reference[2:])})
    angular_controls = []
    for q in (2, 3):
        for t in (.001, .003, .03, .249999, .25, 1.):
            value = angular_contrast(t,q)
            reference = _independent_angular(t,q)
            angular_controls.append({'q': q, 't': t, 'D': value[0], 'D_prime': value[1],
                                     'absolute_errors': [abs(value[i]-reference[i]) for i in (0,1)]})
    cases = []
    for q in (2, 3):
        levels = {name: flux_integrals(q, cfg) for name,cfg in configs.items()}
        finite_differences = []
        for step in (.002, .001, .0005):
            plus = flux_integrals(q, configs['fine'], step)['action']
            minus = flux_integrals(q, configs['fine'], -step)['action']
            finite_differences.append({'eta_step': step, 'action_plus': plus,
                                      'action_minus': minus, 'central_derivative': (plus-minus)/(2*step)})
        richardson = [(4*finite_differences[i+1]['central_derivative']-
                       finite_differences[i]['central_derivative'])/3 for i in (0,1)]
        exact_response = levels['fine']['derivative_integration_by_parts']
        local = -q*q*math.log(q)/(24*math.pi**2)
        repaired = local+(q-q*q)/(24*math.pi**2)
        constant_ir = (q-q*q)*hyperbolic_heat(configs['fine']['ir_start'],128)[0]/(2*math.pi**2)
        cases.append({'q': q, 'configurations': levels,
                      'finite_difference': finite_differences, 'richardson_derivatives': richardson,
                      'response': exact_response,
                      'response_main_to_fine_absolute_change': abs(exact_response-levels['main']['derivative_integration_by_parts']),
                      'action_main_to_fine_absolute_change': abs(levels['fine']['action']-levels['main']['action']),
                      'direct_vs_ibp_absolute_difference': abs(exact_response-levels['fine']['derivative_direct']),
                      'richardson_vs_ibp_absolute_difference': abs(exact_response-richardson[-1]),
                      'uv_coefficients_t3_t4_t5': list(_uv_coefficients(q).values()),
                      'ir_response_constant_D_control': constant_ir,
                      'ir_response_minus_constant_control': levels['fine']['contributions_uv_core_ir'][2][2]-constant_ir,
                      'local_log_action': -q*q*math.log(q)/(96*math.pi**2),
                      'local_log_response': local,
                      'local_plus_minimal_trace_repair_response': repaired,
                      'response_minus_local': exact_response-local,
                      'response_minus_local_plus_trace_repair': exact_response-repaired,
                      'response_to_local_ratio': exact_response/local})
    thresholds = {'heat_relative_error': 2e-10, 'angular_absolute_error': 2e-11,
                  'response_refinement_absolute': 2e-9, 'action_refinement_absolute': 2e-9,
                  'direct_vs_ibp_absolute': 2e-9, 'richardson_vs_ibp_absolute': 2e-9}
    tests = {
        'heat_defining_integral': max(max(r['relative_errors']) for r in heat_controls)<thresholds['heat_relative_error'],
        'angular_direct_spectrum': max(max(r['absolute_errors']) for r in angular_controls)<thresholds['angular_absolute_error'],
        'response_refinement': all(c['response_main_to_fine_absolute_change']<thresholds['response_refinement_absolute'] for c in cases),
        'action_refinement': all(c['action_main_to_fine_absolute_change']<thresholds['action_refinement_absolute'] for c in cases),
        'integration_by_parts': all(c['direct_vs_ibp_absolute_difference']<thresholds['direct_vs_ibp_absolute'] for c in cases),
        'finite_difference': all(c['richardson_vs_ibp_absolute_difference']<thresholds['richardson_vs_ibp_absolute'] for c in cases),
    }
    return {
        'passed': all(tests.values()), 'tests': tests, 'thresholds': thresholds,
        'normalization': 'One minimal massless 4D Dirac field; unit reference radii; dimensionless action density per homogeneous 4-volume, with fixed a^2*r^2=1.',
        'effective_action_convention': 'Real parity-even W_E=-1/2 log det(D^dagger D)=+1/2 integral dt K(t)/t; no determinant-phase claim.',
        'equations': {
            'H': 'H(t)=integral_0^infinity lambda*coth(pi*lambda)*exp(-t*lambda^2) dlambda',
            'A': 'A_q(t)=q+2*sum_{n>=1}(q+2n)*exp(-t*n*(n+q))',
            'contrast': 'D_q A=A_q-q^2*A_1+(q^2-1)*A_0',
            'deformation': 'a=exp(eta), r=exp(-eta); w_q=integral dt/t H(t*exp(-2eta))*D_q A(t*exp(2eta))/(8*pi^2)',
            'response': "w'_q(0)=integral dt[-H'*D+H*D']/(4*pi^2)=-integral dt H'*D/(2*pi^2)",
            'boundary_terms': 'H*D tends to zero as t->0 and t->infinity; D begins at t^3 and tends to q-q^2.',
            'local_comparator': 'local_R=-q^2*log(q)/(24*pi^2); local_R+(q-q^2)/(24*pi^2) adds the specified minimal homogeneous trace-repair comparator.',
        },
        'numerics': {'configurations': configs, 'angular_precision_digits_below_point25': 45,
                     'angular_exponent_cutoff': 100, 'independent_angular_precision_digits': 55,
                     'uv_treatment': 'Audited D series through t^5 and H through t^3, integrated analytically; varying cutoff controls but does not rigorously bound omitted terms.',
                     'ir_treatment': 'Full spectrum retained; x=1/sqrt(t) compactifies the complete infinite tail. Constant-D analytic control is separate.',
                     'error_scope': 'Observed numerical convergence and independent integral comparisons, not interval-certified total error.'},
        'heat_controls': heat_controls, 'angular_controls': angular_controls, 'cases': cases,
        'external_independent_q2_response': {'value': -.03242773217598199,
            'method': 'Separate parent calculation: 35-digit direct angular sums, adaptive Gaussian H moments, log-time integration and independent tail treatment.',
            'absolute_difference': abs(cases[0]['response']+.03242773217598199)},
        'external_independent_q3_response': {'value': -.1031821229596531,
            'method': 'Separate parent calculation: 35-digit direct angular sums, adaptive Gaussian H moments, log-time integration and independent tail treatment.',
            'absolute_difference': abs(cases[1]['response']+.1031821229596531)},
        'scope': [
            'Homogeneous Euclidean H2 x S2 minimal Dirac spectrum with no electric field and standard H2 continuum density.',
            'Finite counterterm-canceling flux contrast at q=2,3, not the absolute determinant at a single flux.',
            'A product-geometry anisotropic shape response, not a Lorentzian quantum state, full stress tensor, mouth matching, or loop Casimir.',
            'Low-flux comparison does not supply a high-flux error bar or rule out other nonlocal/Wilsonian completions.',
            'A difference from the specified local plus minimal trace-repair approximation rejects that approximation as an exact product contrast only.',
        ],
        'source': {'url': 'https://arxiv.org/html/2603.15730v1',
                   'title': 'Exact Path Integral Methods in Supersymmetric AdS2 x S2 Backgrounds',
                   'equations': ['3.48','3.49','3.58','3.59'],
                   'mapping': 'Use minimally coupled spin-1/2 spectrum in Sec.3.2.1, paper monopole g=q/2 and electric e=0; do not use later nonminimal supersymmetric operator.'},
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == '__main__':
    started = time.perf_counter()
    result = dirac_flux_response_checks()
    destination = project_path('dirac_flux_response_results.json')
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'passed': result['passed'], 'tests': result['tests'],
                      'cases': [{k:c[k] for k in ('q','response','response_main_to_fine_absolute_change','direct_vs_ibp_absolute_difference','richardson_vs_ibp_absolute_difference')} for c in result['cases']],
                      'seconds': time.perf_counter()-started}, indent=2))
    if not result['passed']:
        raise SystemExit('Numerical controls failed; diagnostic results were retained.')
