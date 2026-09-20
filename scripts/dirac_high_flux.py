"""All-high-flux pressure contrast on the same fixed H2 x S2 vacuum as increment81.

This matches a vacuum pressure combination, not a varying throat or loop state.
"""
from project_paths import project_path
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from dirac_flux_response import hyperbolic_heat, _log_nodes


WEIGHTS = ((1, 5), (2, -8), (3, 3))
C0 = -(27*math.log(3)-32*math.log(2))/(24*math.pi**2)
C1_H = -1/(12*math.pi**2)
C1_S = -3*float(mp.zeta(3))/(4*math.pi**4)


@lru_cache(maxsize=10000)
def angular_scaled(y, q, digits=40, exponent_cut=100):
    """B_q(y/q)/q and its y derivative, retaining the exact angular spectrum.

    Windows longdouble equals double here: high precision prevents cancellation
    between the three large traces. Recurrences avoid one exp call per mode.
    """
    if type(q) is not int or q < 1 or y <= 0:
        raise ValueError('Require integer q>=1 and y>0')
    with mp.workdps(digits):
        yy, qq = mp.mpf(y), mp.mpf(q)
        value, derivative = mp.mpf(0), mp.mpf(0)
        for p, weight in WEIGHTS:
            count = math.ceil(2*exponent_cut/(y*(p+math.sqrt(p*p+4*exponent_cut/(q*y)))))+2
            factor = mp.exp(-yy*(p+1/qq))
            ratio = mp.exp(-yy*(p+3/qq))
            ratio_step = mp.exp(-2*yy/qq)
            a, ap = mp.mpf(p), mp.mpf(0)
            for n in range(1, count+1):
                term = 2*(p+2*n/qq)*factor
                a += term
                ap -= n*(p+n/qq)*term
                factor *= ratio
                ratio *= ratio_step
            value += weight*a
            derivative += weight*ap
        return float(value), float(derivative)


def uv_integrals(q, cutoff):
    # ponytail: a finite audited series suffices at this cutoff; decrease it
    # and compare rather than introducing a general symbolic-series engine.
    bc = {3: -1/3, 4: -5/(21*q), 5: 1/9-1/(6*q*q)}
    hc = {-1: .5, 0: 1/12, 1: -1/120, 2: 1/504, 3: -1/1440}
    direct = ibp = 0.
    for m, hm in hc.items():
        for n, bn in bc.items():
            integral = hm*bn*q**(-m-1)*cutoff**(m+n)/(m+n)
            direct += (n-m)*integral/(4*math.pi**2)
            ibp -= m*integral/(2*math.pi**2)
    return np.array([direct, ibp])


def response(q, *, lower=.02, upper=40., order=12, width=.75, heat_order=80):
    """Return response/q^2 using direct variation and integration by parts."""
    core = np.zeros(2)
    for y, weight in _log_nodes(lower, upper, width, order):
        b, bp = angular_scaled(y, q)
        h, hp = hyperbolic_heat(y/q, heat_order)
        core += weight*y*np.array([-hp*b+q*h*bp, -2*hp*b])/(4*math.pi**2*q*q)
    uv = uv_integrals(q, lower)
    h_end = hyperbolic_heat(upper/q, heat_order)[0]
    tail_ibp = -h_end/(math.pi**2*q)
    tail = np.array([tail_ibp/2, tail_ibp])
    # Excited angular modes beyond upper are bounded by geometric sums.
    # H decreases, so sup|B+2q| times integral(-H') bounds the IBP tail error.
    b_tail_bound = 0.
    for p, weight in WEIGHTS:
        z = math.exp(-p*upper)
        b_tail_bound += abs(weight)*(2*p*z/(1-z)+4*z/(q*(1-z)**2))
    tail_error_bound = h_end*b_tail_bound/(2*math.pi**2*q)
    total = uv+core+tail
    normalized = float(total[1])
    return {'q': q, 'response_over_q2': normalized,
            'response': normalized*q*q,
            'direct_over_q2': float(total[0]),
            'direct_vs_ibp_over_q2': abs(float(total[0]-total[1])),
            'uv_core_ir_over_q2': [uv.tolist(), core.tolist(), tail.tolist()],
            'excited_tail_error_bound_over_q2': tail_error_bound,
            'curvature_coefficient': (normalized-C0)*q,
            'missing_after_minimal_repair_over_q': (normalized-C0)*q-C1_H,
            'remainder_after_two_term_expansion': (normalized-C0-(C1_H+C1_S)/q)*q*q}


def direct_angular_check(y, q):
    """Independent nonrecursive sums at selected points."""
    with mp.workdps(55):
        yy, qq = mp.mpf(y), mp.mpf(q)
        a, ap = mp.mpf(0), mp.mpf(0)
        for p, weight in WEIGHTS:
            count = math.ceil(2*140/(y*(p+math.sqrt(p*p+4*140/(q*y)))))+2
            terms = [2*(p+2*n/qq)*mp.exp(-yy*n*(p+n/qq)) for n in range(1,count+1)]
            a += weight*(p+mp.fsum(terms))
            ap -= weight*mp.fsum(n*(p+n/qq)*term for n,term in enumerate(terms,1))
        return float(a), float(ap)


def uv_truncation_bounds():
    """Analytic EM8/thermal-series bounds for the actual q=100000 UV pieces.

    The bound is evaluated in high precision, not interval arithmetic. It
    excludes rounding, the main quadrature, and the separate infrared tail.
    """
    with mp.workdps(85):
        q, degree = mp.mpf(100000), 8
        polynomial, remainder = {}, {}
        for m in range(2, degree+1):
            for j in range(m-1):
                n, power = 2*m-j-1, 2*m-2*j
                spectral = sum(c*p**power for p,c in WEIGHTS)
                coefficient = 2*mp.bernoulli(2*m)*(-1)**j*spectral/(q**j*mp.factorial(j)*mp.factorial(power))
                polynomial[n] = polynomial.get(n,0)+coefficient
        prefactor = 4*mp.zeta(2*degree)/(2*mp.pi)**(2*degree)
        for j in range(degree+1):
            r = degree-j
            for k in range(r+1):
                n = 2*degree-j-k-1
                spectral = sum(abs(c)*p**(2*(r-k)) for p,c in WEIGHTS)
                coefficient = prefactor*mp.factorial(2*degree+1)*mp.factorial(r)*(4/q)**k*spectral
                coefficient /= mp.factorial(j)*mp.factorial(2*degree+1-2*j)*q**j*mp.factorial(r-k)
                remainder[n] = remainder.get(n,0)+coefficient
        out = []
        for cut in ('.01','.02'):
            c = mp.mpf(cut)
            bound = mp.mpf(0)
            for n in range(6,2*degree):
                cn = abs(polynomial.get(n,0))+remainder.get(n,0)
                bound += cn*(c**(n-1)/(4*mp.pi**2*(n-1))+c**(n+1)/(240*mp.pi**2*q*q*(n+1)))
            bound += sum(abs(polynomial[n])*c**(n+4)/(n+4) for n in (3,4,5))/(1584*mp.pi**2*q**5)
            out.append({'q':100000, 'cutoff':float(c),
                        'ibp_normalized_truncation_bound':math.nextafter(float(bound),math.inf),
                        'linear_coefficient_truncation_bound':math.nextafter(float(q*bound),math.inf),
                        'scope':'UV truncation only; 85-digit non-interval evaluation of an analytic bound.'})
        return out


def high_flux_checks():
    configs = {
        'main': dict(lower=.02, upper=40., order=12, width=.75, heat_order=80),
        'fine': dict(lower=.01, upper=48., order=16, width=.5, heat_order=128),
    }
    cases = []
    for q in (1, 100, 1000, 10000, 100000):
        main = response(q, **configs['main'])
        fine = response(q, **configs['fine'])
        cases.append({'q': q, 'main': main, 'fine': fine,
                      'normalized_refinement_change': abs(fine['response_over_q2']-main['response_over_q2']),
                      'linear_coefficient_refinement_change': abs(fine['curvature_coefficient']-main['curvature_coefficient'])})
    angular_controls = []
    for q, y in ((1,.03),(100,.03),(100000,.03),(100000,.2),(100000,1.)):
        value = angular_scaled(y,q)
        reference = direct_angular_check(y,q)
        angular_controls.append({'q':q, 'y':y, 'value':value,
                                 'absolute_errors':[abs(a-b) for a,b in zip(value,reference)]})
    old = json.loads(project_path('dirac_flux_response_results.json').read_text())
    small_q_expected = 3*old['cases'][1]['response']-8*old['cases'][0]['response']
    small_q_error = abs(cases[0]['fine']['response']-small_q_expected)
    high = cases[-1]['fine']
    uv_bounds = uv_truncation_bounds()
    tests = {
        'flux_counterterm_cancellation': sum(c for p,c in WEIGHTS)==0 and sum(c*p*p for p,c in WEIGHTS)==0,
        'independent_angular_sums': max(max(c['absolute_errors']) for c in angular_controls)<1e-13,
        'saved_low_flux_crosscheck': small_q_error<2e-10,
        'direct_variation_vs_ibp': max(c['fine']['direct_vs_ibp_over_q2'] for c in cases)<2e-10,
        'normalized_refinement': max(c['normalized_refinement_change'] for c in cases)<2e-10,
        'linear_coefficient_resolved': cases[-1]['linear_coefficient_refinement_change']<2e-7,
        'uv_truncation_resolves_linear_term': uv_bounds[0]['linear_coefficient_truncation_bound']<1e-8,
        'high_flux_asymptotic_agreement': abs(high['curvature_coefficient']-C1_H-C1_S)<2e-5,
    }
    return {'passed':all(tests.values()), 'tests':tests, 'configurations':configs,
            'contrast':'3 W(3q) - 8 W(2q) + 5 W(q), same homogeneous minimal Dirac vacuum',
            'observable':'Fixed-volume shape derivative, a=exp(eta), r=exp(-eta), at eta=0; units hbar*c/L^4.',
            'asymptotic':{'quadratic_coefficient':C0, 'linear_hyperbolic':C1_H,
                          'linear_spherical':C1_S, 'linear_total':C1_H+C1_S},
            'cases':cases, 'angular_controls':angular_controls,
            'uv_truncation_bounds':uv_bounds,
            'low_flux_control':{'from_increment81':small_q_expected, 'absolute_error':small_q_error},
            'scope':['Counterterms independent of flux or proportional to q^2 cancel; all three fluxes are high at q=100000.',
                     'Exact angular spectrum retained; UV series and quadrature have observed convergence, not a certified total error bound.',
                     'Analytic excited-mode IR tail bound is separate from unbounded quadrature/UV truncation errors.',
                     'This homogeneous product vacuum has zero radial null stress; not a returning-loop state or varying-throat source.',
                     'A matching coefficient is not a new physical law or a construction recipe.'],
            'sources':['https://arxiv.org/pdf/2603.15730','https://arxiv.org/pdf/0906.2430'],
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'reused_source_sha256':hashlib.sha256(project_path('dirac_flux_response.py').read_bytes()).hexdigest()}


if __name__ == '__main__':
    result = high_flux_checks()
    project_path('dirac_high_flux_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':result['passed'], 'tests':result['tests'],
                      'asymptotic':result['asymptotic'],
                      'cases':[{k:c['fine'][k] for k in ('q','response','curvature_coefficient','remainder_after_two_term_expansion')} for c in result['cases']]},indent=2))
    if not result['passed']:
        raise SystemExit('One or more declared controls failed; results retained for inspection.')
