"""Exact-arithmetic angular quadrature error for the projected magnetic tail.

The bound compares exact Gauss-Legendre quadrature to the exact finite-J
Galerkin integral. Floating nodes, harmonics, matrix arithmetic and omitted
angular channels remain separate. No exterior or radial solve is performed.
"""
from project_paths import project_path
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.special import eval_gegenbauer

from magnetic_tail_map import _AngularCoupling


INPUTS = ('magnetic_tail_map.py', 'spherical_weyl_leads.py',
          'threshold_return_bound.py', 'threshold_return_bound_results.json')


def rational_bounds(radius=2, separation=1, charge=1, order=64, cutoff=9):
    """Positive rational bounds on C, C', C'' projection errors at r>=R."""
    radius, separation, charge = map(Fraction, (radius, separation, charge))
    if (radius <= separation/2 or separation <= 0 or charge < 0
            or isinstance(order, bool) or isinstance(cutoff, bool)
            or not isinstance(order, int) or not isinstance(cutoff, int)
            or cutoff < 1 or order < cutoff+1):
        raise ValueError('Require R>d/2>0, q>=0, integers Q>=J+1 and J>=1')
    M = order-cutoff-1
    u = separation/(2*radius)
    s, gap, a = u*u, 1-u*u, M+2
    # Tail moments of (j+1)^p, j=M+1,..., after factoring u^(2M+3).
    t0 = 1/gap
    t1 = s/gap**2
    t2 = s*(1+s)/gap**3
    t3 = s*(1+4*s+s*s)/gap**4
    moments = (a*t0+t1, a*a*t0+2*a*t1+t2,
               a**3*t0+3*a*a*t1+3*a*t2+t3)
    factor = charge*u**(2*M+3)
    scalar = factor*moments[0]
    return dict(M=M, maximum_integrand_degree=2*M+2*cutoff+1,
        GL_exact_degree=2*order-1, f_error=scalar,
        C_error=2*scalar/radius,
        C_prime_error=4*factor*moments[1]/radius**2,
        C_second_error=factor*(8*moments[2]+4*moments[1])/radius**3,
        integrated_C_error=2*charge*Fraction(M+2,2*M+3)
                           *u**(2*M+3)/(1-u*u)**2)


def outward(value):
    lo = math.nextafter(float(value), -math.inf)
    hi = math.nextafter(float(value), math.inf)
    assert Fraction.from_float(lo) <= value <= Fraction.from_float(hi)
    return dict(lower=lo, upper=hi,
                exact_numerator=str(value.numerator), exact_denominator=str(value.denominator))


def series_controls():
    ratios, endpoint_errors, derivative_errors = [], [], []
    with mp.workdps(110):
        for radius, M in ((2,2), (2,6), (2,54), (mp.mpf('1.25'),6)):
            u = 1/(2*radius)
            u = mp.mpf(u)  # The selected binary fractions are exact.
            eps = u**(2*M+3)*((M+2)-(M+1)*u*u)/(1-u*u)**2
            for z in map(mp.mpf, ('-1','-.8','0','.3','1')):
                if abs(z)==1:
                    actual = u/(1-u*u)**2
                else:
                    actual = ((z+u)/mp.sqrt(1+2*u*z+u*u)
                              -(z-u)/mp.sqrt(1-2*u*z+u*u))/(2*(1-z*z))
                polynomial = sum(mp.gegenbauer(2*j,mp.mpf('1.5'),z)*u**(2*j+1)/(2*j+1)
                                 for j in range(M+1))
                ratio = abs(actual-polynomial)/eps
                ratios.append(float(ratio))
                if abs(z)==1: endpoint_errors.append(float(abs(ratio-1)))
        # Independent differentiation of the closed positive endpoint remainder.
        for M in (0,2,6,54):
            def envelope(r):
                u = 1/(2*r)
                return 2*u**(2*M+3)*((M+2)-(M+1)*u*u)/(r*(1-u*u)**2)
            data = rational_bounds(order=M+10)
            for derivative, key in ((0,'C_error'),(1,'C_prime_error'),(2,'C_second_error')):
                expected = mp.mpf(data[key].numerator)/data[key].denominator
                actual = abs(mp.diff(envelope,mp.mpf(2),derivative))
                derivative_errors.append(float(abs(actual/expected-1)))
    return dict(maximum_series_error_over_bound=max(ratios),
                maximum_endpoint_relative_error=max(endpoint_errors),
                maximum_radial_derivative_relative_error=max(derivative_errors))


def polynomial_projection_controls():
    records=[]
    for mu in (-2.5,-.5,.5,1.5,2.5):
        J,Q=3,8
        M=Q-J-1
        angulars=[_AngularCoupling(mu=mu,nmax=J,angular_order=n) for n in (Q,40)]
        projected=[]
        for angular in angulars:
            z=angular.cosine
            f=sum(eval_gegenbauer(2*j,1.5,z)*.25**(2*j+1)/(2*j+1) for j in range(M+1))
            projected.append(np.einsum('p,pmn->mn',(1-z*z)*f,angular.kernel))
        records.append(dict(mu=mu,J=J,Q=Q,polynomial_degree=2*M,
            projection_difference=float(np.max(abs(projected[0]-projected[1]))),
            maximum_harmonic_gram_error=max(v.gram_error for v in angulars)))
    return records


def run():
    start=time.perf_counter()
    hashes={p:hashlib.sha256(project_path(p).read_bytes()).hexdigest() for p in INPUTS}
    data=rational_bounds()
    prior=json.loads(project_path('threshold_return_bound_results.json').read_text())
    assert prior['passed'] and prior['parameters']['outer_nmax']==9 and prior['parameters']['R']==2
    series=series_controls();projection=polynomial_projection_controls()
    records={k:outward(v) if isinstance(v,Fraction) else v for k,v in data.items()}
    tests=dict(
        exact_rational_export=all(Fraction.from_float(v['lower']) <= data[k] <= Fraction.from_float(v['upper'])
            for k,v in records.items() if isinstance(v,dict)),
        polynomial_degree=data['maximum_integrand_degree']<=data['GL_exact_degree'],
        independent_potential_series=series['maximum_series_error_over_bound']<=1+1e-13,
        endpoint_tail_equality=series['maximum_endpoint_relative_error']<1e-35,
        radial_derivatives=series['maximum_radial_derivative_relative_error']<1e-90,
        exact_polynomial_projection=max(r['projection_difference'] for r in projection)<1e-12,
        harmonic_gram_normalization=max(r['maximum_harmonic_gram_error'] for r in projection)<1e-12,
        frozen_inputs_preserved=hashes=={p:hashlib.sha256(project_path(p).read_bytes()).hexdigest() for p in INPUTS})
    result=dict(passed=all(tests.values()),tests=tests,seconds=time.perf_counter()-start,
        parameters=dict(q=1,d=1,R=2,J=9,Q=64,mu_scope='Every retained half-odd-integer sector'),
        bounds=records,series_controls=series,polynomial_projection_controls=projection,
        equations=dict(
            potential='u=d/(2r), z=cos(theta): Aphi/(1-z²)=q/2 integral_-u^u (1+2vz+v²)^(-3/2) dv.',
            series='f=q sum_j>=0 C_(2j)^(3/2)(z) u^(2j+1)/(2j+1).',
            remainder='epsilon=q u^(2M+3)[(M+2)-(M+1)u²]/(1-u²)².',
            degree='sin(theta)*Phi_plus_dagger sigma_y Phi_minus is a polynomial of degree<=2J+1. 2M+2J+1<=2Q-1.',
            projection='Exact and positive-weight GL projection each have operator norm<=sup|f-p|; error<=2epsilon/r for C.',
            derivative='Differentiate the uniformly convergent series; radial powers yield the saved positive moment formulas.'),
        sources=['https://dlmf.nist.gov/18.12.E4','https://dlmf.nist.gov/18.14.E4','https://dlmf.nist.gov/3.5.v'],
        scope=[
            'Exact-arithmetic Gauss-Legendre truncation versus exact finite-J angular Galerkin projection, uniformly for r>=2; not omitted J>9.',
            'Rational bound values and outward float endpoints are exact; sampled controls check the implementation, not the all-z proof.',
            'Floating nodes, weights, harmonics and radial/matrix arithmetic are not certified by this tiny truncation bound.',
            'No completed source, supporting geometry, stability or chronology conclusion; no exterior solve is repeated.'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashes)
    target=project_path('tail_quadrature_bound_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,seconds=result['seconds'],
                         bounds={k:v['upper'] if isinstance(v,dict) else v for k,v in records.items()},
                         series_controls=series),indent=2))
    assert result['passed'],tests
    return result


if __name__=='__main__':
    run()
