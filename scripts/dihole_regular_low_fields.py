"""Pole-regular low-order fields of the fixed balanced dihole calibration.

``build_evaluator(context)(R, x)`` returns zeta, zeta_R, zeta_x and
G, G_R, G_x, where x is the local polar cosine and R the local radius.
All R derivatives hold x fixed. The charge-independent potential is the
same a as in dihole_channel_bounds.field_function and PointwiseDirac.

The exact balanced identities used here are
  sqrt(N*P) = 1 + (1-x*x)*zeta,
  a = (1-x)/2 + (1-x*x)*K,       G = K_R.
No evaluation divides by 1-x*x or subtracts nearly equal polar values.
The expressions apply for 0<R<2*k, -1<=x<=1, including both poles.
They support mpmath and its interval context; a point self-check is not a
whole-buffer enclosure. Physical conventions: Emparan hep-th/9906160,
Eqs. 8--16, and the project's fixed charge-included gauge calibration.
"""
from functools import lru_cache
from time import perf_counter

import sympy as sp
from mpmath import mp

from dihole_channel_bounds import calibration, field_function


@lru_cache(maxsize=1)
def _expressions():
    """Remove the polar factor by exact polynomial division, then differentiate."""
    u, v, z, w = sp.symbols('u v z w', real=True)
    # z=M/(k+a_g) rationalizes both a_g**2+M**2=k**2 and balance.
    m = 2*z/(1+z*z)
    ag = (1-z*z)/(1+z*z)
    B = z*(1+z*z)/(1+z)**2
    Q = 2*z*(1+z)**2/(1+z*z)**2
    r, cosine = m+(u+v)/2, (v-u)/2
    delta = ((u+v)**2-4)/4
    bigW = (r*r-ag*ag)**2+delta*ag*ag*(1-cosine*cosine)
    V = (delta+ag*ag*(1-cosine*cosine)
         + 2*B*m*r*ag*(1-cosine*cosine)
         + B*B*bigW*(1-cosine*cosine)/4)
    potential_numerator = 2*m*r*ag+B*bigW/2
    numerator, denominator = sp.together(
        u*(v-u+2)*potential_numerator-Q*(u+v+2)*V
    ).as_numer_denom()
    quotient, remainder = sp.div(sp.Poly(numerator, v), sp.Poly(u+v-2, v))
    if not remainder.is_zero:
        raise ArithmeticError('Balanced potential did not have its polar factor')
    expected_denominator = 128*(1+z)**2*(1+z*z)**4
    if sp.cancel(denominator/expected_denominator) != 1:
        raise ArithmeticError('Unexpected exact potential normalization')
    # w=(v-2)/u has the stable evaluation (u+4*x)/(v+2). Shifting before
    # polynomial evaluation also removes the zero at u=0 symbolically.
    polynomial = sp.cancel(
        quotient.as_expr().subs(v, 2+u*w)/(4*u*z*(1+z*z)**2)
    )
    polynomial = sp.horner(sp.Poly(polynomial, u, w, z).as_expr(), u, w, z)
    v = 2+u*w
    r = m+1+u*(1+w)/2
    x = w+u*(w*w-1)/4
    sine_squared = u*(1-w)*(1-u*(1-w)/4)
    ell = (2*B*m*ag+B*B*(r+m+1)*(r*r+(m+1)**2-2*ag*ag)/4)/(r-m+1)
    ell += B*B*ag*ag*sine_squared/4
    zeta = u*ell/v
    h = 1+(1-x*x)*zeta
    K = u*polynomial/(32*(1+z)**4*(4+u*(1+w))*v*h*(4+u*(w-1)))
    values = (zeta, sp.diff(zeta, u), sp.diff(zeta, w),
              sp.diff(K, u), sp.diff(K, w), sp.diff(K, u, 2),
              sp.diff(K, u, w), sp.diff(K, w, 2))
    return (u, w, z), values


def build_evaluator(context=mp):
    """Compile once. Recompute the linked calibration if context.dps changes."""
    symbols, expressions = _expressions()
    function = sp.lambdify(symbols, expressions, modules=[{'sqrt': context.sqrt}], cse=True)
    parameter_cache = {}

    def evaluate(R, x):
        R, x = context.mpf(R), context.mpf(x)
        precision = context.dps
        if precision not in parameter_cache:
            k, M, ag, _, _ = calibration(context)
            parameter_cache[precision] = k, M/(k+ag)
        k, z = parameter_cache[precision]
        if not (R > 0 and R < 2*k and x >= -1 and x <= 1):
            raise ValueError('Expected 0<R<2*k and -1<=x<=1, including interval endpoints')
        u = R/k
        v = context.sqrt(u*u+4+4*u*x)
        w = (u+4*x)/(v+2)
        Z, Zu, Zw, Ku, Kw, Kuu, Kuw, Kww = function(u, w, z)
        # x=w+u*(w*w-1)/4: these derivatives hold the original x fixed.
        wu, wx = (1-w*w)/(2*v), 2/v
        wuu = -(2*w*wu+u*wu*wu)/v
        wux = -2*(w+u*wu)/(v*v)
        return {
            'zeta': Z,
            'zeta_R': (Zu+Zw*wu)/k,
            'zeta_x': Zw*wx,
            'G': (Ku+Kw*wu)/k,
            'G_R': (Kuu+2*Kuw*wu+Kww*wu*wu+Kw*wuu)/(k*k),
            'G_x': (Kuw*wx+Kww*wu*wx+Kw*wux)/k,
        }

    return evaluate


def self_check():
    """Compare the untouched potential jets, including independent polar limits."""
    started = perf_counter()
    with mp.workdps(65):
        evaluate = build_evaluator(mp)
        original = field_function(mp)
        parameters = calibration(mp)
        maximum = mp.mpf(0)
        cases = 0
        for Rtext, xtext in (('0.025', '-0.7'), ('0.047', '0.2'), ('0.070', '0.8'),
                             ('0.025', '-1'), ('0.025', '1'),
                             ('0.070', '-1'), ('0.070', '1')):
            R, x = mp.mpf(Rtext), mp.mpf(xtext)
            got = evaluate(R, x)
            raw = original(R, x, *parameters)
            h = mp.sqrt(raw[0])
            hR, hx = raw[3]/(2*h), raw[6]/(2*h)
            if abs(x) < 1:
                s2 = 1-x*x
                Z = (h-1)/s2
                expected = (Z, hR/s2, hx/s2+2*x*Z/s2,
                            raw[5]/s2, raw[11]/s2, raw[14]/s2+2*x*raw[5]/s2**2)
            else:
                Z = -x*hx/2
                hRx = raw[12]/(2*h)-raw[3]*raw[6]/(4*h**3)
                hxx = raw[15]/(2*h)-raw[6]**2/(4*h**3)
                G = -x*raw[14]/2
                aRRx = mp.diff(lambda xx: original(R, xx, *parameters)[11], x)
                aRxx = mp.diff(lambda xx: original(R, xx, *parameters)[14], x)
                expected = (Z, -x*hRx/2, -x*(hxx+2*Z)/4,
                            G, -x*aRRx/2, -x*(aRxx+2*G)/4)
            maximum = max(maximum, *(abs(value-reference)/max(1, abs(reference))
                                    for value, reference in zip(got.values(), expected)))
            cases += 1
        from mpmath.ctx_iv import MPIntervalContext
        iv = MPIntervalContext()
        iv.dps = 60
        interval_evaluate = build_evaluator(iv)
        contained = True
        for pole in (-1, 1):
            interval_x = iv.mpf([-1, '-0.99'] if pole == -1 else ['0.99', 1])
            interval_values = interval_evaluate(iv.mpf(['0.049', '0.051']), interval_x)
            point_values = evaluate(mp.mpf('0.05'), mp.mpf(pole))
            for name, value in point_values.items():
                interval = interval_values[name]
                contained &= bool(interval.a <= iv.mpf(str(value)) <= interval.b)
        assert maximum < mp.mpf('1e-55'), maximum
        assert contained, 'Pole-containing intervals missed a checked point'
        return {'passed': True, 'original_jet_cases': cases,
                'maximum_scaled_error': float(maximum),
                'both_polar_interval_containments': bool(contained),
                'seconds': perf_counter()-started}


if __name__ == '__main__':
    import json
    print(json.dumps(self_check(), indent=2))
