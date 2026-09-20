"""Proper-clock null maps for prescribed timelike mouths in a flat exterior.

This checks transport needed by the driven JT test. It supplies no mouth
equation of motion, four-dimensional matching, or preparation mechanism.
"""
from project_paths import project_path
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


def proper_ray_jet(velocity, acceleration, jerk, receiver):
    """Return E', E'', E''' for receiver proper time -> emitter proper time.

    A stays at x=0; B follows x=d(t). Evaluate A's jet at the B emission
    event t (A receives at t+d(t)); B's jet at its reception event t.
    """
    if receiver not in ('A', 'B') or not math.isfinite(velocity) or abs(velocity)>=1:
        raise ValueError('Choose receiver A/B and a finite timelike velocity')
    if not all(math.isfinite(x) for x in (acceleration,jerk)):
        raise ValueError('Acceleration and jerk must be finite')
    v,a,j=velocity,acceleration,jerk
    q=math.sqrt(1-v*v)
    # Derivatives of B proper time with respect to the exterior lab time.
    q1=-v*a/q
    q2=-(a*a+v*j)/q-v*v*a*a/q**3
    f1,f2,f3,g1,g2,g3=(q,q1,q2,1+v,a,j) if receiver=='A' else (1-v,-a,-j,q,q1,q2)
    return (f1/g1,(f2*g1-f1*g2)/g1**3,
            (f3*g1*g1-f1*g3*g1-3*f2*g1*g2+3*f1*g2*g2)/g1**5)


def prescribed_exterior(amplitude, delay=.7, duration=4.):
    """Return exact ray maps for the compact displacement used in increment78.

    All callable inputs use lab time. The emission-map derivatives relate
    the two proper clocks, not the lab clocks. Motion is externally prescribed.
    """
    if not all(math.isfinite(x) for x in (amplitude,delay,duration)) or delay<=0 or duration<=0 or amplitude<=-1:
        raise ValueError('Require finite parameters, positive D/T, and separated mouths')
    polynomial=np.polynomial.Polynomial([0,0,0,0,256,-1024,1536,-1024,256])
    jets=[polynomial.deriv(k) for k in range(4)]
    extrema=[float(z.real) for z in jets[2].roots() if abs(z.imag)<1e-12 and 0<z.real<1]
    peak_speed=max(abs(amplitude*delay*jets[1](z)/duration) for z in [0.,1.,*extrema])
    if peak_speed>=1:raise ValueError('The prescribed trajectory is not timelike')
    def motion(t):
        if not 0<t<duration:return (0.,0.,0.,0.)
        return tuple(float(amplitude*delay*p(t/duration)/duration**k) for k,p in enumerate(jets))
    def proper_lapse(t):return math.sqrt(1-motion(t)[1]**2)
    # Integrating the clock difference avoids loss of a small proper-time lag.
    clock=solve_ivp(lambda t,y:(proper_lapse(t)-1,), (0.,duration), (0.,),
                    method='DOP853',rtol=2e-12,atol=2e-16,max_step=duration/80,dense_output=True)
    assert clock.success
    lag=-float(clock.y[0,-1])
    def proper_time(t):
        return t if t<=0 else t-lag if t>=duration else t+float(clock.sol(t)[0])
    minimum=delay*(1+min(amplitude,0.));maximum=delay*(1+max(amplitude,0.))
    def emission(t,side):
        if side not in (0,1):raise ValueError('Receiver side must be 0 (A) or 1 (B)')
        if side==0:
            guess=t-delay
            if motion(guess)[0]==0:emitter=guess
            else:
                margin=8*np.finfo(float).eps*(1+abs(t)+maximum)
                emitter=brentq(lambda e:e+delay+motion(e)[0]-t,
                               t-maximum-margin,t-minimum+margin,xtol=5e-14)
            _,v,a,j=motion(emitter)
            jet=proper_ray_jet(v,a,j,'A')
        else:
            h,v,a,j=motion(t);emitter=t-delay-h
            jet=proper_ray_jet(v,a,j,'B')
        schwarzian=jet[2]/jet[0]-1.5*(jet[1]/jet[0])**2
        assert emitter<t and jet[0]>0
        return float(emitter),float(jet[0]),float(schwarzian)
    return dict(motion=motion,proper_time=proper_time,proper_lapse=proper_lapse,
                emission=emission,final_proper_clock_lag=lag,maximum_speed=peak_speed,
                delay=delay,duration=duration,amplitude=amplitude)


def prescribed_exterior_checks():
    """Check emission roots, independent clock quadrature, and the static limit."""
    rows=[]
    for amplitude in (0.,.001,.05):
        exterior=prescribed_exterior(amplitude)
        roots=[];clocks=[]
        for t in np.linspace(-1,8,49):
            e,_,_=exterior['emission'](float(t),0)
            roots.append(abs(e+.7+exterior['motion'](e)[0]-t))
            e,_,_=exterior['emission'](float(t),1)
            roots.append(abs(t-e-.7-exterior['motion'](float(t))[0]))
            reference=t+quad(lambda u:exterior['proper_lapse'](u)-1,0,max(0.,min(4.,t)),
                             epsabs=3e-15,epsrel=2e-12)[0]
            clocks.append(abs(exterior['proper_time'](float(t))-reference))
            if amplitude==0:
                assert exterior['emission'](float(t),0)==(t-.7,1.,0.)
                assert exterior['emission'](float(t),1)==(t-.7,1.,0.)
                assert exterior['proper_time'](float(t))==t
        assert max(roots)<3e-13 and max(clocks)<3e-13
        rows.append(dict(amplitude=amplitude,maximum_null_ray_residual=max(roots),
            maximum_clock_quadrature_difference=max(clocks),maximum_speed=exterior['maximum_speed'],
            final_proper_clock_lag=exterior['final_proper_clock_lag']))
    return dict(all_declared_checks_pass=True,cases=rows)


def moving_energy_identity_checks():
    """Algebraic controls for the improved-charge and moving-field balance."""
    import sympy as sp
    v,a,ad,B,qi,qo=sp.symbols('v a a_proper_prime B q_in q_out', real=True)
    su,sv=-ad-a*a/2,ad-a*a/2
    pin,pout=qi+su/3,qo+sv/3
    improved=B-2*a/3
    # Boundary A has already cancelled its adjacent fixed-end field flux.
    field_and_A=-(1+v)*pin+(1-v)*pout
    boosted_derivative=v*a*improved+qi-qo-2*ad/3
    drive=v*(a*B-qi-qo-a*a/3)
    balance=sp.simplify(field_and_A+boosted_derivative-drive)
    primitive=sp.simplify((1-v)*(qo+sv/3)-(1-v)*(qo+a*a/6)-(1-v)*(ad-a*a)/3)
    passive=sp.simplify(drive.subs({B:2*a/3,qi:-su/3,qo:-sv/3}))
    h,hd,h2,hd2,om=sp.symbols('h h_delayed h_second h_delayed_second Omega_squared')
    old_effective=(h2+hd2+om*(h+hd))/3
    new_total=old_effective-om*hd/3+(h2-hd2)/3-2*h2/3
    first_order=sp.simplify(new_total-om*h/3)
    residuals=[balance,primitive,passive,first_order]
    assert residuals==[0,0,0,0]
    return dict(all_declared_checks_pass=True,
        residuals=dict(zip(('boosted_balance','integrated_by_parts_field_energy',
                           'passive_vacuum_drive','first_order_vacuum_work'),map(str,residuals))),
        reference_choice='B_i=-b0*(K_i-1/2); no material rest-energy constant supplied',
        scope='Exact identities of the specified reduced model. The flat bulk energy is physical within that model; the improved boundary charge is not an identified material rest mass or complete apparatus energy.')


def moving_exterior_map_checks():
    """Independent exact-map controls for the first-order motion source."""
    D,T=.7,4.
    omega2=(2*math.pi/(math.pi+D))**2
    q0=-omega2/6
    bump=np.polynomial.Polynomial([0,0,0,0,256,-1024,1536,-1024,256])
    derivatives=[bump.deriv(k) for k in range(4)]
    def displacement(t,amplitude):
        if not 0<t<T:return np.zeros(4)
        return amplitude*D*np.array([p(t/T)/T**k for k,p in enumerate(derivatives)])
    samples=np.linspace(.02,T-.02,99)
    cases=[];rapidity_error=0.;vacuum_error=0.;minimum_slope=1.
    for amplitude in (.008,.004,.002,.001):
        map_errors=[];source_errors=[];peak_velocity=0.;peak_source=0.
        for t in samples:
            h,v,a,j=displacement(t,amplitude)
            sigma=t+quad(lambda u:math.sqrt(1-displacement(u,amplitude)[1]**2)-1,
                         0,t,epsabs=2e-15,epsrel=2e-12)[0]
            peak_velocity=max(peak_velocity,abs(v))
            eta_dot=a/(1-v*v)
            eta_ddot=j/(1-v*v)+2*v*a*a/(1-v*v)**2
            for receiver in ('A','B'):
                p1,p2,p3=proper_ray_jet(v,a,j,receiver)
                schwarzian=p3/p1-1.5*(p2/p1)**2
                if receiver=='A':
                    receiving=t+D+h;emitting=sigma;argument=receiving-D
                    rapidity_schwarzian=(-eta_ddot/(1+v)**2+eta_dot*a/(1+v)**3
                                        -.5*eta_dot**2/(1+v)**2)
                else:
                    receiving=sigma;emitting=t-D-h;argument=receiving
                    rapidity_schwarzian=-(eta_ddot+v*a*eta_dot/(1-v*v)+.5*eta_dot**2)/(1-v*v)
                hl,vl,_,jl=displacement(argument,amplitude)
                map_errors.append(abs(emitting-(receiving-D-hl)))
                exact_source=(p1*p1-1)*q0-schwarzian/3
                source_errors.append(abs(exact_source-(jl+omega2*vl)/3))
                peak_source=max(peak_source,abs(exact_source))
                rapidity_error=max(rapidity_error,abs(schwarzian-rapidity_schwarzian),
                                   abs(p1-math.sqrt((1-v)/(1+v))))
                minimum_slope=min(minimum_slope,p1)
                # Passive rechart of Minkowski vacuum: normal ordering alone
                # is not covariant stress. Compute the geometric term separately.
                sigma_prime=p2/(2*p1)
                sigma_second=(p3/p1-(p2/p1)**2)/2
                geometric=2*(sigma_second-sigma_prime*sigma_prime)/3
                vacuum_error=max(vacuum_error,abs(-schwarzian/3+geometric))
        cases.append(dict(amplitude=amplitude,maximum_sampled_speed=peak_velocity,
            maximum_map_linearization_error=max(map_errors),maximum_source_linearization_error=max(source_errors),
            maximum_exact_transformed_source_change=peak_source,
            final_proper_clock_lag=-quad(lambda u:math.sqrt(1-displacement(u,amplitude)[1]**2)-1,
                                       0,T,epsabs=2e-15,epsrel=2e-12)[0]))
    for coarse,fine in zip(cases,cases[1:]):
        for key in ('maximum_map_linearization_error','maximum_source_linearization_error','final_proper_clock_lag'):
            ratio=coarse[key]/fine[key]
            assert 3.7<ratio<4.3,(key,ratio)
    for v in (-.3,0.,.3):
        expected=math.sqrt((1-v)/(1+v))
        for receiver in ('A','B'):
            jet=proper_ray_jet(v,0.,0.,receiver)
            assert abs(jet[0]-expected)<1e-15 and jet[1:]==(0.,0.)
    assert rapidity_error<2e-14 and vacuum_error<2e-15 and minimum_slope>0
    return dict(model='Exact flat-exterior null transport for prescribed A stationary and B timelike; no gravitational closure',
        parameters=dict(initial_separation=D,motion_duration=T,normalized_vacuum_stress=q0),
        equations=dict(proper_clock='sigma(t)=integral_0^t sqrt(1-d_dot²)dt',
            maps='E_A(t+d(t))=sigma(t); E_B(sigma(t))=t-d(t). Both E_prime=exp(-atanh(d_dot)), at different physical events.',
            normal_ordered_transport='Q_received=E_prime²*Q_emitted(E)-Schwarzian(E)/3; this is not by itself covariant observer stress.',
            linear_source='d=D+h: i_A(s)=[h_tripleprime+Omega²*h_prime](s-D)/3; i_B(s)=[h_tripleprime+Omega²*h_prime](s)/3.',
            vacuum_control='For a passive flat null rechart, geometric 2*(sigma_second-sigma_prime²)/3 cancels -Schwarzian(E)/3.'),
        cases=cases,verification=dict(all_declared_checks_pass=True,
            rapidity_and_parametric_jet_maximum_difference=rapidity_error,
            passive_Minkowski_vacuum_covariant_stress_maximum_error=vacuum_error,
            inertial_affine_map_checks_pass=True,second_order_linearization_scaling_pass=True,
            minimum_sampled_map_derivative=minimum_slope),
        scope='Exact transport on an assumed flat one-dimensional exterior with specified transparent proper-clock sewing. The driven feedback uses only its first order in amplitude. Source-chart terms alone are not radiation; net production and motion cost require higher-order/global dynamics. Not a physical MMP magnetic loop or backwards-time experiment.',
        sources=['https://arxiv.org/pdf/math-ph/0412028: normal-ordered stress transformation, Eqs2.7-2.9',
                 'https://arxiv.org/html/1807.04726v3: geometric anomaly must be retained, AppendixF'])


if __name__=='__main__':
    result=moving_exterior_map_checks()
    project_path('moving_exterior_map_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['verification'],indent=2))
