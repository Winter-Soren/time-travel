"""Leading gapped-mode bulk action, unequal-radius matching and source variation.

The global returning-field state and dynamical stability are not supplied here.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad

from dirac_flux_response import hyperbolic_heat, _log_nodes
from dirac_high_flux import angular_scaled, C0, C1_H, C1_S


D_GAP = -C1_S/8


def shape_integrals(q, eta, lower=.01, upper=48., order=16, width=.5, heat_order=128):
    """Action/q² and two shape-response/q² integrals at fixed a²r²=1."""
    s = math.exp(-4*eta)
    uv = np.zeros(3)
    for m, hm in {-1:.5, 0:1/12, 1:-1/120, 2:1/504, 3:-1/1440}.items():
        for n, bn in {3:-1/3, 4:-5/(21*q), 5:1/9-1/(6*q*q)}.items():
            term = hm*bn*s**m*q**(-m-1)*lower**(m+n)/(m+n)
            uv += term*np.array([1/8, (n-m)/4, -m/2])/math.pi**2
    core = np.zeros(3)
    for y, weight in _log_nodes(lower, upper, width, order):
        b, bp = angular_scaled(y,q)
        h, hp = hyperbolic_heat(s*y/q,heat_order)
        core += weight*np.array([h*b/(8*math.pi**2*q),
                                 y*(-s*hp*b+q*h*bp)/(4*math.pi**2*q*q),
                                 -y*s*hp*b/(2*math.pi**2*q*q)])
    u0 = s*upper/q
    xmax = 1/math.sqrt(u0)
    # The constant angular zero-mode tail is integrated to infinity using
    # x=1/sqrt(u). The kernel is regular at x=0; split its numerical branches.
    def tail_integrand(x):
        return 1/math.sqrt(math.pi) if x==0 else 2*hyperbolic_heat(1/(x*x),heat_order)[0]/x
    stop = min(1.,xmax)
    tail_heat = quad(tail_integrand,0,stop,epsabs=1e-12,epsrel=1e-12)[0]
    if xmax>1:
        tail_heat += quad(lambda z:2*hyperbolic_heat(math.exp(-2*z),heat_order)[0],
                          0,math.log(xmax),epsabs=1e-11,epsrel=1e-12)[0]
    tail_response = -hyperbolic_heat(u0,heat_order)[0]/(math.pi**2*q)
    tail = np.array([-tail_heat/(4*math.pi**2*q), tail_response/2, tail_response])
    total = uv+core+tail
    predicted = C1_H+C1_S/s
    coefficient = (float(total[2])-C0/s)*q
    return {'q':q, 'eta':eta, 'a_over_r':math.exp(2*eta),
            'action_over_q2':float(total[0]), 'response_over_q2':float(total[2]),
            'direct_response_over_q2':float(total[1]),
            'direct_vs_ibp':abs(float(total[1]-total[2])),
            'linear_coefficient':coefficient, 'predicted_linear_coefficient':predicted,
            'linear_coefficient_residual':coefficient-predicted,
            'hat_curvature':-2*s, 'first_gap_squared':q+1,
            'curvature_squared_logq_over_q':4*s*s*math.log(q)/q,
            'uv_core_ir': [uv.tolist(),core.tolist(),tail.tolist()]}


def source_identities():
    r,q,d,Q,k,S,A,e,w = sp.symbols('r q d Q k S A e w',positive=True)
    U=4*sp.pi*d*q/r**2
    J=-sp.diff(U,r)
    rho=U/(4*sp.pi*r**2); pr=-rho; pt=J/(8*sp.pi*r)
    E=Q+2*e*sp.log(S)
    prod_equation=S-E-k*(1-2*e/S)
    product_A=S*S/(S-2*e)
    product_Q=S-2*e*sp.log(S)-k*(1-2*e/S)
    c=k/S
    full_flare=c*((1-c)*w-(E-2*e)/S**2)/(1-3*c)
    reduced_flare=k/(S-3*k)*(w-(S-2*e)/S**2)
    # Use the radial constraint E=S*(1-k*w) only for this central identity.
    constrained_flare=full_flare.subs(Q,S*(1-k*w)-2*e*sp.log(S))
    residuals={
        'ward':-sp.diff(U,r)-J,
        'four_dimensional_trace':-rho+pr+2*pt,
        'radial_nec':rho+pr,
        'static_conservation':sp.diff(pr,r)+2*(pr-pt)/r,
        'combined_weyl_trace':-2*U+r*J,
        'product_tt':(1/S-E/S**2-k/(S*A)).subs(A,product_A).subs(Q,product_Q),
        'product_angular':(1/A-(E-2*e)/S**2-k/(S*A)).subs(A,product_A).subs(Q,product_Q),
        'classical_log_reference':prod_equation.subs({S:1,Q:1,k:0}),
        'zero_log_reference':prod_equation.subs({e:0,S:Q+k}),
        'central_flare_constraint':constrained_flare-reduced_flare,
    }
    return {name:str(sp.simplify(value)) for name,value in residuals.items()}


def reference_geometry():
    with mp.workdps(60):
        alpha=mp.mpf('1e-9'); q=mp.mpf(100000)
        k=alpha/3; epsilon=q*alpha/6
        delta=3*mp.zeta(3)*alpha/mp.pi**2
        cases=[]
        for name,shift in [('minimal_completion',mp.mpf(0)),('matched_gapped_potential',delta)]:
            Q=1+shift
            f=lambda S:S-Q-2*epsilon*mp.log(S)-k*(1-2*epsilon/S)
            S=mp.findroot(f,(mp.mpf(1),1+2*(delta+k)))
            A=S*S/(S-2*epsilon)
            E=Q+2*epsilon*mp.log(S)
            eqs=[1/S-E/S**2-k/(S*A),1/A-(E-2*epsilon)/S**2-k/(S*A)]
            cases.append({'name':name,'sphere_radius_squared':mp.nstr(S,45),
                          'ads_radius_squared':mp.nstr(A,45),
                          'sphere_radius':mp.nstr(mp.sqrt(S),45),
                          'ads_radius':mp.nstr(mp.sqrt(A),45),
                          'reference_state_omega2':mp.nstr(1/A,45),
                          'strict_first_order':{
                              'sphere_radius_squared':mp.nstr(1+shift+k,45),
                              'ads_radius_squared':mp.nstr(1+shift+k+2*epsilon,45),
                              'reference_state_omega2':mp.nstr(1-shift-k-2*epsilon,45)},
                          'maximum_equation_residual':float(max(abs(v) for v in eqs))})
        return {'alpha':float(alpha),'epsilon':float(epsilon),'k':float(k),
                'delta_energy_coefficient':float(delta),'cases':cases,
                'scope':'Resummed algebraic product reference of a truncated one-loop bulk action; includes increment80 magnetic logarithm once. Printed precision verifies equations, not physical accuracy. No global loop state or stability result.',
                'ordering_caveat':'Strict first order counts delta, k and epsilon as one-loop quantities at fixed q. Terms quadratic in them are not higher-loop predictions. Here epsilon^2 is comparable to k; finely tuned opening thresholds require omitted terms to be bounded.',
                'epsilon_squared':float(epsilon**2)}


def matched_source_checks():
    main=dict(lower=.02,upper=40.,order=12,width=.75,heat_order=80)
    cases=[]
    for q in (1000,100000):
        for eta in (-.3,-.15,0.,.15,.3):
            coarse=shape_integrals(q,eta,**main)
            fine=shape_integrals(q,eta)
            cases.append({'q':q,'eta':eta,'fine':fine,
                          'response_refinement':abs(fine['response_over_q2']-coarse['response_over_q2']),
                          'linear_coefficient_refinement':abs(fine['linear_coefficient']-coarse['linear_coefficient']),
                          'action_refinement':abs(fine['action_over_q2']-coarse['action_over_q2'])})
    differences=[]
    q=100000
    for eta in (-.3,.3):
        values=[]
        for step in (.002,.001):
            plus=shape_integrals(q,eta+step)['action_over_q2']
            minus=shape_integrals(q,eta-step)['action_over_q2']
            values.append((plus-minus)/(2*step))
        richardson=(4*values[1]-values[0])/3
        exact=next(c['fine']['response_over_q2'] for c in cases if c['q']==q and c['eta']==eta)
        differences.append({'eta':eta,'steps':[.002,.001],'central_derivatives':values,
                            'richardson':richardson,'absolute_error':abs(richardson-exact)})
    residuals=source_identities()
    geometry=reference_geometry()
    prior=json.loads(project_path('dirac_high_flux_results.json').read_text())
    old=prior['cases'][-1]['fine']['response_over_q2']
    same=next(c['fine']['response_over_q2'] for c in cases if c['q']==100000 and c['eta']==0.)
    tests={
        'action_variation_and_geometry_identities':all(v=='0' for v in residuals.values()),
        'equal_radius_saved_control':abs(old-same)<1e-13,
        'direct_variation_vs_ibp':max(c['fine']['direct_vs_ibp'] for c in cases)<1e-11,
        'quadrature_refinement':max(c['response_refinement'] for c in cases)<1e-11,
        'action_refinement':max(c['action_refinement'] for c in cases)<1e-11,
        'independent_action_derivative':max(d['absolute_error'] for d in differences)<2e-11,
        'unequal_radius_high_flux_matching':max(abs(c['fine']['linear_coefficient_residual']) for c in cases if c['q']==100000)<2e-6,
        'product_equations':max(c['maximum_equation_residual'] for c in geometry['cases'])<1e-50,
    }
    return {'passed':all(tests.values()),'tests':tests,'d_gap':D_GAP,
            'action_E':'4*pi*d_gap*q*integral sqrt(h)/R^2',
            'action_L':'-4*pi*d_gap*q*integral sqrt(-h)/R^2',
            'uplifted_stress':'(rho, p_radial, p_angular)=d_gap*q/R^4*(1,-1,1)',
            'shape_prediction':'R(eta)=C0*q^2*exp(4eta)+q*(C1_H+C1_S*exp(4eta))+o(q)',
            'cases':cases,'action_derivative_controls':differences,
            'symbolic_residuals':residuals,'product_reference':geometry,
            'regime':'Common four-dimensional renormalization; lowest modes kept once as q W_D[h/R^2;state], with zero added CFT cosmological term. Fixed smooth spherical bulk; unexcited gapped tower.',
            'limitations':['Only the leading O(q) gapped bulk action is matched; higher curvature terms are lower order for fixed smooth geometry.',
                           'The gap criterion alone does not bound all omitted terms for q-dependent geometries; curvature-squared and derivative bounds are needed.',
                           'Euler terms require their boundary completion. Physical interfaces and the returning-path state are not determined.',
                           'The added potential has zero radial null stress. A shifted product reference is not an opening, stability or formation proof.',
                           'Photon/graviton loops and remaining subleading one-loop terms are excluded.'],
            'sources':['https://arxiv.org/pdf/2603.15730','https://arxiv.org/pdf/1803.06948','https://arxiv.org/html/1807.04726v3'],
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result=matched_source_checks()
    project_path('dirac_matched_source_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':result['passed'],'tests':result['tests'],
                      'high_flux_cases':[{k:c['fine'][k] for k in ('eta','a_over_r','linear_coefficient','predicted_linear_coefficient','linear_coefficient_residual')} for c in result['cases'] if c['q']==100000],
                      'reference':result['product_reference']},indent=2))
    if not result['passed']:
        raise SystemExit('A source-matching control failed; diagnostic results retained.')
