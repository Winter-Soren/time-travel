"""Corrected spherical bulk source with a prescribed weak exterior return path.

This closes only the spherical average. It does not solve the exterior
Einstein-Maxwell problem, prepare the NS state, or establish stability.
"""
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from polyakov_response import _magnetic_closure_paths


def central_radius(alpha, w, epsilon, delta):
    with mp.workdps(50):
        a, frequency, e, d = map(mp.mpf, map(str, (alpha, w, epsilon, delta)))
        u = mp.findroot(lambda u: (1+u)**2*(1-a*frequency/3)-1-d-4*e*mp.log1p(u),
                        (mp.mpf(0), a))
        return float(u)


def collar(w, alpha=1e-9, epsilon=None, delta=None, refinement=1, outer=100., *, return_profile=False):
    if epsilon is None: epsilon=100000*alpha/6
    if delta is None: delta=float(3*mp.zeta(3)/mp.pi**2)*alpha
    if not (alpha>0 and 0<=epsilon<.1 and w>0 and outer>1):
        raise ValueError('Require a positive weak-coupling branch and outer radius > 1')
    k=alpha/3; offset=central_radius(alpha,w,epsilon,delta); r0=1+offset

    def fields(state):
        y,v,n,a,x,z=state
        u=offset+alpha*y; r=1+u; b=alpha*v
        E=1+delta+4*epsilon*np.log1p(u); W=w*np.exp(-2*n)
        D=r*r-k; B=r*a-b
        if np.any(r*r<=3*k): raise ValueError('Left the regular source branch R²>3k')
        ap=D/(r*r-3*k)*((E-2*epsilon+k*B*B)/r**4-a*a-2*a*b/r
                                  -k*(r*r+k)/(r*r*D)*W)
        vp=a*v+r/(3*D)*(W-ap)
        # Keep the O(alpha) central radius offset separate from the unit radius.
        C=-u*(2+u)+delta+4*epsilon*np.log1p(u)+r*r*b*b+2*r**3*b*a+k*(B*B+W*r*r)
        return r,b,a,ap,alpha*vp,C/r**2

    def rhs(ell,state):
        r,b,a,ap,bp,C=fields(state)
        return [state[1],bp/alpha,a,ap,math.exp(-state[2]),1/r]

    def end(ell,state): return 1+offset+alpha*state[0]-outer
    end.terminal=True; end.direction=1
    sol=solve_ivp(rhs,(0.,outer+100.),np.zeros(6),method='DOP853',
                  rtol=2e-11/refinement,atol=2e-12/refinement,
                  max_step=.35/refinement,events=end,dense_output=True)
    if not sol.success or len(sol.t_events[0])!=1:
        raise RuntimeError('Corrected collar did not reach the requested outer radius: '+sol.message)
    anchor=sol.y[:,-1]; r,b,a,ap,bp,C=fields(anchor)
    clock=math.exp(anchor[2])/b
    # A pure local magnetic exterior provides only an isotropic-coordinate anchor.
    # The omitted returning-field tail is not silently continued to infinity.
    mass=r*(1-b*b)/2+(1+delta+4*epsilon*(math.log(r)+1))/(2*r)
    def isotropic_integrand(t):
        if t==0: return mass
        f=1-2*mass*t+(1+delta+4*epsilon*(1-math.log(t)))*t*t
        root=math.sqrt(f)
        return (1-f)/(root*(1+root)*t)
    integral=quad(isotropic_integrand,0,1/r,epsabs=2e-13,epsrel=2e-12)[0]
    rho_anchor=r*math.exp(-integral)
    cuts={}
    for radius in (10.,30.,100.):
        if radius>outer: continue
        ell=sol.t[-1] if radius==outer else brentq(lambda t:fields(sol.sol(t))[0]-radius,0,sol.t[-1],xtol=1e-12)
        state=sol.sol(ell); R,b,a,ap,bp,C=fields(state)
        cuts[str(radius)]={'optical_x':float(state[4]),'lapse':math.exp(state[2]),
            'normal_radius_derivative':b,'local_osculating_clock':math.exp(state[2])/b,
            'isotropic_radius_from_collar':rho_anchor*math.exp(state[5]-anchor[5])}
    center=fields(np.zeros(6))
    sample=np.linspace(0,sol.t[-1],241); samples=sol.sol(sample)
    R,b,a,ap,bp,C=fields(samples)
    W=w*np.exp(-2*samples[2]); B=R*a-b; Bd=b*a+R*ap-bp
    E=1+delta+4*epsilon*np.log1p(offset+alpha*samples[0])
    rho=E/R**4+k/R**4*(2*R*Bd+B*B-W*R*R)
    pr=-E/R**4-k/R**4*(B*B+W*R*R)
    pt=(E-2*epsilon)/R**4+k/R**4*(R*Bd+B*B)
    einstein=np.array([(1-b*b-2*R*bp)/R**2,(-1+b*b+2*R*b*a)/R**2,
                      ap+a*a+(a*b+bp)/R])
    curvature=-2*(R*Bd+B*B)
    record={'clock':clock,'cuts':cuts,'central_radius_offset':offset,
            'central_flare_over_alpha':center[4]/(r0*alpha),
            'central_lapse_curvature':center[3], 'mass_anchor':mass,
            'max_constraint':float(np.max(abs(C))),
            'max_constraint_near_throat_over_alpha':float(np.max(abs(C[R<1.01]))/alpha),
            'max_einstein_residual':float(np.max(abs(einstein-np.array([rho,pr,pt])))),
            'max_trace_residual':float(np.max(abs(-rho+pr+2*pt+4*epsilon/R**4))),
            'max_hat_curvature':float(np.max(abs(curvature))),
            'curvature_squared_logq_over_q':float(np.max(curvature**2)*math.log(100000)/100000),
            'proper_half_length':float(sol.t[-1]),'integration_evaluations':sol.nfev,
            'parameters':dict(alpha=alpha,epsilon=epsilon,delta=delta,w=w)}
    if return_profile:
        return record, {'solution':sol,'central_radius_offset':offset,'alpha':alpha}
    return record


def closure_checks():
    prior=json.loads(Path(__file__).with_name('polyakov_closure_results.json').read_text())
    distance=prior['parameters']['separation_over_re']; alpha=1e-9
    old=prior['quadrature_refinement'][-1]['roots'][-1]
    records=[]
    for order,refinement in ((64,1),(128,2)):
        nu,weights,f,pieces,path_checks=_magnetic_closure_paths(order,distance)
        background=lru_cache(maxsize=None)(lambda w:collar(float(w),refinement=refinement))
        def state_mean(w,radius=100.,path_scale=1.):
            bg=background(w); cut=bg['cuts'][str(radius)]; piece=pieces[str(radius)]
            external=path_scale*distance*f-2*piece['isotropic']+bg['mass_anchor']*piece['shapiro']
            period=2*cut['optical_x']+external/bg['clock']
            if np.any(external<=0) or np.any(period<=0): raise ValueError('Nonpositive modeled return path')
            return float(weights@((2*math.pi/period)**2)),period,external
        roots=[]
        for radius in (10.,30.,100.):
            residual=lambda w:w-state_mean(w,radius)[0]
            bracket=[float(residual(1.1)),float(residual(3.5))]
            assert bracket[0]<0<bracket[1]
            w=brentq(residual,1.1,3.5,xtol=2e-12,rtol=2e-12)
            bg=background(w); mean,period,external=state_mean(w,radius)
            cut=bg['cuts'][str(radius)]; piece=pieces[str(radius)]
            rho_error=abs(cut['isotropic_radius_from_collar']-piece['isotropic'])
            radius_bound=4*rho_error*(1+4*bg['mass_anchor']/(piece['isotropic']-rho_error))/bg['clock']
            axis=distance-2*piece['isotropic']+4*bg['mass_anchor']*math.log((distance-piece['isotropic'])/piece['isotropic'])
            roots.append({'cutoff':radius,'omega_squared':w,'root_residual':w-mean,
                'bracket_residuals':bracket,'background':bg,
                'throat_time':2*cut['optical_x'],'shortest_exterior_time':axis/bg['clock'],
                'throat_to_exterior_ratio':2*cut['optical_x']/(axis/bg['clock']),
                'isotropic_radius_difference':rho_error,'radius_map_period_bound':radius_bound,
                'channel_variance':float(weights@(((2*math.pi/period)**2-mean)**2)),
                'periods':period.tolist(),'nu':nu.tolist(),'weights':weights.tolist()})
        root=roots[-1]; w=root['omega_squared']; bg=root['background']
        derivatives=[]
        for h in (1e-4,5e-5):
            derivatives.append((state_mean(w+h)[0]-state_mean(w-h)[0])/(2*h))
        conditioning=1/(1-derivatives[-1])
        sensitivity=[]
        for scale in (.999,1.001):
            value=brentq(lambda trial:trial-state_mean(trial,path_scale=scale)[0],1.1,3.5,
                         xtol=2e-12,rtol=2e-12)
            sensitivity.append({'path_scale':scale,'omega_squared':value,
                                'central_flare_over_alpha':background(value)['central_flare_over_alpha']})
        records.append({'order':order,'refinement':refinement,'roots':roots,
            'path_checks':path_checks,'fixed_point_derivatives':derivatives,
            'conditioning_factor':conditioning,'path_scale_sensitivity':sensitivity,
            'cutoff_root_span':float(np.ptp([r['omega_squared'] for r in roots]))})
    coarse,fine=records
    root=fine['roots'][-1]; bg=root['background']; w=root['omega_squared']
    refinement=max(abs(a['omega_squared']-b['omega_squared']) for a,b in zip(coarse['roots'],fine['roots']))
    all_bg=[r['background'] for record in records for r in record['roots']]
    from dirac_matched_source import reference_geometry
    omega_ref=float(reference_geometry()['cases'][-1]['reference_state_omega2'])
    central_expected=(w-(1/(1+bg['central_radius_offset'])**2
        -2*bg['parameters']['epsilon']/(1+bg['central_radius_offset'])**4))/(3*((1+bg['central_radius_offset'])**2-alpha))
    tests={
        'closed_return_state':max(abs(r['root_residual']) for record in records for r in record['roots'])<1e-10,
        'quadrature_and_ode_refinement':refinement<2e-9,
        'matching_cutoff_refinement':fine['cutoff_root_span']<2e-9,
        'einstein_and_trace':max(b['max_einstein_residual'] for b in all_bg)<1e-10 and max(b['max_trace_residual'] for b in all_bg)<2e-14,
        'near_throat_constraint':max(b['max_constraint_near_throat_over_alpha'] for b in all_bg)<1e-3,
        'central_flare_identity':abs(central_expected-bg['central_flare_over_alpha'])<1e-13,
        'derivative_resolution':abs(fine['fixed_point_derivatives'][0]-fine['fixed_point_derivatives'][1])<2e-7,
        'positive_opening_and_lengths':all(b['central_flare_over_alpha']>0 for b in all_bg) and all(min(r['periods'])>0 for record in records for r in record['roots']),
    }
    tests={name:bool(value) for name,value in tests.items()}
    return {'passed':all(tests.values()),'tests':tests,
        'parameters':{'alpha':alpha,'q':100000,'distance_over_re':distance,
                      'distance_convention':prior['parameters']['separation_convention']},
        'records':records,'root_refinement':refinement,
        'old_bare_polyakov_root':old['omega_squared'],'corrected_root':w,
        'change_from_old_source':w-old['omega_squared'],'reference_omega_squared':omega_ref,
        'opening_margin':w-omega_ref,
        'equations':{'state':'w=<[2pi/L_nu]^2>, uniform full flux measure dnu/2; NS combined spin/Wilson phase.',
            'period':'L_nu=2*x_b+[d*f(nu)-2*(Rb-1)+M_anchor*Shapiro_nu]/A; A=N100/b100.',
            'source':'q W_D[h/R²] plus magnetic log and d*q/R4 potential; all varied including radius.',
            'constraint':'C=-1+b²+2Rab+E/R²+k[(a-b/R)²+w/N²]; C_l=0 on the ODE.',
            'isotropic_anchor':'Pure local magnetic exterior f=1-2M/R+[Q+4eps(logR+1)]/R², used only to estimate the cap-coordinate discrepancy.'},
        'scope':'Self-consistent spherical average within the declared truncated bulk action and prescribed point-magnetic/weak-exterior paths. Not a complete exterior match or dynamical stability result.',
        'limitations':['Corrected nonlinear Maxwell field paths and nonspherical metric/current remain unsolved; the exterior path is prescribed.',
            'The NS combined phase is assumed, not dynamically prepared. Matching cuts are transparent, not reflecting walls.',
            'One-loop resummation is not higher-loop accuracy; lower-order curvature/boundary terms and higher loops remain omitted.',
            'The 0.1 percent path rescaling is a sensitivity probe, not a physical uncertainty bound or solved field deformation.',
            'Fixed-point iteration convergence is not physical stability. No creation mechanism, material recipe or past-time passage is established.'],
        'sources':['https://arxiv.org/html/1807.04726v3#S5.SS5','https://arxiv.org/abs/2603.15730'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result=closure_checks()
    Path(__file__).with_name('dirac_return_closure_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','tests','corrected_root','change_from_old_source','opening_margin','root_refinement')},indent=2))
    print(json.dumps({k:result['records'][-1][k] for k in ('fixed_point_derivatives','conditioning_factor','path_scale_sensitivity','cutoff_root_span')},indent=2))
    if not result['passed']: raise SystemExit('A closure control failed; diagnostics retained.')
