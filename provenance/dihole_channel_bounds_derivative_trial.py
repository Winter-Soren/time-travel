"""Uniform regular-band moments and conservative collar channel residual bounds.

Directed geometry enclosures and analytic moments replace selected-point
asymptotics. No exterior/source trace, spectral gap or atom transit is bounded.
"""
import json
import math
from time import perf_counter

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext

from balanced_dihole import symbolic_fields
from dihole_local_source import local_first
from dihole_slice_channel import local_fields, sha
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME='dihole_channel_bounds'
FIELDS=('NP','W','a')
JETS=('value','R','x','RR','Rx','xx')


def field_function(context):
    """Stable factored expressions, differentiated before interval evaluation."""
    R,x,k,M,ag,B,Q=sp.symbols('R x k M ag B Q',real=True)
    d=sp.sqrt(R**2+4*k*k+4*k*R*x)
    r,c=M+(R+d)/2,(d-R)/(2*k)
    S=(d-R-2*k*x)/(2*k*k); delta=(R+d+2*k*x)/2
    w=(r*r-ag*ag)**2+R*R*delta*ag*ag*S
    h=2*M*r*ag+B*w/2
    z=ag*ag+2*B*M*r*ag+B*B*w/4
    V=delta+S*z; sigma=r*r-ag*ag*c*c
    values=((V/d)**2,(sigma/(R*d))**2,S*h/(2*Q*V))
    expressions=[]
    for op in (lambda f:f,lambda f:sp.diff(f,R),lambda f:sp.diff(f,x),
               lambda f:sp.diff(f,R,2),lambda f:sp.diff(f,R,x),lambda f:sp.diff(f,x,2)):
        expressions.extend(op(f) for f in values)
    return sp.lambdify((R,x,k,M,ag,B,Q),expressions,modules=[{'sqrt':context.sqrt}],cse=True)


def calibration(iv):
    k,Q=iv.mpf(1)/2,iv.mpf(1)/5
    M=(iv.sqrt(k*k+4*Q*k)-k)/2
    ag=iv.sqrt(k*k-M*M);B=M/((k+ag)*(M+k))
    return k,M,ag,B,Q


def lower(value): return float(np.nextafter(float(value.a),-np.inf))
def upper(value): return float(np.nextafter(float(value.b),np.inf))
def magnitude(value): return max(abs(lower(value)),abs(upper(value)))


def evaluate(function,parameters,R,theta,iv):
    raw=function(R,iv.cos(theta),*parameters)
    return {field:{jet:raw[3*j+i] for j,jet in enumerate(JETS)} for i,field in enumerate(FIELDS)}


def cover(iv,function,parameters,theta_lo,theta_hi,radial=16,angular=32):
    bounds={field:{jet:[math.inf,-math.inf] for jet in JETS} for field in FIELDS}
    extra={name:[math.inf,-math.inf] for name in ('a_theta','L_theta_theta','L_theta_theta_theta','L_s_theta_theta','log_W_R','log_W_theta')}
    nu=iv.mpf([iv.mpf(17)/20,iv.mpf(19)/20])
    for i in range(radial):
        R=iv.mpf([iv.mpf(1)/40+iv.mpf(i)/(40*radial),iv.mpf(1)/40+iv.mpf(i+1)/(40*radial)])
        for j in range(angular):
            theta=iv.mpf([theta_lo+(theta_hi-theta_lo)*j/angular,theta_lo+(theta_hi-theta_lo)*(j+1)/angular])
            data=evaluate(function,parameters,R,theta,iv)
            for field in FIELDS:
                for jet in JETS:
                    b=bounds[field][jet];v=data[field][jet]
                    b[0]=min(b[0],lower(v));b[1]=max(b[1],upper(v))
            sine,cosine=iv.sin(theta),iv.cos(theta)
            f,a=data['NP'],data['a'];W=data['W']
            at=-sine*a['x'];ft=-sine*f['x']
            derived={'a_theta':at,'log_W_R':W['R']/W['value'],'log_W_theta':-sine*W['x']/W['value']}
            if sine.a>0:
                att=-cosine*a['x']+sine*sine*a['xx']
                ftt=-cosine*f['x']+sine*sine*f['xx']
                art=-sine*a['Rx'];frt=-sine*f['Rx']
                g=nu-a['value'];inverse=1/sine;it=-cosine/sine**2;itt=1/sine+2*cosine**2/sine**3
                first=ft*g-f['value']*at
                derived['L_theta_theta']=first*inverse+f['value']*g*it
                derived['L_theta_theta_theta']=(ftt*g-2*ft*at-f['value']*att)*inverse+2*first*it+f['value']*g*itt
                # |R'|<=1. The symmetric interval includes all actual rates.
                partial=(frt*g-ft*a['R']-f['R']*at-f['value']*art)*inverse+(f['R']*g-f['value']*a['R'])*it
                derived['L_s_theta_theta']=iv.mpf([-magnitude(partial),magnitude(partial)])
            for name,value in derived.items():
                b=extra[name];b[0]=min(b[0],lower(value));b[1]=max(b[1],upper(value))
    # Endpoint-singular individual L derivatives are never used globally.
    for name in extra:
        if not all(math.isfinite(v) for v in extra[name]):extra[name]=None
    return dict(radial_tiles=radial,angular_tiles=angular,fields=bounds,derived=extra)


def edge_cover(iv,function,parameters,theta):
    lo,hi=math.inf,-math.inf
    for i in range(32):
        R=iv.mpf([iv.mpf(1)/40+iv.mpf(i)/1280,iv.mpf(1)/40+iv.mpf(i+1)/1280])
        value=evaluate(function,parameters,R,theta,iv)['a']['value']
        lo=min(lo,lower(value));hi=max(hi,upper(value))
    return [lo,hi]


def point_controls():
    function=field_function(np);original=symbolic_fields()[0];p=exact_target_parameters()
    errors=[]
    for R,theta in ((.026,.4),(.037,2.5),(.049,2.9)):
        got=np.asarray(function(R,np.cos(theta),p['k'],p['M'],p['a'],p['B'],p['Q'])).reshape(6,3)
        def reference(radius,angle):
            local=local_fields(original,p,radius,angle)
            _,_,angular=local_first(original,p,radius,angle)
            value=np.array([local['NP'],local['W'],local['a']])
            dr=np.array([local['NP']*local['log_NP_r'],local['W']*local['log_W_r'],local['ar']])
            dx=-np.array([local['NP']*(angular[0]+angular[1]),local['W']*(angular[1]-angular[0]),local['at']])/np.sin(angle)
            return value,dr,dx
        zero=reference(R,theta);tiny=1e-25j
        radial=reference(R+tiny,theta);angle=reference(R,theta+tiny)
        expected=np.array([*zero,np.imag(radial[1])/tiny.imag,-np.imag(angle[1])/tiny.imag/np.sin(theta),-np.imag(angle[2])/tiny.imag/np.sin(theta)])
        errors.append(float(np.max(abs(got-expected)/np.maximum(1,abs(expected)))))
    return errors


def constants(iv,global_cover,center,plateau,envelope):
    abs_bound=lambda item:iv.mpf(max(abs(item[0]),abs(item[1])))
    g,c,p=global_cover['fields'],center['fields'],plateau['derived']
    wmin=iv.mpf(g['W']['value'][0]);wmax=iv.mpf(g['W']['value'][1])
    center_wmin=iv.mpf(c['W']['value'][0]);NPmax=iv.mpf(g['NP']['value'][1])
    sine=iv.sin(iv.mpf(277)/100)
    atmin=iv.mpf(center['derived']['a_theta'][0]);atmax=abs_bound(center['derived']['a_theta'])
    drift=abs_bound(c['a']['R'])/atmin
    att=abs_bound(c['a']['x'])+abs_bound(c['a']['xx'])
    acceleration=(abs_bound(c['a']['RR'])+160*abs_bound(c['a']['R'])+2*abs_bound(c['a']['Rx'])*drift+att*drift**2)/atmin
    tilt=drift/20;tilt_s=drift+acceleration/20
    lammin=1/iv.sqrt(1+tilt**2);lam_s=tilt*tilt_s
    # Analytic global step bounds: S'<=8, |S''|<120, with collar width1/20.
    F,Fs,Ct=iv.mpf(4),iv.mpf(1200),iv.mpf(40)
    beta=tilt;beta_s=tilt_s;beta_dot=F*tilt
    beta_s_dot=F*tilt_s+Fs*tilt+2*F*tilt**2*tilt_s
    lam_dot=F*tilt**2
    lam_s_dot=2*F*tilt*tilt_s+Fs*tilt**2+3*F*tilt**3*tilt_s
    one_minus_lam=tilt**2/2
    M=iv.mpf(c['NP']['value'][1])*atmax/(2*sine)
    qmin=iv.mpf(1000);radius=iv.mpf(1)/50;clearance=iv.mpf(1)/25
    length=iv.mpf(min(lower(radius*iv.sqrt(qmin)),lower(1/iv.sqrt(2*M))))
    Zlower=2*center_wmin*length*iv.exp(-2*M*length**2)
    gamma=iv.mpf(envelope['log_envelope_coefficient']);gaussian=2*gamma*lammin
    moments=[]
    for j in range(1,5):
        factorial=math.prod(range(1,2*j,2))
        integral=factorial*iv.sqrt(iv.pi)/(2**j*gaussian**j*iv.sqrt(gaussian))
        moments.append(wmax/Zlower*integral)
    C1,C2,_,C4=moments
    L2,L3,Ls2=(abs_bound(p[n]) for n in ('L_theta_theta','L_theta_theta_theta','L_s_theta_theta'))
    # The polar singularities cancel before bounding: NP_R=a_R=0 at both poles.
    Lsmax=2*(abs_bound(g['NP']['Rx'])+NPmax*abs_bound(g['a']['Rx']))
    Vmax=NPmax/iv.sin(iv.mpf(3)/20);Lmax=iv.pi*Vmax
    Aout=lam_s*Lmax+Lsmax+drift*(Vmax+Ct*one_minus_lam*Lmax)
    Adotout=lam_s_dot*Lmax+lam_dot*Lsmax+drift*lam_dot*(Vmax+Ct*Lmax)
    max_upper=lambda *values:iv.mpf(max(upper(v) for v in values))
    A2=max_upper((lam_s*L2+Ls2+drift*L3)/2,Aout/clearance**2)
    Adot2=max_upper((lam_s_dot*L2+lam_dot*(Ls2+drift*L3))/2,Adotout/clearance**2)
    J2=lam_dot*max_upper(L2/2,Lmax/clearance**2)
    d0=abs_bound(global_cover['derived']['log_W_R'])+Ct*drift+drift*abs_bound(global_cover['derived']['log_W_theta'])
    b=A2*C1+d0/2
    bdot=Adot2*C1+2*J2*A2*(C2+C1**2)+2*d0*J2*C1
    H=J2*(iv.sqrt(C2)+C1);eta_dot=beta_dot/2
    # Combined angular coefficient vanishes on both C=0 and C=1 regions.
    # Interpolation error bounds for sin(C beta), cos(C beta) retain O(t^2).
    E1=tilt**2/8+tilt**3/4
    E1dot=F*(iv.mpf(17)/8*tilt**3+2*tilt**2)
    matrix=iv.sqrt(1+tilt**2)
    E=40*(E1*Vmax+matrix*Ct*one_minus_lam*Lmax)
    Edot=40*(E1dot*Vmax+Ct*Lmax*(lam_dot*matrix+one_minus_lam*F*tilt*(1+matrix/2)))
    K2=(A2+E/clearance**2)/wmin
    K2dot=((2*F+eta_dot)*A2+Adot2+(F*E+Edot)/clearance**2)/wmin
    D=(beta_s/2+40*Ct*beta/2)/wmin
    Ddot=(F*beta_s+beta_s_dot/2+beta_s*beta_dot/4+Fs+40*(F*Ct*beta/2+Ct*beta_dot/2+Ct*beta*beta_dot/4))/wmin
    B0=1/wmin;B0dot=(2*F+eta_dot)/wmin
    R0=K2*iv.sqrt(C2)+D+b*B0
    R0dot=K2dot*iv.sqrt(C2)+Ddot+bdot*B0+b*B0dot+J2*K2*(iv.sqrt(C4)+C1*iv.sqrt(C2))+H*(D+b*B0)+2*(H+eta_dot)*R0
    R1dot=B0dot+(3*H+2*eta_dot)*B0
    table={'W_lower':wmin,'W_upper':wmax,'center_W_lower':center_wmin,'NP_upper':NPmax,'a_theta_lower':atmin,
        'theta0_s_bound':drift,'theta0_ss_bound':acceleration,'tilt_bound':tilt,'tilt_s_bound':tilt_s,'lambda_lower':lammin,
        'local_lower_log_M':M,'Z_coefficient_lower':Zlower,'gaussian_exponent_lower':gaussian,
        'L2_bound':L2,'L3_bound':L3,'Ls2_bound':Ls2,'Ls_global_bound':Lsmax,'L_on_window_bound':Lmax,
        'A2_bound':A2,'Adot2_bound':Adot2,'J2_bound':J2,'normalization_b_bound':b,'normalization_bdot_bound':bdot,
        'log_weight_transport_bound':d0,'profile_metric_derivative_bound':H,'spin_metric_derivative_bound':eta_dot,
        'K2_bound':K2,'K2dot_bound':K2dot,'D_bound':D,'Ddot_bound':Ddot,'B0_bound':B0,'B0dot_bound':B0dot,
        'r0_bound_D':R0,'dot_r0_bound_D':R0dot,'r1_bound':B0,'dot_r1_bound':R1dot}
    table={name:[lower(value),upper(value)] for name,value in table.items()}
    return dict(constants=table,even_moment_coefficients=[upper(x) for x in moments],
        q_min=1000,local_radius=.02,plateau_clearance=.04,step_bounds=dict(S_first=8,S_second=120,F=4,F_s=1200,C_theta=40,R_second=160))


def run():
    started=perf_counter()
    path=project_path('dihole_aligned_channel_results.json');previous=json.loads(path.read_text(encoding='utf-8'))
    inputs=dict(previous['input_sha256']);inputs.update({'dihole_aligned_channel.py':previous['source_sha256'],path.name:sha(path)})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact prior inputs required.')
    iv=MPIntervalContext();iv.dps=60;parameters=calibration(iv);function=field_function(iv)
    edges=[edge_cover(iv,function,parameters,t) for t in (iv.mpf(11)/5,iv.mpf(11)/4)]
    regions={}
    for name,lo,hi in (('global',iv.mpf(0),iv.pi),('center',iv.mpf(109)/50,iv.mpf(277)/100),('plateau',iv.mpf(7)/20,iv.pi-iv.mpf(7)/20)):
        regions[name]=cover(iv,function,parameters,lo,hi)
        print(json.dumps({'region_complete':name,'elapsed_seconds':perf_counter()-started}),flush=True)
    bound=constants(iv,regions['global'],regions['center'],regions['plateau'],previous['anchored_log_envelope'])
    point=point_controls()
    c=bound['constants'];coeff=bound['even_moment_coefficients']
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        interior_flux_band=edges[0][1]<.85 and edges[1][0]>.95,
        window_clearance=bool(iv.mpf(11)/4+iv.mpf(1)/25<iv.pi-iv.mpf(7)/20),
        positive_weights=c['W_lower'][0]>0 and c['center_W_lower'][0]>0,
        positive_angular_field=c['a_theta_lower'][0]>0,
        positive_width=c['lambda_lower'][0]>0,
        positive_normalization=c['Z_coefficient_lower'][0]>0,
        finite_positive_moments=all(math.isfinite(x) and x>0 for x in coeff),
        independent_field_jets=max(point)<1e-8,
        finite_residual_bounds=all(math.isfinite(c[n][1]) and c[n][1]>0 for n in ('r0_bound_D','dot_r0_bound_D','r1_bound','dot_r1_bound')),
        saved_high_flux_states_enclosed=all(x['r0_norm_D']<=c['r0_bound_D'][1] and x['dot_r0_norm_D']<=c['dot_r0_bound_D'][1] and x['r1_norm']<=c['r1_bound'][1] and x['dot_r1_norm']<=c['dot_r1_bound'][1] for x in previous['records'] if x['q']>=1000))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,edge_potential_enclosures=edges,regions=regions,uniform_bound=bound,
        point_jet_relative_errors=point,band=dict(nu=[.85,.95],R_over_D=[.025,.05],theta_center=[2.2,2.75]),
        scope=['Uniform collar trial residual bounds for q>=1000 and regular half-integer mu=q nu in the declared band.',
            'Bounds include pole-restoring transitions and normalization/metric derivatives. They are deliberately conservative and are not exact residual values.',
            'Baseline geometry and first paired metric derivative only, not a finite-epsilon stability interval.',
            'No finite-q spectral gap, exterior/null/open channels, global source trace, supported geometry or atom propagation is established.'],
        full_source_bound=None,atomic_transit=None,source_sha256=sha(__file__),input_sha256=inputs,elapsed_seconds=perf_counter()-started)
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:raise ArithmeticError('Uniform channel bound checks failed; results retained.')
    print(json.dumps({'passed':result['passed'],'tests':tests,'uniform_bound':bound},indent=2))


if __name__=='__main__':run()
