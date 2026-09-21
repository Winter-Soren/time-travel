"""Metric-inserted finite collar trace, with full radial graph control.

Fixed and moving effective filters are distinguished. Neither is a vacuum
stress calculation, and the global realization must preserve axial sectors.
"""
import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

from balanced_dihole import symbolic_fields
from dihole_aligned_channel import geometry, norm, trial
from dihole_channel_bounds import MPIntervalContext, upper
from dihole_collar_resolvent import bounds as collar_bounds
from dihole_local_source import local_first
from dihole_projected_gap import bounds as gap_bounds
from dihole_slice_channel import sha
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME='dihole_metric_resolvent'


def bounds(iv,previous,window,q):
    base=collar_bounds(iv,previous,q);gap=gap_bounds(iv,previous,window,q)
    c=previous['uniform_bound']['constants'];u=previous['uniform_bound']
    val=lambda key:iv.mpf(c[key][1])
    absolute=lambda key:iv.mpf(max(abs(x) for x in previous['regions']['global']['derived'][key]))
    w=iv.mpf(c['W_lower'][0]);W=val('W_upper');F=iv.mpf(4);Fs=iv.mpf(1200)
    T,Ts,drift=val('tilt_bound'),val('tilt_s_bound'),val('theta0_s_bound')
    Lwt,Lws=absolute('log_W_theta'),absolute('log_W_R');Ct=iv.mpf(40)
    C1,C2,C3=map(iv.mpf,u['even_moment_coefficients'][:3])
    A2,Ad2,J2,d0=(val(k) for k in ('A2_bound','Adot2_bound','J2_bound','log_weight_transport_bound'))
    R1=iv.mpf(base['r1_upper']);v=iv.mpf(base['velocity_lower']);vs=iv.mpf(base['velocity_derivative_upper_D'])
    M=iv.mpf(gap['profile_metric_derivative_norm_upper']);Lf=iv.mpf(base['f_theta_lipschitz_upper'])
    # Normalization removes the component of dot(chi) parallel to chi.
    dotr1=2*F*R1+2*M*R1+2*M/w
    alpha=2*F+2*M*R1/v
    Ldf=F*T*T*(2*Ct+Lwt)/w
    G=(Ts*T+Lws)/w+drift*Lf
    dotG=(F*Ts*T+T*(Fs*T+F*Ts+2*F*T*T*Ts)+F*T*T*Lws)/w+drift*Ldf
    K5=iv.sqrt(C2*C3)+2*C1*iv.sqrt(C1*C2)+iv.sqrt(C1)*C2+4*C1*C1*iv.sqrt(C1)
    pairs=2*(Ad2*Lf+A2*Ldf)*iv.sqrt(C1*C2)+d0*Ldf*iv.sqrt(C1)
    triples=4*A2*Lf*J2*K5+4*d0*Lf*J2*(iv.sqrt(C1*C2)+3*C1*iv.sqrt(C1))
    dotvs=2*F*vs+2*Fs/w+dotG+2*J2*iv.sqrt(C2)*G+(pairs+triples)/iv.sqrt(q)
    # Export before composition so every downstream bound uses declared inputs.
    dotr1,alpha,dotvs=map(lambda x:iv.mpf(upper(x)),(dotr1,alpha,dotvs))
    beta=(dotvs+alpha*vs)/2
    ad=dotr1/v;bd=dotr1*vs/(2*v)+val('dot_r0_bound_D')+800*(2*F+M)/w
    fields=previous['regions']['global']['fields']
    rx=lambda name:iv.mpf(max(abs(x) for x in fields[name]['Rx']))
    amin,amax=map(iv.mpf,fields['a']['value'])
    nua=iv.mpf(max(upper(iv.mpf(19)/20-amin),upper(amax-iv.mpf(17)/20)))
    magnetic=40*(rx('NP')*nua+val('NP_upper')*rx('a'))
    if upper(magnetic)>=31:raise ArithmeticError('Pole-safe magnetic bound not enclosed by 31q.')
    Cg=400+31*iv.mpf(q)
    rows=[]
    for b in base['frequency_rows']:
        y=iv.mpf(b['y_D']);e=iv.mpf(b['resolvent_error_upper_over_D']);hs2=iv.mpf(b['filter_HS_norm_squared_upper'])
        CH=F+iv.sqrt((F*iv.sqrt(W)+Fs/(iv.sqrt(w)*y))**2+Cg*F*F/(w*y*y))/iv.sqrt(w)
        Ch=iv.sqrt(alpha**2+(beta/y)**2);CT=iv.sqrt(ad**2+(bd/y)**2)
        derivative=(CH+Ch)*e+CT/y
        inserted=(derivative+M*e)*hs2
        moving=inserted+iv.sqrt(iv.mpf(11))*alpha*e*hs2
        rows.append(dict(y_D=b['y_D'],filter_Lambda_D=b['filter_Lambda_D'],
            full_metric_graph_upper=upper(CH),reduced_metric_graph_upper=upper(Ch),defect_metric_graph_upper=upper(CT),
            resolvent_derivative_upper_over_D=upper(derivative),scaled_resolvent_derivative_upper=upper(y*derivative),
            fixed_filter_derivative_and_insertion_error_upper_over_D=upper(inserted),
            moving_filter_derivative_error_upper_over_D=upper(moving),
            band_insertion_error_upper_over_D=upper(inserted*base['band']['count'])))
    # At fixed Lambda this positive 1/y envelope cannot be integrated to infinity.
    tail=(F*(1+iv.sqrt(W/w))+alpha+M)*iv.mpf(base['graph_h_coefficient_upper'])+ad
    return dict(q=q,band=base['band'],profile_metric_norm_upper=upper(M),dot_r1_upper=upper(dotr1),
        relative_velocity_derivative_upper=upper(alpha),velocity_mixed_derivative_upper_D=upper(dotvs),
        reduced_graph_zero_upper_D=upper(beta),defect_graph_h_upper=upper(ad),defect_graph_zero_upper_D=upper(bd),
        magnetic_radial_coefficient_upper=upper(magnetic),
        fixed_filter_error_envelope_tail_coefficient_per_HS_squared=upper(tail),frequency_rows=rows)


def field_controls(iv,previous,window):
    function=symbolic_fields()[0];parameters=exact_target_parameters();records=[];tiny=1e-25
    for q in (1005,10005):
        limit=bounds(iv,previous,window,q)
        for u,nu in ((.0125,.85),(.025,.9),(.0375,.95)):
            mu=math.floor(q*nu)+.5;data=geometry(function,parameters,u,mu/q,1024,80)
            state=trial(data,q);varied=trial(data,q,1j*tiny)
            dotchi=np.imag(varied['chi'])/tiny;dotr1=np.imag(varied['r1'])/tiny
            dotv=np.imag(varied['velocity'])/tiny;dotvs=np.imag(varied['velocity_s'])/tiny
            F,Fs=data['mask'][:2];W=data['W'];C,Ct=data['window'],data['window_theta']
            _,_,angular=local_first(function,parameters,data['radius'],data['theta'])
            Lwt=angular[1]-angular[0];drift=data['limit']['peak_drift_D'];zeta=C*drift
            d=data['log_W_s']+Ct*drift+zeta*Lwt
            def transported(epsilon):
                st=trial(data,q,epsilon);p=st['probability'];factor=np.exp(epsilon*F)
                t=factor*data['tilt'];ts=factor*(data['tilt_s']+epsilon*Fs*data['tilt'])
                beta=C*np.arctan(t);bt=Ct*np.arctan(t);bs=C*ts/(1+t*t)
                lam=1/np.sqrt(1+t*t);ls=-t*ts/(1+t*t)**1.5
                f=np.cos(beta)/W;ft=(-np.sin(beta)*bt-np.cos(beta)*Lwt)/W
                fs=(-np.sin(beta)*bs-np.cos(beta)*data['log_W_s'])/W
                P=C*ls*data['log_profile_per_q']+(1+C*(lam-1))*data['log_profile_s_per_q']
                A=P+zeta*st['log_t']/q
                cov=lambda a,b:p@((a-p@a)*(b-p@b))
                return factor**2*(p@(fs+zeta*ft)+2*q*cov(A,f)+cov(d,f)+2*epsilon*Fs*(p@f))
            sz=np.array([[1.,0.],[0.,-1.]])
            independent=2*F*state['r1']-(dotv-2*F*state['velocity'])*state['chi']+(sz@dotchi)/W-state['velocity']*dotchi
            records.append(dict(q=q,u=u,mu=mu,dot_r1=norm(data,dotr1),relative_velocity_derivative=float(abs(dotv/state['velocity'])),
                velocity_mixed_derivative=float(dotvs),dot_r1_identity_error=norm(data,dotr1-independent),
                velocity_transport_derivative_error=float(abs(np.imag(transported(1j*tiny))/tiny-dotvs)),
                bounds_enclose=bool(norm(data,dotr1)<=limit['dot_r1_upper'] and abs(dotv/state['velocity'])<=limit['relative_velocity_derivative_upper'] and abs(dotvs)<=limit['velocity_mixed_derivative_upper_D'])))
    return records


def matrix_controls():
    rng=np.random.default_rng(168);raw=rng.normal(size=(7,7))+1j*rng.normal(size=(7,7))
    H=(raw+raw.conj().T)/2;J=np.linalg.qr(raw)[0][:,:3];Q=np.eye(7)-J@J.conj().T
    dJ=Q@raw@J/10;K=dJ@J.conj().T-J@dJ.conj().T
    F=np.diag(np.linspace(-.3,.7,7))
    # The angular term must commute with the scalar radial mask.
    Ha=np.diag(np.linspace(-1.,1.,7));Hr=H-Ha;dH=F@Hr+Hr@F+F@Ha
    h=J.conj().T@H@J;dh=dJ.conj().T@H@J+J.conj().T@dH@J+J.conj().T@H@dJ
    T=H@J-J@h;dT=dH@J+H@dJ-dJ@h-J@dh;X=expm(-h@h/2)
    records=[]
    for y in (.3,1.,3.):
        GH=np.linalg.inv(H-1j*y*np.eye(7));Gh=np.linalg.inv(h-1j*y*np.eye(3))
        E=GH@J-J@Gh;Em=GH.conj().T@J-J@Gh.conj().T
        dE=GH@dH@GH@T@Gh-GH@dT@Gh+GH@T@Gh@dh@Gh
        insert=J.conj().T@GH@dH@GH@J-Gh@dh@Gh
        identity=-J.conj().T@dE+Em.conj().T@dJ
        def varied(epsilon):
            c=expm(epsilon*F);he=c@Hr@c+c@Ha;je=expm(epsilon*K)@J;small=je.conj().T@he@je
            ge=np.linalg.inv(he-1j*y*np.eye(7));gs=np.linalg.inv(small-1j*y*np.eye(3))
            ee=ge@je-je@gs;xe=expm(-small@small/2)
            return ee,np.trace(X.conj().T@je.conj().T@ee@X),np.trace(xe.conj().T@je.conj().T@ee@xe)
        step=1e-5;plus,minus=varied(step),varied(-step)
        dX=(expm(-(h+step*dh)@(h+step*dh)/2)-expm(-(h-step*dh)@(h-step*dh)/2))/(2*step)
        fixed=np.trace(X.conj().T@(dJ.conj().T@E+J.conj().T@dE)@X)
        moving=fixed+np.trace(dX.conj().T@J.conj().T@E@X+X.conj().T@J.conj().T@E@dX)
        size=lambda a:np.linalg.norm(a,2)
        e=size(GH)*size(T@Gh);cap=(size(dH@GH)+size(dh@Gh)+size(dJ))*e+size(dT@Gh)/y
        records.append(dict(y=y,parallel_phase_error=size(J.conj().T@dJ),resolvent_derivative_error=size((plus[0]-minus[0])/(2*step)-dE),
            insertion_identity_error=size(insert-identity),fixed_trace_derivative_error=float(abs((plus[1]-minus[1])/(2*step)-fixed)),
            moving_trace_derivative_error=float(abs((plus[2]-minus[2])/(2*step)-moving)),
            insertion_bound_holds=bool(abs(np.trace(X.conj().T@insert@X))<=cap*np.linalg.norm(X,'fro')**2)))
    return records


def gaussian_controls():
    weight=lambda x:math.exp(-x*x)/math.sqrt(math.pi)
    upper_coefficient=quad(lambda x:(1+x*x)**2*weight(x),-np.inf,np.inf)[0]
    scaling_coefficient=quad(lambda x:(x*x-1)**2*weight(x),-np.inf,np.inf)[0]
    return dict(general_squared_coefficient=upper_coefficient,uniform_velocity_squared_coefficient=scaling_coefficient)


def run():
    path=project_path('dihole_collar_resolvent_results.json');previous_collar=json.loads(path.read_text())
    inputs=dict(previous_collar['input_sha256']);inputs.update({path.name:sha(path),'dihole_collar_resolvent.py':previous_collar['source_sha256']})
    previous=json.loads(project_path('dihole_channel_bounds_results.json').read_text())
    gap=json.loads(project_path('dihole_projected_gap_results.json').read_text());window=gap['window_cover']
    if not previous_collar['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing prior inputs required.')
    iv=MPIntervalContext();iv.dps=60
    rows=[bounds(iv,previous,window,x['q']) for x in previous_collar['records']]
    fields=field_controls(iv,previous,window);matrices=matrix_controls();gaussian=gaussian_controls()
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        pole_safe_magnetic_enclosure=all(x['magnetic_radial_coefficient_upper']<31 for x in rows),
        actual_momentum_derivative_identity=max(x['dot_r1_identity_error'] for x in fields)<1e-12,
        actual_velocity_mixed_identity=max(x['velocity_transport_derivative_error'] for x in fields)<1e-8,
        actual_metric_quantities_enclosed=all(x['bounds_enclose'] for x in fields),
        parallel_phase=max(x['parallel_phase_error'] for x in matrices)<1e-14,
        finite_matrix_resolvent_derivative=max(x['resolvent_derivative_error'] for x in matrices)<2e-8,
        metric_insertion_identity=max(x['insertion_identity_error'] for x in matrices)<1e-12,
        fixed_trace_derivative=max(x['fixed_trace_derivative_error'] for x in matrices)<2e-8,
        moving_trace_derivative=max(x['moving_trace_derivative_error'] for x in matrices)<2e-8,
        finite_matrix_trace_bounds=all(x['insertion_bound_holds'] for x in matrices),
        gaussian_derivative_coefficients=abs(gaussian['general_squared_coefficient']-11/4)<1e-12 and abs(gaussian['uniform_velocity_squared_coefficient']-3/4)<1e-12,
        finite_exported_bounds=all(math.isfinite(v) and v>0 for x in rows for row in x['frequency_rows'] for v in row.values()))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=rows,field_controls=fields,matrix_controls=matrices,gaussian_controls=gaussian,
        source_sha256=sha(__file__),input_sha256=inputs,total_source=None,source_sign=None,atom_transit=None,
        current_error_envelope_absolutely_frequency_integrable=False,
        scope=['Baseline and first paired metric variation on the actual collar, with fixed y,Lambda,chi,eta; no finite-epsilon stability theorem.',
            'Global self-adjoint realization must preserve the declared axial sector; the broader exterior-mixing class of167 is not covered by this new graph estimate.',
            'Full radial graph control uses F itself as a compact multiplier. Polar cancellations retained. No global spectral gap assumed.',
            'Finite metric-inserted resolvent trace difference is bounded, not measured; effective h variation includes the moving normalized trial.',
            'Fixed-filter and varying-effective-Gaussian derivative bounds are separate. Neither is the unfiltered renormalized vacuum stress.',
            'Atomic-reference entries evaluate bounding formulas and exact sector counts only; no huge-flux wave function, supported opening, atom survival or clock shift.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Metric resolvent controls failed; evidence retained.')
    plot(result)
    print(json.dumps({'passed':result['passed'],'tests':tests,'atomic_bound':rows[-1]},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    qs=[x['q'] for x in result['records']];rows=[x['frequency_rows'][1] for x in result['records']]
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    axes[0].loglog(qs,[r['full_metric_graph_upper'] for r in rows],'o-',label='Full geometry response')
    axes[0].loglog(qs,[r['reduced_metric_graph_upper'] for r in rows],'s-',label='Reduced channel response')
    axes[0].set(title='Metric-response graph bounds',ylabel='Dimensionless upper bound');axes[0].legend()
    axes[1].loglog(qs,[r['fixed_filter_derivative_and_insertion_error_upper_over_D'] for r in rows],'o-',label='Metric insertion / fixed filter')
    axes[1].loglog(qs,[r['moving_filter_derivative_error_upper_over_D'] for r in rows],'s-',label='Also varying effective filter')
    axes[1].set(title='Finite selected trace: error bounds',ylabel='Per-sector bound / D');axes[1].legend()
    for ax in axes:ax.set_xlabel('Integer flux q (evaluated bounding formulas)');ax.grid(alpha=.2)
    fig.suptitle('One-atom project: geometry variation is now included')
    fig.supxlabel('y D = Lambda D = sqrt(q). Axial sectors preserved. No stress value, stable opening or atom passage established.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
