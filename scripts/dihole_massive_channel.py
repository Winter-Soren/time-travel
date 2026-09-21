"""Proper-mass channel defects at fixed regulator and nonzero frequency.

The mass multiplies the actual lapse. These finite selected comparisons do
not remove the regulator or give a renormalized supporting source.
"""
import json
import math

import numpy as np
from scipy.linalg import expm

from balanced_dihole import symbolic_fields
from dihole_aligned_channel import geometry, norm, trial
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_local_source import local_first
from dihole_slice_channel import sha
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME='dihole_massive_channel'
COEFFICIENTS=(1,-3,3,-1)


def lapse_bounds(iv,previous,window,q):
    c=previous['uniform_bound']['constants'];region=previous['regions']['global']
    val=lambda k:iv.mpf(c[k][1])
    absolute=lambda xs:iv.mpf(max(abs(x) for x in xs))
    npmin,npmax=map(iv.mpf,region['fields']['NP']['value'])
    w=iv.mpf(c['W_lower'][0]);W=val('W_upper')
    nt=absolute(region['fields']['NP']['x'])/npmin+absolute(region['derived']['log_W_theta'])
    ns=absolute(region['fields']['NP']['R'])/npmin+absolute(region['derived']['log_W_R'])
    nmin=iv.sqrt(npmin/W);nmax=iv.sqrt(npmax/w);LN=nmax*nt/2;Ns=nmax*ns/2
    C1,C2,C3=map(iv.mpf,previous['uniform_bound']['even_moment_coefficients'][:3])
    K6=C3+2*iv.sqrt(C1*C2*C3)+3*C1*C2+4*C1*C1*iv.sqrt(C2)+6*C1**3
    Mwindow=absolute(window['derived']['L_theta_theta'])/2
    F=iv.mpf(4);T=val('tilt_bound');J2=F*T*T*Mwindow
    RN=LN*iv.sqrt(C1/q)
    dotRN=F*RN+LN/iv.sqrt(q)*iv.sqrt(J2*J2*K6+F*F*T*T*C1/4)
    averaged_s=Ns+val('theta0_s_bound')*LN+LN*iv.sqrt(C1/q)*(2*val('A2_bound')*iv.sqrt(C2)+val('log_weight_transport_bound'))
    CL=iv.sqrt(npmax*W)*(ns/2+40)
    return dict(lapse_lower=lower(nmin),lapse_upper=upper(nmax),lapse_theta_lipschitz_upper=upper(LN),
        lapse_s_upper_D=upper(Ns),mass_defect_upper=upper(RN),metric_mass_defect_upper=upper(dotRN),
        averaged_lapse_s_upper_D=upper(averaged_s),metric_log_quadratic_upper=upper(J2),
        sixth_moment_combination_upper=upper(K6),mass_radial_graph_coefficient_upper_D=upper(CL))


def bounds(iv,base,metric,lapse,constants):
    q=base['q'];reference=base['frequency_rows'][1];y=iv.mpf(reference['y_D'])
    w=iv.mpf(constants['W_lower'][0]);W=iv.mpf(constants['W_upper'][1])
    a=iv.mpf(base['graph_h_coefficient_upper']);b=iv.mpf(base['graph_zero_coefficient_upper_D'])
    ad=iv.mpf(metric['defect_graph_h_upper']);bd=iv.mpf(metric['defect_graph_zero_upper_D'])
    alpha=iv.mpf(metric['relative_velocity_derivative_upper']);beta=iv.mpf(metric['reduced_graph_zero_upper_D'])
    Mprofile=iv.mpf(metric['profile_metric_norm_upper']);n=iv.mpf(lapse['lapse_upper'])
    RN=iv.mpf(lapse['mass_defect_upper']);dotRN=iv.mpf(lapse['metric_mass_defect_upper'])
    CL=iv.mpf(lapse['mass_radial_graph_coefficient_upper_D']);dotn=4*n+2*Mprofile*RN
    hs2=2*iv.mpf(reference['filter_HS_norm_squared_upper']);cases=[]
    for ratio in (.1,1.,10.):
        # M D is the declared binary float; sqrt(j)*M is then interval evaluated.
        regulator_value=ratio*reference['y_D'];regulator=iv.mpf(regulator_value);rows=[]
        for j in range(4):
            mass=iv.sqrt(iv.mpf(j))*regulator
            bm=b+mass*(a*n+RN);bdm=bd+mass*(ad*n+dotRN)
            betam=beta+mass*(alpha*n+dotn)
            CT=iv.sqrt(a*a+(bm/y)**2);dotCT=iv.sqrt(ad*ad+(bdm/y)**2)
            CH=4+iv.sqrt((4*iv.sqrt(W)+1200/(iv.sqrt(w)*y))**2+(400+31*q+mass*CL)*16/(w*y*y))/iv.sqrt(w)
            Ch=iv.sqrt(alpha*alpha+(betam/y)**2);e=CT/y
            insertion=hs2*((Mprofile+CH+Ch)*e+dotCT/y)
            rows.append(dict(j=j,absolute_mass_D_interval=[max(0.,lower(mass)),upper(mass)],
                scaled_resolvent_error_upper=upper(CT),resolvent_error_upper_over_D=upper(e),
                full_metric_graph_upper=upper(CH),reduced_metric_graph_upper=upper(Ch),
                defect_metric_graph_upper=upper(dotCT),
                finite_metric_insertion_error_upper_over_D=upper(insertion)))
        # Only a finite-frequency trace combination: not the sign source.
        combined=sum(iv.mpf(abs(c))*iv.mpf(row['finite_metric_insertion_error_upper_over_D']) for c,row in zip(COEFFICIENTS,rows))/2
        cases.append(dict(regulator_M_D=regulator_value,nominal_M_over_y=ratio,masses=rows,
            finite_signed_mass_trace_combination_error_upper_over_D=upper(combined)))
    return dict(q=q,y_D=reference['y_D'],filter_Lambda_D=reference['filter_Lambda_D'],
        four_component_filter_HS_squared_upper=upper(hs2),averaged_lapse_metric_upper=upper(dotn),
        lapse=lapse,cases=cases)


def field_controls(limits):
    function=symbolic_fields()[0];parameters=exact_target_parameters();records=[];tiny=1e-25
    for q in (1005,10005):
        limit=limits[q]
        for u,nu in ((.0125,.85),(.025,.9),(.0375,.95)):
            mu=math.floor(q*nu)+.5;data=geometry(function,parameters,u,mu/q,1024,80)
            st=trial(data,q);varied=trial(data,q,1j*tiny);F=data['mask'][0]
            logs,radial,angular=local_first(function,parameters,data['radius'],data['theta'])
            N=np.exp(logs[0]);Ns=N*radial[0]*data['radius_rate'];Nt=N*angular[0]
            chi=st['chi'];dotchi=np.imag(varied['chi'])/tiny;p=st['probability'];n=p@N
            r=(N-n)*chi;dotN=F*N;variedN=np.exp(1j*tiny*F)*N
            variedn=varied['probability']@variedN;dotn=float(np.imag(variedn)/tiny)
            dr=np.imag((variedN-variedn)*varied['chi'])/tiny
            inner=lambda left,right:data['weights']@np.sum(left*right,axis=0)
            direct=F*n+2*inner(dotchi,r);reflected=F*r+(N-n)*dotchi
            ns=p@Ns+2*inner(st['chi_s'],r)
            Lwt=data['log_W_s'];drift=data['limit']['peak_drift_D']
            zeta=data['window']*drift;C=data['window'];Ct=data['window_theta']
            t,ts=data['tilt'],data['tilt_s'];lam=1/np.sqrt(1+t*t);ls=-t*ts/(1+t*t)**1.5
            A=C*ls*data['log_profile_per_q']+(1+C*(lam-1))*data['log_profile_s_per_q']+zeta*st['log_t']/q
            d=Lwt+Ct*drift+zeta*(angular[1]-angular[0])
            cov=lambda a,b:p@((a-p@a)*(b-p@b))
            transported=p@(Ns+zeta*Nt)+2*q*cov(A,N)+cov(d,N)
            records.append(dict(q=q,u=u,mu=mu,averaged_lapse=float(n),mass_defect=norm(data,r),
                metric_mass_defect=norm(data,dr),averaged_lapse_s_D=float(ns),
                normalization_derivative_error=abs(dotn-direct),mass_derivative_identity_error=norm(data,dr-((dotN-dotn)*chi+(N-n)*dotchi)),
                reflection_norm_error=abs(norm(data,dr)-norm(data,reflected)),
                transported_lapse_derivative_error=abs(float(ns-transported)),
                bounds_enclose=bool(limit['lapse_lower']<=n<=limit['lapse_upper'] and norm(data,r)<=limit['mass_defect_upper'] and norm(data,dr)<=limit['metric_mass_defect_upper'] and abs(ns)<=limit['averaged_lapse_s_upper_D'])))
    return records


def matrix_controls():
    rng=np.random.default_rng(170);raw=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
    H=(raw+raw.conj().T)/2;U=np.linalg.qr(raw)[0][:,:2];V=np.diag(np.linspace(-.2,.6,5))
    N=np.diag(np.linspace(.4,1.7,5));f=.7;dU=(np.eye(5)-U@U.conj().T)@raw@U/10
    K=dU@U.conj().T-U@dU.conj().T
    h=U.conj().T@H@U;n=U.conj().T@N@U
    dh=dU.conj().T@H@U+U.conj().T@V@U+U.conj().T@H@dU
    dn=f*n+dU.conj().T@N@U+U.conj().T@N@dU
    t1=np.array([[0.,1.],[1.,0.]]);t3=np.diag([1.,-1.]);I=np.eye(2)
    U4=np.kron(I,U);dU4=np.kron(I,dU);X=np.kron(I,expm(-h@h/2));rows=[]
    R=H@U-U@h;r=N@U-U@n
    for mass in (0.,.7,3.):
        Hm=np.kron(t3,H)+mass*np.kron(t1,N);hm=np.kron(t3,h)+mass*np.kron(t1,n)
        Vm=np.kron(t3,V)+mass*f*np.kron(t1,N);vm=np.kron(t3,dh)+mass*np.kron(t1,dn)
        Rm=Hm@U4-U4@hm;dRm=Vm@U4+Hm@dU4-dU4@hm-U4@vm
        y=.8;G=np.linalg.inv(Hm-1j*y*np.eye(10));g=np.linalg.inv(hm-1j*y*np.eye(4))
        E=G@U4-U4@g;Em=G.conj().T@U4-U4@g.conj().T
        dE=G@Vm@G@Rm@g-G@dRm@g+G@Rm@g@vm@g
        inserted=U4.conj().T@G@Vm@G@U4-g@vm@g
        def varied(epsilon):
            u=expm(epsilon*K)@U;big=np.kron(t3,H+epsilon*V)+mass*np.exp(epsilon*f)*np.kron(t1,N)
            small=np.kron(t3,u.conj().T@(H+epsilon*V)@u)+mass*np.exp(epsilon*f)*np.kron(t1,u.conj().T@N@u)
            embed=np.kron(I,u)
            return np.linalg.solve(big-1j*y*np.eye(10),embed)-embed@np.linalg.inv(small-1j*y*np.eye(4))
        step=1e-5;size=lambda a:float(np.linalg.norm(a,2))
        cap=((size(Vm@G)+size(vm@g)+size(dU4))*size(Rm@g)+size(dRm@g))/y
        chirality=np.kron(t3,np.eye(5));negative=np.kron(t3,H)-mass*np.kron(t1,N)
        negativeV=np.kron(t3,V)-mass*f*np.kron(t1,N)
        def source(big,variation):
            values,vectors=np.linalg.eigh(big);sign=(vectors*np.sign(values))@vectors.conj().T
            return float(np.real(np.trace(X.conj().T@U4.conj().T@variation@sign@U4@X)))
        rows.append(dict(mass=mass,projection_error=size(Rm-(np.kron(t3,R)+mass*np.kron(t1,r))),
            square_compression_error=size(U4.conj().T@Hm@Hm@U4-hm@hm-Rm.conj().T@Rm),
            square_compression_shortcut_mismatch=size(Rm.conj().T@Rm),
            resolvent_identity_error=size(E+G@Rm@g),metric_derivative_error=size((varied(step)-varied(-step))/(2*step)-dE),
            metric_insertion_identity_error=size(inserted+U4.conj().T@dE-Em.conj().T@dU4),
            finite_trace_bound_holds=bool(abs(np.trace(X.conj().T@inserted@X))<=cap*np.linalg.norm(X,'fro')**2),
            signed_mass_unitary_error=size(negative-chirality@Hm@chirality),signed_mass_source_difference=abs(source(Hm,Vm)-source(negative,negativeV))))
    P=raw.conj().T@raw+.2*np.eye(5);Pv=U.conj().T@P@U;identities=[]
    for matrix in (P,Pv):
        eye=np.eye(len(matrix));M=.9;eta=.4;factors=[matrix+(eta+j*M*M)*eye for j in range(4)]
        lhs=sum(c*np.linalg.inv(a) for c,a in zip(COEFFICIENTS,factors));product=eye
        for a in factors:product=product@a
        rhs=6*M**6*np.linalg.inv(product);identities.append(float(np.linalg.norm(lhs-rhs,2)))
    coefficient=-3+3*math.sqrt(2)-math.sqrt(3);growth=[]
    for M in (1.,10.,100.):
        selected=0.
        for j,c in enumerate(COEFFICIENTS):
            for sign in (-1,1):
                Hm=sign*math.sqrt(j)*M*t1
                values,vectors=np.linalg.eigh(Hm);signH=(vectors*np.sign(values))@vectors.conj().T
                selected+=c*float(np.trace(Hm@signH))/4
        growth.append(dict(M=M,signed_pair_average=selected,analytic_error=abs(selected-coefficient*M)))
    return dict(records=rows,common_positive_operator_pv_errors=identities,
        finite_selection_mass_growth_coefficient=coefficient,
        finite_selection_mass_growth=growth,
        polynomial_pv_moments=[sum(c*j**k for j,c in enumerate(COEFFICIENTS)) for k in range(3)])


def run():
    path=project_path('dihole_paired_frequency_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_paired_frequency.py':previous['source_sha256']})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing prior inputs required.')
    collar=json.loads(project_path('dihole_collar_resolvent_results.json').read_text())
    metric=json.loads(project_path('dihole_metric_resolvent_results.json').read_text())
    fields=json.loads(project_path('dihole_channel_bounds_results.json').read_text())
    window=json.loads(project_path('dihole_projected_gap_results.json').read_text())['window_cover']
    iv=MPIntervalContext();iv.dps=60
    limits={b['q']:lapse_bounds(iv,fields,window,b['q']) for b in collar['records']}
    records=[bounds(iv,b,m,limits[b['q']],fields['uniform_bound']['constants']) for b,m in zip(collar['records'],metric['records'],strict=True)]
    profiles=field_controls(limits);matrices=matrix_controls();mr=matrices['records']
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        actual_lapse_bounds=all(r['bounds_enclose'] for r in profiles),
        normalized_lapse_metric_identity=max(r['normalization_derivative_error'] for r in profiles)<1e-12,
        actual_mass_metric_identity=max(r['mass_derivative_identity_error'] for r in profiles)<1e-12,
        normalization_reflection_identity=max(r['reflection_norm_error'] for r in profiles)<1e-12,
        lapse_transport_identity=max(r['transported_lapse_derivative_error'] for r in profiles)<1e-10,
        massive_projection_identity=max(r['projection_error'] for r in mr)<1e-12,
        massive_square_compression=max(r['square_compression_error'] for r in mr)<1e-12 and min(r['square_compression_shortcut_mismatch'] for r in mr)>1,
        massive_resolvent_identity=max(r['resolvent_identity_error'] for r in mr)<1e-12,
        massive_metric_derivative=max(r['metric_derivative_error'] for r in mr)<2e-8,
        massive_metric_insertion=max(r['metric_insertion_identity_error'] for r in mr)<1e-12,
        finite_trace_enclosures=all(r['finite_trace_bound_holds'] for r in mr),
        signed_mass_evenness=max(r['signed_mass_source_difference'] for r in mr)<1e-12 and max(r['signed_mass_unitary_error'] for r in mr)<1e-12,
        common_operator_pv=max(matrices['common_positive_operator_pv_errors'])<1e-12 and matrices['polynomial_pv_moments']==[0,0,0],
        finite_selection_not_uniform_in_mass=abs(matrices['finite_selection_mass_growth_coefficient'])>.4 and max(r['analytic_error'] for r in matrices['finite_selection_mass_growth'])<1e-11,
        finite_exported_bounds=all(math.isfinite(c['finite_signed_mass_trace_combination_error_upper_over_D']) and c['finite_signed_mass_trace_combination_error_upper_over_D']>0 for r in records for c in r['cases']))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=records,field_controls=profiles,matrix_controls=matrices,
        source_sha256=sha(__file__),input_sha256=inputs,total_renormalized_source=None,source_sign=None,atom_transit=None,
        scope=['Fixed finite proper masses +/-sqrt(j)*M and nonzero y; Gaussian spatial filter frozen and mass independent. Both Weyl blocks, one collar, one axial sector.',
            'The scalar reduced lapse has a chosen bounded positive complete-line continuation inside the exported collar bounds, with its metric variation fixed outside the collar. No spatial derivative bound on that continuation is used.',
            'The full exterior lapse may be unbounded. Full metric graph estimates use lapse data only where the compact mask is supported; global axial self-adjoint realization assumed.',
            'Finite resolvent-insertion error and weighted finite-mass error are not sign-source values. No frequency integral, regulator removal, omitted-mode bound or physical counterterms evaluated.',
            'The exact common-P rational identity survives common form compression; it does not identify a squared compressed Hamiltonian with the compressed full square.',
            'Atomic-flux entries evaluate bounding formulas only. No supported opening, payload propagation, stability, future clock advantage or new law established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Massive channel checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'],tests=tests,atomic_bound=records[-1]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    qs=[r['q'] for r in result['records']]
    for key,label,marker in (('mass_defect_upper','Lapse leakage','o'),('metric_mass_defect_upper','Its geometry response','s')):
        axes[0].loglog(qs,[r['lapse'][key] for r in result['records']],marker+'-',label=label)
    axes[0].set(title='Supporting field: mass-channel defects',xlabel='Integer flux q (bounding formulas)',ylabel='Dimensionless upper bound');axes[0].legend()
    atomic=result['records'][-1]
    for j in (0,1,3):
        axes[1].loglog([c['nominal_M_over_y'] for c in atomic['cases']],
            [c['masses'][j]['scaled_resolvent_error_upper'] for c in atomic['cases']],'o-',label=f'|m| = sqrt({j}) M')
    axes[1].set(title='Atomic reference: fixed-mass comparison',xlabel='Regulator scale M / frequency y',ylabel='Error bound / universal resolvent scale');axes[1].legend()
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('One-atom target: the proper mass now includes the varying lapse')
    fig.supxlabel('Auxiliary regulator fields, not the hydrogen payload. Finite M and y only. Supporting source and passage remain unknown.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
