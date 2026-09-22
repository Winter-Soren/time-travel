"""Lapse-weighted collar channel and a common proper-P regulator filter.

Only first-order domains are used. The finite proper metric insertion is
bounded with a frozen selection; its infrared shift remains positive.
"""
import itertools
import json
import math

import numpy as np

from balanced_dihole import symbolic_fields
from dihole_aligned_channel import geometry, norm, trial
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_collar_resolvent import cutoff
from dihole_local_source import local_first
from dihole_massive_channel import COEFFICIENTS
from dihole_slice_channel import sha
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME='dihole_proper_channel'


def graph_resolvent(iv,c,d,g,m):
    """Maximize c||aG u||+d||G u|| on its spectral norm circle."""
    if min(c.a,d.a,g.a)<0 or m.a<=0:raise ValueError('Nonnegative graph data and positive mass required.')
    interior=iv.sqrt(c*c+(d/m)**2);boundary=(c*g+d)/iv.sqrt(g*g+m*m)
    condition=c*m*m-d*g
    if condition.a>=0:return interior
    if condition.b<=0:return boundary
    return iv.mpf(max(upper(interior),upper(boundary)))


def product_bound(iv,c,d,g,eta,M):
    if eta.a<=0 or M.a<=0:raise ValueError('Positive proper infrared shift and regulator required.')
    masses=[iv.sqrt(eta+j*M*M) for j in range(4)];scale=iv.sqrt(6)*M**3
    C=[graph_resolvent(iv,c,d,g,m) for m in masses]
    orders=[]
    for order in itertools.permutations(range(4)):
        error=iv.mpf(0);reduced=iv.mpf(1)
        for j in order:
            error=(error+C[j]*reduced)/masses[j]
            reduced/=iv.sqrt(g*g+masses[j]**2)
        orders.append(dict(order=list(order),error_upper=upper(scale*error),
            graph_error_upper=upper(scale*error*masses[order[-1]])))
    best=min(orders,key=lambda r:r['error_upper']);worst=max(orders,key=lambda r:r['error_upper'])
    reduced=scale;universal=scale
    for m in masses:
        reduced/=iv.sqrt(g*g+m*m);universal/=m
    action=[]
    for omit in range(4):
        value=scale
        for j,m in enumerate(masses):
            if j!=omit:value/=iv.sqrt(g*g+m*m)
        action.append(upper(value))
    return dict(best=best,worst=worst,best_graph=min(orders,key=lambda r:r['graph_error_upper']),
        reduced_action_norm_upper=min(action),
        reduced_norm_upper=upper(reduced),universal_norm_upper=upper(universal),
        graph_resolvent_bounds=[upper(x) for x in C],
        mass_intervals=[[lower(m),upper(m)] for m in masses])


def bounds(iv,base,massive,constants):
    lapse=massive['lapse'];nmin=iv.mpf(lapse['lapse_lower']);nmax=iv.mpf(lapse['lapse_upper'])
    RN=iv.mpf(lapse['mass_defect_upper']);old_a=iv.mpf(base['graph_h_coefficient_upper'])
    c_value=upper(old_a*iv.sqrt(nmax/nmin)+RN/nmin)
    d_value=upper(iv.mpf(base['graph_zero_coefficient_upper_D'])/nmin)
    kappa=massive['y_D'];g_value=lower(iv.mpf(kappa)/nmax)
    c,d,g=map(iv.mpf,(c_value,d_value,g_value));hs2=iv.mpf(massive['four_component_filter_HS_squared_upper'])
    w=iv.mpf(constants['W_lower'][0]);W=iv.mpf(constants['W_upper'][1]);joint=iv.sqrt(1+W/w)
    ca_value=upper(4*iv.sqrt(nmax/nmin)*joint)
    da_value=upper((1200*joint+600+4*iv.sqrt(400+31*iv.mpf(base['q'])))/(w*nmin))
    ca,da=map(iv.mpf,(ca_value,da_value))
    cases=[]
    for fraction in (.0001,.01,1.):
        eta_value=fraction*kappa*kappa;eta=iv.mpf(eta_value)
        for ratio in (.1,1.,10.):
            M_value=ratio*kappa;M=iv.mpf(M_value);b=product_bound(iv,c,d,g,eta,M)
            error=iv.mpf(b['best']['error_upper']);reduced=iv.mpf(b['reduced_norm_upper'])
            # J is contractive, not isometric. Compare ||T_A J X|| and ||J T_a X||.
            trace=error*(2*reduced+error)*hs2
            graph_error=iv.mpf(b['best_graph']['graph_error_upper'])
            approximate_action=(1+c)*iv.mpf(b['reduced_action_norm_upper'])+d*reduced
            metric_error=hs2*(graph_error*(ca*approximate_action+da*reduced)
                +(approximate_action+graph_error)*(ca*graph_error+da*error))
            source_upper=hs2*(approximate_action+graph_error)*(ca*(approximate_action+graph_error)+da*(reduced+error))
            scale=iv.sqrt(6)*M**3
            for j in range(4):scale/=iv.sqrt(eta+j*M*M)
            cases.append(dict(eta_D2=eta_value,regulator_M_D=M_value,nominal_eta_over_kappa_squared=fraction,
                nominal_M_over_kappa=ratio,best_order=b['best']['order'],worst_order=b['worst']['order'],
                common_filter_error_upper_over_D=b['best']['error_upper'],worst_order_error_upper_over_D=b['worst']['error_upper'],
                first_graph_filter_error_upper=b['best_graph']['graph_error_upper'],best_graph_order=b['best_graph']['order'],
                scaled_common_filter_error_upper=upper(iv.mpf(kappa)*error),
                error_over_universal_filter_scale_upper=upper(error/scale),
                reduced_filter_norm_upper_over_D=b['reduced_norm_upper'],universal_filter_norm_upper_over_D=b['universal_norm_upper'],
                finite_positive_trace_error_upper_over_D2=upper(trace),
                reduced_filter_action_norm_upper=b['reduced_action_norm_upper'],
                approximate_full_action_norm_upper=upper(approximate_action),
                finite_proper_metric_insertion_error_upper=upper(metric_error),
                finite_proper_metric_insertion_absolute_upper=upper(source_upper),
                graph_resolvent_bounds=b['graph_resolvent_bounds'],absolute_mass_D_intervals=b['mass_intervals']))
    return dict(q=base['q'],kappa_D=kappa,proper_graph_a_upper=c_value,proper_graph_zero_upper_D=d_value,
        reduced_frequency_gap_lower_D=g_value,filter_HS_squared_upper=float(hs2.mid),
        proper_metric_graph_a_upper=ca_value,proper_metric_graph_zero_upper_D=da_value,cases=cases)


def field_controls():
    function=symbolic_fields()[0];parameters=exact_target_parameters();rows=[];sz=np.diag([1.,-1.])
    for q in (1005,10005):
        kappa=math.sqrt(q)
        for u,nu in ((.0125,.85),(.025,.9),(.0375,.95)):
            mu=math.floor(q*nu)+.5;data=geometry(function,parameters,u,mu/q,1024,80);state=trial(data,q)
            logs,radial,_=local_first(function,parameters,data['radius'],data['theta'])
            N=np.exp(logs[0]);Ns=N*radial[0]*data['radius_rate'];U=state['chi'];p=state['probability'];n=p@N
            inner=lambda a,b:data['weights']@np.sum(a.conj()*b,axis=0)
            rN=(N-n)*U;ns=p@Ns+2*inner(state['chi_s'],rN)
            V=np.sqrt(N/n)*U;v,vs=state['velocity'],state['velocity_s'];chi,chit=cutoff(u)
            for k in (0.,1.,10.):
                derivative=1j*k-ns/(2*n)
                direct=(1j*kappa*U-1j*(state['radial_action']+state['angular_action']+derivative*(sz@U)/data['W']))/np.sqrt(N*n)
                a=(1j*kappa-1j*(v*derivative+vs/2))/n
                R=direct-V*a
                predicted=(-1j*(state['r1']*derivative+state['r0'])-rN*a)/np.sqrt(N*n)
                T0=-1j*(chi*(state['r1']*derivative+state['r0'])+chit*(sz@U)/data['W'])/np.sqrt(n)
                T=T0/np.sqrt(N)-chi*rN*a/np.sqrt(N*n)
                direct_cut=chi*direct-1j*chit*(sz@U)/(data['W']*np.sqrt(N*n))-chi*V*a
                F,Fs=data['mask'][:2];rad=state['radial_action']+derivative*(sz@U)/data['W']
                dotA=(-1j*kappa*F*U-1j*F*rad-1j*Fs*(sz@U)/(2*data['W']))/np.sqrt(N*n)
                def varied(epsilon):
                    return (np.exp(-epsilon*F)*1j*kappa*U
                        -1j*np.exp(epsilon*F)*(rad+epsilon*Fs*(sz@U)/(2*data['W']))
                        -1j*state['angular_action'])/np.sqrt(N*n)
                step=1e-5;metric_error=norm(data,(varied(step)-varied(-step))/(2*step)-dotA)/max(1.,norm(data,dotA))
                size=norm(data,direct)**2
                rows.append(dict(q=q,u=u,mu=mu,k_D=k,isometry_error=abs(inner(V,V)-1),
                    compression_error=abs(inner(V,direct)-a),proper_defect_identity_error=norm(data,R-predicted),
                    cutoff_defect_identity_error=norm(data,T-direct_cut),
                    cutoff_omission_norm=norm(data,T-chi*R),
                    proper_metric_derivative_relative_error=metric_error,
                    positive_form_relative_error=abs(size-abs(a)**2-norm(data,R)**2)/max(1.,size)))
    return [{k:float(v) if isinstance(v,(np.floating,np.complexfloating)) else v for k,v in r.items()} for r in rows]


def matrix_controls():
    rng=np.random.default_rng(171);raw=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
    H=(raw+raw.conj().T)/2;U=np.linalg.qr(raw)[0][:,:2];N=np.diag(np.linspace(.4,1.8,5))
    n=U.conj().T@N@U;nv,nz=np.linalg.eigh(n);root=(nz*np.sqrt(nv))@nz.conj().T;inverse=(nz/nv**.5)@nz.conj().T
    t1=np.array([[0.,1.],[1.,0.]]);t3=np.diag([1.,-1.]);I=np.eye(2)
    U4=np.kron(I,U);V=np.kron(I,np.sqrt(N)@U@inverse);B=np.kron(I,np.diag(1/np.sqrt(np.diag(N))))
    H0=np.kron(t3,H);beta=np.kron(t1,np.eye(5));h0=U4.conj().T@H0@U4
    n4=np.kron(I,n);root4=np.kron(I,root);inverse4=np.kron(I,inverse);N4=np.kron(I,N)
    C=np.diag([.3,1.,.3,1.]);J=V@C;X=np.diag([.7,1.,.8,.9]);size=lambda a:float(np.linalg.norm(a,2))
    iv=MPIntervalContext();iv.dps=60;rows=[]
    for kappa in (.2,1.,3.):
        A=B@beta@(1j*kappa*np.eye(10)+H0)@B;a=inverse4@np.kron(t1,I)@(1j*kappa*np.eye(4)+h0)@inverse4
        R=A@V-V@a;R0=H0@U4-U4@h0;rN=N4@U4-U4@n4
        predicted=B@beta@(R0-rN@np.linalg.solve(n4,1j*kappa*np.eye(4)+h0))@inverse4
        T=A@J-J@a;gap=kappa/max(nv);actualgap=min(abs(np.linalg.eigvalsh(1j*a)))
        P=A.conj().T@A;p=a.conj().T@a;form=V.conj().T@P@V
        for eta,M in ((.01,.7),(.1,1.)):
            Fa=np.sqrt(6)*M**3*np.eye(4,dtype=complex);FA=np.sqrt(6)*M**3*np.eye(10,dtype=complex)
            rational=np.zeros((10,10),dtype=complex)
            for j,cj in enumerate(COEFFICIENTS):
                m=math.sqrt(eta+j*M*M)
                FA=np.linalg.solve(A+m*np.eye(10),FA);Fa=np.linalg.solve(a+m*np.eye(4),Fa)
                rational+=cj*np.linalg.inv(P+(eta+j*M*M)*np.eye(10))
            bound=product_bound(iv,iv.mpf(0),iv.mpf(np.nextafter(size(T),np.inf)),iv.mpf(np.nextafter(gap,0)),iv.mpf(eta),iv.mpf(M))
            error=bound['best']['error_upper'];actual=size(FA@J-J@Fa)
            graph_actual=size(A@(FA@J-J@Fa));graph_bound=bound['best_graph']['graph_error_upper']
            fullnorm=np.linalg.norm(FA@J@X,'fro')**2;approx=np.linalg.norm(J@Fa@X,'fro')**2
            tracebound=error*(2*bound['reduced_norm_upper']+error)*np.linalg.norm(X,'fro')**2
            F=np.kron(I,np.diag(np.linspace(-.3,.5,5)));dH=F@H0+H0@F
            dotA=B@beta@dH@B-(F@A+A@F)/2;da=size(dotA)
            u=FA@J@X;v=J@Fa@X;source=lambda z:-float(np.real(np.trace((A@z).conj().T@dotA@z)))
            Lv=bound['reduced_action_norm_upper']+size(T)*bound['reduced_norm_upper']
            metricbound=np.linalg.norm(X,'fro')**2*da*(graph_bound*bound['reduced_norm_upper']+(Lv+graph_bound)*error)
            def varied_source(epsilon):
                exp=np.diag(np.exp(epsilon*np.diag(F)));be=np.diag(np.exp(-epsilon*np.diag(F)/2))@B
                ae=be@beta@(1j*kappa*np.eye(10)+exp@H0@exp)@be
                return -.5*np.linalg.norm(ae@u,'fro')**2
            step=1e-5;metric_fd=(varied_source(step)-varied_source(-step))/(2*step)
            rows.append(dict(kappa=kappa,eta=eta,M=M,skew_error=size(A+A.conj().T),isometry_error=size(V.conj().T@V-np.eye(4)),
                weighted_compression_error=size(V.conj().T@A@V-a),weighted_defect_error=size(R-predicted),
                positive_form_error=size(form-p-R.conj().T@R),square_shortcut_mismatch=size(form-p),
                cutoff_form_error=size(J.conj().T@P@J-(a@C).conj().T@(a@C)-(R@C).conj().T@(R@C)),
                cutoff_isometry_mismatch=size(J.conj().T@J-np.eye(4)),reduced_gap_margin=float(actualgap-gap),
                common_filter_identity_error=size(FA.conj().T@FA-rational),actual_filter_error=actual,filter_error_bound=error,
                actual_first_graph_error=graph_actual,first_graph_error_bound=graph_bound,
                finite_trace_difference=float(abs(fullnorm-approx)),finite_trace_error_bound=float(tracebound),
                finite_metric_insertion_difference=abs(source(u)-source(v)),finite_metric_insertion_error_bound=float(metricbound),
                frozen_metric_derivative_relative_error=abs(metric_fd-source(u))/max(1.,abs(source(u))),
                metric_insertion_enclosed=bool(abs(source(u)-source(v))<=metricbound),
                filter_enclosed=bool(actual<=error),first_graph_enclosed=bool(graph_actual<=graph_bound),trace_enclosed=bool(abs(fullnorm-approx)<=tracebound)))
    return rows


def graph_controls(iv):
    rows=[]
    for c,d,g,m in ((.1,2.,1.,.2),(.1,2.,1.,10.),(1.,.1,3.,1.),(0.,2.,1.,.5),(1.,0.,0.,2.)):
        bound=upper(graph_resolvent(iv,*map(iv.mpf,(c,d,g,m))))
        rmax=1/math.hypot(g,m);r=np.linspace(0,rmax,10001)
        values=c*np.sqrt(np.maximum(0.,1-m*m*r*r))+d*r
        rows.append(dict(c=c,d=d,g=g,m=m,sampled_max=float(max(values)),bound=bound,enclosed=bool(max(values)<=bound*(1+1e-14))))
    return rows


def run():
    path=project_path('dihole_massive_channel_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_massive_channel.py':previous['source_sha256']})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing prior inputs required.')
    base=json.loads(project_path('dihole_collar_resolvent_results.json').read_text());iv=MPIntervalContext();iv.dps=60
    constants=json.loads(project_path('dihole_channel_bounds_results.json').read_text())['uniform_bound']['constants']
    records=[bounds(iv,b,m,constants) for b,m in zip(base['records'],previous['records'],strict=True)]
    fields=field_controls();matrices=matrix_controls();graphs=graph_controls(iv)
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        field_isometry=max(r['isometry_error'] for r in fields)<1e-12,
        field_compression=max(r['compression_error'] for r in fields)<1e-10,
        field_proper_defect=max(r['proper_defect_identity_error'] for r in fields)<1e-10,
        field_cutoff_defect=max(r['cutoff_defect_identity_error'] for r in fields)<1e-10 and max(r['cutoff_omission_norm'] for r in fields)>1,
        field_positive_form=max(r['positive_form_relative_error'] for r in fields)<1e-12,
        actual_proper_metric_derivative=max(r['proper_metric_derivative_relative_error'] for r in fields)<1e-8,
        matrix_skew_isometry=max(max(r['skew_error'],r['isometry_error']) for r in matrices)<1e-12,
        matrix_weighted_channel=max(max(r['weighted_compression_error'],r['weighted_defect_error']) for r in matrices)<1e-12,
        matrix_positive_form=max(max(r['positive_form_error'],r['cutoff_form_error']) for r in matrices)<1e-12,
        rejects_compression_shortcuts=min(r['square_shortcut_mismatch'] for r in matrices)>1 and min(r['cutoff_isometry_mismatch'] for r in matrices)>.5,
        auxiliary_frequency_gap=min(r['reduced_gap_margin'] for r in matrices)>-1e-12,
        common_filter_identity=max(r['common_filter_identity_error'] for r in matrices)<1e-10,
        product_error_bounds=all(r['filter_enclosed'] and r['first_graph_enclosed'] and r['trace_enclosed'] for r in matrices),
        frozen_proper_metric_insertion=all(r['metric_insertion_enclosed'] for r in matrices) and max(r['frozen_metric_derivative_relative_error'] for r in matrices)<1e-8,
        spectral_graph_bounds=all(r['enclosed'] for r in graphs),
        finite_exported_bounds=all(math.isfinite(c['finite_positive_trace_error_upper_over_D2']) and c['finite_positive_trace_error_upper_over_D2']>0 for r in records for c in r['cases']))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=records,field_controls=fields,matrix_controls=matrices,graph_controls=graphs,
        source_sha256=sha(__file__),input_sha256=inputs,physical_source=None,source_sign=None,atom_transit=None,
        scope=['The lapse-weighted V is isometric on the uncut collar. J=chi V is contractive; cutoff commutators and J in the output norm remain explicit.',
            'The auxiliary complete-line n,v continuation stays positive and bounded with complete optical coordinate. Its frequency gap is not a global full-operator gap.',
            'The common first-order product squared equals the proper-P rational filter. The approximation J T_a X is not the squared compressed P form or a bare reduced trace.',
            'Only fixed nonzero kappa, finite regulator M and positive eta; no frequency integration, eta removal, regulator removal or omitted insertion bound.',
            'The finite selected proper metric insertion is -Re<Au,dotA u> with u=T_A JX. Its input JX is frozen; varying the channel or filter would be a different derivative.',
            'Its absolute value and comparison error are bounded, not measured. No full renormalized source, new material, stability, atom propagation or clock difference is established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Proper-channel controls failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'],tests=tests,atomic_record=records[-1]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained');atomic=result['records'][-1]
    for fraction in (.0001,.01,1.):
        rows=[r for r in atomic['cases'] if r['nominal_eta_over_kappa_squared']==fraction]
        axes[0].semilogx([r['nominal_M_over_kappa'] for r in rows],[r['error_over_universal_filter_scale_upper'] for r in rows],'o-',label=f'eta / kappa^2 = {fraction:g}')
    axes[0].set(title='Atomic reference: whole regulator product',xlabel='Regulator M / frequency kappa',ylabel='Error bound / universal filter norm scale');axes[0].legend()
    rows=[r['cases'][4] for r in result['records']]
    axes[1].loglog([r['q'] for r in result['records']],[r['finite_proper_metric_insertion_error_upper'] for r in rows],'o-',color='#187c80')
    axes[1].set(title='Finite proper metric insertion: error',xlabel='Integer flux q (bounding formulas)',ylabel='Dimensionless error bound at one frequency')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('The actual collar channel now keeps the common regulator cancellation')
    fig.supxlabel('Right: M = kappa, eta = 0.01 kappa^2. Fixed positive shift and frozen selection. Full supporting energy and atom transit remain unknown.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
