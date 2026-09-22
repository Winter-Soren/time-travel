"""Remove the proper infrared shift in a compact, fixed-frequency insertion.

Zero shift means a localized resolvent limit, not a global inverse or a
whole-space positive-filter trace. The physical regulator M stays fixed.
"""
import json
import math

import numpy as np
from mpmath import mp
from scipy.integrate import quad

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_proper_channel import graph_resolvent
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_local_zero_shift'


def components(iv,c,d,g,eta,M):
    if eta.a<0 or M.a<=0 or g.a<=0:raise ValueError('Require eta>=0, M>0 and a positive auxiliary frequency gap.')
    masses=[iv.sqrt(eta+j*M*M) for j in range(4)];scale=iv.sqrt(6)*M**3;residues=[]
    for j,m in enumerate(masses):
        value=scale
        for k,n in enumerate(masses):
            if k!=j:value*=(n+m)/((k-j)*M*M)
        residues.append(value)
    C=[c+d/g if m.b==0 else graph_resolvent(iv,c,d,g,m) for m in masses]
    r=scale
    for m in masses:r/=iv.sqrt(g*g+m*m)
    actions=[]
    for omit in range(4):
        value=scale
        for j,m in enumerate(masses):
            if j!=omit:value/=iv.sqrt(g*g+m*m)
        actions.append(upper(value))
    # The singular mass is applied last. The final A-graph estimate never
    # divides by m0, so its endpoint and whole-interval envelope stay finite.
    error=iv.mpf(0);reduced=iv.mpf(1)
    for j in (1,2,3):
        error=(error+C[j]*reduced)/masses[j]
        reduced/=iv.sqrt(g*g+masses[j]**2)
    graph=scale*(error+C[0]*reduced)
    return dict(masses=masses,residues=residues,C=C,r=r,a=iv.mpf(min(actions)),graph=graph)


def local_bound(iv,c,d,g,eta,M,local_inverse,ca,da,hs2,uniform=False):
    high=components(iv,c,d,g,eta,M)
    low=components(iv,c,d,g,iv.mpf(0),M) if uniform else high
    weights=[abs(x) for x in high['residues']];C=low['C'];masses=high['masses']
    summed=sum(w*z for w,z in zip(weights,C));e0=local_inverse*summed
    cancelled=c*low['a']+d*low['r']+local_inverse*sum(w*m*z for w,m,z in zip(weights,masses,C))
    e1=iv.mpf(min(upper(summed),upper(cancelled),upper(low['graph'])))
    L=(1+c)*low['a']+d*low['r'];r=low['r']
    insertion_error=hs2*(e1*(ca*L+da*r)+(L+e1)*(ca*e1+da*e0))
    absolute=hs2*(L+e1)*(ca*(L+e1)+da*(r+e0))
    return dict(local_filter_error_upper_over_D=upper(e0),local_first_graph_error_upper=upper(e1),
        first_graph_triangle_upper=upper(summed),first_graph_cancellation_upper=upper(cancelled),
        first_graph_product_upper=upper(low['graph']),reduced_filter_norm_upper_over_D=upper(r),
        reduced_action_norm_upper=upper(low['a']),approximate_full_action_norm_upper=upper(L),
        finite_local_metric_insertion_error_upper=upper(insertion_error),
        finite_local_metric_insertion_absolute_upper=upper(absolute),
        residue_intervals=[[lower(x),upper(x)] for x in high['residues']],
        graph_resolvent_upper=[upper(x) for x in C])


def frequency_tail(iv,c,d,N,M,ca,da,hs2,cutoffs,eta_cap=0):
    """Positive inverse-frequency polynomials; fixed M and frozen filter."""
    endpoint=components(iv,c,d,iv.mpf(1),iv.mpf(eta_cap),M)
    W0=sum(abs(x) for x in endpoint['residues'])
    W1=sum(abs(x)*m for x,m in zip(endpoint['residues'],endpoint['masses']))
    A3=iv.sqrt(6)*M**3*N**3;A4=A3*N;zero=iv.mpf(0)
    e0=[zero,N*W0*c,N*N*W0*d,zero,zero]
    e1=[zero,N*W1*c,N*N*W1*d,c*A3,d*A4]
    trial=[zero,zero,zero,(1+c)*A3,d*A4];r=[zero,zero,zero,zero,A4]
    u1=[a+b for a,b in zip(trial,e1)];u0=[a+b for a,b in zip(r,e0)]
    absolute={};error={}
    for p in range(2,9):
        pairs=[(j,p-j) for j in range(5) if 0<=p-j<5]
        absolute[p]=iv.mpf(upper(hs2*sum(u1[j]*(ca*u1[k]+da*u0[k]) for j,k in pairs)))
        error[p]=iv.mpf(upper(hs2*sum(e1[j]*(ca*trial[k]+da*r[k])+u1[j]*(ca*e1[k]+da*e0[k]) for j,k in pairs)))
    rows=[]
    for value in cutoffs:
        K=iv.mpf(value)
        if K.a<=0:raise ValueError('The high-frequency tail requires a strictly positive lower cutoff.')
        total=sum(coefficient/((p-1)*K**(p-1)) for p,coefficient in absolute.items())
        defect=sum(coefficient/((p-1)*K**(p-1)) for p,coefficient in error.items())
        rows.append(dict(cutoff_K_D=value,positive_frequency_absolute_integral_upper_D=upper(total),
            two_sided_frequency_measure_absolute_upper_D=upper(total/iv.pi),
            two_sided_frequency_measure_comparison_error_upper_D=upper(defect/iv.pi)))
    return dict(eta_cap_D2=float(iv.mpf(eta_cap).mid),absolute_coefficient_upper_D_power={p:float(x.mid) for p,x in absolute.items()},
        comparison_error_coefficient_upper_D_power={p:float(x.mid) for p,x in error.items()},
        coefficient_units='Coefficient of |kappa|^-p has units D^-p; the insertion bound is dimensionless.',
        tail_bounds=rows)


def bounds(iv,previous,lapse):
    c=iv.mpf(previous['proper_graph_a_upper']);d=iv.mpf(previous['proper_graph_zero_upper_D'])
    ca=iv.mpf(previous['proper_metric_graph_a_upper']);da=iv.mpf(previous['proper_metric_graph_zero_upper_D'])
    hs2=iv.mpf(previous['filter_HS_squared_upper']);N=iv.mpf(lapse['lapse_upper'])
    M_value=previous['kappa_D'];M=iv.mpf(M_value);eta_max=.01*M_value*M_value;frequencies=[]
    for fraction in (.1,1.,10.):
        kappa_value=fraction*M_value;kappa=iv.mpf(kappa_value);g_value=lower(kappa/N)
        g=iv.mpf(g_value);local=iv.mpf(upper(N/kappa));samples=[]
        for eta_ratio in (0.,1e-8,1e-4,.01):
            eta_value=eta_ratio*M_value*M_value
            row=local_bound(iv,c,d,g,iv.mpf(eta_value),M,local,ca,da,hs2)
            samples.append(dict(eta_D2=eta_value,nominal_eta_over_M_squared=eta_ratio,**row))
        envelope=local_bound(iv,c,d,g,iv.mpf(eta_max),M,local,ca,da,hs2,uniform=True)
        frequencies.append(dict(kappa_D=kappa_value,nominal_kappa_over_fixed_M=fraction,
            auxiliary_gap_lower_D=g_value,compact_inverse_upper_over_D=upper(local),
            eta_max_D2=eta_max,uniform_eta_interval=envelope,samples=samples))
    return dict(q=previous['q'],fixed_regulator_M_D=M_value,fixed_filter_Lambda_D=M_value,
        filter_HS_squared_upper=previous['filter_HS_squared_upper'],frequencies=frequencies,
        zero_shift_frequency_tail=frequency_tail(iv,c,d,N,M,ca,da,hs2,[f*M_value for f in (.1,1.,10.,100.)]),
        uniform_eta_frequency_tail=frequency_tail(iv,c,d,N,M,ca,da,hs2,[f*M_value for f in (.1,1.,10.,100.)],eta_max))


def residue_controls():
    mp.dps=90;rows=[]
    for eta,M in ((0,1),(.01,.7),(1e12,1)):
        e,v=mp.mpf(eta),mp.mpf(M);m=[mp.sqrt(e+j*v*v) for j in range(4)];scale=mp.sqrt(6)*v**3
        stable=[scale*mp.fprod((m[k]+m[j])/((k-j)*v*v) for k in range(4) if k!=j) for j in range(4)]
        direct=[scale/mp.fprod(m[k]-m[j] for k in range(4) if k!=j) for j in range(4)]
        z=mp.mpc('.3','.8');product=scale/mp.fprod(z+x for x in m)
        fraction=sum(w/(z+x) for w,x in zip(stable,m))
        moments=[sum(w*x**k for w,x in zip(stable,m))+(scale if k==3 else 0) for k in range(4)]
        rows.append(dict(eta=eta,M=M,coefficient_relative_error=float(max(abs(a-b)/max(1,abs(b)) for a,b in zip(stable,direct))),
            rational_relative_error=float(abs(fraction-product)/abs(product)),scaled_moment_error=float(max(abs(x)/max(1,scale) for x in moments))))
    return rows


def matrix_controls(iv):
    H=np.diag([-.7,.1,.8,1.3,2.])+np.diag([.5,.4,.3,.2],1)+np.diag([.5,.4,.3,.2],-1)
    N=np.diag([.4,.9,1.5,10.,100.]);U=np.array([[1.],[.6],[0.],[0.],[0.]])/math.sqrt(1.36)
    n=float((U.T@N@U)[0,0]);t1=np.array([[0.,1.],[1.,0.]]);t3=np.diag([1.,-1.]);I=np.eye(2)
    H0=np.kron(t3,H);beta=np.kron(t1,np.eye(5));B=np.kron(I,np.diag(1/np.sqrt(np.diag(N))))
    V=np.kron(I,np.sqrt(N)@U/math.sqrt(n));J=.8*V;X=np.diag([.8,1.]);S=np.kron(I,np.diag([1.,1.,1.,0.,0.]))
    F=np.kron(I,np.diag([.2,.5,0.,0.,0.]));size=lambda x:float(np.linalg.norm(x,2));rows=[]
    for kappa in (.2,1.,3.):
        A=B@beta@(1j*kappa*np.eye(10)+H0)@B;a=V.conj().T@A@V;T=A@J-J@a
        d=np.nextafter(size(T),np.inf);g=np.nextafter(kappa/n,0);local=np.nextafter(1.5/kappa,np.inf)
        dotA=B@beta@(F@H0+H0@F)@B-(F@A+A@F)/2;da=np.nextafter(size(dotA),np.inf)
        fixed=list(map(iv.mpf,(0.,d,g)));tail=list(map(iv.mpf,(1.,local,0.,da,np.linalg.norm(X,'fro')**2)))
        envelope=local_bound(iv,*fixed,iv.mpf(.01),*tail,uniform=True)
        for eta in (0.,1e-8,1e-4,.01):
            b=local_bound(iv,*fixed,iv.mpf(eta),*tail)
            masses=[math.sqrt(eta+j) for j in range(4)];alpha=math.sqrt(6)
            residues=[alpha*math.prod((masses[k]+masses[j])/(k-j) for k in range(4) if k!=j) for j in range(4)]
            FA=alpha*np.eye(10,dtype=complex);Fa=alpha*np.eye(2,dtype=complex);partial=np.zeros_like(FA)
            graphterm=np.zeros((10,2),dtype=complex);greens=[]
            for m,w in zip(masses,residues):
                G=np.linalg.inv(A+m*np.eye(10));small=np.linalg.inv(a+m*np.eye(2))
                FA=G@FA;Fa=small@Fa;partial+=w*G;graphterm+=w*m*G@T@small
                Hm=H0+m*np.kron(t1,N)
                identity=np.kron(I,np.sqrt(N))@np.linalg.inv(Hm+1j*kappa*np.eye(10))@beta@np.kron(I,np.sqrt(N))
                greens.append(size(S@(G-identity)@S))
            E=FA@J-J@Fa;actual0=size(S@E);actual1=size(S@A@E)
            u=FA@J@X;v=J@Fa@X;source=lambda z:-float(np.real(np.trace((A@z).conj().T@dotA@z)))
            error=abs(source(u)-source(v))
            rows.append(dict(kappa=kappa,eta=eta,defect_support_error=size((np.eye(10)-S)@T),
                metric_support_error=size(dotA-S@dotA@S),compact_green_identity_error=max(greens),
                partial_fraction_relative_error=size(partial-FA)/max(1.,size(FA)),
                graph_cancellation_error=size(A@E+T@Fa-graphterm),local_output_error=actual0,local_graph_error=actual1,
                metric_insertion_difference=error,metric_insertion_value=source(u),
                endpoint_bounds_hold=bool(actual0<=b['local_filter_error_upper_over_D'] and actual1<=b['local_first_graph_error_upper'] and error<=b['finite_local_metric_insertion_error_upper']),
                uniform_bounds_hold=bool(actual0<=envelope['local_filter_error_upper_over_D'] and actual1<=envelope['local_first_graph_error_upper'] and error<=envelope['finite_local_metric_insertion_error_upper'])))
    return rows


def uniform_controls(iv):
    c,d,g,M,local,ca,da,hs2=map(iv.mpf,(.2,.7,1.3,.9,.4,2.,3.,1.7));cap=iv.mpf(.02)
    envelope=local_bound(iv,c,d,g,cap,M,local,ca,da,hs2,uniform=True);ratios=[]
    for eta in np.linspace(0,.02,51):
        b=local_bound(iv,c,d,g,iv.mpf(float(eta)),M,local,ca,da,hs2)
        ratios.append(max(b[key]/envelope[key] for key in ('local_filter_error_upper_over_D','local_first_graph_error_upper','finite_local_metric_insertion_error_upper')))
    rejected=[]
    for eta,mass,gap in ((-1.,1.,1.),(0.,0.,1.),(0.,1.,0.)):
        try:components(iv,c,d,iv.mpf(gap),iv.mpf(eta),iv.mpf(mass))
        except ValueError:rejected.append(True)
        else:rejected.append(False)
    return dict(sampled_uniform_bound_ratios=ratios,invalid_parameters_rejected=rejected,
        unbounded_lapse_example=[dict(dimension=n,global_zero_shift_inverse_norm=float(n),fixed_four_site_inverse_norm=4.) for n in (10,100,1000)])


def tail_controls(iv):
    parameters=list(map(iv.mpf,(.2,.7,1.1,.9,2.,3.,1.7)))
    polynomial=frequency_tail(iv,*parameters,[.5,2.,10.]);errors=[]
    coefficients=polynomial['absolute_coefficient_upper_D_power']
    for row in polynomial['tail_bounds']:
        K=row['cutoff_K_D'];scale=row['positive_frequency_absolute_integral_upper_D']
        # Rescale the independent integral to a unit lower limit.
        integral=quad(lambda z:sum(a/(K*z)**int(p) for p,a in coefficients.items())*K/scale,1,np.inf,epsabs=1e-12)[0]
        errors.append(abs(integral-1))
    toy=[]
    for K in (.2,1.,10.):
        M=.7;hs2=1.3
        bound=frequency_tail(iv,*map(iv.mpf,(0.,0.,1.,M,1.,0.,hs2)),[K])['tail_bounds'][0]['positive_frequency_absolute_integral_upper_D']
        # A=i*kappa*sigma1, dotA=-A, J=I: exact selected insertion.
        actual=quad(lambda k:hs2*6*M**6/((k*k+M*M)*(k*k+2*M*M)*(k*k+3*M*M)),K,np.inf,epsabs=1e-12)[0]
        toy.append(dict(K=K,exact_model_integral=actual,polynomial_bound=bound,enclosed=bool(actual<=bound)))
    return dict(polynomial_integral_relative_errors=errors,commuting_model_tails=toy)


def run():
    path=project_path('dihole_proper_channel_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_proper_channel.py':previous['source_sha256']})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing prior inputs required.')
    lapse=json.loads(project_path('dihole_massive_channel_results.json').read_text())
    iv=MPIntervalContext();iv.dps=60
    records=[bounds(iv,p,m['lapse']) for p,m in zip(previous['records'],lapse['records'],strict=True)]
    residues=residue_controls();matrices=matrix_controls(iv);uniform=uniform_controls(iv);tails=tail_controls(iv)
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        stable_residues=max(r['coefficient_relative_error'] for r in residues)<1e-60,
        rational_residue_identity=max(r['rational_relative_error'] for r in residues)<1e-40,
        regulator_moments=max(r['scaled_moment_error'] for r in residues)<1e-40,
        compact_defect_and_metric=max(max(r['defect_support_error'],r['metric_support_error']) for r in matrices)<1e-12,
        localized_green_identity=max(r['compact_green_identity_error'] for r in matrices)<1e-12,
        matrix_common_product=max(r['partial_fraction_relative_error'] for r in matrices)<1e-12,
        graph_cancellation=max(r['graph_cancellation_error'] for r in matrices)<1e-12,
        pointwise_local_enclosures=all(r['endpoint_bounds_hold'] for r in matrices),
        whole_interval_local_enclosures=all(r['uniform_bounds_hold'] for r in matrices) and max(uniform['sampled_uniform_bound_ratios'])<=1+1e-14,
        zero_shift_included=sum(r['eta']==0 for r in matrices)==3,
        invalid_parameters_rejected=all(uniform['invalid_parameters_rejected']),
        finite_zero_shift_insertions=all(math.isfinite(f['samples'][0]['finite_local_metric_insertion_absolute_upper']) for r in records for f in r['frequencies']),
        integrated_tail_polynomial=max(tails['polynomial_integral_relative_errors'])<1e-10,
        exact_commuting_tail_enclosed=all(x['enclosed'] for x in tails['commuting_model_tails']),
        high_frequency_tails_decrease=all(all(a['two_sided_frequency_measure_absolute_upper_D']>b['two_sided_frequency_measure_absolute_upper_D'] for a,b in zip(r['zero_shift_frequency_tail']['tail_bounds'],r['zero_shift_frequency_tail']['tail_bounds'][1:])) for r in records),
        uniform_frequency_tail_covers_endpoint=all(all(e['two_sided_frequency_measure_absolute_upper_D']<=u['two_sided_frequency_measure_absolute_upper_D'] for e,u in zip(r['zero_shift_frequency_tail']['tail_bounds'],r['uniform_eta_frequency_tail']['tail_bounds'])) for r in records),
        fixed_regulator_across_frequencies=all(all(abs(f['kappa_D']/r['fixed_regulator_M_D']-f['nominal_kappa_over_fixed_M'])<1e-14 for f in r['frequencies']) for r in records))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=records,residue_controls=residues,matrix_controls=matrices,uniform_controls=uniform,tail_controls=tails,
        source_sha256=sha(__file__),input_sha256=inputs,physical_source=None,source_sign=None,atom_transit=None,
        global_zero_shift_positive_filter_trace=None,quantitative_eta_convergence_rate=None,
        scope=['Compact output and defect supports lie in the certified collars; no global upper lapse bound or full proper frequency gap.',
            'At eta=0 the massless resolvent denotes the localized strong/HS limit from the common Hamiltonian core, not a globally bounded A inverse.',
            'The whole eta interval is bounded by proven monotonic coefficient/graph envelopes. Sampled controls do not replace that argument.',
            'The m0-last product controls the first-action error without dividing by m0. It bounds the same commuting product as the localized partial-fraction error.',
            'Regulator M and Gaussian scale are fixed across kappa. A positive inverse-frequency polynomial bounds the absolute two-sided high-frequency tail with measure d kappa/(2 pi), uniformly over the declared eta interval. This permits eta removal under the integral for |kappa|>=K>0; no small-kappa control, omitted selections or regulator removal.',
            'Only a frozen finite proper metric insertion; no stable source, obtainable material, atomic propagation, clock gain or new law established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Local zero-shift checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'],tests=tests,atomic_endpoint=records[-1]['frequencies'][1]['samples'][0]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    atomic=result['records'][-1];central=atomic['frequencies'][1];rows=central['samples']
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    axes[0].plot([r['nominal_eta_over_M_squared'] for r in rows],[central['kappa_D']*r['local_filter_error_upper_over_D'] for r in rows],'o-',color='#187c80')
    axes[0].set_xscale('symlog',linthresh=1e-8)
    axes[0].set(title='The localized bound remains finite at zero',xlabel='Artificial shift eta / fixed M^2',ylabel='Frequency times local filter error bound')
    tail=atomic['zero_shift_frequency_tail']['tail_bounds'];M=atomic['fixed_regulator_M_D']
    axes[1].loglog([r['cutoff_K_D']/M for r in tail],
        [r['two_sided_frequency_measure_absolute_upper_D'] for r in tail],'o-',color='#7840aa')
    axes[1].set(title='Absolute high-frequency tail at zero shift',xlabel='Tail starts at |kappa| / fixed M',ylabel='Tail bound times D, including d kappa / (2 pi)')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Atomic reference: remove the artificial shift from the compact insertion')
    fig.supxlabel('Fixed regulator and spatial filter. The low-frequency contribution, full supporting source and atom passage remain unresolved.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
