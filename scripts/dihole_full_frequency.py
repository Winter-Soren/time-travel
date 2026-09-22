"""Full-field, fixed-regulator frequency bound with an integrable tail.

Local frequency coercivity preserves the uncut exterior. The first-order
source form avoids replacing a value/derivative cross term by two derivatives.
"""
import json
import math

import numpy as np
from mpmath import mp

from dihole_buffer_smoothing import full_bound, weighted_graph
from dihole_channel_bounds import MPIntervalContext, upper
from dihole_slice_channel import sha
from dirac_green_mass_limit import benchmark
from project_paths import project_path


NAME='dihole_full_frequency'


def multiply(left, right):
    answer={}
    for p,a in left.items():
        for q,b in right.items():
            answer[p+q]=answer.get(p+q,0)+a*b
    return answer


def tail_polynomial(iv, geometry, comparison, constants, q, M_value, eta_cap):
    """All coefficients are positive outward floats reused as exact inputs."""
    point=lambda value:iv.mpf(upper(value))
    get=lambda key:iv.mpf(geometry[key])
    M=iv.mpf(M_value);nu=iv.mpf(constants['inverse_lapse_upper'])
    cw=iv.mpf(constants['weight_upper_D'])
    G=weighted_graph(iv,M,iv.mpf([0,eta_cap]),iv.mpf(0),nu,get('N_upper'),cw)
    cq=point(get('scalar_curvature_negative_upper_D2')/4+q*get('magnetic_unit_flux_upper_D2')+get('log_lapse_gradient_squared_upper_D2')/4)
    cq2=point(get('scalar_curvature_negative_upper_D2')/4+q*get('magnetic_unit_flux_upper_D2')+get('log_lapse_gradient_squared_upper_D2')/2)
    h=point(2*get('N_upper'));start=upper(h*iv.sqrt(cq2))
    g=point(3200/get('P_lower'))
    lap=point((iv.mpf(19200000)+3200*get('radial_laplacian_drift_upper_D'))/get('P_lower')**2)
    D=point(lap+2*g*(2*g+iv.sqrt(cq)))
    u0=point(h**3*(G[3]+3*g*G[2]+3*g*g*G[1]+g**3*G[0]))
    u1=point(h*h*(G[3]+2*g*G[2]+g*g*G[1]))
    u2=point(h*(G[3]+g*G[2]))
    H0={1:u2,2:point(2*g*u1),3:point(D*u0)}
    H1={0:G[3],1:point(2*g*u2),2:point(D*u1)}
    c=iv.mpf(comparison['scalar_eigenvalue_factor_lower']);L=iv.mpf(9)/100;r0=iv.mpf(1)/40
    S={1:point(4*L*r0*r0*h/c**iv.mpf('1.5')),
       2:point(4*L*r0*h*h/c),3:point(2*L*h**3/iv.sqrt(c))}
    fN=point(4/get('N_lower'));fP=point(4/get('P_lower'));P=get('P_upper')
    radial=point(P*(2*g+iv.sqrt(cq))+get('radial_connection_upper_D'))
    zero=multiply(H0,H0);cross=multiply(H0,H1)
    inside={p:fN*v for p,v in zero.items()}
    for p,v in cross.items():inside[p+1]=inside.get(p+1,0)+fP*P*v
    for p,v in zero.items():inside[p+1]=inside.get(p+1,0)+fP*radial*v
    coefficients={p:upper(v) for p,v in multiply(S,inside).items()}
    export=lambda poly:{str(p):float(v.mid) for p,v in poly.items()}
    return dict(minimum_cutoff_K_D=start,coefficient_upper_D_power=coefficients,
        Cq_upper_D2=float(cq.mid),Cq2_upper_D2=float(cq2.mid),inverse_gap_coefficient=float(h.mid),
        cutoff_gradient_upper_D=float(g.mid),cutoff_laplacian_upper_D2=float(lap.mid),
        first_order_coefficient_N=float(fN.mid),first_order_coefficient_P=float(fP.mid),
        radial_HS_conversion_upper_D=float(radial.mid),P_upper=geometry['P_upper'],
        local_graph_upper=[float(x.mid) for x in G],H0_polynomial=export(H0),H1_polynomial=export(H1),
        heat_polynomial=export(S),iteration_norm_coefficients=[float(x.mid) for x in (u0,u1,u2,G[3])],
        cutoff_supports_over_D=[[-.005,.055],[-.0075,.0575],[-.010,.060],[-.0125,.0625],[-.015,.065],[-.0175,.0675]],
        auxiliary_domain_over_D=[-.02,.07],auxiliary_inverse_shift=0.)


def tail_integral(iv, data, cutoff):
    if not math.isfinite(cutoff) or cutoff<data['minimum_cutoff_K_D']:
        raise ValueError('Full tail requires the proved strictly positive local frequency floor.')
    K=iv.mpf(cutoff)
    return upper(sum(iv.mpf(v)/((int(p)-1)*K**(int(p)-1))
                     for p,v in data['coefficient_upper_D_power'].items())/iv.pi)


def sharp_low(iv, geometry, comparison, constants, q, M, cutoff, eta_cap):
    selection=dict(q=q,fixed_regulator_M_D=M,selected_split=dict(split_K_D=cutoff))
    old=full_bound(iv,geometry,comparison,constants,selection,eta_cap)
    get=lambda key:iv.mpf(geometry[key])
    P=get('P_upper');g=iv.mpf(old['cutoff_gradient_upper_D']);cq=iv.mpf(old['Cq_upper_D2'])
    nu=iv.mpf(constants['inverse_lapse_upper']);radius=iv.mpf(M)/(2*nu)
    legs=[]
    for name,rho in (('real',iv.mpf(0)),('complex_disk',radius)):
        x0,xA=map(iv.mpf,old[name]['localized_HS_upper'])
        xs=P*(xA+(2*g+iv.sqrt(cq)+rho*nu)*x0)+get('radial_connection_upper_D')*x0
        legs.append((iv.mpf(upper(x0)),iv.mpf(upper(xs))))
    (x0,xs),(x0disk,xsdisk)=legs
    x0prime=x0disk/radius;xsprime=xsdisk/radius
    fN=4/get('N_lower');fP=4/get('P_lower');K=iv.mpf(cutoff)
    B0=fP*x0*xs
    slope=fN*x0*x0+2*fN*K*x0*x0prime+fP*(x0prime*xs+x0*xsprime)
    return dict(cutoff_K_D=cutoff,paired_low_upper_D=upper(B0/2+K*slope/iv.pi),
        zero_source_trace_upper_D=upper(B0),source_lipschitz_upper=upper(slope),
        scalar_HS_upper=upper(x0),radial_HS_upper_D=upper(xs),
        scalar_HS_derivative_upper_over_D=upper(x0prime),radial_HS_derivative_upper=upper(xsprime),
        scalar_disk_HS_upper=upper(x0disk),radial_disk_HS_upper_D=upper(xsdisk),
        old_combined_jet_low_upper_D=old['full_paired_low_integral_absolute_upper_D'],
        reused_buffer_estimate=old)


def complete_bound(iv, geometry, comparison, constants, reference, eta_cap):
    q=reference['q'];M=reference['regulator_M_D']
    polynomial=tail_polynomial(iv,geometry,comparison,constants,q,M,eta_cap)
    rows=[];best_low=None;best=None
    # ponytail: a finite split search supplies a bound, not a global optimum.
    for j in range(17):
        K=polynomial['minimum_cutoff_K_D']*2**(j/2)
        low=sharp_low(iv,geometry,comparison,constants,q,M,K,eta_cap)
        high=tail_integral(iv,polynomial,K)
        total=upper(iv.mpf(low['paired_low_upper_D'])+iv.mpf(high))
        row=dict(cutoff_K_D=K,low_upper_D=low['paired_low_upper_D'],high_upper_D=high,complete_upper_D=total)
        rows.append(row)
        if best is None or total<best['complete_upper_D']:best,best_low=row,low
    previous=sharp_low(iv,geometry,comparison,constants,q,M,reference['split_K_D'],eta_cap)
    return dict(q=q,regulator_M_D=M,eta_ratio_interval=[0.,eta_cap],tail=polynomial,
        selected=best,selected_low_detail=best_low,split_grid=rows,previous_split_low=previous,
        full_frequency_paired_absolute_upper_D=best['complete_upper_D'])


def insertion_controls():
    rows=[];I=np.eye(6);rng=np.random.default_rng(176)
    raw=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6));C=(raw+raw.conj().T)/2
    derivative=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    F=np.diag([1.,.4,.2]*2)
    norm=lambda x:float(np.linalg.norm(x,'fro'))
    trace_norm=lambda x:float(np.linalg.svd(x,compute_uv=False).sum())
    for largest in (10.,1000.):
        for kappa in (0.,.1,2.):
            _,beta,N,_,_,A,_,_=benchmark(largest,kappa)
            W=np.diag(np.sqrt(np.diag(N)));gen=1j*beta@np.linalg.inv(N);T=F@beta@np.linalg.inv(N)
            fN=np.linalg.norm(T,2);fP=np.linalg.norm(C,2)
            for M in (.7,2.):
                masses=M*np.sqrt(np.arange(1,4)+.01);alpha=np.sqrt(6)*M**3
                def values(shift):
                    z=kappa+shift;R=[np.linalg.inv(A+shift*gen+m*I) for m in masses]
                    X=alpha*np.linalg.multi_dot(R)@W;Y=derivative@X
                    B=beta@(-1j*z*X.conj().T@T@X+(X.conj().T@C@Y-Y.conj().T@C@X)/2)
                    return R,X,Y,B
                R,X,Y,B=values(0.)
                Xp=-alpha*sum(np.linalg.multi_dot(R[:j]+[R[j],gen]+R[j:]) for j in range(3))@W
                Yp=derivative@Xp
                Bp=beta@(-1j*X.conj().T@T@X-1j*kappa*(Xp.conj().T@T@X+X.conj().T@T@Xp)
                    +(Xp.conj().T@C@Y+X.conj().T@C@Yp-Yp.conj().T@C@X-Y.conj().T@C@Xp)/2)
                step=1e-5;finite=(values(step)[3]-values(-step)[3])/(2*step)
                ceiling=abs(kappa)*fN*norm(X)**2+fP*norm(X)*norm(Y)
                derivative_ceiling=fN*norm(X)**2+2*fN*abs(kappa)*norm(X)*norm(Xp)+fP*(norm(Xp)*norm(Y)+norm(X)*norm(Yp))
                rows.append(dict(largest_lapse=largest,kappa=kappa,M=M,
                    trace_bound_ratio=trace_norm(B)/ceiling,derivative_trace_bound_ratio=trace_norm(Bp)/derivative_ceiling,
                    derivative_relative_error=norm(finite-Bp)/max(1,norm(Bp))))
    return rows


def polynomial_controls(records):
    mp.dps=90;rows=[]
    for record in records:
        data=record['tail'];K=mp.mpf(record['selected']['cutoff_K_D'])
        poly=lambda key,t:mp.fsum(mp.mpf(v)*t**int(p) for p,v in data[key].items())
        def direct(t):
            S=poly('heat_polynomial',t);H0=poly('H0_polynomial',t);H1=poly('H1_polynomial',t)
            return S*(mp.mpf(data['first_order_coefficient_N'])*H0*H0
                +mp.mpf(data['first_order_coefficient_P'])*t*H0*(mp.mpf(data['P_upper'])*H1+mp.mpf(data['radial_HS_conversion_upper_D'])*H0))
        expanded=lambda t:poly('coefficient_upper_D_power',t)
        discrepancies=[float((expanded(t)-direct(t))/max(mp.mpf(1),direct(t))) for t in (1/K,1/(2*K),1/(10*K))]
        analytic=mp.fsum(mp.mpf(v)/((int(p)-1)*K**(int(p)-1)) for p,v in data['coefficient_upper_D_power'].items())/mp.pi
        numerical=mp.quad(lambda x:direct(x/K)*K/(x*x) if x else 0,[0,1])/mp.pi
        rows.append(dict(q=record['q'],eta_cap=record['eta_ratio_interval'][1],polynomial_margin_min=min(discrepancies),
                         integral_relative_difference=float(abs(analytic-numerical)/max(1,analytic)),
                         integrated_upper_covers=mp.mpf(record['selected']['high_upper_D'])>=analytic,
                         tail_integral_encloses_unexpanded=analytic>=numerical))
    return rows


def run():
    path=project_path('dihole_buffer_smoothing_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_buffer_smoothing.py':previous['source_sha256']})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing buffered-smoothing inputs required.')
    old=json.loads(project_path('dihole_source_completion_results.json').read_text())
    paired=json.loads(project_path('dihole_paired_low_frequency_results.json').read_text())
    constants={**paired['global_lapse'],**old['insertion_factor'],'weight_upper_D':old['global_weight']['global_commutator_upper_D']}
    geometry=previous['geometry'];comparison=previous['scalar_comparison'];iv=MPIntervalContext();iv.dps=60
    records=[complete_bound(iv,geometry,comparison,constants,row,row['eta_ratio_interval'][1]) for row in previous['records']]
    insertion=insertion_controls();polynomials=polynomial_controls(records)
    rejected=False
    try:tail_integral(iv,records[0]['tail'],0.)
    except ValueError:rejected=True
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        positive_frequency_guard=rejected,
        tail_starts_at_inverse_cubic=all(min(map(int,r['tail']['coefficient_upper_D_power']))==3 for r in records),
        tail_coefficients_nonnegative=all(v>0 for r in records for v in r['tail']['coefficient_upper_D_power'].values()),
        positive_local_frequency_floor=all(r['selected']['cutoff_K_D']>=r['tail']['minimum_cutoff_K_D']>0 for r in records),
        nested_cutoffs_fit=all(r['tail']['cutoff_supports_over_D'][-1][0]>-.02 and r['tail']['cutoff_supports_over_D'][-1][1]<.07 for r in records),
        tails_decrease=all(all(a['high_upper_D']>b['high_upper_D'] for a,b in zip(r['split_grid'],r['split_grid'][1:])) for r in records),
        first_order_trace=max(r['trace_bound_ratio'] for r in insertion)<=1+1e-12,
        differentiated_first_order_trace=max(r['derivative_trace_bound_ratio'] for r in insertion)<=1+1e-12,
        source_derivative=max(r['derivative_relative_error'] for r in insertion)<2e-7,
        polynomial_expansion=all(r['polynomial_margin_min']>=0 for r in polynomials),
        independent_tail_integral=max(r['integral_relative_difference'] for r in polynomials)<1e-13,
        outward_tail_integral=all(r['integrated_upper_covers'] and r['tail_integral_encloses_unexpanded'] for r in polynomials),
        sharp_low_improves_same_split=all(r['previous_split_low']['paired_low_upper_D']<r['previous_split_low']['old_combined_jet_low_upper_D'] for r in records),
        finite_full_frequency_bounds=all(math.isfinite(r['full_frequency_paired_absolute_upper_D']) and r['full_frequency_paired_absolute_upper_D']>0 for r in records),
        uniform_shift_interval=all(records[i+1]['full_frequency_paired_absolute_upper_D']>=records[i]['full_frequency_paired_absolute_upper_D'] for i in (0,2)))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=records,insertion_controls=insertion,polynomial_controls=polynomials,
        source_sha256=sha(__file__),input_sha256=inputs,renormalized_source=None,source_sign=None,atom_transit=None,
        scope=['All angular states, both collars and the full exterior at fixed common regulator; no Gaussian or transverse selection.',
               'Low frequencies are paired before taking absolute value. High frequencies are absolutely integrable on each sign.',
               'This finite-regulator bound is not a renormalized source value or sign, a material construction or atomic passage.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Full-frequency controls failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,tests=tests,atomic=[dict(eta=r['eta_ratio_interval'],selected=r['selected'],previous_low=r['previous_split_low']['paired_low_upper_D']) for r in records[-2:]]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-2];M=row['regulator_M_D'];data=row['split_grid']
    fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    x=[r['cutoff_K_D']/M for r in data]
    for key,label in (('low_upper_D','Paired low frequencies'),('high_upper_D','Absolute high-frequency tail'),('complete_upper_D','Complete frequency bound')):
        axes[0].loglog(x,[r[key] for r in data],label=label)
    axes[0].set(xlabel='Frequency split / fixed regulator M',ylabel='Upper bound times D',title='The complete frequency integral is bounded')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    previous=row['previous_split_low']
    axes[1].bar(['Previous combined\nvalue/derivative bound','First-order source\nretained explicitly'],[previous['old_combined_jet_low_upper_D'],previous['paired_low_upper_D']],color=['#888888','#2c7899'])
    axes[1].set(yscale='log',ylabel='Paired low-frequency ceiling times D',title='Same atomic parameters and frequency split')
    fig.supxlabel('Fixed-regulator upper bounds, not a physical energy prediction, stress sign, or atom-transmission result.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
