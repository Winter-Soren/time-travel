"""Finite nested-buffer heat bounds, retaining the unchanged global operator.

The terminal bound is the existing global filter. All recurrences stop at a
finite depth within its plateau; no exterior boundary condition is imposed.
"""
import json

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext
from dihole_heat_defect import unpack
from dihole_heat_source_graph import calculate
from dihole_local_transport import packed
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_nested_heat_bound'
DEPTHS = (0,1,2,4,6,8)


def smaller(left,right):
    """Select a proved ceiling; never recover an error by subtraction."""
    return left if left.b <= right.b else right


def bootstrap(iv,rows,start,depth):
    if not isinstance(start,int) or not isinstance(depth,int) or start < 0 or depth < 0:
        raise ValueError('Nonnegative integer padding and finite depth required.')
    if depth >= len(rows) or (start+depth+1)*iv.mpf(2)**-32 >= iv.mpf(1)/400:
        raise ValueError('Finite rows and filter-plateau coverage required.')
    gamma = lambda n:iv.sqrt(1+40*n*iv.mpf(2)**-32)
    bound = [unpack(iv,x) for x in rows[depth]['direct_heat_HS_upper']]
    if depth == 0:
        error = [gamma(start)*unpack(iv,rows[0][k]) for k in
                 ('value_HS_error_upper','Dirac_HS_error_upper')]
        approx = [gamma(start)*unpack(iv,rows[0][k]) for k in
                  ('frozen_value_HS_upper','frozen_Dirac_HS_upper')]
        return error,[smaller(x,q+f) for x,q,f in zip(bound,approx,error)]
    for j in range(depth-1,-1,-1):
        row=rows[j]; g=gamma(start+j)
        # This bounds the INTEGRAL of defect norms, not a triangle heat error.
        eta=unpack(iv,row['early_defect_operator_integral_upper'])
        late=[unpack(iv,row['HS_error_terms'][k]) for k in ('value_late','Dirac_late')]
        previous=[g*unpack(iv,row[k]) for k in ('value_HS_error_upper','Dirac_HS_error_upper')]
        error=[smaller(g*l+eta*b,old) for l,b,old in zip(late,bound,previous)]
        approx=[g*unpack(iv,row[k]) for k in ('frozen_value_HS_upper','frozen_Dirac_HS_upper')]
        bound=[smaller(unpack(iv,d),q+f) for d,q,f in zip(row['direct_heat_HS_upper'],approx,error)]
    return error,bound


def comparison(iv,rows,T,depth):
    core,X=bootstrap(iv,rows,0,depth)
    buffer,_=bootstrap(iv,rows,2,depth)
    a,b=(unpack(iv,rows[0]['geometry'][k]) for k in ('insertion_graph_a','insertion_graph_b'))
    Q1=iv.sqrt(1+80*iv.mpf(2)**-32)*unpack(iv,rows[0]['frozen_Dirac_HS_upper'])
    terms=[2*T*buffer[1]*(a*X[1]+b*X[0]),
           2*T*Q1*(a*buffer[1]+b*buffer[0]),
           2*T*core[1]*(a*core[1]+b*core[0])]
    norm=iv.mpf(rows[0]['q'])**2
    return dict(depth=depth,core_error_upper=[packed(iv,v) for v in core],
        buffered_error_upper=[packed(iv,v) for v in buffer],
        core_heat_upper=[packed(iv,v) for v in X],
        signed_reference_error_terms_per_q2_upper=[packed(iv,v/norm) for v in terms],
        signed_reference_error_per_q2_upper=packed(iv,sum(terms)/norm),
        absolute_heat_response_per_q2_upper=packed(iv,2*T*X[1]*(a*X[1]+b*X[0])/norm),
        normalized_Dirac_error_to_frozen_ceiling=packed(iv,core[1]/unpack(iv,rows[0]['frozen_Dirac_HS_upper'])))


def controls():
    """Independent exterior-coupled matrix semigroup, both derivative orders."""
    mp.dps=80
    P=mp.diag([3,4,5,6,7])
    for j in range(4):P[j,j+1]=P[j+1,j]=mp.mpf('.4')
    A=mp.diag([P[j,j] for j in range(5)])
    lam,U=mp.eigsy(P); root=U*mp.diag([mp.sqrt(v) for v in lam])*U.T
    heat=lambda t:U*mp.diag([mp.exp(-t*v) for v in lam])*U.T
    projection=lambda n:mp.diag([int(j<min(n+1,5)) for j in range(5)])
    delta=mp.norm(P-A)
    ratios=[]; support=[]; leakage=[]
    for t in map(mp.mpf,('.01','.1','.3')):
        for depth in (1,2,4):
            for order in (0,1):
                Da=mp.eye(5) if order==0 else root
                bound=mp.norm(Da*heat(t/2**depth))
                for j in range(depth-1,-1,-1):
                    time=t/2**j;Pi=projection(j)
                    Q=mp.diag([mp.exp(-time*A[k,k]) for k in range(5)])*Pi
                    defect=(P-A)*Q
                    support.append(mp.norm((mp.eye(5)-projection(j+1))*defect))
                    late=delta*(time/2 if order==0 else mp.sqrt(time/mp.e))
                    error=late+time*delta*bound/2
                    actual_error=mp.norm(Da*(heat(time)*Pi-Q))
                    ratios.append(actual_error/error)
                    bound=min(mp.norm(Da*heat(time)),mp.norm(Da*Q)+error)
                    ratios.append(mp.norm(Da*heat(time)*Pi)/bound)
                leakage.append(mp.norm((mp.eye(5)-projection(0))*heat(t)*projection(0)))
    # Independent positive recurrence expansion checks the finite terminal term.
    eta=[mp.mpf('.2'),mp.mpf('.1'),mp.mpf('.05')]
    local=list(map(mp.mpf,(2,3,4)));terminal=mp.mpf(9)
    recur=terminal
    for a,b in reversed(list(zip(local,eta))):recur=a+b*recur
    unrolled=local[0]+eta[0]*local[1]+eta[0]*eta[1]*local[2]+mp.fprod(eta)*terminal
    return dict(matrix_bound_max_ratio=str(max(ratios)),defect_support_error=str(max(support)),
        exterior_leakage_min=str(min(leakage)),unrolled_terminal_error=str(abs(recur-unrolled)))


def run():
    iv=MPIntervalContext();iv.dps=120
    path=project_path('dihole_signed_kernel_results.json');prior191=json.loads(path.read_text())
    inputs={**prior191['input_sha256'],path.name:sha(path),'dihole_signed_kernel.py':prior191['source_sha256']}
    if not prior191['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing191 inputs required.')
    names=('dihole_frozen_magnetic','dihole_regular_jets','dihole_positive_frame',
        'dihole_local_transport','dihole_buffer_smoothing','dihole_regulator_heat_tail','dihole_full_frequency')
    frozen,prior,frame,local,buffer,heat,frequency=(json.loads(project_path(n+'_results.json').read_text()) for n in names)
    records=[]
    for row,hr in zip(frozen['records'],heat['records']):
        if row['q']!=hr['q']:raise ArithmeticError('Flux order mismatch.')
        fr=next(x for x in frequency['records'] if x['q']==row['q'] and x['eta_ratio_interval']==[0.,0.])
        steps=[calculate(iv,frozen,row,hr,fr,prior,frame,local,buffer['geometry'],iv.mpf(2)**-j)
               for j in range(4+max(DEPTHS)+1)]
        for exponent in (0,1,2,4):
            fraction=iv.mpf(2)**-exponent;T=iv.mpf(hr['cutoff_s0_over_D_squared'])*fraction
            trials=[comparison(iv,steps[exponent:],T,depth) for depth in DEPTHS]
            best=min(trials,key=lambda r:unpack(iv,r['signed_reference_error_per_q2_upper']).b)
            old=unpack(iv,trials[0]['signed_reference_error_per_q2_upper'])
            new=unpack(iv,best['signed_reference_error_per_q2_upper'])
            records.append(dict(q=row['q'],heat_time_fraction=float(fraction.a),trials=trials,
                selected_depth=best['depth'],selected_reference_error_per_q2_upper=best['signed_reference_error_per_q2_upper'],
                reference_error_ceiling_improvement_factor=packed(iv,old/new),
                strict_improvement=bool(new.b<old.a),signed_full_heat_interval=None,
                renormalized_source=None,atom_transit=None))
    control=controls()
    invalid=0
    for start,depth in ((-1,0),(0,-1),(0,99),(0,1.5),(2**24,0)):
        try:bootstrap(iv,steps,start,depth)
        except ValueError:invalid+=1
    atomic=[r for r in records if r['q']==prior191['q']]
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        nested_matrix_bounds=mp.mpf(control['matrix_bound_max_ratio'])<=1,
        exact_defect_support=mp.mpf(control['defect_support_error'])==0,
        exterior_propagation_retained=mp.mpf(control['exterior_leakage_min'])>0,
        finite_terminal_recurrence=mp.mpf(control['unrolled_terminal_error'])<mp.mpf('1e-70'),
        unsupported_nesting_rejected=invalid==5,
        source_output_plateau=(3+max(DEPTHS))*iv.mpf(2)**-32<iv.mpf(1)/400,
        atomic_strict_improvement=all(r['strict_improvement'] for r in atomic),
        no_quadrature_or_source_claim=all(r['signed_full_heat_interval'] is None
            and r['renormalized_source'] is None and r['atom_transit'] is None for r in records))
    # Convert interval comparisons to JSON booleans; retain exact dyadic ceilings.
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        maximum_depth=max(DEPTHS),maximum_source_output_padding_over_D=packed(iv,(3+max(DEPTHS))*iv.mpf(2)**-32),
        prior_numerical_reference_per_q2=prior191['records'][-1]['control_variate_reference_per_q2'],
        numerical_reference_error_certified=False,source_sha256=sha(__file__),input_sha256=inputs,
        scope=['Finite nesting at half heat time; original global-filter bound terminates every chain.',
            'All exterior inputs and the unchanged full semigroup are retained; source projections are never differentiated.',
            'Error radii concern the exact191 spatial-cutoff signed reference, not its un-enclosed quadrature.',
            'The time parameter is Euclidean heat time, not transit time or a clock shift.',
            'No renormalized supporting source, stable opening, atom transport or new physical law.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Nested heat controls failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,atomic=[dict(fraction=r['heat_time_fraction'],
        depth=r['selected_depth'],radius_log=r['selected_reference_error_per_q2_upper']['log10_upper']) for r in atomic]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=max((r for r in result['records'] if r['heat_time_fraction']==1),key=lambda r:r['q'])
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3),layout='constrained')
    trials=row['trials']
    axes[0].plot([r['depth'] for r in trials],
        [r['signed_reference_error_per_q2_upper']['log10_upper'] for r in trials],'o-',color='#247886')
    axes[0].set(xlabel='Finite nesting depth',ylabel='log10 of normalized error ceiling',
        title='The coarse smoothing bound is substantially reduced')
    reference=float(mp.log10(abs(mp.mpf(result['prior_numerical_reference_per_q2']))))
    axes[1].bar(['Numerical reference','Remaining error ceiling'],
        [reference,row['selected_reference_error_per_q2_upper']['log10_upper']],color=['#247886','#b05740'])
    axes[1].axhline(0,color='#888888',linewidth=.7)
    axes[1].set(ylabel='log10 of normalized magnitude',title='The full response sign is still unresolved')
    fig.supxlabel('Same physical field and one-atom target. These are upper bounds, not observed errors.\nNo supported opening, atomic passage or operational time shift is established.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
