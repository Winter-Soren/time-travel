"""Retain ordered-time simplex factors in the unchanged physical transport.

Historical transport sources stay frozen. The common physical-operator and
defect helpers are reused; only the parallel-transport majorant changes.
"""
import copy
import json
import math

from mpmath import mp

from dihole_centered_gauge import saved_normal
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import apply_operator, assemble, physical_coefficients, unpack
from dihole_local_transport import compose, dirac_sequences, dprod, packed, pmul
from dihole_scaled_graph import record_bound
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_ordered_transport'


def ordered_majorant(K):
    """Coefficients of exp(sum_{j>=1} K_j z^j), through the supplied order."""
    U = [K[0]*0+1]
    for n in range(1,len(K)):
        U.append(sum(j*K[j]*U[n-j] for j in range(1,n+1))/n)
    return U


def determinant_sqrt(iv, delta_det, n):
    """Positive detJ in [1/16,81/16]; delta_det contains 16 times its jets."""
    binomial=iv.mpf(1);coefficients=[iv.mpf(9)/4]
    for k in range(1,n+1):
        binomial*= (iv.mpf(1)/2-(k-1))/k
        coefficients.append(abs(binomial)/4)
    return compose(coefficients,delta_det,n)


def transport_bounds(iv, prior, record, normal, real):
    q=record['q'];_,_,(frame,invframe,gauge,V,Vi)=dirac_sequences(iv,prior,q)
    R=normal['radius'];Z=normal['Z'];E=normal['exp'];Y=normal['inverse'];n=16
    spin=dprod(dprod(dprod(frame,invframe),invframe),frame[1:])
    omega=[768*x+abs(q)*y for x,y in zip(spin,gauge)]
    omega_taylor=[2**(k+1)*omega[k]/math.factorial(k) for k in range(n+1)]
    path=[0]+[R*Z[k] for k in range(1,n+1)]
    K=pmul(compose(omega_taylor,path,n),[2*R]+path[1:],n)
    # Real-path unitarity removes K0. Ordered integration retains 1/m!
    # for m insertions without commuting any connection matrices.
    U=ordered_majorant(K)
    J=[iv.mpf(3)/2]+[(k+1)*Z[k+1] for k in range(1,n+1)]
    determinant=[1]+[0]*n
    for _ in range(4):determinant=pmul(determinant,J,n)
    delta_det=[0]+[24*x*16 for x in determinant[1:]]
    inverse_root=compose([4*iv.mpf(math.comb(2*k,k))/4**k for k in range(n+1)],delta_det,n)
    vicomposed=compose([2**k*Vi[k]/math.factorial(k) for k in range(n+1)],path,n)
    u0w=[iv.sqrt(iv.mpf(real['volume_upper']))*x for x in pmul(pmul(vicomposed,inverse_root,n),U,n)]
    u0=compose(u0w,[0]+[Y[k]/R for k in range(1,n+1)],n)
    inverse0=iv.sqrt(iv.mpf(real['volume_upper'])/iv.mpf(real['volume_lower']))*iv.mpf(9)/4
    u0inv=[inverse0]
    for k in range(1,n+1):u0inv.append(inverse0*sum(u0[j]*u0inv[k-j] for j in range(1,k+1)))
    inverse_legacy=u0inv
    # u0^-1=V(x)/V(y)*sqrt(det Dexp)*U*. Adjoint derivatives have the
    # same unitary bounds; do not invert their majorant as a geometric series.
    vcomposed=compose([2**k*V[k]/math.factorial(k) for k in range(n+1)],path,n)
    inverse_w=[x/iv.sqrt(iv.mpf(real['volume_lower'])) for x in
        pmul(pmul(vcomposed,determinant_sqrt(iv,delta_det,n),n),U,n)]
    inverse_direct=compose(inverse_w,[0]+[Y[k]/R for k in range(1,n+1)],n)
    u0inv=[a if a.b<=b.b else b for a,b in zip(inverse_legacy,inverse_direct)]
    physical=physical_coefficients(iv,prior,record)
    X=compose(E,Y,n);X[0]=0
    coefficients=[u0]
    for j in range(1,5):
        order=n-2*j
        applied=apply_operator(physical,coefficients[-1],order)
        integrand=compose(pmul(u0inv,applied,order),X,order)
        coefficients.append([x/j for x in pmul(u0,integrand,order)])
    return dict(q=q,
        convention='Physical-volume heat coefficients: ||D_x^n uj||/n! in coordinate Euclidean Frechet norm, source fixed. Same gauge and physical operator as184.',
        u_normalized_frechet_upper=[[packed(iv,x) for x in row] for row in coefficients],
        derivative_orders=[len(row)-1 for row in coefficients],
        u0_inverse_normalized_frechet_upper=[packed(iv,x) for x in u0inv],
        u0_inverse_geometric_comparison_upper=[packed(iv,x) for x in inverse_legacy],
        u0_inverse_direct_upper=[packed(iv,x) for x in inverse_direct],
        parallel_transport_scaled_normalized_frechet_upper=[packed(iv,x) for x in U],
        uniform_generator_scaled_normalized_frechet_upper=[packed(iv,x) for x in K],
        unitary_transport_used=True,ordered_simplex_factors_used=True,heat_defect_B4=None)


def controls(iv):
    mp.dps=160
    K=[iv.mpf(0),iv.mpf(7)]+[iv.mpf(0)]*15
    linear=ordered_majorant(K)
    linear_error=max(abs(mp.mpf(x._mpi_[1])/(mp.mpf(7)**n/math.factorial(n))-1) for n,x in enumerate(linear))
    # Independent finite sum of powers/m!, rather than the recurrence.
    k=[mp.mpf(0),mp.mpf(2),mp.mpf(3),mp.mpf(5)]+[mp.mpf(0)]*9
    recurrence=ordered_majorant(k);power=[mp.mpf(1)]+[mp.mpf(0)]*12;series=power.copy()
    for m in range(1,13):
        power=pmul(power,k,12)
        series=[a+b/math.factorial(m) for a,b in zip(series,power)]
    series_error=max(abs(a-b)/max(1,abs(b)) for a,b in zip(recurrence,series))
    # Independent polynomial-time solution of U'=(w Ax+t*w^2 Az)U.
    # Ax and Az are anti-Hermitian and do not commute. Integrate each
    # Taylor coefficient exactly in time via polynomial coefficients.
    Ax=1j*mp.matrix([[0,1],[1,0]]);Az=1j*mp.diag([1,-1])
    solutions=[{0:mp.eye(2)}];ratios=[]
    majorant=ordered_majorant([mp.mpf(0),mp.mpf(1),mp.mpf(1)]+[mp.mpf(0)]*8)
    for n in range(1,11):
        forcing={}
        for power,value in solutions[n-1].items():forcing[power]=Ax*value
        if n>=2:
            for power,value in solutions[n-2].items():
                forcing[power+1]=forcing.get(power+1,mp.zeros(2))+Az*value
        solution={p+1:a/(p+1) for p,a in forcing.items()};solutions.append(solution)
        endpoint=sum(solution.values(),mp.zeros(2))
        ratios.append(float(max(mp.svd(endpoint,compute_uv=False))/majorant[n]))
    # Known inverse heat-leading factor: density exp(w/7), positive
    # determinant (1+w/10)(1+w/20), and unitary exp(7iw sigma_x).
    center=mp.mpf('.2');order=10
    det=lambda w:(1+w/10)*(1+w/20)
    delta=[iv.mpf(0)]+[16*iv.mpf(str(abs(mp.diff(det,center,k))/math.factorial(k))) for k in range(1,order+1)]
    sqrt_bound=determinant_sqrt(iv,delta,order)
    density=[iv.exp(iv.mpf(1)/35)/(iv.mpf(7)**k*math.factorial(k)) for k in range(order+1)]
    phase=ordered_majorant([iv.mpf(0),iv.mpf(7)]+[iv.mpf(0)]*(order-1))
    inverse_bound=pmul(pmul(density,sqrt_bound,order),phase,order)
    def exact_inverse(w):
        return mp.exp(w/7)*mp.sqrt(det(w))*mp.matrix([[mp.cos(7*w),-1j*mp.sin(7*w)],[-1j*mp.sin(7*w),mp.cos(7*w)]])
    inverse_ratios=[]
    for n in range(order+1):
        derivative=mp.matrix([[mp.diff(lambda w:exact_inverse(w)[i,j],center,n)/math.factorial(n) for j in range(2)] for i in range(2)])
        inverse_ratios.append(float(max(mp.svd(derivative,compute_uv=False))/mp.mpf(inverse_bound[n]._mpi_[0])))
    return dict(constant_linear_relative_error=float(linear_error),
        exponential_series_relative_error=float(series_error),
        noncommuting_generator_norm=float(mp.norm(Ax*Az-Az*Ax)),
        noncommuting_time_ordered_max_ratio=max(ratios),noncommuting_orders=10,
        direct_inverse_max_ratio=max(inverse_ratios),direct_inverse_orders=10,
        pure_constant_generator_gives_no_parameter_derivatives=all(x==0 for x in ordered_majorant([iv.mpf(100)]+[iv.mpf(0)]*16)[1:]))


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_positive_frame_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_positive_frame.py':previous['source_sha256']}
    trial_path=project_path('dihole_ordered_transport_simplex_trial.json')
    trial_source=project_path('provenance/dihole_ordered_transport_simplex_trial.py')
    trial=json.loads(trial_path.read_text())
    if not trial['passed'] or trial['source_sha256']!=sha(trial_source):
        raise ArithmeticError('The simplex-only pilot must remain intact.')
    for name in (trial_path.name,'provenance/dihole_ordered_transport_simplex_trial.py',
                 'dihole_ordered_transport_simplex_trial_progress.png','dihole_ordered_transport_simplex_trial_progress.pdf'):
        inputs[name]=sha(project_path(name))
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing184 inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    local=json.loads(project_path('dihole_local_transport_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    heat=json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    updated=copy.deepcopy(prior);fields=updated['jets']['total_order_derivative_upper']
    for key in ('coframe','inverse_coframe'):fields[key]=previous['frame_bounds'][key]
    fields['gauge_per_abs_q']=previous['gauge']['gauge_per_abs_q_derivative_upper']
    normal=saved_normal(iv,local);records=[];parallel_ok=[];heat_ok=[]
    for old,mr,hr in zip(previous['records'],magnetic['records'],heat['records']):
        if not old['q']==mr['q']==hr['q']:raise ArithmeticError('Flux records do not align.')
        operator=old['operator'];point=old['point_evaluation']
        transport=transport_bounds(iv,updated,operator,normal,local['real_geometry'])
        for a,b in zip(transport['parallel_transport_scaled_normalized_frechet_upper'],old['heat_transport']['parallel_transport_scaled_normalized_frechet_upper']):
            parallel_ok.append(bool(unpack(iv,a).b<=unpack(iv,b).b))
        for row,oldrow in zip(transport['u_normalized_frechet_upper'],old['heat_transport']['u_normalized_frechet_upper']):
            heat_ok.extend(bool(unpack(iv,a).b<=unpack(iv,b).b) for a,b in zip(row,oldrow))
        defect=assemble(iv,updated,local,operator,transport,point,hr,mr)
        scaled=record_bound(iv,updated,operator,defect,mr)
        baseline=old['scaled_bound']['heat_defect_B4_uniform_upper']
        ratio=iv.log(unpack(iv,baseline)/unpack(iv,scaled['heat_defect_B4_uniform_upper']),10)
        records.append(dict(q=old['q'],operator=operator,point_evaluation=point,heat_transport=transport,
            heat_defect=defect,scaled_bound=scaled,previous_B4_upper=baseline,
            ordered_ceiling_log10_change_interval=[lower(ratio),upper(ratio)]))
    control=controls(iv)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        constant_linear_control=control['constant_linear_relative_error']<1e-95,
        exponential_series_control=control['exponential_series_relative_error']<1e-150,
        noncommuting_control=control['noncommuting_generator_norm']>0,
        ordered_time_control=control['noncommuting_time_ordered_max_ratio']<=1,
        constant_generator_removed=control['pure_constant_generator_gives_no_parameter_derivatives'],
        direct_inverse_control=control['direct_inverse_max_ratio']<1,
        simplex_trial_preserved=trial['passed'] and trial['source_sha256']==sha(trial_source),
        parallel_bounds_do_not_increase=all(parallel_ok),heat_bounds_do_not_increase=all(heat_ok),
        derivative_schedule=all(r['heat_transport']['derivative_orders']==[16,14,12,10,8] for r in records),
        operator_and_point_unchanged=all(r['operator']==o['operator'] and r['point_evaluation']==o['point_evaluation'] for r,o in zip(records,previous['records'])),
        full_bounds_improve=all(r['ordered_ceiling_log10_change_interval'][0]>0 for r in records),
        physical_source_not_claimed=all(r['scaled_bound']['physical_source'] is None and r['scaled_bound']['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,bound_encoding=previous['bound_encoding'],
        scope=['Ordered-time m-fold integrals retain their 1/m! factors. The inverse leading coefficient retains the adjoint unitary transport explicitly.',
               'The finite-jet argument does not require analyticity or commuting matrices; zero-parameter-order evolution is absorbed unitarily.',
               'Physical geometry, field, gauge, operator, point bound, regulator, normal domain and one-atom target remain unchanged.',
               'Full downstream heat coefficients and defects are recomputed. These ceilings do not establish source, stability, atom transit or clock gain.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Ordered transport checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,records=[dict(q=r['q'],
        B4_log10=r['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper'],
        improvement=r['ordered_ceiling_log10_change_interval'],
        Cstar_log10=r['scaled_bound']['exact_low_order_remainder_Cstar_upper']['log10_upper']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-1]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['Previous transport estimate','Ordered transport and inverse'],
        [row['previous_B4_upper']['log10_upper'],row['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper']],
        color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of B4 ceiling, D = 1',title='Retaining the order of quantum-field propagation')
    fig.supxlabel('Same one-atom target and physical field. Conservative ceilings, not actual errors.\nSupport, stability, atom passage and clock gain remain unproved.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
