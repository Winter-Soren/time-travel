"""Retain the exact radial-gauge cancellation along the same geodesics.

This changes a generator bound, not the magnetic field or its gauge.
Historical sources remain frozen; shared coefficient/defect helpers are reused.
"""
import copy
import json
import math

from mpmath import mp

from dihole_centered_gauge import saved_normal
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import apply_operator, assemble, physical_coefficients, unpack
from dihole_local_transport import compose, dirac_sequences, dprod, packed, pmul
from dihole_ordered_transport import determinant_sqrt, ordered_majorant
from dihole_scaled_graph import record_bound
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_radial_cancellation'


def radial_magnetic(a, velocity, acceleration, n):
    old=pmul(a,velocity,n)
    new=[x/2 for x in pmul(a,acceleration,n)]
    return old,new,[x if x.b<=y.b else y for x,y in zip(old,new)]


def generator_bounds(iv, prior, q, normal, local):
    _,_,(frame,invframe,gauge,_,_)=dirac_sequences(iv,prior,q)
    R=normal['radius'];Z=normal['Z'];n=16
    path=[0]+[R*Z[k] for k in range(1,n+1)]
    velocity=[2*R]+path[1:]
    gamma=local['normal_geometry']['christoffel_component_upper']
    C=[2**(k+3)*unpack(iv,gamma[k])/math.factorial(k) for k in range(n+1)]
    acceleration=pmul(compose(C,path,n),pmul(velocity,velocity,n),n)
    spin=dprod(dprod(dprod(frame,invframe),invframe),frame[1:])
    spin_taylor=[2**(k+1)*768*spin[k]/math.factorial(k) for k in range(n+1)]
    magnetic_taylor=[2**(k+1)*abs(q)*gauge[k]/math.factorial(k) for k in range(n+1)]
    spin_generator=pmul(compose(spin_taylor,path,n),velocity,n)
    a=compose(magnetic_taylor,path,n)
    old,new,selected=radial_magnetic(a,velocity,acceleration,n)
    combined=[s+a for s,a in zip(spin_generator,selected)]
    return combined,dict(
        spin_generator_upper=[packed(iv,x) for x in spin_generator],
        magnetic_generator_original_upper=[packed(iv,x) for x in old],
        magnetic_generator_radial_upper=[packed(iv,x) for x in new],
        magnetic_generator_selected_upper=[packed(iv,x) for x in selected],
        acceleration_scaled_normalized_frechet_upper=[packed(iv,x) for x in acceleration],
        magnetic_zero_ratio_interval=[lower(new[0]/old[0]),upper(new[0]/old[0])],
        expected_zero_ratio_interval=[lower(C[0]*R),upper(C[0]*R)],
        improved_magnetic_orders=[k for k,(a,b) in enumerate(zip(old,new)) if b.b<a.a])


def transport_from_generator(iv, prior, record, normal, real, K):
    """Same185 heat construction, taking the explicitly bounded generator."""
    _,_,(_,_,_,V,Vi)=dirac_sequences(iv,prior,record['q'])
    R=normal['radius'];Z=normal['Z'];E=normal['exp'];Y=normal['inverse'];n=16
    path=[0]+[R*Z[k] for k in range(1,n+1)]
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
    inverse_legacy=[inverse0]
    for k in range(1,n+1):inverse_legacy.append(inverse0*sum(u0[j]*inverse_legacy[k-j] for j in range(1,k+1)))
    vcomposed=compose([2**k*V[k]/math.factorial(k) for k in range(n+1)],path,n)
    inverse_w=[x/iv.sqrt(iv.mpf(real['volume_lower'])) for x in
        pmul(pmul(vcomposed,determinant_sqrt(iv,delta_det,n),n),U,n)]
    inverse_direct=compose(inverse_w,[0]+[Y[k]/R for k in range(1,n+1)],n)
    u0inv=[a if a.b<=b.b else b for a,b in zip(inverse_legacy,inverse_direct)]
    physical=physical_coefficients(iv,prior,record)
    X=compose(E,Y,n);X[0]=0;coefficients=[u0]
    for j in range(1,5):
        order=n-2*j;applied=apply_operator(physical,coefficients[-1],order)
        integrand=compose(pmul(u0inv,applied,order),X,order)
        coefficients.append([x/j for x in pmul(u0,integrand,order)])
    return dict(q=record['q'],convention='Same physical-volume heat coefficients and fixed-source normalized Frechet derivative convention as185.',
        u_normalized_frechet_upper=[[packed(iv,x) for x in row] for row in coefficients],
        derivative_orders=[len(row)-1 for row in coefficients],
        u0_inverse_normalized_frechet_upper=[packed(iv,x) for x in u0inv],
        u0_inverse_geometric_comparison_upper=[packed(iv,x) for x in inverse_legacy],
        u0_inverse_direct_upper=[packed(iv,x) for x in inverse_direct],
        parallel_transport_scaled_normalized_frechet_upper=[packed(iv,x) for x in U],
        uniform_generator_scaled_normalized_frechet_upper=[packed(iv,x) for x in K],
        unitary_transport_used=True,ordered_simplex_factors_used=True,radial_cancellation_used=True,heat_defect_B4=None)


def controls(iv):
    mp.dps=110;B=mp.mpf(5)
    # A curved geodesic in the flat metric pulled back by (x,y)->(x,y-x^2).
    # Its sole nonzero Christoffel symbol is Gamma^y_xx=-2.
    path=lambda t,w:mp.matrix([w*t,w*w*t*t])
    velocity=lambda t,w:mp.matrix([w,2*w*w*t])
    acceleration=lambda t,w:mp.matrix([0,2*w*w])
    A=lambda x:mp.matrix([-B*x[1]/2,B*x[0]/2])
    errors=[];residual=[];integral_errors=[];derivative_errors=[];ratios=[]
    for w in map(mp.mpf,('-.3','0','.2')):
        for t in map(mp.mpf,('0','.01','.5','1')):
            x=path(t,w);v=velocity(t,w);a=A(x)
            W=v-x/t if t else mp.zeros(2,1)
            exact_W=mp.matrix([0,w*w*t])
            errors.extend((abs((a.T*x)[0]),abs((a.T*v)[0]-(a.T*W)[0]),mp.norm(W-exact_W)))
            residual.append(mp.norm(acceleration(t,w)+mp.matrix([0,-2*v[0]**2])))
            if t:
                integrated=mp.matrix([mp.quad(lambda s:s*acceleration(s,w)[i],[0,t])/t for i in range(2)])
                integral_errors.append(mp.norm(W-integrated))
            raw=lambda z:(A(path(t,z)).T*velocity(t,z))[0]
            canceled=lambda z:(A(path(t,z)).T*mp.matrix([0,z*z*t]))[0]
            for n in range(7):
                derivative_errors.append(abs(mp.diff(raw,w,n)-mp.diff(canceled,w,n)))
            for n in range(3):
                numerator=abs(mp.diff(lambda z:z*z*t,w,n))
                denominator=abs(mp.diff(lambda z:2*z*z,w,n))/2
                if denominator:ratios.append(float(numerator/denominator))
        phase=mp.quad(lambda t:(A(path(t,w)).T*velocity(t,w))[0],[0,1])
        errors.append(abs(phase-B*w**3/6))
    z=iv.mpf(0);o=iv.mpf(1)
    old,new,selected=radial_magnetic([o,o,o],[o,o,o],[z,z,z],2)
    return dict(radial_identity_error=float(max(errors)),geodesic_residual=float(max(residual)),
        acceleration_integral_error=float(max(integral_errors)),parameter_derivative_error=float(max(derivative_errors)),
        half_acceleration_max_ratio=max(ratios),flat_straight_generator_zero=all(x==0 for x in new+selected),
        old_flat_control_positive=all(x>0 for x in old),exact_geodesic_parameter_orders=6)


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_ordered_transport_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_ordered_transport.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing185 inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    frame=json.loads(project_path('dihole_positive_frame_results.json').read_text())
    local=json.loads(project_path('dihole_local_transport_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    heat=json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    updated=copy.deepcopy(prior);fields=updated['jets']['total_order_derivative_upper']
    for key in ('coframe','inverse_coframe'):fields[key]=frame['frame_bounds'][key]
    fields['gauge_per_abs_q']=frame['gauge']['gauge_per_abs_q_derivative_upper']
    normal=saved_normal(iv,local);records=[];nonincrease=[];baseline_agreement=[]
    for old,mr,hr in zip(previous['records'],magnetic['records'],heat['records']):
        if not old['q']==mr['q']==hr['q']:raise ArithmeticError('Flux records do not align.')
        operator=old['operator'];point=old['point_evaluation']
        proposed,generator=generator_bounds(iv,updated,old['q'],normal,local)
        baseline=[unpack(iv,x) for x in old['heat_transport']['uniform_generator_scaled_normalized_frechet_upper']]
        K=[a if a.b<=b.b else b for a,b in zip(proposed,baseline)]
        transport=transport_from_generator(iv,updated,operator,normal,local['real_geometry'],K)
        for key in ('parallel_transport_scaled_normalized_frechet_upper','u0_inverse_normalized_frechet_upper'):
            nonincrease.extend(bool(unpack(iv,a).b<=unpack(iv,b).b) for a,b in zip(transport[key],old['heat_transport'][key]))
        for row,oldrow in zip(transport['u_normalized_frechet_upper'],old['heat_transport']['u_normalized_frechet_upper']):
            nonincrease.extend(bool(unpack(iv,a).b<=unpack(iv,b).b) for a,b in zip(row,oldrow))
        # A historical-generator replay verifies the factored heat assembly.
        replay=transport_from_generator(iv,updated,operator,normal,local['real_geometry'],baseline)
        for key in ('u_normalized_frechet_upper','u0_inverse_normalized_frechet_upper'):
            baseline_agreement.append(replay[key]==old['heat_transport'][key])
        defect=assemble(iv,updated,local,operator,transport,point,hr,mr)
        scaled=record_bound(iv,updated,operator,defect,mr)
        baseline_B=old['scaled_bound']['heat_defect_B4_uniform_upper']
        ratio=iv.log(unpack(iv,baseline_B)/unpack(iv,scaled['heat_defect_B4_uniform_upper']),10)
        records.append(dict(q=old['q'],generator=generator,operator=operator,point_evaluation=point,
            heat_transport=transport,heat_defect=defect,scaled_bound=scaled,previous_B4_upper=baseline_B,
            radial_ceiling_log10_change_interval=[lower(ratio),upper(ratio)]))
    control=controls(iv)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        radial_identity=control['radial_identity_error']<1e-90,
        curved_geodesic_control=control['geodesic_residual']<1e-90,
        acceleration_integral=control['acceleration_integral_error']<1e-90,
        parameter_derivatives=control['parameter_derivative_error']<1e-90,
        half_acceleration_factor=control['half_acceleration_max_ratio']<=1,
        flat_straight_control=control['flat_straight_generator_zero'] and control['old_flat_control_positive'],
        historical_generator_replay=all(baseline_agreement),
        transport_bounds_do_not_increase=all(nonincrease),
        zeroth_magnetic_improves=all(r['generator']['magnetic_zero_ratio_interval'][1]<1 for r in records),
        derivative_schedule=all(r['heat_transport']['derivative_orders']==[16,14,12,10,8] for r in records),
        operator_and_point_unchanged=all(r['operator']==o['operator'] and r['point_evaluation']==o['point_evaluation'] for r,o in zip(records,previous['records'])),
        full_bounds_improve=all(r['radial_ceiling_log10_change_interval'][0]>0 for r in records),
        physical_source_not_claimed=all(r['scaled_bound']['physical_source'] is None and r['scaled_bound']['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,bound_encoding=previous['bound_encoding'],
        scope=['The existing source-centred magnetic radial condition removes the radial part of each geodesic velocity exactly.',
               'Acceleration bounds replace that canceled part in normalized generator jets through16, minimized against the old bounds; spin remains included.',
               'Geometry, gauge, operator, point bounds, source variation, normal domain, regulator and one-atom target are unchanged.',
               'The full heat/defect estimates are recomputed. No source sign, stable opening, atom transmission or clock advantage is established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Radial cancellation checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,records=[dict(q=r['q'],
        improved_magnetic_orders=r['generator']['improved_magnetic_orders'],
        B4_log10=r['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper'],
        improvement=r['radial_ceiling_log10_change_interval'],
        Cstar_log10=r['scaled_bound']['exact_low_order_remainder_Cstar_upper']['log10_upper']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-1]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['Previous propagation bound','Radial cancellation retained'],
        [row['previous_B4_upper']['log10_upper'],row['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper']],
        color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of B4 ceiling, D = 1',title='An exact cancellation tightens the field estimate')
    fig.supxlabel('Same physical field and one-atom target. Conservative ceilings, not actual errors.\nSupport, stability, atom passage and clock gain remain unproved.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
