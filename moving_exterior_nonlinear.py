"""Finite prescribed motion with retarded JT/CFT proper-clock feedback.

The sewing and flat exterior are assumptions, not a four-dimensional solution.
Eight clock variables evolve in lab time; three passive integrals audit energy.
"""
import hashlib
import json
import math
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from moving_exterior_maps import prescribed_exterior, prescribed_exterior_checks, moving_energy_identity_checks


def solve_moving_exterior(amplitude, maximum_step=.05, loops=4, rtol=1e-11, atol=1e-13):
    exterior=prescribed_exterior(amplitude)
    delay,duration=exterior['delay'],exterior['duration']
    length=math.pi+delay;omega2=(2*math.pi/length)**2
    baseline=(1-omega2)/6;q0=-omega2/6;b0=math.pi*(omega2-1)/6
    final=loops*length
    segments,ends=[],[];end_clocks=[[],[]]
    diagnostics=dict(rejected_segments=0,minimum_external_lab_gap=math.inf,
        minimum_throat_lab_gap=math.inf,maximum_incoming_formula_difference=0.)

    class HistoryUnavailable(RuntimeError):
        pass

    def history(t):
        t=float(t)
        if t<=0:return np.array([t,1.,0.,.5,t,1.,0.,.5,q0*t,q0*t,0.])
        if not ends or t>ends[-1]:
            raise HistoryUnavailable(('lab history not completed',t,ends[-1:]))
        index=int(np.searchsorted(ends,t,side='left'))
        return segments[index].sol(t)

    @lru_cache(maxsize=32768)
    def inverse(target,side):
        target=float(target)
        if target<=0:return target
        clocks=end_clocks[side]
        if not clocks or target>clocks[-1]:
            raise HistoryUnavailable(('clock target beyond completed lab history',target,side,clocks[-1:]))
        index=int(np.searchsorted(clocks,target,side='left'))
        start=0. if index==0 else ends[index-1]
        return float(brentq(lambda t:segments[index].sol(t)[4*side]-target,
                           start,ends[index],xtol=5e-14))

    @lru_cache(maxsize=32768)
    def quantum_history(t,side):
        return incoming(float(t),history(float(t)),side)

    def incoming(t,y,side):
        if t<=0:return baseline
        i,j=4*side,4*(1-side)
        emitter,eprime,schwarzian=exterior['emission'](float(t),side)
        diagnostics['minimum_external_lab_gap']=min(diagnostics['minimum_external_lab_gap'],t-emitter)
        old=history(emitter)
        entry=inverse(float(old[j]-math.pi),side)
        gap=emitter-entry
        if gap<=0:raise RuntimeError(('nonretarded throat query; formulation unavailable',emitter,entry))
        diagnostics['minimum_throat_lab_gap']=min(diagnostics['minimum_throat_lab_gap'],gap)
        past=quantum_history(entry,side)
        own_s=y[i+3]-.5*y[i+1]**2;old_s=old[j+3]-.5*old[j+1]**2
        ratio=(eprime*old[j+1]/y[i+1])**2
        anomaly=(own_s-eprime**2*old_s-schwarzian)/(3*y[i+1]**2)
        # Equivalent deviation form preserves the exact static vacuum fixed point.
        value=baseline+ratio*(past-baseline)+(ratio-1)*(baseline-1/6)+anomaly
        literal=1/6+ratio*(past-1/6)+anomaly
        diagnostics['maximum_incoming_formula_difference']=max(
            diagnostics['maximum_incoming_formula_difference'],abs(value-literal))
        return float(value)

    def fluxes(t,y):
        if min(y[1],y[5])<=0:raise RuntimeError('Proper-clock monotonicity lost')
        entries=[inverse(float(y[4*side]-math.pi),1-side) for side in (0,1)]
        gaps=[t-entry for entry in entries]
        if min(gaps)<=0:raise RuntimeError('Retarded outgoing throat history unavailable')
        diagnostics['minimum_throat_lab_gap']=min(diagnostics['minimum_throat_lab_gap'],*gaps)
        jin=np.array([incoming(t,y,side) for side in (0,1)])
        jout=np.array([quantum_history(entries[side],1-side) for side in (0,1)])
        velocity=y[[1,5]];schwarzian=y[[3,7]]-.5*velocity**2
        qin=velocity**2*(jin-1/6)-schwarzian/3
        qout=velocity**2*(jout-1/6)-schwarzian/3
        return jin,jout,qin,qout,entries

    def motion_energy(t,y,flux=None):
        if flux is None:flux=fluxes(t,y)
        h,velocity,acceleration,jerk=exterior['motion'](float(t))
        gamma=1/math.sqrt(1-velocity**2)
        ap=gamma**3*acceleration
        ap_prime=gamma**4*jerk+3*gamma**6*velocity*acceleration**2
        su,sv=-ap_prime-.5*ap**2,ap_prime-.5*ap**2
        boundary=-b0*(y[[3,7]]-.5)
        power=velocity*(ap*boundary[1]-flux[2][1]-flux[3][1]-ap**2/3)
        return dict(h=h,velocity=velocity,gamma=gamma,proper_acceleration=ap,
                    S_U=su,S_V=sv,k=math.sqrt((1-velocity)/(1+velocity)),
                    boundary=boundary,power=power)

    def rhs(t,y):
        flux=fluxes(float(t),y);energy=motion_energy(float(t),y,flux)
        derivative=np.zeros(11)
        for side in (0,1):
            i=4*side;v,a,k=y[i+1:i+4]
            lapse=1. if side==0 else exterior['proper_lapse'](float(t))
            derivative[i:i+4]=lapse*np.array([v,a,v*(k-.5*v*v)+1.5*a*a/v,
                                               v*v*(flux[1][side]-flux[0][side])/b0])
        derivative[8]=flux[3][0]
        derivative[9]=(1-energy['velocity'])*(flux[3][1]+energy['proper_acceleration']**2/6)
        derivative[10]=energy['power']
        return derivative

    start=0.;initial=history(0.);evaluations=0
    while start<final:
        segment_width=min(delay/2,final-start)
        while True:
            stop=min(start+segment_width,final)
            try:
                segment=solve_ivp(rhs,(start,stop),initial,method='DOP853',dense_output=True,
                                  rtol=rtol,atol=atol,max_step=maximum_step)
                if not segment.success:raise RuntimeError(segment.message)
                break
            except HistoryUnavailable:
                diagnostics['rejected_segments']+=1;segment_width/=2
                if segment_width<1e-10:raise
        evaluations+=segment.nfev
        segments.append(segment);ends.append(stop)
        initial=segment.y[:,-1];start=stop
        for side in (0,1):
            clock=float(initial[4*side])
            if end_clocks[side] and clock<=end_clocks[side][-1]:
                raise RuntimeError('Completed clock history is not monotone')
            end_clocks[side].append(clock)
    diagnostics.update(segments=len(segments),solver_evaluations=evaluations,
        quantum_history_cache=quantum_history.cache_info()._asdict(),inverse_cache=inverse.cache_info()._asdict())
    return dict(history=history,inverse=inverse,fluxes=fluxes,motion_energy=motion_energy,
        exterior=exterior,diagnostics=diagnostics,parameters=dict(amplitude=amplitude,delay=delay,
            duration=duration,length=length,omega_squared=omega2,baseline_j=baseline,q0=q0,b0=b0,
            final_lab_time=final,maximum_step=maximum_step,rtol=rtol,atol=atol))


def nonlinear_solution_record(solution):
    p=solution['parameters'];final=p['final_lab_time'];delay=p['delay'];q0=p['q0']
    history,inverse=solution['history'],solution['inverse'];exterior=solution['exterior']
    times=np.unique(np.r_[np.linspace(0,final,401),delay,p['duration'],p['duration']+delay])
    values=np.array([history(float(t)) for t in times])
    fluxes=[solution['fluxes'](float(t),y) for t,y in zip(times,values)]
    energies=[solution['motion_energy'](float(t),y,f) for t,y,f in zip(times,values,fluxes)]
    lapses=np.array([exterior['proper_lapse'](float(t)) for t in times])
    proper=np.column_stack((times,[exterior['proper_time'](float(t)) for t in times]))
    lab_rates=values[:,[1,5]]*np.column_stack((np.ones(len(times)),lapses))
    total=[];bulk=[];boundary=[];improvement=[]
    for t,y,e in zip(times,values,energies):
        emitter=exterior['emission'](float(t),0)[0]
        past=history(emitter);older=history(float(t)-delay-e['h'])
        ep=solution['motion_energy'](emitter,past)
        endpoint=(e['k']*e['proper_acceleration']-ep['k']*ep['proper_acceleration'])/3
        outside=y[8]-older[8]+y[9]-past[9]+endpoint
        charges=e['boundary'][0]+e['gamma']*(e['boundary'][1]-2*e['proper_acceleration']/3)
        bulk.append(float(outside));boundary.append(float(charges));improvement.append(float(endpoint))
        total.append(float(outside-2*delay*q0+charges))
    total=np.asarray(total);residual=total-values[:,10]

    launches=np.unique(np.r_[np.linspace(0,final-math.pi-1,161),0.,.2,.5])
    probes=[];returns=np.full((len(launches),2),np.nan);one_way=np.full_like(returns,np.nan)
    end=history(final)
    for n,launch in enumerate(launches):
        initial=history(float(launch))
        for side in (0,1):
            target=float(initial[4*side]+math.pi);other=1-side
            exit_time=inverse(target,other) if target<=end[4*other] else None
            arrival=None
            if exit_time is not None:
                if side==0:arrival=exit_time+delay+exterior['motion'](exit_time)[0]
                else:
                    if final-delay-exterior['motion'](final)[0]>=exit_time:
                        arrival=brentq(lambda t:t-delay-exterior['motion'](t)[0]-exit_time,
                                       exit_time,final,xtol=5e-14)
                if arrival is not None and arrival>final:arrival=None
                one_way[n,side]=exit_time-launch
            if arrival is not None:returns[n,side]=arrival-launch
            probes.append(dict(launch_lab_time=float(launch),launch_side=side,
                throat_exit_lab_time=exit_time,return_arrival_lab_time=arrival,
                complete_return_delay=None if arrival is None else float(arrival-launch),
                complete_within_recorded_interval=arrival is not None))
    censored=int(np.count_nonzero(~np.isfinite(returns)))
    def column_stat(table,operation):
        return [None if not np.any(np.isfinite(column)) else float(operation(column[np.isfinite(column)]))
                for column in table.T]
    def nullable(column):
        return [float(value) if math.isfinite(value) else None for value in column]
    # Independent direct integrals use the equivalent jerk-containing null flux.
    quad_error=0.;work_quad_error=0.
    for t in (.731,2.113,4.0,5.317):
        y=history(t);e=solution['motion_energy'](t,y);emitter=exterior['emission'](t,0)[0]
        def a_flux(u):return float(solution['fluxes'](u,history(u))[3][0])
        def b_flux(u):
            state=history(u);f=solution['fluxes'](u,state);en=solution['motion_energy'](u,state,f)
            return (1-en['velocity'])*(f[3][1]+en['S_V']/3)
        cuts_a=[c for c in (0.,p['duration']) if t-delay-e['h']<c<t]
        cuts_b=[c for c in (0.,p['duration']) if emitter<c<t]
        direct=quad(a_flux,t-delay-e['h'],t,points=cuts_a,epsabs=2e-10,epsrel=2e-10,limit=120)[0]
        direct+=quad(b_flux,emitter,t,points=cuts_b,epsabs=2e-10,epsrel=2e-10,limit=120)[0]
        ep=solution['motion_energy'](emitter,history(emitter))
        primitive=y[8]-history(t-delay-e['h'])[8]+y[9]-history(emitter)[9]+(e['k']*e['proper_acceleration']-ep['k']*ep['proper_acceleration'])/3
        quad_error=max(quad_error,abs(direct-primitive))
        work=quad(lambda u:solution['motion_energy'](u,history(u))['power'],0,min(t,p['duration']),
                  epsabs=2e-11,epsrel=2e-10,limit=120)[0]
        work_quad_error=max(work_quad_error,abs(work-y[10]))
    jin=np.array([f[0] for f in fluxes]);jout=np.array([f[1] for f in fluxes])
    return dict(parameters=p,history_diagnostics=solution['diagnostics'],
        minimum_sampled_proper_clock_rate=float(values[:,[1,5]].min()),
        minimum_sampled_lab_clock_rate=float(lab_rates.min()),
        maximum_sampled_lab_clock_rate_change=float(np.max(abs(lab_rates-1))),
        maximum_sampled_clock_shift_from_lab=float(np.max(abs(values[:,[0,4]]-times[:,None]))),
        maximum_sampled_source_change=float(np.max(abs(jin-p['baseline_j']))),
        censored_return_probes=censored,
        minimum_sampled_one_way_time=column_stat(one_way,np.min),
        minimum_sampled_complete_return=column_stat(returns,np.min),
        maximum_sampled_complete_return=column_stat(returns,np.max),
        minimum_sampled_return_change=column_stat(returns-p['length'],np.min),
        maximum_sampled_return_change=column_stat(returns-p['length'],np.max),
        final_proper_clock_lag=exterior['final_proper_clock_lag'],
        reduced_drive_at_motion_end=float(history(p['duration'])[10]),
        reduced_drive_motion_end_to_final_difference=float(values[-1,10]-history(p['duration'])[10]),
        final_reduced_drive=float(values[-1,10]),maximum_reduced_balance_residual=float(np.max(abs(residual))),
        direct_bulk_quadrature_vs_improved_primitive_error=float(quad_error),
        direct_drive_quadrature_error=float(work_quad_error),
        first_order_improved_energy_remainder=float(np.max(abs(total-p['omega_squared']*
            np.array([exterior['motion'](float(t))[0] for t in times])/3))),
        selected_probes=[row for row in probes if row['launch_lab_time'] in (0.,.2,.5)],
        probes=probes,
        plot=dict(lab_time=times.tolist(),clock_shift_A=(values[:,0]-times).tolist(),
            clock_shift_B=(values[:,4]-times).tolist(),proper_clock_response_A=(values[:,0]-proper[:,0]).tolist(),
            proper_clock_response_B=(values[:,4]-proper[:,1]).tolist(),
            proper_clock_rate_A=values[:,1].tolist(),proper_clock_rate_B=values[:,5].tolist(),
            lab_clock_rate_A=lab_rates[:,0].tolist(),lab_clock_rate_B=lab_rates[:,1].tolist(),
            K_A=values[:,3].tolist(),K_B=values[:,7].tolist(),incoming_j_A=jin[:,0].tolist(),
            incoming_j_B=jin[:,1].tolist(),outgoing_j_A=jout[:,0].tolist(),outgoing_j_B=jout[:,1].tolist(),
            displacement=[exterior['motion'](float(t))[0] for t in times],
            return_launch_lab_time=launches.tolist(),return_A=nullable(returns[:,0]),return_B=nullable(returns[:,1]),
            reduced_bulk_energy=bulk,improved_boundary_energy=boundary,
            improved_total_energy_change=total.tolist(),reduced_drive=values[:,10].tolist(),
            energy_balance_residual=residual.tolist(),null_flux_endpoint_improvement=improvement))


def linear_comparison(record, maximum_step=.025):
    from mmp_pulse import _distributed_clock_response
    p=record['parameters'];amplitude=p['amplitude'];exterior=prescribed_exterior(amplitude)
    def motion_array(t):
        a=np.asarray(t);jets=np.array([exterior['motion'](float(x)) for x in a.ravel()]).T
        return tuple(jets[k].reshape(a.shape) for k in range(4))
    history,_,data=_distributed_clock_response([p['delay']],[1.],maximum_step,
        p['final_lab_time']+math.pi/2,lambda t:0.,exterior_motion=motion_array)
    def clocks(t):
        now=history(np.asarray(t)+math.pi/2);old=history(np.asarray(t)-math.pi/2)
        return ((now[..., [6,7]]-old[..., [7,6]])/data['b0'],
                (now[..., [0,3]]-old[..., [3,0]])/data['b0'])
    times=np.asarray(record['plot']['lab_time']);shift,rate=clocks(times)
    launch=np.asarray(record['plot']['return_launch_lab_time']);before=clocks(launch)[0]
    after=clocks(launch+math.pi)[0]
    geometric=np.array([[exterior['motion'](float(s+math.pi))[0],
                         exterior['motion'](float(s+math.pi+p['delay']))[0]] for s in launch])
    returns=p['length']+before-after[:,::-1]+geometric
    actual_shift=np.column_stack([record['plot']['clock_shift_A'],record['plot']['clock_shift_B']])
    actual_rates=np.column_stack([record['plot']['lab_clock_rate_A'],record['plot']['lab_clock_rate_B']])
    actual_returns=np.asarray([record['plot']['return_A'],record['plot']['return_B']],dtype=float).T
    completed=np.isfinite(actual_returns)
    return dict(maximum_step=maximum_step,
        maximum_clock_shift_remainder=float(np.max(abs(actual_shift-shift))),
        maximum_lab_clock_rate_remainder=float(np.max(abs(actual_rates-1-rate))),
        maximum_return_remainder=None if not np.any(completed) else float(np.max(abs(actual_returns[completed]-returns[completed]))),
        prediction=dict(clock_shift=shift.tolist(),clock_rate=(1+rate).tolist(),return_time=returns.tolist()),
        comparison='Actual nonlinear clocks and completed null returns versus a separately integrated linear model at identical lab times; no interpolation of the saved increment78 plot.')


def moving_exterior_nonlinear_checks(progress=None):
    begin=time.perf_counter();records=[]
    for amplitude,step in ((0.,.05),(.001,.05),(.0005,.05),(.05,.05),(.05,.025)):
        case_start=time.perf_counter();solution=solve_moving_exterior(amplitude,step)
        row=nonlinear_solution_record(solution)
        if amplitude in (.001,.0005):row['linear_comparison']=linear_comparison(row)
        row['runtime_seconds']=time.perf_counter()-case_start
        records.append(row)
        if progress is not None:progress(row)
    static=records[0];sp=static['plot']
    static_error=max(max(abs(np.asarray(sp[key]))) for key in ('clock_shift_A','clock_shift_B'))
    static_error=max(static_error,max(abs(np.asarray(sp['return_A'])-static['parameters']['length'])),
                     max(abs(np.asarray(sp['return_B'])-static['parameters']['length'])))
    static_error=max(static_error,static['maximum_sampled_source_change'],
        *(float(np.max(abs(np.asarray(sp[key])-reference))) for key,reference in
          (('proper_clock_rate_A',1.),('proper_clock_rate_B',1.),('K_A',.5),('K_B',.5))))
    ratios={key:records[1]['linear_comparison'][key]/records[2]['linear_comparison'][key]
            for key in ('maximum_clock_shift_remainder','maximum_lab_clock_rate_remainder','maximum_return_remainder')}
    ratios['first_order_improved_energy_remainder']=records[1]['first_order_improved_energy_remainder']/records[2]['first_order_improved_energy_remainder']
    keys=('clock_shift_A','clock_shift_B','lab_clock_rate_A','lab_clock_rate_B','K_A','K_B',
          'incoming_j_A','incoming_j_B','return_A','return_B','improved_total_energy_change','reduced_drive')
    refinement={key:float(np.max(abs(np.asarray(records[3]['plot'][key])-np.asarray(records[4]['plot'][key])))) for key in keys}
    # Resolve the linear comparison independently of the finite-motion tolerance.
    linear_fine=linear_comparison(records[2],maximum_step=.0125)
    linear_errors={key:float(np.max(abs(np.asarray(records[2]['linear_comparison']['prediction'][key])-
        np.asarray(linear_fine['prediction'][key])))) for key in ('clock_shift','clock_rate','return_time')}
    accepted=static_error<2e-10 and all(3.5<ratio<4.5 for ratio in ratios.values())
    accepted=accepted and max(refinement.values())<2e-8 and max(linear_errors.values())<2e-9
    accepted=accepted and all(row['censored_return_probes']==0 and row['maximum_reduced_balance_residual']<2e-8 and
        row['direct_bulk_quadrature_vs_improved_primitive_error']<2e-8 and
        row['direct_drive_quadrature_error']<2e-8 for row in records)
    files=('moving_exterior_maps.py','mmp_pulse.py')
    result=dict(model='Finite prescribed-motion JT/CFT feedback using proper clocks and lab-time retarded history',
        dependency_sha256={f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files},
        cases=records,exact_map_controls=prescribed_exterior_checks(),energy_identity_controls=moving_energy_identity_checks(),
        verification=dict(all_declared_checks_pass=bool(accepted),static_maximum_error=float(static_error),
            amplitude_halving_remainder_ratios=ratios,moderate_timestep_differences=refinement,
            linear_reference_timestep_differences=linear_errors,maximum_step_threshold=2e-8,
            linear_reference_threshold=2e-9,total_runtime_seconds=time.perf_counter()-begin),
        equations=dict(clock_state='c,v,a,K use proper derivatives; common lab RHS multiplies each four-state block by g_A=1 or g_B=sqrt(1-d_dot squared).',
            incoming='Jin_i=1/6+(Eprime*v_j(emitter)/v_i)^2*(Jin_i(entry)-1/6)+(S_i-Eprime^2*S_j(emitter)-Schwarzian(E))/(3*v_i^2), entry=inverse_c_i(c_j(emitter)-pi).',
            outgoing='Jout_i=Jin_other(inverse_c_other(c_i(t)-pi)). All queried histories are completed lab history; no proper-time ordering is assumed between different mouths.',
            return_A='Exit B solves c_B(exit)=c_A(launch)+pi; lab return=exit+d(exit).',
            return_B='Exit A solves c_A(exit)=c_B(launch)+pi; lab return solves arrival-d(arrival)=exit.',
            energy='Eext-2D Q0+B_A+gamma*(B_B-2a_proper/3)=W. B_i=-b0*(K_i-.5). Wdot=d_dot*(a_proper*B_B-qin_B-qout_B-a_proper squared/3).',
            bulk='F_A prime=qout_A; F_B prime=(1-d_dot)*(qout_B+a_proper squared/6). Eext=F_A(t)-F_A(t-d(t))+F_B(t)-F_B(eA)+(k(t)*a_proper(t)-k(eA)*a_proper(eA))/3, k=exp(-rapidity).'),
        scope='Dimensionless finite-motion result under assumed transparent proper-clock sewing, a flat one-dimensional exterior and externally prescribed timelike motion. Physical clock states are eight dimensional; three passive integrals audit an independently derived reduced balance. No injected packet. This is not physical 4D matching, a mouth equation of motion, apparatus work, preparation or a stability theorem. The inertial mouth charge reference is M0=0; an arbitrary rest energy adds M0*(gamma-1) to energy and M0*gamma_dot to the drive and is not fixed here. This reference drops out of completed-cycle net reduced transfer because both endpoint Lorentz factors equal one; it does not establish apparatus work. Missing endpoints would be censored, not zero delays or permanent trapping. Retarded-history step rejection is numerical bookkeeping, not chronology protection.')
    return result


if __name__=='__main__':
    target=Path(__file__).with_name('moving_exterior_nonlinear_results.json');checkpoint=[]
    def save_progress(row):
        checkpoint.append(row)
        target.write_text(json.dumps(dict(status='RUNNING',completed_cases=checkpoint),indent=2)+'\n',encoding='utf-8')
        print(json.dumps(dict(amplitude=row['parameters']['amplitude'],step=row['parameters']['maximum_step'],
            seconds=row['runtime_seconds'],ledger=row['maximum_reduced_balance_residual'],
            linear=row.get('linear_comparison',{}).get('maximum_return_remainder'))),flush=True)
    result=moving_exterior_nonlinear_checks(save_progress)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['verification'],indent=2),flush=True)
    if not result['verification']['all_declared_checks_pass']:raise SystemExit(1)
