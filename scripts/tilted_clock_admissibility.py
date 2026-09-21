"""Causal orientation and coupled frozen principal checks on saved149 states.

No background is evolved. Principal health is not global perturbation stability.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline

from expanding_field_completion import coefficients
from nonlinear_packet_constraints import background
from project_paths import project_path


def principal(P,X,theta):
    c=X;U=X*(3*X-1)/2
    gate=U-c*(P-9*c)*theta*theta/(P+3*c)
    return gate,P*c/((P+3*c)*gate)


def direct_block(P,X,K):
    """Quadratic ADM Hessian in a transverse spatial gauge; retain tensor mixing."""
    c=X;theta=np.trace(K);U=X*(3*X-1)/2
    maps=np.zeros((6,3,3))
    maps[0]=np.diag([1.,1.,0.])
    maps[1]=np.diag([.5,-.5,0.])
    maps[2,0,1]=maps[2,1,0]=.5
    maps[3,2,2]=-1.
    maps[4,0,2]=maps[4,2,0]=-.5
    maps[5,1,2]=maps[5,2,1]=-.5
    tr=np.trace(maps,axis1=1,axis2=2)
    matrix=np.zeros((7,7))
    matrix[:6,:6]=P/2*np.einsum('iab,jab->ij',maps,maps)-(P/2+c)*np.outer(tr,tr)
    matrix[:6,6]=matrix[6,:6]=(-P*np.einsum('ab,iab->i',K,maps)+(P+6*c)*theta*tr)/2
    matrix[6,6]=U+P/2*(np.sum(K*K)-theta*theta)-6*c*theta*theta
    keep=[0,1,2,6];shift=[3,4,5]
    reduced=matrix[np.ix_(keep,keep)]-matrix[np.ix_(keep,shift)]@np.linalg.solve(
        matrix[np.ix_(shift,shift)],matrix[np.ix_(shift,keep)])
    A0=reduced[:3,:3];h=reduced[:3,3];D=reduced[3,3]
    schur=D-h@np.linalg.solve(A0,h)
    kinetic=A0-np.outer(h,h)/D
    e=np.array([1.,0.,0.])
    gyro=P/D*(np.outer(e,h)-np.outer(h,e))
    stiffness=P*P/D*np.outer(e,e)
    evolution=np.block([[np.zeros((3,3)),np.eye(3)],
        [-np.linalg.solve(kinetic,stiffness),-np.linalg.solve(kinetic,gyro)]])
    eigen=np.linalg.eigvals(evolution)
    leading=eigen[np.argmax(abs(eigen))]
    return dict(schur=float(schur),minimum_kinetic_eigenvalue=float(np.min(np.linalg.eigvalsh(kinetic))),
        k4_coefficient=float(-(leading*leading).real),leading_real_part=float(leading.real),D=float(D))


def run():
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    prior=json.loads(project_path('tilted_clock_throat_results.json').read_text())
    diagnostic=json.loads(project_path('tilted_clock_diagnostics_results.json').read_text())
    for stem,saved in [('tilted_clock_throat',prior),('tilted_clock_diagnostics',diagnostic)]:
        if (sha(project_path(stem+'.py'))!=saved['source_sha256']
                or sha(project_path(stem+'_state.npz'))!=saved['state_sha256']):
            raise ArithmeticError('Saved-state provenance failed: '+stem)
    if not prior['passed'] or not diagnostic['near_horizon_checks_passed']:
        raise ArithmeticError('Required local checks have not passed.')
    P=prior['parameters']['P'];profiles={};summary={};block_errors=[];block_min=[];block_re=[]
    clock={};orientation={}
    with np.load(project_path('tilted_clock_throat_state.npz')) as local, np.load(project_path('tilted_clock_diagnostics_state.npz')) as near:
        for name,row in prior['branches'].items():
            for label,state in [('static',local),('horizon',near)]:
                y=state[name+'_y'];X=state[name+'_X'];theta=state[name+'_theta']
                if label=='static':
                    r,v,eta,n,w,chi=y
                    Kperp=-w*v/r;coordinate=state[name+'_l']
                else:
                    r,R,k,W,chi,Q=y
                    Kperp=-W*R/r;coordinate=state[name+'_x']
                Krad=theta-2*Kperp
                gate,a4=principal(P,X,theta)
                key=name+'_'+label
                profiles[key]=dict(coordinate=coordinate,gate=gate,a4=a4)
                summary[key]=dict(minimum_A_star=float(np.min(gate)),maximum_A_star=float(np.max(gate)),
                    minimum_k4_coefficient=float(np.min(a4)),maximum_k4_coefficient=float(np.max(a4)),
                    samples=len(X),scope='Frozen preferred-clock principal system only; lower orders and global stability unresolved.')
                indices=sorted(set([0,len(X)//2,len(X)-1,int(np.argmin(gate)),int(np.argmax(a4))]))
                for i in indices:
                    for angle in np.linspace(0,np.pi/2,7):
                        ca,sa=np.cos(angle),np.sin(angle)
                        rotation=np.array([[-sa,ca,0.],[0.,0.,1.],[ca,sa,0.]])
                        K=rotation@np.diag([Krad[i],Kperp[i],Kperp[i]])@rotation.T
                        out=direct_block(P,float(X[i]),K)
                        block_errors.append([abs(out['schur']/gate[i]-1),abs(out['k4_coefficient']/a4[i]-1)])
                        block_min.append(out['minimum_kinetic_eigenvalue']);block_re.append(abs(out['leading_real_part']))
            r,v,eta,n,w,chi=local[name+'_y'];l=local[name+'_l']
            N=np.exp(eta);integrand=1/N
            optical=float(CubicSpline(l,integrand).integrate(0,l[-1]))
            simpson_value=float(simpson(integrand,x=l))
            clock[name]=dict(half_patch_optical_time=optical,quadrature_difference=abs(optical-simpson_value),
                symmetric_patch_light_crossing=2*optical,same_station_reflected_round_trip=4*optical,
                scope='Killing-coordinate time differences Delta t in units 1/Lambda, not station proper time or generally preferred-clock time. Conditional ideal metric-null test between symmetric N=.1 stations; no reflector or apparatus constructed.')
            y=near[name+'_y'];x=near[name+'_x'];W=y[3];Q=y[5];A=Q-W*W
            plus=-W+np.sqrt(Q);minus=-W-np.sqrt(Q)
            outside=A<0;inside=A>0
            orientation[name]=dict(right_exterior_plus_max=float(np.max(plus[outside])),
                right_exterior_minus_max=float(np.max(minus[outside])),
                inside_plus_min=float(np.min(plus[inside])),inside_minus_max=float(np.max(minus[inside])),
                minimum_dphi_along_normalized_future_null=float(np.min(np.sqrt(2*near[name+'_X']))),
                scope='q>0,W>0: both right exterior future metric-null directions point toward the static patch.')
        figure_x=near['flow_reversal_x'].copy();figure_y=near['flow_reversal_y'].copy()
    base=background();old=coefficients(.01,base['H'])
    _,ds_a4=principal(P,base['X'],3*base['H'])
    desitter_error=abs(ds_a4/old['dispersion_A_times_cutoff2']-1)
    errors=np.max(block_errors,axis=0)
    checks=dict(saved_provenance=True,all_sampled_principal_gates_positive=all(row['minimum_A_star']>0 for row in summary.values()),
        directional_schur_identity=bool(errors[0]<1e-8),coupled_k4_identity=bool(errors[1]<1e-8),
        sampled_full_kinetic_matrices_positive=min(block_min)>0,
        high_frequency_eigenvalue_real_part=max(block_re)<1e-8,
        desitter_limit=desitter_error<1e-12,
        right_horizon_incoming=all(row['right_exterior_plus_max']<0 and row['right_exterior_minus_max']<0 for row in orientation.values()),
        static_patch_two_null_directions=all(row['inside_plus_min']>0 and row['inside_minus_max']<0 for row in orientation.values()),
        light_round_trip_positive=all(row['same_station_reflected_round_trip']>0 for row in clock.values()),
        optical_quadrature_consistent=all(row['quadrature_difference']<1e-4 for row in clock.values()))
    result=dict(passed=all(checks.values()),checks=checks,principal_profiles=summary,
        direct_matrix_checks=dict(count=len(block_errors),maximum_relative_schur_error=float(errors[0]),
            maximum_relative_k4_error=float(errors[1]),minimum_kinetic_eigenvalue=min(block_min),
            maximum_high_frequency_real_part=max(block_re),desitter_relative_error=desitter_error),
        horizon_orientation=orientation,conditional_static_clock_tests=clock,
        causal_decision='A smooth real single-valued scalar with X>0 forbids closed future metric-causal curves wherever those assumptions hold. Local data do not establish a complete global extension; dispersive field signals require their own causal analysis.',
        stability_decision='Positive coupled frozen kinetic and k4 terms at saved samples; no lower-order or global stability conclusion.',
        source_sha256=sha(__file__),inputs_sha256={name:sha(project_path(name)) for name in
            ('tilted_clock_throat.py','tilted_clock_throat_results.json','tilted_clock_throat_state.npz',
             'tilted_clock_diagnostics.py','tilted_clock_diagnostics_results.json','tilted_clock_diagnostics_state.npz',
             'expanding_field_completion.py','expanding_field_completion_results.json','nonlinear_packet_constraints.py','project_paths.py')},
        calculation_type='Saved-state analysis and local algebraic matrices only; no background or perturbation evolution.')
    project_path('tilted_clock_admissibility_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Admissibility checks failed; evidence preserved.')
    plot(profiles,figure_x,figure_y,diagnostic['near_horizon']['flow_reversal']['horizon'])
    print(json.dumps(result,indent=2))
    return result


def plot(profiles,x,y,horizon):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,4.3),layout='constrained')
    for name in ('through_flow','flow_reversal'):
        p=profiles[name+'_static'];label=name.replace('_',' ')
        axes[0].plot(p['coordinate'],p['gate'],label=label)
        axes[1].plot(p['coordinate'],p['a4'],label=label)
    axes[0].set(xlabel='Proper distance from throat (1/Lambda)',ylabel='A*',title='Frozen kinetic gate stays positive')
    axes[1].set(xlabel='Proper distance from throat (1/Lambda)',ylabel='Quartic dispersion coefficient (Lambda^-2)',title='Highest-frequency term stays positive')
    axes[0].legend(fontsize=8)
    ax=axes[2];W,Q=y[3],y[5]
    points=np.array([.2,1.4,2.6,3.8,4.6])
    wp=np.interp(points,x,W);qp=np.interp(points,x,Q)
    for row,sign in [(1,1),(0,-1)]:
        direction=np.sign(-wp+sign*np.sqrt(qp))
        ax.quiver(points,np.full(points.shape,row),direction*.35,np.zeros(points.shape),
            angles='xy',scale_units='xy',scale=1,color='#205e80' if sign==1 else '#b85432',width=.009,zorder=3)
    ax.axvline(horizon,color='black',ls='--',lw=1)
    ax.axvspan(horizon,5,color='#eeeeee')
    ax.text(horizon+.1,1.45,'horizon',fontsize=8)
    ax.text(2.6,.52,'Both future directions point inward',fontsize=8,ha='center')
    ax.set(xlim=(0,5),ylim=(-.4,1.65),yticks=[0,1],yticklabels=['Future ray -','Future ray +'],
        xlabel='Right-side coordinate (direction arrows)',title='Selected horizon allows entry')
    for ax in axes[:2]: ax.grid(alpha=.2)
    fig.suptitle('Local short-wavelength health does not create a route into the past',fontsize=13)
    fig.supxlabel('A smooth real clock with timelike gradient advances along every future metric-causal path. Global stability remains unproved.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('tilted_clock_admissibility_progress.'+extension),dpi=160)
    plt.close(fig)


if __name__=='__main__':
    run()
