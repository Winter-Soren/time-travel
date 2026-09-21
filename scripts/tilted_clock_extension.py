"""Continue saved tilted-clock patches in coordinates regular at A=0.

Beyond A=0 the x coordinate is temporal; static-frame density labels no longer
apply. Finite continuation does not establish global completeness or stability.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from project_paths import project_path
from tilted_clock_throat import rhs as proper_rhs


def fields(y,q2,vacuum):
    r,R,k,W,chi,Q=y
    X=q2/(2*Q);A=Q-W*W
    F=(X-1)**2/4-vacuum;V=chi*chi/(4*X)
    D=2*(X*(X-1)/2-V)
    pr=F+V;null=D*A/Q
    return dict(X=X,A=A,D=D,pr=pr,pt=pr+D*W*W/Q,
                rho_formal=null-pr,null_formal=null,theta=-chi/(2*X))


def rhs(x,y,q2,P,vacuum):
    r,R,k,W,chi,Q=y
    f=fields(y,q2,vacuum)
    Wx=chi/(2*f['X'])-2*R*W/r
    return np.array([R,-r*f['D']/(2*P*Q),
        (f['pt']+f['null_formal']/2)/P-2*k*R/r,
        Wx,f['D']*W/Q,2*k+2*W*Wx])


def transformed(y):
    r,v,eta,n,w,chi=y;N=np.exp(eta)
    return np.array([r,v/N,N*n,N*w,chi,N*N*(1+w*w)])


def continue_patch(y0,q2,r0,P,vacuum,refined=False):
    def horizon(x,y): return fields(y,q2,vacuum)['A']
    horizon.terminal=False;horizon.direction=0
    def radius_limit(x,y): return y[0]-.1*r0
    def lower_X(x,y): return fields(y,q2,vacuum)['X']-.1
    def upper_X(x,y): return 3.-fields(y,q2,vacuum)['X']
    def derivative_limit(x,y): return .3-abs(fields(y,q2,vacuum)['theta'])
    events=[horizon,radius_limit,lower_X,upper_X,derivative_limit]
    for event in events[1:]:
        event.terminal=True;event.direction=-1
    sol=solve_ivp(lambda x,y:rhs(x,y,q2,P,vacuum),(0.,4000.),y0,
        method='DOP853',rtol=2e-12 if refined else 2e-10,
        atol=2e-14 if refined else 2e-12,max_step=.05 if refined else .1,
        dense_output=True,events=events)
    if not sol.success:
        raise ArithmeticError(sol.message)
    reason=next((event.__name__ for event,times in zip(events[1:],sol.t_events[1:]) if len(times)),
                'finite_coordinate_limit')
    return sol,reason


def plot(results,states):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,4.1),layout='constrained')
    for name,state in states.items():
        x=state['x'];y=state['y'];row=results[name]
        label=name.replace('_',' ')
        axes[0].plot(x,state['A'],label=label)
        axes[1].plot(x,y[0]/row['initial_throat_radius'],label=label)
        axes[2].plot(x,state['theta'],label=label)
        for hx in row['horizon_positions']:
            for ax in axes: ax.axvline(hx,color='grey',alpha=.25,lw=.7)
    axes[0].axhline(0,color='black',lw=.7)
    axes[0].set(ylabel='A (static where positive)',title='Crossing the coordinate horizon')
    axes[1].set(ylabel='Areal radius / original throat radius',title='Extent of the continued solution')
    axes[2].set(ylabel='Expansion of the field flow / Lambda',title='Selected derivative stop')
    axes[2].axhline(.3,color='black',lw=.6,ls=':')
    axes[2].axhline(-.3,color='black',lw=.6,ls=':')
    for ax in axes:
        ax.set_xlabel('Quasiglobal coordinate from saved endpoint')
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Tilted-field continuation through a horizon: global behavior remains unresolved',fontsize=12)
    fig.supxlabel('Where A<0, the horizontal coordinate is temporal. These curves are not an embedding diagram or travel trajectory.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('tilted_clock_extension_progress.'+extension),dpi=160)
    plt.close(fig)


def run():
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    prior_path=project_path('tilted_clock_throat_results.json')
    prior=json.loads(prior_path.read_text(encoding='utf-8'))
    state_path=project_path('tilted_clock_throat_state.npz')
    if (not prior['passed'] or sha(state_path)!=prior['state_sha256']
            or sha(project_path('tilted_clock_throat.py'))!=prior['source_sha256']):
        raise ArithmeticError('Saved throat source/state provenance failed.')
    P=prior['parameters']['P'];vacuum=prior['parameters']['vacuum']
    results={};states={};checks=dict(saved_provenance=True)
    with np.load(state_path) as saved:
        for name,old in prior['branches'].items():
            q2=old['q_squared'];r0=old['throat_radius_times_cutoff']
            p=saved[name+'_y'][:,-1];initial=transformed(p)
            # Compare to the original proper-distance equations before crossing.
            dp=proper_rhs(0,p,q2,P,vacuum);r,v,eta,n,w,chi=p;N=np.exp(eta)
            change=np.array([dp[0],(dp[1]-v*dp[2])/N,N*(dp[3]+n*dp[2]),
                N*(dp[4]+w*dp[2]),dp[5],2*N*N*(dp[2]*(1+w*w)+w*dp[4])])/N
            transform_error=float(np.max(abs(change-rhs(0,initial,q2,P,vacuum))))
            sol,reason=continue_patch(initial,q2,r0,P,vacuum)
            ref,ref_reason=continue_patch(initial,q2,r0,P,vacuum,refined=True)
            x=np.linspace(0,sol.t[-1],1601);y=sol.sol(x);f=fields(y,q2,vacuum)
            r,R,k,W,chi,Q=y;A=f['A']
            dy=rhs(x,y,q2,P,vacuum)
            constraint=f['pr']-P*((A*R*R-1)/(r*r)+2*k*R/r)
            s=np.sqrt(2*f['X']);q=np.sqrt(q2)
            # Both saved branches have W>0 here: the advanced chart is regular.
            phi_x=-s*s/(s*W+q)
            scalar_norm=-q*phi_x-A*phi_x*phi_x/2
            kretschmann=4*(dy[2]**2+2*(k*R/r)**2+2*((k*R+A*dy[1])/r)**2+((1-A*R*R)/r**2)**2)
            horizon_rows=[]
            for hx,hy in zip(sol.t_events[0],sol.y_events[0]):
                hf=fields(hy,q2,vacuum)
                horizon_rows.append(dict(x=float(hx),r=float(hy[0]),X=float(hf['X']),
                    theta=float(hf['theta']),A_derivative=float(2*hy[2]),
                    mixed_stress_eigenvalues=[float(-hf['rho_formal']),float(hf['pr']),float(hf['pt'])],
                    advanced_radial_stress=float(hf['D']/hy[5])))
            common=np.linspace(0,min(sol.t[-1],ref.t[-1]),1601)
            errors=np.max(abs(sol.sol(common)-ref.sol(common)),axis=1)
            row=dict(initial_throat_radius=r0,stop_reason=reason,coordinate_endpoint=float(x[-1]),
                horizon_positions=list(map(float,sol.t_events[0])),horizons=horizon_rows,
                final_radius_over_throat=float(r[-1]/r0),final_A=float(A[-1]),
                minimum_X=float(np.min(f['X'])),maximum_X=float(np.max(f['X'])),
                maximum_abs_theta=float(np.max(abs(f['theta']))),
                maximum_Kretschmann_fourth_root=float(np.max(kretschmann)**.25),
                maximum_constraint_residual=float(np.max(abs(constraint))),
                maximum_scalar_norm_residual=float(np.max(abs(scalar_norm-f['X']))),
                proper_coordinate_transform_error=transform_error,
                refinement_maximum_errors=dict(zip(('r','r_x','A_x_over_two','W','chi','Q'),map(float,errors))),
                refinement_endpoint_difference=float(abs(sol.t[-1]-ref.t[-1])),
                rhs_evaluations=sol.nfev,
                scope='Finite continuation across the chart horizon. No global asymptotic state, future completeness, stability or time travel follows.')
            results[name]=row;states[name]=dict(x=x,y=y,constraint=constraint,kretschmann=kretschmann,**f)
            checks[name+'_coordinate_equivalence']=transform_error<1e-9
            checks[name+'_constraint']=row['maximum_constraint_residual']<1e-7
            checks[name+'_regular_scalar_chart']=row['maximum_scalar_norm_residual']<1e-10 and bool(np.min(s*W+q)>0)
            checks[name+'_crossed_horizon']=len(horizon_rows)>0 and bool(np.min(A)<0)
            checks[name+'_refinement']=bool(np.max(errors)<1e-6 and row['refinement_endpoint_difference']<1e-6 and reason==ref_reason)
            print(name,json.dumps(row),flush=True)
    out_state=project_path('tilted_clock_extension_state.npz')
    np.savez_compressed(out_state,**{name+'_'+key:value for name,state in states.items() for key,value in state.items()})
    result=dict(complete=True,passed=all(checks.values()),checks=checks,branches=results,
        source_sha256=sha(__file__),state_sha256=sha(out_state),
        inputs_sha256={name:sha(project_path(name)) for name in
            ('tilted_clock_throat.py','tilted_clock_throat_results.json','tilted_clock_throat_state.npz','project_paths.py')},
        parameters=dict(P=P,vacuum=vacuum,alpha=1.,minimum_radius_fraction=.1,minimum_X=.1,maximum_X=3.,maximum_abs_theta=.3),
        calculation_type='New horizon-regular equations initialized from saved149 endpoints; no rerun of the proper-distance integration.')
    project_path('tilted_clock_extension_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Extension check failed; result and states preserved.')
    plot(results,states)
    print(json.dumps(dict(passed=result['passed'],checks=checks),indent=2))
    return result


if __name__=='__main__':
    run()
