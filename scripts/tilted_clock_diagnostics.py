"""Retain the failed long extension; check its horizon and derivative scales.

The short check is separate evidence, never a replacement for the failed test.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from project_paths import project_path
from tilted_clock_extension import fields,rhs,transformed


def run():
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    parent=json.loads(project_path('tilted_clock_extension_results.json').read_text())
    initial=json.loads(project_path('tilted_clock_throat_results.json').read_text())
    for name,result in [('tilted_clock_extension',parent),('tilted_clock_throat',initial)]:
        assert sha(project_path(name+'.py'))==result['source_sha256']
        assert sha(project_path(name+'_state.npz'))==result['state_sha256']
    P=initial['parameters']['P'];vac=initial['parameters']['vacuum']
    z=2*vac/(P+1+np.sqrt((P+1)**2+12*vac))
    exact_X=1+2*z;exact_chi=2*exact_X*np.sqrt(z);radius=1/np.sqrt(z)
    residual=[]
    for x in (0.,radius/2,radius,1.2*radius):
        state=np.array([radius,0.,-z*x,np.sqrt(z)*x,exact_chi,1.])
        residual.append(float(np.max(abs(rhs(x,state,2*exact_X,P,vac)-[0,0,-z,np.sqrt(z),0,0]))))
    exact=dict(radius=radius,X=exact_X,chi=exact_chi,z=z,
        maximum_residual=max(residual),scope='Exact constant-radius, zero-null-stress control; no flaring throat.')
    derivative={};local={};states={}
    with np.load(project_path('tilted_clock_extension_state.npz')) as saved, np.load(project_path('tilted_clock_throat_state.npz')) as seed:
        for name in initial['branches']:
            q2=initial['branches'][name]['q_squared']
            x=saved[name+'_x'];y=saved[name+'_y'];f=fields(y,q2,vac);dy=rhs(x,y,q2,P,vac)
            Xx=-f['X']*dy[5]/y[5]
            tx=-dy[4]/(2*f['X'])+y[4]*Xx/(2*f['X']**2)
            values=dict(abs_DuTheta_over_Lambda2=abs(y[3]*tx),
                abs_DeTheta_over_Lambda2=abs(np.sqrt(y[5])*tx),
                abs_Du_logX_over_Lambda=abs(y[3]*Xx/f['X']),
                acceleration_over_Lambda=abs(np.sqrt(y[5])*Xx/(2*f['X'])))
            derivative[name]={}
            for label,value in values.items():
                index=np.flatnonzero(value>1)
                derivative[name][label]=dict(maximum_saved=float(np.max(value)),
                    first_saved_sample_above_one=float(x[index[0]]) if index.size else None)
            y0=transformed(seed[name+'_y'][:,-1])
            def horizon(t,u): return fields(u,q2,vac)['A']
            horizon.terminal=False
            sols=[]
            for rtol,atol,step in [(2e-10,2e-12,.02),(2e-12,2e-14,.01)]:
                sol=solve_ivp(lambda t,u:rhs(t,u,q2,P,vac),(0.,5.),y0,method='DOP853',
                    rtol=rtol,atol=atol,max_step=step,dense_output=True,events=horizon)
                if not sol.success: raise ArithmeticError(sol.message)
                sols.append(sol)
            x=np.linspace(0,5,1001);y=sols[1].sol(x);f=fields(y,q2,vac)
            error=float(np.max(abs(sols[0].sol(x)-y)))
            herror=float(abs(sols[0].t_events[0][0]-sols[1].t_events[0][0]))
            r,R,k,W,chi,Q=y
            constraint=f['pr']-P*((f['A']*R*R-1)/(r*r)+2*k*R/r)
            dy=rhs(x,y,q2,P,vac)
            Xx=-f['X']*dy[5]/Q
            tx=-dy[4]/(2*f['X'])+chi*Xx/(2*f['X']**2)
            row=dict(interval=[0.,5.],maximum_refinement_error=error,horizon_position_error=herror,
                horizon=float(sols[1].t_events[0][0]),maximum_constraint_residual=float(np.max(abs(constraint))),
                maximum_abs_DuTheta=float(np.max(abs(W*tx))),
                maximum_abs_Du_logX=float(np.max(abs(W*Xx/f['X']))),
                minimum_X=float(np.min(f['X'])),maximum_X=float(np.max(f['X'])),
                passed=bool(error<1e-8 and herror<1e-8 and np.max(abs(constraint))<1e-8),
                scope='Separate finite near-horizon check. The longer integration failure remains in force.')
            local[name]=row;states[name]=dict(x=x,y=y,**f)
    result=dict(near_horizon_checks_passed=all(row['passed'] for row in local.values()),
        parent_long_extension_passed=parent['passed'],global_medium_established=False,
        exact_constant_radius_control=exact,near_horizon=local,saved_derivative_diagnostics=derivative,
        derivative_scope='Sampled nominal-scale indicators, not suprema or the actual EFT cutoff. The long reversing branch is unconverged at the requested tolerance.',
        source_sha256=sha(__file__),inputs_sha256={name:sha(project_path(name)) for name in
            ('tilted_clock_extension.py','tilted_clock_extension_results.json','tilted_clock_extension_state.npz',
             'tilted_clock_throat.py','tilted_clock_throat_results.json','tilted_clock_throat_state.npz','project_paths.py')})
    if not result['near_horizon_checks_passed'] or exact['maximum_residual']>1e-12:
        raise ArithmeticError('Separate diagnostic check failed.')
    state_path=project_path('tilted_clock_diagnostics_state.npz')
    np.savez_compressed(state_path,**{name+'_'+key:value for name,state in states.items() for key,value in state.items()})
    result['state_sha256']=sha(state_path)
    project_path('tilted_clock_diagnostics_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11.5,4),layout='constrained')
    for name,state in states.items():
        x=state['x'];label=name.replace('_',' ')
        axes[0].plot(x,state['A'],label=label)
        axes[1].plot(x,state['X'],label=label)
        axes[2].plot(x,state['theta'],label=label)
        for ax in axes: ax.axvline(local[name]['horizon'],color='gray',lw=.6,alpha=.4)
    axes[0].axhline(0,color='black',lw=.7)
    axes[0].set(ylabel='A',title='Regular chart crosses A = 0')
    axes[1].set(ylabel='Scalar norm X',title='Timelike field remains finite')
    axes[2].set(ylabel='Field expansion / Lambda',title='Finite local derivative')
    for ax in axes:
        ax.set_xlabel('Coordinate from the saved static endpoint')
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Near-horizon checks pass; the long reversing continuation misses its accuracy target',fontsize=12)
    fig.supxlabel('A < 0 changes the coordinate roles. This is a local solution of a hypothetical action, not a verified travel route.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('tilted_clock_diagnostics_progress.'+extension),dpi=160)
    plt.close(fig)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    run()
