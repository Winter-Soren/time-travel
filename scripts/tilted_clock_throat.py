"""Zero-current stationary tilted-clock equations and finite throat continuations.

Lambda=1. A local throat and a finite continuation are not a global wormhole.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from nonlinear_packet_constraints import background
from project_paths import project_path


def matter(y, q2, vacuum, alpha=1.):
    r, velocity, eta, n, w, chi = y
    X=q2*np.exp(-2*eta)/(2*(1+w*w))
    F=(X-1)**2/4-vacuum
    V=chi*chi/(4*alpha*X)
    B=((X-1)/2-V/X)/(1+w*w)
    null=2*X*B
    pr=F+V
    return dict(X=X,F=F,V=V,B=B,null=null,pr=pr,rho=null-pr,
                pt=pr+null*w*w,theta=-chi/(2*alpha*X))


def rhs(l, y, q2, P, vacuum):
    r,v,eta,n,w,chi=y
    m=matter(y,q2,vacuum)
    acceleration=(1-v*v-m['rho']*r*r/P)/(2*r)
    return np.array([v,acceleration,n,
        m['pt']/P-n*n-n*v/r-acceleration/r,
        chi/(2*m['X'])-(n+2*v/r)*w,m['null']*w])


def integrate(seed, P, vacuum, refined=False):
    X0,w0,chi0=seed
    q2=2*X0*(1+w0*w0)
    pr0=(X0-1)**2/4-vacuum+chi0*chi0/(4*X0)
    if pr0>=0:
        raise ValueError('A regular static throat requires negative radial pressure.')
    r0=np.sqrt(-P/pr0)
    y0=np.array([r0,0.,0.,0.,w0,chi0])
    if matter(y0,q2,vacuum)['null']>=0:
        raise ValueError('Initial data do not flare outward.')
    # ponytail: finite static patches only; a horizon needs regular coordinates
    # and new boundary conditions, not extrapolation past this lapse cutoff.
    def lapse_cut(l,y): return y[2]-np.log(.1)
    def radius_cut(l,y): return y[0]-.1*r0
    def low_X(l,y): return matter(y,q2,vacuum)['X']-.1
    def high_X(l,y): return 3.-matter(y,q2,vacuum)['X']
    def derivative_cut(l,y): return .3-abs(matter(y,q2,vacuum)['theta'])
    events=[lapse_cut,radius_cut,low_X,high_X,derivative_cut]
    for event in events:
        event.terminal=True
        event.direction=-1
    sol=solve_ivp(lambda l,y:rhs(l,y,q2,P,vacuum),(0.,400.),y0,method='DOP853',
        rtol=2e-12 if refined else 2e-10,atol=2e-14 if refined else 2e-12,
        max_step=.1 if refined else .2,dense_output=True,events=events)
    if not sol.success:
        raise ArithmeticError(sol.message)
    reason=next((event.__name__ for event,times in zip(events,sol.t_events) if len(times)),
                'finite_distance_limit')
    return sol,q2,r0,reason


def inspect(sol,q2,r0,P,vacuum):
    l=np.linspace(0,sol.t[-1],1601)
    y=sol.sol(l)
    r,v,eta,n,w,chi=y
    m=matter(y,q2,vacuum)
    dy=np.stack([rhs(t,u,q2,P,vacuum) for t,u in zip(l,y.T)],axis=1)
    constraint=m['pr']-P*((v*v-1)/(r*r)+2*n*v/r)
    # Independent original covariant stress, before zero-current elimination.
    gamma=np.sqrt(1+w*w);s=np.sqrt(2*m['X'])
    Yprime=-dy[5]/2
    DuY=-w*Yprime
    iso=m['F']+m['X']*m['theta']**2+2*DuY
    uu=2*m['X']*((m['X']-1)/2-m['theta']**2)-2*DuY
    rho_original=-iso+uu*gamma**2
    pr_original=iso+uu*w*w+4*w*Yprime
    flux_original=uu*gamma*w+2*gamma*Yprime
    current_original=-s*((m['X']-1)/2-m['theta']**2)*w-2*(1+w*w)*Yprime/s
    # Conservation uses independent differentiation of F(X)+chi²/(4X).
    Xprime=-2*m['X']*(n+w*dy[4]/(1+w*w))
    prprime=((m['X']-1)/2-m['V']/m['X'])*Xprime+chi*dy[5]/(2*m['X'])
    conservation=prprime+m['null']*n-2*(m['pt']-m['pr'])*v/r
    finite=(np.isfinite(y).all() and all(np.isfinite(a).all() for a in m.values()))
    summary=dict(end_distance_times_cutoff=float(l[-1]),nfev=sol.nfev,
        throat_radius_times_cutoff=float(r0),q_squared=q2,
        initial={key:float(value[0]) for key,value in m.items()},
        final={key:float(value[-1]) for key,value in m.items()},
        final_radius_over_throat=float(r[-1]/r0),maximum_radius_over_throat=float(np.max(r)/r0),
        final_lapse=float(np.exp(eta[-1])),minimum_lapse=float(np.exp(np.min(eta))),
        minimum_X=float(np.min(m['X'])),maximum_X=float(np.max(m['X'])),
        maximum_abs_theta=float(np.max(abs(m['theta']))),
        minimum_radius_derivative=float(np.min(v)),maximum_radius_derivative=float(np.max(v)),
        noncentral_tilt_sign_changes=int(np.count_nonzero(w[1:-1]*w[2:]<0)),
        throat_flare=-float(m['null'][0])*r0*r0/P,
        maximum_constraint_residual=float(np.max(abs(constraint))),
        maximum_original_current_residual=float(np.max(abs(current_original))),
        maximum_original_flux_residual=float(np.max(abs(flux_original))),
        maximum_covariant_stress_residual=float(max(np.max(abs(rho_original-m['rho'])),
                                                   np.max(abs(pr_original-m['pr'])),
                                                   np.max(abs(iso-m['pt'])))),
        maximum_conservation_residual=float(np.max(abs(conservation))),finite=bool(finite),
        scope='Sampled finite static patch. No global ends, perturbation stability, EFT completion or chronology change established.')
    state=dict(l=l,y=y,constraint=constraint,**m)
    return summary,state


def plot(result,states):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')
    for name,state in states.items():
        l=state['l'];r,v,eta,n,w,chi=state['y']
        label=name.replace('_',' ')
        axes[0,0].plot(l,r/r[0],label=label)
        axes[0,1].plot(l,np.exp(eta),label=label)
        axes[1,0].plot(l,state['null'],label=label)
        axes[1,1].plot(l,w,label=label)
    axes[0,0].set(ylabel='Areal radius / initial throat radius',title='Does the opening keep widening?')
    axes[0,1].set(ylabel='Static clock rate N',title='Finite static patch only')
    axes[1,0].set(ylabel='Radial null stress / Lambda^4',title='Negative stress alone is insufficient')
    axes[1,1].set(ylabel='Radial tilt w',title='Flow reversal avoids one conditional obstruction')
    for ax in axes.flat:
        ax.set_xlabel('Proper distance from throat (1/Lambda)')
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    axes[1,0].axhline(0,color='black',lw=.6)
    fig.suptitle('Self-consistent local tilted-clock throats: finite continuations',fontsize=13)
    fig.supxlabel('The throat is imposed as initial data. Formation, complete traversability and stability are not proved.',fontsize=9)
    for extension in ('png','pdf'):
        fig.savefig(project_path('tilted_clock_throat_progress.'+extension),dpi=160)
    plt.close(fig)


def run():
    base=background();P,vacuum=base['P'],base['rho']
    seeds=dict(through_flow=(.8,1.,0.),flow_reversal=(1.,0.,.02))
    results={};states={};checks={}
    for name,seed in seeds.items():
        sol,q2,r0,reason=integrate(seed,P,vacuum)
        row,state=inspect(sol,q2,r0,P,vacuum)
        refined,_,_,refined_reason=integrate(seed,P,vacuum,refined=True)
        common=np.linspace(0,min(sol.t[-1],refined.t[-1]),1601)
        errors=np.max(abs(sol.sol(common)-refined.sol(common)),axis=1)
        row.update(seed=dict(X=seed[0],w=seed[1],chi=seed[2]),stop_reason=reason,
            refinement_maximum_absolute_errors=dict(zip(('r','r_prime','log_N','log_N_prime','w','chi'),map(float,errors))),
            refinement_endpoint_difference=float(abs(sol.t[-1]-refined.t[-1])))
        checks[name+'_finite']=row['finite']
        checks[name+'_constraint']=row['maximum_constraint_residual']<1e-7
        checks[name+'_current_and_flux']=max(row['maximum_original_current_residual'],row['maximum_original_flux_residual'])<1e-10
        checks[name+'_covariant_stress']=row['maximum_covariant_stress_residual']<1e-10
        checks[name+'_conservation']=row['maximum_conservation_residual']<1e-10
        checks[name+'_refinement']=bool(np.max(errors)<1e-6 and row['refinement_endpoint_difference']<1e-6 and reason==refined_reason)
        results[name]=row;states[name]=state
        print(name,reason,json.dumps(row),flush=True)
    checks['flow_reversal_regular_center']=bool(states['flow_reversal']['y'][4,0]==0 and states['flow_reversal']['y'][5,0]!=0)
    checks['both_initial_throats_flare']=all(row['throat_flare']>0 for row in results.values())
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    state_path=project_path('tilted_clock_throat_state.npz')
    np.savez_compressed(state_path,**{name+'_'+key:value for name,state in states.items() for key,value in state.items()})
    result=dict(complete=True,passed=all(checks.values()),checks=checks,
        parameters=dict(P=P,alpha=1.,vacuum=vacuum,units='Lambda=1',maximum_proper_distance=400.,
            stops=dict(minimum_lapse=.1,minimum_radius_fraction=.1,minimum_X=.1,maximum_X=3.,maximum_abs_theta=.3)),
        branches=results,source_sha256=sha(__file__),state_sha256=sha(state_path),
        inputs_sha256={name:sha(project_path(name)) for name in
            ('nonlinear_packet_constraints.py','expanding_field_completion.py','expanding_field_completion_results.json','project_paths.py')},
        scope='First-order stationary zero-current field plus Einstein equations. A finite local throat continuation is not a global solution or stable supporting material.')
    project_path('tilted_clock_throat_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Tilted-clock checks failed; state and results preserved.')
    plot(result,states)
    print(json.dumps(dict(passed=result['passed'],checks=checks),indent=2))
    return result


if __name__=='__main__':
    run()
