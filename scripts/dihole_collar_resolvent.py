"""Full-operator collar resolvent comparison and a filtered trace error.

Uses the actual radial trial residual with an indispensable smooth cutoff.
The finite selected observable is not a stress tensor or total source.
"""
import json
import math

import numpy as np
from scipy.integrate import quad

from balanced_dihole import symbolic_fields
from dihole_aligned_channel import geometry, norm, trial
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_local_source import local_first
from dihole_slice_channel import sha
from dihole_stretched_collar import exact_target_parameters, smooth_step
from project_paths import project_path


NAME='dihole_collar_resolvent'


def cutoff(s,inner=False):
    left,right,width=(.015,.035,.005) if inner else (.005,.045,.01)
    a,da=smooth_step((np.asarray(s)-left)/width)
    b,db=smooth_step((right-np.asarray(s))/width)
    return a*b,(da*b-a*db)/width


def sector_count(q):
    if type(q) is not int or q<1000:raise ValueError('Require an integer q>=1000.')
    lo=(17*q-10+19)//20;hi=(19*q-10)//20
    return dict(first_n=lo,last_n=hi,count=hi-lo+1,mu_equals='n+1/2')


def bounds(iv,previous,q):
    band=sector_count(q);c=previous['uniform_bound']['constants'];u=previous['uniform_bound']
    value=lambda key:iv.mpf(c[key][1])
    wmin=iv.mpf(c['W_lower'][0]);wmax=value('W_upper')
    T=value('tilt_bound');Ts=value('tilt_s_bound');drift=value('theta0_s_bound')
    Ct=iv.mpf(u['step_bounds']['C_theta'])
    absolute=lambda key:iv.mpf(max(abs(x) for x in previous['regions']['global']['derived'][key]))
    Ltheta,Ls=absolute('log_W_theta'),absolute('log_W_R')
    C1,C2=(iv.mpf(x) for x in u['even_moment_coefficients'][:2])
    Lf=(Ct*T*T+Ltheta)/wmin
    r1=iv.sqrt(Lf*Lf*C1/q+T*T/wmin**2)
    vmin=1/(iv.sqrt(1+T*T)*wmax)
    velocity_s=(Ts*T+Ls)/wmin+drift*Lf+Lf*iv.sqrt(C1)*(2*value('A2_bound')*iv.sqrt(C2)+value('log_weight_transport_bound'))/iv.sqrt(q)
    # Compose graph constants from the exported bounds, not hidden tighter values.
    r1_upper=upper(r1);v_lower=lower(vmin);vs_upper=upper(velocity_s)
    a=iv.mpf(r1_upper)/iv.mpf(v_lower)
    # Outer cutoff ramps are disjoint; S'<=8 and ramp width=1/100.
    chi1=iv.mpf(800)
    b=iv.mpf(r1_upper)*iv.mpf(vs_upper)/(2*iv.mpf(v_lower))+value('r0_bound_D')+chi1/wmin
    a_upper=upper(a);b_upper=upper(b)
    optical_volume=iv.mpf(1)/50/iv.mpf(v_lower)
    rows=[]
    for frequency_factor in (.1,1.,10.):
        y_value=float(frequency_factor*math.sqrt(q));y=iv.mpf(y_value)
        lambda_value=float(math.sqrt(q));lam=iv.mpf(lambda_value)
        scaled=iv.sqrt(iv.mpf(a_upper)**2+(iv.mpf(b_upper)/y)**2)
        error=scaled/y;hs2=lam*optical_volume/(2*iv.sqrt(iv.pi));trace=error*hs2
        rows.append(dict(y_D=y_value,filter_Lambda_D=lambda_value,
            scaled_resolvent_error_upper=upper(scaled),resolvent_error_upper_over_D=upper(error),
            filter_HS_norm_squared_upper=upper(hs2),per_sector_trace_error_upper_over_D=upper(trace),
            band_trace_error_upper_over_D=upper(trace*band['count'])))
    return dict(q=q,band=band,r1_upper=r1_upper,velocity_lower=v_lower,velocity_derivative_upper_D=vs_upper,
        f_theta_lipschitz_upper=upper(Lf),graph_h_coefficient_upper=a_upper,graph_zero_coefficient_upper_D=b_upper,
        outer_cutoff_derivative_upper_D=800,inner_optical_volume_upper_over_D=upper(optical_volume),frequency_rows=rows)


def field_controls(iv,previous):
    function=symbolic_fields()[0];parameters=exact_target_parameters();records=[]
    for q in (1005,10005):
        bound=bounds(iv,previous,q)
        for u,nu in ((.0125,.85),(.025,.9),(.0375,.95)):
            mu=math.floor(q*nu)+.5;data=geometry(function,parameters,u,mu/q,1024,80)
            state=trial(data,q);p=state['probability'];W=data['W'];C=data['window'];Ct=data['window_theta']
            t=data['tilt'];ts=data['tilt_s'];drift=data['limit']['peak_drift_D']
            beta=C*np.arctan(t);betas=C*ts/(1+t*t);betat=Ct*np.arctan(t)
            _,_,angular=local_first(function,parameters,data['radius'],data['theta'])
            logWt=angular[1]-angular[0];f=np.cos(beta)/W
            f_s=(-np.sin(beta)*betas-np.cos(beta)*data['log_W_s'])/W
            f_t=(-np.sin(beta)*betat-np.cos(beta)*logWt)/W
            lam=1/np.sqrt(1+t*t);lams=-t*ts/(1+t*t)**1.5
            P=C*lams*data['log_profile_per_q']+(1+C*(lam-1))*data['log_profile_s_per_q']
            zeta=C*drift;A=P+zeta*state['log_t']/q
            d=data['log_W_s']+Ct*drift+zeta*logWt
            cov=lambda a,b:float(p@(a*b)-(p@a)*(p@b))
            transported=float(p@(f_s+zeta*f_t)+2*q*cov(A,f)+cov(d,f))
            r1=norm(data,state['r1']);variance=float(p@((f-p@f)**2)+p@(np.sin(beta)**2/W**2))
            records.append(dict(q=q,u=u,mu=mu,r1=r1,velocity=float(state['velocity']),velocity_derivative=float(state['velocity_s']),
                r1_variance_identity_error=abs(r1*r1-variance),transport_identity_error=abs(transported-state['velocity_s']),
                bounds_enclose=bool(r1<=bound['r1_upper'] and state['velocity']>=bound['velocity_lower'] and abs(state['velocity_s'])<=bound['velocity_derivative_upper_D'] and max(abs(f_t))<=bound['f_theta_lipschitz_upper'])))
    return records


def analytic_controls():
    points=np.linspace(0,.05,1001);chi,chit=cutoff(points);eta,etat=cutoff(points,True)
    gaussian=[]
    for lam in (.1,1.,10.):
        spectral=quad(lambda p:math.exp(-p*p/lam**2)/(2*math.pi),-np.inf,np.inf,epsabs=1e-12)[0]
        spatial=quad(lambda x:lam**2*math.exp(-lam**2*x*x)/(2*math.pi),-np.inf,np.inf,epsabs=1e-12)[0]
        exact=lam/(2*math.sqrt(math.pi));gaussian.extend((abs(spectral/exact-1),abs(spatial/exact-1)))
    spectral_errors=[]
    for y in (.1,1.,10.):
        eigenvalues=np.linspace(-50,50,501);weights=np.exp(-eigenvalues**2/100);weights/=weights.sum()
        a,b=.17,3.2
        hg2=float(weights@(eigenvalues**2/(eigenvalues**2+y*y)))
        g2=float(weights@(1/(eigenvalues**2+y*y)))
        spectral_errors.append(abs(hg2+y*y*g2-1))
        if a*math.sqrt(hg2)+b*math.sqrt(g2)>math.sqrt(a*a+(b/y)**2)*(1+1e-14):raise ArithmeticError('Spectral graph combination failed.')
    counts=[]
    for q in (1000,1005,10005):
        brute=sum(17*q<=20*n+10<=19*q for n in range(q))
        counts.append(brute==sector_count(q)['count'])
    return dict(cutoff_overlap_error=float(max(abs(chi*chi*eta-eta))),cutoff_derivative_max_sample=float(max(abs(chit))),
        inner_eta_squared_integral=float(np.trapezoid(eta*eta,points)),gaussian_HS_relative_errors=gaussian,
        spectral_graph_identity_errors=spectral_errors,exact_sector_counts=counts)


def run():
    path=project_path('dihole_projected_gap_results.json');previous_gap=json.loads(path.read_text())
    inputs=dict(previous_gap['input_sha256']);inputs.update({path.name:sha(path),'dihole_projected_gap.py':previous_gap['source_sha256']})
    previous=json.loads(project_path('dihole_channel_bounds_results.json').read_text())
    atom=json.loads(project_path('atomic_aperture_budget_results.json').read_text())
    if not previous_gap['passed'] or not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact source/input provenance required.')
    q_continuous=atom['geometry']['continuous_flux_requirement_at_illustrative_alpha'];q_integer=int(q_continuous)
    iv=MPIntervalContext();iv.dps=60
    rows=[bounds(iv,previous,q) for q in (1000,1005,10005,10**8,10**12,10**16,q_integer)]
    fields=field_controls(iv,previous);controls=analytic_controls()
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        actual_r1_variance_identity=max(x['r1_variance_identity_error'] for x in fields)<1e-16,
        actual_velocity_transport_identity=max(x['transport_identity_error'] for x in fields)<2e-9,
        actual_field_quantities_enclosed=all(x['bounds_enclose'] for x in fields),
        positive_velocity=all(x['velocity_lower']>0 for x in rows),
        improved_momentum_residual=all(x['r1_upper']<previous['uniform_bound']['constants']['r1_bound'][1] for x in rows),
        nested_cutoffs=controls['cutoff_overlap_error']==0 and controls['cutoff_derivative_max_sample']<=800,
        gaussian_trace_normalization=max(controls['gaussian_HS_relative_errors'])<1e-11,
        exact_spectral_graph_identity=max(controls['spectral_graph_identity_errors'])<1e-14,
        exact_integer_band_count=all(controls['exact_sector_counts']),
        finite_trace_bounds=all(math.isfinite(y['band_trace_error_upper_over_D']) and y['band_trace_error_upper_over_D']>0 for x in rows for y in x['frequency_rows']),
        atomic_scaled_resolvent_error_below_one=all(y['scaled_resolvent_error_upper']<1 for y in rows[-1]['frequency_rows']))
    tests={key:bool(value) for key,value in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=rows,field_controls=fields,analytic_controls=controls,
        continuous_atomic_flux_estimate=q_continuous,chosen_integer_near_estimate=q_integer,
        source_sha256=sha(__file__),input_sha256=inputs,total_source=None,metric_trace_derivative_bound=None,atom_transit=None,
        cutoff_definitions=dict(outer_support=[.005,.045],outer_ramp_width=.01,inner_support=[.015,.035],inner_ramp_width=.005),
        scope=['Uniform collar radial comparison for q integer>=1000 and all half-integer labels in [.85q,.95q], not a frozen-carrier calculation.',
            'H is a self-adjoint global realization agreeing with the collar and containing its regular compact core. h is a positive complete-line velocity extension; reflecting/circle completions are not included.',
            'Scaled resolvent error is y*||G_H J-J G_h||, normalized to the universal 1/y scale, not an actual relative error or probability.',
            'Finite observable Tr[(JX)*G_H(JX)-X*G_h X], X=eta exp(-h^2/(2Lambda^2)); trace-class error bound, not its measured value or a stress trace.',
            'Integer near atomic estimate is used only to evaluate bounds and exact sector count; no huge-flux wave function, supporting geometry or atomic trajectory is simulated.',
            'Metric derivative, frequency/filter removal, omitted channels and renormalization still prevent a total supporting-source conclusion.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Collar resolvent controls failed; evidence retained.')
    plot(result,previous,iv)
    print(json.dumps({'passed':result['passed'],'tests':tests,'atomic_bound':rows[-1]},indent=2))


def plot(result,previous,iv):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    qs=[int(q) for q in np.geomspace(1000,1e25,160)];rows=[bounds(iv,previous,q) for q in qs]
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    axes[0].loglog(qs,[x['frequency_rows'][1]['scaled_resolvent_error_upper'] for x in rows],lw=2,color='#187c80')
    axes[0].axhline(1,color='gray',ls='--');axes[0].set(title='Full radial operator: localized channel comparison',ylabel='Upper bound on y times the resolvent error')
    axes[1].loglog(qs,[x['frequency_rows'][1]['per_sector_trace_error_upper_over_D'] for x in rows],lw=2,color='#7840aa')
    axes[1].set(title='Selected Gaussian-filtered observable',ylabel='Per-sector trace-error bound / D')
    for ax in axes:ax.set_xlabel('Integer flux q (bounding formulas)');ax.grid(alpha=.2)
    fig.suptitle('Radial coupling retained, with a full-operator comparison bound')
    fig.supxlabel('y D = filter width times D = sqrt(q). Universal norm scale, not a passage probability or supporting-stress error.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
