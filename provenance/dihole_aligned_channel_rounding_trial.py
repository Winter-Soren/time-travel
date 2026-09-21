"""Pole-regular field-following trial in the actual collar Dirac operator.

The full finite-q residual and its paired metric derivative are measured.
This is one supporting-field state per flux, not atom passage or a trace bound.
"""
import json
import math

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from balanced_dihole import symbolic_fields
from dihole_local_source import mask_jets
from dihole_slice_channel import channel, local_fields, localization, prepare, profile_logs, sha
from dihole_stretched_collar import exact_target_parameters, radius_profile, smooth_step
from project_paths import project_path


NAME = 'dihole_aligned_channel'
FLUXES = (5, 25, 105, 1005, 10005)
SX = np.array([[0., 1.], [1., 0.]])
ISY = np.array([[0., 1.], [-1., 0.]])
SZ = np.diag([1., -1.])


def polar_window(theta, edge=.15, width=.2):
    """C-infinity window: exact frozen profile and spin near both poles."""
    left, dl = smooth_step((theta-edge)/width)
    right, dr = smooth_step((np.pi-edge-theta)/width)
    return left*right, (dl*right-left*dr)/width


def geometry(function, p, u, nu, order=1024, inner=80):
    data = prepare(function, p, u, nu, order, inner)
    limit = localization(function, p, u, nu)
    theta = limit['theta_peak_limit']
    if polar_window(np.array([theta]))[0][0] != 1:
        raise ValueError('The field-line center must lie inside the window plateau.')
    radius, rate, second = radius_profile(u, order=128)
    fields = local_fields(function, p, radius, theta)
    tiny = 1e-25j
    rr = np.imag(local_fields(function, p, radius+tiny, theta)['ar'])/tiny.imag
    rt = np.imag(local_fields(function, p, radius, theta+tiny)['ar'])/tiny.imag
    tt = np.imag(local_fields(function, p, radius, theta+tiny)['at'])/tiny.imag
    drift = limit['peak_drift_D']
    acceleration = -(rr*rate**2+fields['ar']*second+2*rt*rate*drift+tt*drift**2)/fields['at']
    data['tilt'], data['tilt_s'] = radius*drift, rate*drift+radius*acceleration
    # Center L before a theta-dependent width change: no arbitrary log constant.
    center = profile_logs(function, p, u, nu, np.array([theta]), inner)
    data['log_center'] = center['log_profile_per_q'][0]
    data['log_profile_per_q'] -= center['log_profile_per_q'][0]
    data['log_profile_s_per_q'] -= center['log_profile_s_per_q'][0]
    data['window'], data['window_theta'] = polar_window(data['theta'])
    data['limit'] = limit
    data['mask'] = mask_jets(u) if 0 < u < .05 else np.zeros(6)
    return data


def trial(data, q, epsilon=0.):
    """Analytic epsilon family in one fixed canonical Hamiltonian space."""
    mu = q*data['nu']
    if type(q) is not int or q < 1 or not (.5 <= mu <= q-.5) or not math.isclose(mu%1, .5, abs_tol=1e-11):
        raise ValueError('Require a regular half-integer zero-mode sector at integer flux.')
    mask, mask_s = data['mask'][:2]
    stretch = np.exp(epsilon*mask)
    tilt = stretch*data['tilt']
    tilt_s = stretch*(data['tilt_s']+epsilon*mask_s*data['tilt'])
    beta0 = np.arctan(tilt)
    beta0_s = tilt_s/(1+tilt**2)
    lam0 = 1/np.sqrt(1+tilt**2)
    lam0_s = -tilt*tilt_s/(1+tilt**2)**1.5
    window, window_t = data['window'], data['window_theta']
    lam, lam_s, lam_t = 1+window*(lam0-1), window*lam0_s, window_t*(lam0-1)
    beta, beta_s, beta_t = window*beta0, window*beta0_s, window_t*beta0
    spin = np.array([np.cos(beta/2), np.sin(beta/2)])
    perpendicular = np.array([-spin[1], spin[0]])/2
    spin_s, spin_t = perpendicular*beta_s, perpendicular*beta_t
    log_profile = q*lam*data['log_profile_per_q']
    amplitude = np.exp(log_profile-np.max(np.real(log_profile)))
    weights, W = data['weights'], data['W']
    profile = amplitude/np.sqrt(weights@(W*amplitude**2))
    probability = weights*W*profile**2
    log_s = q*(lam_s*data['log_profile_per_q']+lam*data['log_profile_s_per_q'])
    ell = log_s-probability@(log_s+data['log_W_s']/2)
    log_t = q*(lam_t*data['log_profile_per_q']+lam*data['angular_log_derivative_per_q'])
    velocity = stretch**2*(probability@(np.cos(beta)/W))
    velocity_s = stretch**2*(probability@((2*(ell+epsilon*mask_s)*np.cos(beta)-beta_s*np.sin(beta))/W))
    r1_coefficient = stretch**2*(SZ@spin)/W-velocity*spin
    radial = stretch**2*(SZ@((ell+epsilon*mask_s)*spin+spin_s))/W
    angular = stretch*(SX@(log_t*spin+spin_t)+q*data['angular_log_derivative_per_q']*(ISY@spin))/(data['radius']*W)
    r0_coefficient = radial+angular-velocity_s*spin/2
    scalar = np.sqrt(W)*profile
    return dict(chi=scalar*spin, chi_s=scalar*((ell+data['log_W_s']/2)*spin+spin_s),
        r1=scalar*r1_coefficient, r0=scalar*r0_coefficient, probability=probability,
        velocity=velocity, velocity_s=velocity_s, window_mass=np.sum(probability[window < 1]),
        profile=profile, ell=ell, spin=spin, spin_s=spin_s, spin_t=spin_t, log_t=log_t,
        radial_action=scalar*radial, angular_action=scalar*angular)


def norm(data, vector):
    return float(np.sqrt(data['weights']@np.sum(abs(vector)**2, axis=0)))


def measure(data, q):
    state = trial(data, q)
    varied = trial(data, q, 1e-25j)
    frozen, _ = channel(data, q)
    row = dict(q=q, mu=q*data['nu'], u_over_D=data['u'], order=data['order'],
        inner_order=data['inner_order'], tilt=float(data['tilt']), tilt_s_D=float(data['tilt_s']),
        normalization_error=abs(float(np.sum(state['probability']))-1),
        projection_error=max(abs(float(data['weights']@np.sum(state['chi']*state[k], axis=0))) for k in ('r0','r1')),
        velocity=float(state['velocity']), velocity_s_D=float(state['velocity_s']),
        r0_norm_D=norm(data,state['r0']), r1_norm=norm(data,state['r1']),
        dot_r0_norm_D=norm(data,np.imag(varied['r0'])/1e-25), dot_r1_norm=norm(data,np.imag(varied['r1'])/1e-25),
        frozen_r0_norm_D=frozen['r0_norm_D'], frozen_r1_norm=frozen['r1_norm'],
        frozen_dot_r0_norm_D=frozen['dot_r0_norm_D'], frozen_dot_r1_norm=frozen['dot_r1_norm'],
        polar_window_probability=float(state['window_mass']))
    row['r0_to_frozen_ratio'] = row['r0_norm_D']/frozen['r0_norm_D']
    row['dot_r0_to_frozen_ratio'] = row['dot_r0_norm_D']/frozen['dot_r0_norm_D']
    row['local_gap_ratio'] = row['r0_norm_D']/(math.sqrt(q)*data['limit']['leading_transverse_gap_over_sqrt_q'])
    # For a local carrier exp(ik s), real residual vectors have no cross term.
    row['local_carriers'] = [dict(k_D=k,
        residual_D=math.hypot(row['r0_norm_D'],k*row['r1_norm']),
        frozen_residual_D=math.hypot(row['frozen_r0_norm_D'],k*row['frozen_r1_norm']),
        metric_derivative_D=math.hypot(row['dot_r0_norm_D'],k*row['dot_r1_norm']),
        frozen_metric_derivative_D=math.hypot(row['frozen_dot_r0_norm_D'],k*row['frozen_dot_r1_norm'])) for k in (0.,1.,10.,100.)]
    return row


def controls(function, p, data, q):
    controls = []
    for epsilon in (0., .03):
        state = trial(data,q,epsilon)
        for step in (1e-5,5e-6):
            pairs = [trial(geometry(function,p,data['u']+sign*step,data['nu'],data['order'],data['inner_order']),q,epsilon) for sign in (-1,1)]
            derivative = (pairs[1]['chi']-pairs[0]['chi'])/(2*step)
            speed = (pairs[1]['velocity']-pairs[0]['velocity'])/(2*step)
            controls.append(dict(epsilon=epsilon, step=step,
                chi_s_relative_error=norm(data,derivative-state['chi_s'])/norm(data,state['chi_s']),
                velocity_s_relative_error=float(abs(speed/state['velocity_s']-1))))
    dot = {key:np.imag(trial(data,q,1e-25j)[key])/1e-25 for key in ('r0','r1')}
    metric = []
    for step in (1e-4,5e-5):
        states = [trial(data,q,sign*step) for sign in (-1,1)]
        metric.append(dict(step=step, **{key+'_relative_error':norm(data,(states[1][key]-states[0][key])/(2*step)-dot[key])/norm(data,dot[key]) for key in dot}))
    # Direct fixed-space operator action on chi, independently differentiated in s.
    state = trial(data,q,.03)
    step=5e-6
    states=[trial(geometry(function,p,data['u']+sign*step,data['nu'],data['order'],data['inner_order']),q,.03) for sign in (-1,1)]
    # W^-1/2 chi = g eta; its derivative contains the entire normalization/spin change.
    gs=((states[1]['profile']*states[1]['spin'])-(states[0]['profile']*states[0]['spin']))/(2*step)
    F, Fs=data['mask'][:2]; stretch=np.exp(.03*F)
    radial=stretch**2*(SZ@(gs+.03*Fs*state['profile']*state['spin']))/np.sqrt(data['W'])
    direct=radial+state['angular_action']-state['velocity_s']*state['chi']/2
    # Angular action from independently shifted profiles; hold the norm fixed.
    tilt=stretch*data['tilt']; lam0=1/np.sqrt(1+tilt**2); beta0=np.arctan(tilt)
    lam=1+data['window']*(lam0-1)
    normalization=np.sqrt(data['weights']@(data['W']*np.exp(2*q*lam*data['log_profile_per_q'])))
    angular_checks=[]
    for h in (2e-6,1e-6):
        step=h*np.minimum(1,np.minimum(data['theta'],np.pi-data['theta']))
        shifted=[]
        for sign in (-1,1):
            theta=data['theta']+sign*step
            window,_=polar_window(theta)
            logs=profile_logs(function,p,data['u'],data['nu'],theta,data['inner_order'])
            log_profile=logs['log_profile_per_q']-data['log_center']
            profile=np.exp(q*(1+window*(lam0-1))*log_profile)/normalization
            beta=window*beta0
            shifted.append(profile*np.array([np.cos(beta/2),np.sin(beta/2)]))
        derivative=(shifted[1]-shifted[0])/(2*step)
        exact=state['profile']*(state['log_t']*state['spin']+state['spin_t'])
        error=stretch*(SX@(derivative-exact))/(data['radius']*np.sqrt(data['W']))
        angular_checks.append(dict(step=h,action_absolute_error=norm(data,error)))
    return dict(radial_finite_differences=controls,metric_finite_differences=metric,
        direct_radial_action_absolute_error=norm(data,direct-state['r0']),angular_finite_differences=angular_checks)


def anchored_log_bound(inputs):
    """A global unnormalized envelope from the saved directed collar bounds."""
    path=project_path('dihole_strong_field_geometry_results.json')
    previous=json.loads(path.read_text(encoding='utf-8'))
    hashes=dict(previous['input_sha256'])
    hashes.update({'dihole_strong_field_geometry.py':previous['source_sha256'],path.name:sha(path)})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in hashes.items()):
        raise ArithmeticError('Exact collar interval provenance required.')
    inputs.update(hashes)
    bounds=previous['interval_zero_exclusion']['records'][-1]
    if not bounds['passed']: raise ArithmeticError('Positive closed-collar bounds required.')
    iv=MPIntervalContext();iv.dps=70
    max_Rd=iv.mpf(1)/20*iv.mpf(21)/20
    n=(iv.mpf(bounds['potential_denominator_lower'])/max_Rd)**2
    m=iv.mpf(-bounds['A_mu_enclosure'][1])
    coefficient=n*m/4
    down=lambda x:float(np.nextafter(float(x.a),-np.inf))
    return dict(NP_lower=down(n),a_theta_over_sin_theta_lower=down(m),log_envelope_coefficient=down(coefficient),
        identity='NP=V^2/(R distance)^2; L(theta0)=0, L_theta=NP(nu-a)/sin(theta).',
        bound='L(theta) <= -log_envelope_coefficient*(theta-theta0)^2 on the entire closed collar and 0<theta<pi for any interior flux center.',
        proof='a_theta>=m sin(theta), NP>=n; concavity of sin gives abs(cos(theta)-cos(theta0))>=sin(theta)*abs(theta-theta0)/2. Integrate L_theta on either side.',
        scope='Unnormalized logarithmic envelope only. A normalized moment bound also needs a local lower envelope and weight bounds; no trial-residual or source bound follows alone.')


def run():
    previous_path=project_path('dihole_slice_channel_results.json')
    previous=json.loads(previous_path.read_text(encoding='utf-8'))
    inputs=dict(previous['input_sha256'])
    inputs.update({'dihole_slice_channel.py':previous['source_sha256'],previous_path.name:sha(previous_path)})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact previous source/input provenance required.')
    function,p=symbolic_fields()[0],exact_target_parameters()
    envelope=anchored_log_bound(inputs)
    rows=[]; finest={}
    for order,inner in ((512,48),(1024,80)):
        for u in (.0125,.025,.0375):
            data=geometry(function,p,u,.9,order,inner)
            rows.extend(measure(data,q) for q in FLUXES)
            if order==1024: finest[u]=data
    metrics=('r0_norm_D','r1_norm','dot_r0_norm_D','dot_r1_norm','velocity')
    refinement={k:max(abs(a[k]/b[k]-1) for a,b in zip(rows[:15],rows[15:])) for k in metrics}
    checked=controls(function,p,finest[.0125],105)
    central=geometry(function,p,0.,.5,512,80)
    aligned=trial(central,1); frozen,_=channel(central,1)
    center_errors=[abs(norm(central,aligned['r1'])-frozen['r1_norm']),norm(central,aligned['r0']),abs(aligned['velocity']-frozen['projected_velocity'])]
    pole=np.array([.01,.1,np.pi-.1,np.pi-.01])
    pole_window=polar_window(pole)
    last=[r for r in rows[15:] if r['q']==FLUXES[-1]]
    tests=dict(exact_current_input_hashes=all(sha(project_path(n))==h for n,h in inputs.items()),
        canonical_normalization=max(r['normalization_error'] for r in rows)<2e-15,
        projection_orthogonality=max(r['projection_error'] for r in rows)<2e-11,
        angular_and_inner_quadrature=max(refinement.values())<2e-6,
        normalized_radial_derivative=max(r['chi_s_relative_error'] for r in checked['radial_finite_differences'] if r['step']==5e-6)<2e-5,
        speed_radial_derivative=max(r['velocity_s_relative_error'] for r in checked['radial_finite_differences'] if r['step']==5e-6)<2e-5,
        common_metric_derivative=max(checked['metric_finite_differences'][-1][k+'_relative_error'] for k in ('r0','r1'))<2e-5,
        direct_action=checked['direct_radial_action_absolute_error']<2e-7,
        independent_angular_action=checked['angular_finite_differences'][-1]['action_absolute_error']<2e-7,
        exact_regular_poles=all(np.array_equal(x,np.zeros_like(x)) for x in pole_window),
        zero_tilt_cylinder_recovered=max(center_errors)<3e-13,
        positive_analytic_log_envelope=envelope['log_envelope_coefficient']>0 and all(np.max(data['log_profile_per_q']+envelope['log_envelope_coefficient']*(data['theta']-data['limit']['theta_peak_limit'])**2)<1e-12 for data in finest.values()),
        high_flux_residual_improved=max(r['r0_to_frozen_ratio'] for r in last)<1,
        high_flux_metric_residual_improved=max(r['dot_r0_to_frozen_ratio'] for r in last)<1,
        tested_high_flux_carriers_improved=all(c['residual_D']<c['frozen_residual_D'] and c['metric_derivative_D']<c['frozen_metric_derivative_D'] for r in last for c in r['local_carriers']))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=rows,quadrature_refinement=refinement,
        independent_controls=checked,central_cylinder_errors=[float(x) for x in center_errors],anchored_log_envelope=envelope,
        trial='g=Z^-1/2 exp(q lambda(theta,s)L), L centered at a(theta0)=nu; lambda=1+C(theta)(1/sqrt(1+t^2)-1), beta=C(theta)atan(t), t=R theta0_s; chi=sqrt(W)g(cos(beta/2),sin(beta/2)).',
        polar_window='C=S((theta-.15)/.2)S((pi-.15-theta)/.2), S is the existing C-infinity step. Exact frozen profile and spin on polar neighborhoods; all C_theta terms retained.',
        common_metric_family='H_e=exp(eF) Hrad exp(eF)+exp(eF) Hang; t_e=exp(eF)t. Both U_e and h_e are differentiated, including F_s.',
        scope=['One regular collar sector mu/q=.9; no sector sum, exterior propagation, or supporting-source trace.',
            'The complete finite-q trial residual and its metric derivative are evaluated at three slices, including pole-restoring derivatives.',
            'Numerical refinement and derivative checks are sensitivity tests, not interval enclosures or uniform atomic-flux remainder bounds.',
            'Local carriers use exp(ik s), not normalizable packets; a general envelope requires both coefficients. The derivative coefficient r1 can increase and the small-flux trial can be worse.',
            'Local gap ratios use the previous oscillator approximation; no finite-q/global gap is proved.',
            'No stable opening, source material, hydrogen survival, actual atom passage, chronology change or future clock advantage is established.'],
        full_source_bound=None,atomic_transit=None,source_sha256=sha(__file__),input_sha256=inputs)
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']: raise ArithmeticError('Aligned trial controls failed; results retained.')
    plot(result)
    print(json.dumps({k:result[k] for k in ('passed','tests','quadrature_refinement','independent_controls')},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.6),layout='constrained')
    for u in (.0125,.025,.0375):
        rows=[r for r in result['records'] if r['order']==1024 and r['u_over_D']==u]
        q=[r['q'] for r in rows]
        for ax,key in zip(axes,('r0_norm_D','dot_r0_norm_D')):
            line,=ax.loglog(q,[r[key] for r in rows],'o-',label=f'u / D = {u:.4f}')
            ax.loglog(q,[r['frozen_'+key] for r in rows],'--',color=line.get_color(),alpha=.6)
    axes[0].set(title='Residual at local carrier k = 0',ylabel='Norm of r0 times D')
    axes[1].set(title='Metric derivative at k = 0',ylabel='Norm of the variation of r0 times D')
    for ax in axes:
        ax.set_xlabel('Integer flux q');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Field-following trial with regular poles: solid; frozen slices: dashed')
    fig.supxlabel('General transport also needs r1. One state per flux; no supporting-source bound or atom passage.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
