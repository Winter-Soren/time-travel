"""Transfer the actual slice gap to the aligned complement, with its variation.

Frozen radial carriers are included. No global resolvent or stress trace is
claimed; the moving-complement derivative is distinguished from its embedding.
"""
import json
import math

import numpy as np
from scipy.linalg import expm, null_space

from balanced_dihole import symbolic_fields
from dihole_aligned_channel import geometry, norm, trial
from dihole_channel_bounds import MPIntervalContext, calibration, cover, field_function, lower, upper
from dihole_slice_channel import channel, sha
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


NAME='dihole_projected_gap'


def bounds(iv,previous,window,q):
    if not math.isfinite(q) or q<1000:raise ValueError('This bound requires finite q>=1000.')
    c=previous['uniform_bound']['constants'];u=previous['uniform_bound'];g=previous['frozen_slice_gap']
    tilt=iv.mpf(c['tilt_bound'][1]);F=iv.mpf(u['step_bounds']['F'])
    Mlocal=iv.mpf(c['local_lower_log_M'][1])
    Mwindow=iv.mpf(max(abs(x) for x in window['derived']['L_theta_theta']))/2
    C2=iv.mpf(u['even_moment_coefficients'][1]);radius=iv.mpf(1)/50
    exponent=iv.mpf(c['gaussian_exponent_lower'][0])*q*radius**2/2
    # For |delta|>=radius, split the Gaussian exponent in half before integrating.
    tail=iv.mpf(min(1.,upper(2**iv.mpf('2.5')*iv.exp(-exponent))))
    stat=iv.sqrt(C2*(Mlocal**2+(Mwindow**2-Mlocal**2)*tail))
    width=tilt**2/2;spin=iv.atan2(tilt,iv.mpf(1))/2
    distance=iv.sqrt((width*stat)**2+spin**2)
    # Compose after exporting the upper angle bound, as downstream users do.
    distance_upper=upper(distance);d=iv.mpf(distance_upper)
    retained_lower=lower(1-d**2);retained=iv.mpf(retained_lower)
    metric=iv.sqrt((F*tilt**2*stat)**2+(F*tilt/2)**2)
    gap=iv.sqrt(2*q*iv.mpf(g['kappa_lower'])-iv.mpf(g['log_W_curvature_upper']))/iv.mpf(g['Rmax_Wmax_upper'])
    delta=gap*retained;wmin=iv.mpf(c['W_lower'][0])
    if retained.a<=0:raise ArithmeticError('Projector distance bound does not transfer a gap.')
    tan=d/iv.sqrt(retained);Mdot=iv.mpf(upper(metric));base=iv.mpf(lower(delta))
    carriers=[]
    for fraction in (0.,.25,.5,.9):
        # Declared k range, not an actual longitudinal wave packet or evolution.
        k_limit=max(0.,lower(iv.mpf(fraction)*base*wmin));k=iv.mpf(k_limit)
        size=k/wmin;alpha=size/base;separation=base-size
        beta=(tan+alpha)/(1-alpha)
        for frequency_ratio in (0.,1.,10.):
            omega_value=max(0.,lower(iv.mpf(frequency_ratio)*base));omega=iv.mpf(omega_value)
            denominator=iv.sqrt(separation**2+omega**2)
            derivative=(F+2*beta*Mdot)/denominator+F*size/denominator**2
            embedded=derivative+2*Mdot/denominator
            carriers.append(dict(momentum_fraction=fraction,k_D_limit=k_limit,omega_D=omega_value,
                gap_lower_D=lower(separation),resolvent_norm_upper_over_D=upper(1/denominator),
                kato_resolvent_derivative_upper_over_D=upper(derivative),embedded_resolvent_derivative_upper_over_D=upper(embedded)))
    return dict(q=q,tail_fourth_moment_fraction_upper=upper(tail),tail_exponent_lower=lower(exponent),
        window_quadratic_coefficient_upper=upper(Mwindow),log_statistic_rms_upper=upper(stat),
        projector_distance_upper=distance_upper,gap_retained_fraction_lower=retained_lower,
        profile_metric_derivative_norm_upper=upper(metric),aligned_complement_gap_lower_D=lower(delta),
        carriers=carriers)


def profile_controls(iv,previous,window):
    function=symbolic_fields()[0];parameters=exact_target_parameters();rows=[]
    for q in (1005,10005):
        bound=bounds(iv,previous,window,q)
        for u,nu in ((.0125,.85),(.025,.9),(.0375,.95)):
            mu=math.floor(q*nu)+.5
            records=[]
            for order in (512,1024):
                data=geometry(function,parameters,u,mu/q,order,80)
                state=trial(data,q);_,frozen=channel(data,q)
                zero=np.array([np.sqrt(data['W'])*frozen['profile'],np.zeros_like(data['theta'])])
                overlap=float(data['weights']@np.sum(zero*state['chi'],axis=0))
                distance=norm(data,zero-overlap*state['chi'])
                derivative=norm(data,np.imag(trial(data,q,1e-25j)['chi'])/1e-25)
                stat=float(np.sqrt(state['probability']@((q*data['window']*data['log_profile_per_q'])**2)))
                records.append(dict(order=order,projector_distance=distance,metric_derivative_norm=derivative,log_statistic_rms=stat))
            rows.append(dict(q=q,u=u,mu=mu,records=records,
                bounds_enclose=all(x['projector_distance']<=bound['projector_distance_upper'] and x['metric_derivative_norm']<=bound['profile_metric_derivative_norm_upper'] and x['log_statistic_rms']<=bound['log_statistic_rms_upper'] for x in records),
                refinement_absolute_change=max(abs(records[0][key]-records[1][key]) for key in ('projector_distance','metric_derivative_norm','log_statistic_rms'))))
    return rows


def compression_controls():
    """Indefinite matrices test the compression, connection and embedding signs."""
    rng=np.random.default_rng(166);raw=rng.normal(size=(7,7))+1j*rng.normal(size=(7,7))
    unitary=np.linalg.qr(raw)[0];H=(unitary*np.arange(-3.,4.))@unitary.conj().T
    frozen=unitary[:,3];B=(raw+raw.conj().T)/2;B/=np.linalg.norm(B,2)
    rows=[]
    for angle in (.1,.6):
        a=np.cos(angle)*frozen+np.sin(angle)*unitary[:,2]
        Q=null_space(a.conj()[None,:]);projection=Q@Q.conj().T
        d=projection@(unitary[:,1]+.3*unitary[:,5]);d/=np.linalg.norm(d)
        K=np.outer(d,a.conj())-np.outer(a,d.conj());F=.7
        C=Q.conj().T@H@Q;b=Q.conj().T@H@a;v=Q.conj().T@frozen
        identity=np.linalg.solve(C,b)+v/np.vdot(a,frozen)
        for k in (0.,.1):
            Hk=H+k*B;Ck=Q.conj().T@Hk@Q;bk=Q.conj().T@Hk@a;dq=Q.conj().T@d
            derivative=F*Ck+F*k*(Q.conj().T@B@Q)-np.outer(bk,dq.conj())-np.outer(dq,bk.conj())
            def varied(epsilon):
                rotation=expm(epsilon*K)
                return Q.conj().T@rotation.conj().T@(np.exp(F*epsilon)*H+np.exp(2*F*epsilon)*k*B)@rotation@Q
            h=1e-6;fd=(varied(h)-varied(-h))/(2*h)
            z=.4j;G=np.linalg.inv(Ck-z*np.eye(len(Ck)));dotG=-G@derivative@G
            fp=(np.linalg.inv(varied(h)-z*np.eye(len(Ck)))-np.linalg.inv(varied(-h)-z*np.eye(len(Ck))))/(2*h)
            embedded=Q@G@Q.conj().T;embedded_dot=K@embedded-embedded@K+Q@dotG@Q.conj().T
            def full(e):
                U=expm(e*K);return U@Q@np.linalg.inv(varied(e)-z*np.eye(len(Ck)))@Q.conj().T@U.conj().T
            fq=(full(h)-full(-h))/(2*h)
            rows.append(dict(angle=angle,k=k,minimum_abs_eigenvalue=float(min(abs(np.linalg.eigvalsh(Ck)))),
                theoretical_gap_lower=float(np.cos(angle)**2-abs(k)),
                inverse_mixing_identity_error=float(np.linalg.norm(identity)),
                connection_derivative_error=float(np.linalg.norm(fd-derivative,2)),
                resolvent_derivative_error=float(np.linalg.norm(fp-dotG,2)),
                embedded_derivative_error=float(np.linalg.norm(fq-embedded_dot,2))))
    return rows


def run():
    path=project_path('dihole_channel_bounds_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_channel_bounds.py':previous['source_sha256']})
    atom_path=project_path('atomic_aperture_budget_results.json');atom=json.loads(atom_path.read_text())
    inputs.update(atom['input_sha256']);inputs.update({atom_path.name:sha(atom_path),'atomic_aperture_budget.py':atom['source_sha256']})
    if not previous['passed'] or not atom['passed'] or not all(sha(project_path(n))==v for n,v in inputs.items()):raise ArithmeticError('Exact input provenance required.')
    iv=MPIntervalContext();iv.dps=60
    window=cover(iv,field_function(iv),calibration(iv),iv.mpf(3)/20,iv.pi-iv.mpf(3)/20)
    q_atom=atom['geometry']['continuous_flux_requirement_at_illustrative_alpha']
    rows=[bounds(iv,previous,window,q) for q in (1000,1005,10005,100000,1000000,10000000,100000000,q_atom)]
    profiles=profile_controls(iv,previous,window);matrices=compression_controls()
    pilot_path=project_path('dihole_projected_gap_parameter_trial.json');pilot=json.loads(pilot_path.read_text())
    pilot_source=project_path('provenance/dihole_projected_gap_parameter_trial.py')
    if sha(pilot_source)!=pilot['source_sha256']:raise ArithmeticError('Exact parameter-reporting pilot required.')
    tests=dict(exact_inputs=all(sha(project_path(n))==v for n,v in inputs.items()),
        finite_window_curvature=all(math.isfinite(v) for v in window['derived']['L_theta_theta']),
        nontrivial_projector_bound=all(0<x['projector_distance_upper']<1 for x in rows),
        positive_compressed_gap=all(x['aligned_complement_gap_lower_D']>0 for x in rows),
        original_profile_controls=all(x['bounds_enclose'] for x in profiles),
        profile_quadrature_sensitivity=max(x['refinement_absolute_change'] for x in profiles)<1e-7,
        indefinite_compression_gap=all(x['minimum_abs_eigenvalue']>=x['theoretical_gap_lower']-1e-12 for x in matrices),
        inverse_mixing_identity=max(x['inverse_mixing_identity_error'] for x in matrices)<1e-12,
        kato_connection_derivative=max(x['connection_derivative_error'] for x in matrices)<1e-8,
        resolvent_derivative=max(x['resolvent_derivative_error'] for x in matrices)<1e-8,
        embedded_projector_derivative=max(x['embedded_derivative_error'] for x in matrices)<1e-8,
        exact_zero_carrier_parameters=all(y['k_D_limit']==0 for x in rows for y in x['carriers'] if y['momentum_fraction']==0),
        all_declared_carrier_gaps_positive=all(y['gap_lower_D']>0 for x in rows for y in x['carriers']))
    result=dict(passed=all(tests.values()),tests=tests,records=rows,window_cover=window,profile_controls=profiles,matrix_controls=matrices,
        parameter_reporting_correction=dict(source='provenance/'+pilot_source.name,source_sha256=sha(pilot_source),results=pilot_path.name,results_sha256=sha(pilot_path),
            finding='Pilot exported rounded-up bounds on k and omega, rather than explicit chosen parameters; even exact zero appeared as a subnormal upper bound. Now choose a rounded-down momentum limit and explicit binary frequency before interval evaluation, and compose the gap from the exported projector-distance bound.',
            scope='Makes the declared all-|k|<=k_limit range and fixed frequency directly reproducible. No test tolerance or physical scope changed.'),
        atomic_flux_is_continuous_estimate=True,source_sha256=sha(__file__),input_sha256=inputs,
        source_trace_bound=None,global_resolvent_bound=None,atom_transit=None,
        scope=['Same collar and regular nu band as165, for q>=1000; moment interpolation covers the normalized path between frozen and aligned profiles.',
            'Compressed frozen transverse gap and imaginary-frequency resolvent with bounded radial carrier k; no s derivatives or exterior propagation are included in this frozen symbol.',
            'First paired metric derivative in Kato-identified complement spaces and separate embedded derivative including both outer projectors.',
            'A whole-space resolvent comparison needs a collar cutoff, its commutator, a self-adjoint longitudinal extension and source trace-ideal estimates. No global spectral gap follows.',
            'No total supporting stress, stable opening, new physical law, atomic survival or time shift is established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Projected gap checks failed; evidence retained.')
    plot(result,previous,iv)
    print(json.dumps({'passed':result['passed'],'tests':tests,'atomic_bounds':rows[-1]},indent=2))


def plot(result,previous,iv):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    q=np.geomspace(1000,1e9,140);rows=[bounds(iv,previous,result['window_cover'],float(x)) for x in q]
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    axes[0].semilogx(q,[x['gap_retained_fraction_lower'] for x in rows],color='#187c80',lw=2)
    axes[0].set(title='Gap transferred to the trial complement',ylabel='Guaranteed fraction of the frozen gap',ylim=(0,1.05))
    axes[1].loglog(q,[x['profile_metric_derivative_norm_upper'] for x in rows],color='#7840aa',lw=2)
    axes[1].set(title='Trial change under the paired metric variation',ylabel='Upper bound on the normalized-state derivative')
    for ax in axes:ax.set_xlabel('Flux q (bounding formulas)');ax.grid(alpha=.2)
    fig.suptitle('Supporting-field reduction: projector mismatch and its variation bounded')
    fig.supxlabel('Declared regular collar band. Frozen transverse operators; no global source, stable opening or hydrogen transit.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
