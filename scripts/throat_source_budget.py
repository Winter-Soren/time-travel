"""Compare saved nonlinear source states with explicit static throat requirements.

No evolution is run. Local algebraic matches are not global wormhole solutions.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from nonlinear_packet_evolution import RadialModel
from project_paths import project_path


def throat(P, radius, flare):
    if not (P>0 and radius>0 and flare>0):
        raise ValueError('Require positive Planck coefficient, radius and flare.')
    scale=P/radius**2
    return dict(radius_times_cutoff=radius,flare=flare,density=(1-flare)*scale,
        radial_pressure=-scale,constant_redshift_tangential_pressure=flare*scale/2,
        radial_null_stress=-flare*scale,proper_radius_second_derivative=flare/(2*radius),
        isotropic_required_redshift_curvature=-(1+flare/2)/radius**2)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.3),layout='constrained')
    available=result['saved_source']['simultaneous_negative_null_magnitude']
    required=-result['benchmarks'][0]['radial_null_stress']
    axes[0].barh(['Strongest saved source','Ellis throat at packet width'],[available,required],color=['#205e80','#b85432'])
    axes[0].set(xscale='log',xlabel='Required / available negative null stress magnitude (cutoff^4)',
                title=f"Specified throat: {required/available:,.0f} times stronger source needed")
    source=result['saved_source'];target=result['benchmarks'][0]
    x=np.arange(3);width=.34
    axes[1].bar(x-width/2,[source['total_density'],source['pressure'],source['pressure']],width,label='Saved source at its peak')
    axes[1].bar(x+width/2,[target['density'],target['radial_pressure'],target['constant_redshift_tangential_pressure']],width,label='Ellis throat requirement')
    axes[1].axhline(0,color='black',lw=.7)
    axes[1].set(xticks=x,xticklabels=['Density','Radial pressure','Tangential pressure'],
                ylabel='Stress / cutoff^4',title='The full stress pattern also differs')
    axes[1].legend(fontsize=8)
    fig.suptitle('Saved nonlinear source versus one explicit throat benchmark',fontsize=13)
    fig.supxlabel('Static finite-redshift comparison only. A local match or a shallow flare is not a constructed two-ended geometry.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('throat_source_budget_progress.'+extension),dpi=160)
    plt.close(fig)


def run():
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    prior_path=project_path('nonlinear_packet_evolution_results.json')
    prior=json.loads(prior_path.read_text(encoding='utf-8'))
    state_path=project_path('nonlinear_packet_evolution_state.npz')
    if (not prior['passed'] or sha(state_path)!=prior['state_sha256']
            or sha(project_path('nonlinear_packet_evolution.py'))!=prior['source_sha256']):
        raise ArithmeticError('Saved nonlinear evolution provenance failed.')
    model=RadialModel(prior['parameters']['modes'],prior['parameters']['xmax'])
    b=model.base;P,R,H=b['P'],b['radius'],b['H']
    with np.load(state_path) as state:
        strength=np.maximum(0,np.minimum(-state['null_plus'],-state['null_minus']))
        index=np.unravel_index(np.argmax(strength),strength.shape)
        available=float(strength[index])
        if not available>0:
            raise ArithmeticError('Saved states contain no simultaneous radial negative null stress.')
        mean=(state['null_plus']+state['null_minus'])/2
        rho=3*P*H*H+state['density_excess']
        pressure=mean-rho
        source=dict(simultaneous_negative_null_magnitude=available,time_Ht=float(state['tau'][index[0]]),
            radius_over_initial_width=float(state['x'][index[1]]),total_density=float(rho[index]),
            pressure=float(pressure[index]),flux=float(state['flux'][index]),
            background_density=3*P*H*H,scope='Maximum over saved samples, not a rigorous continuum supremum.')
        z=state['coefficients'][:,:model.modes]@model.B.T
        areal=R*state['x'][None,:]*np.exp(state['tau'][:,None]+z)
        local_flare=strength*areal*areal/P
        local_index=np.unravel_index(np.argmax(local_flare),local_flare.shape)
        same_radius=dict(largest_saved_nec_flare_budget=float(local_flare[local_index]),
            time_Ht=float(state['tau'][local_index[0]]),radius_over_initial_width=float(state['x'][local_index[1]]),
            areal_radius_times_cutoff=float(areal[local_index]),
            scope='Necessary static NEC magnitude comparison at the same areal radius; no throat or full stress match.')
        fields=[model.fields(t,v) for t,v in zip(state['tau'],state['coefficients'])]
        ricci=(rho-3*pressure)/P
        scales=dict(H_over_cutoff=H,inverse_initial_width_over_cutoff=1/R,
            maximum_retained_local_wavenumber_over_cutoff=float(np.max(model.k[-1]/(R*np.exp(state['tau'][:,None]+z)))),
            maximum_trace_extrinsic_curvature_over_cutoff=max(float(np.max(abs(3*H+f['dk']))) for f in fields),
            maximum_sqrt_abs_Ricci_scalar_over_cutoff=float(np.sqrt(np.max(abs(ricci)))),
            minimum_X=min(float(np.min(b['X']*(1+f['dy']))) for f in fields),
            maximum_X=max(float(np.max(b['X']*(1+f['dy']))) for f in fields),
            true_strong_coupling_cutoff_established=False,
            scope='Derivative indicators relative to the nominal Lambda parameter; not bounds on all omitted operators or all curvature components.')
    benchmarks=[]
    for name,flare in [('Ellis, constant redshift',2.),('Unit flare',1.),('Weak flare',1e-4)]:
        row=throat(P,R,flare);row.update(name=name,required_to_saved_null_ratio=-row['radial_null_stress']/available)
        benchmarks.append(row)
    rmatch=math.sqrt(-P/source['pressure'])
    delta=(source['total_density']+source['pressure'])/source['pressure']
    matched=throat(P,rmatch,delta)
    matched['scope']='Algebraic pointwise matching of rho and pr at a hypothetical static throat. The saved center has zero areal radius, not rmatch; time evolution, topology, flux away from center and global matching remain unsolved.'
    # Stationary clock-aligned subclass: beta=Theta=0, so the higher-derivative
    # source vanishes. These bounds do not cover the evolving source above.
    vacuum=b['rho'];v=4*vacuum-2*math.sqrt(vacuum*(4*vacuum-1));X=1-v
    static=dict(minimum_total_density=vacuum-1/3,largest_negative_null_magnitude=.25,
        maximum_local_flare=X*(1-X)/(vacuum-(X-1)**2/4),X_at_maximum_flare=X,
        radius_at_maximum_flare_times_cutoff=math.sqrt(P/(vacuum-(X-1)**2/4)),
        time_kinetic_factor_at_maximum=(3*X-1)/2,gradient_factor_at_maximum=(X-1)/2,
        scope='Only strictly stationary, clock-aligned beta=Theta=0 configurations of this quartic action and this fixed vacuum term; no global static solution or stability claim.')
    # NIST CODATA2022 central values. Reduced Planck mass = Planck mass/sqrt(8pi).
    lp,tp,planck_GeV=1.616255e-35,5.391247e-44,1.220890e19
    length_unit=math.sqrt(8*math.pi*P)*lp
    time_unit=math.sqrt(8*math.pi*P)*tp
    physical=dict(cutoff_over_reduced_planck=1/math.sqrt(P),cutoff_GeV=planck_GeV/math.sqrt(8*math.pi*P),
        inverse_cutoff_length_m=length_unit,inverse_cutoff_time_s=time_unit,
        initial_packet_width_m=R*length_unit,simulated_duration_s=.5/H*time_unit,
        pressure_matched_radius_m=rmatch*length_unit,
        constants_source='https://physics.nist.gov/cuu/pdf/all.pdf',relative_planck_unit_uncertainty=1.1e-5,
        scope='Conversion for the actual chosen cutoff/Planck ratio. Changing Lambda changes P and the dimensionless solution; this is not a scalable material blueprint.')
    checks=dict(saved_provenance=True,peak_flux_free=abs(source['flux'])<1e-18,
        ellis_tensor_identity=math.isclose(benchmarks[0]['density'],benchmarks[0]['radial_pressure'])
                              and math.isclose(benchmarks[0]['constant_redshift_tangential_pressure'],-benchmarks[0]['density']),
        pointwise_density_match=math.isclose(matched['density'],source['total_density'],rel_tol=1e-12),
        pointwise_pressure_match=math.isclose(matched['radial_pressure'],source['pressure'],rel_tol=1e-12),
        pointwise_nec_match=math.isclose(matched['radial_null_stress'],-available,rel_tol=1e-10),
        saved_source_insufficient_for_ellis=benchmarks[0]['required_to_saved_null_ratio']>1e5,
        static_subclass_ellis_density_obstruction=static['minimum_total_density']>0,
        static_flare_extremum=abs(v*v-8*vacuum*v+4*vacuum)<1e-14 and 1/3<X<1,
        conversion_consistency=math.isclose(length_unit/time_unit,299792458.,rel_tol=3e-7),
        original_width_and_time_are_microscopic=physical['initial_packet_width_m']<1e-29 and physical['simulated_duration_s']<1e-37)
    result=dict(complete=True,passed=all(checks.values()),checks=checks,saved_source=source,benchmarks=benchmarks,
        largest_flare_from_saved_strength_at_packet_width=available*R*R/P,
        ellis_radius_if_strength_could_be_relocated=math.sqrt(2*P/available),
        same_radius_nec_comparison=same_radius,pointwise_isotropic_match=matched,
        stationary_clock_aligned_subclass=static,nominal_EFT_indicators=scales,physical_units=physical,
        source_sha256=sha(__file__),inputs_sha256={name:sha(project_path(name)) for name in
            ('nonlinear_packet_evolution.py','nonlinear_packet_evolution_results.json','nonlinear_packet_evolution_state.npz',
             'nonlinear_packet_constraints.py','expanding_field_completion.py','expanding_field_completion_results.json','project_paths.py')},
        calculation_type='Analysis and field reconstruction of saved states only; no evolution or refinement rerun.',
        decision='Reject direct use of the saved packet as support for the specified Ellis throat at the packet width. Very weak static flare can pass local algebra, but no global two-ended or dynamic throat solution is supplied.')
    project_path('throat_source_budget_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Source-budget checks failed; evidence saved.')
    plot(result)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    run()
