"""Complete N=4 finite-model infrared remainder, using the frozen104 proof.

Only eight new zero-frequency exterior systems are solved. The saved m=0
record is reused exactly. Absolute mouth n labels, the shifted free outer
labels, both physical Weyls and the declared reflected endpoint connection
are kept separately for every sector. No degeneracy multiplier substitutes
for the individual exterior responses.

Bounds use the actual pair tail to infinity at fixed J=9/GL64. They are
conditional finite-DG analytic enclosures evaluated in ordinary floating
arithmetic, not interval-certified or continuum/infinite-angular results.

Run: .venv\\Scripts\\python.exe -B -W error complete_threshold_bound.py
"""
from project_paths import project_path, source_matches
import copy
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import block_diag

from complete_angular_source import CompleteSector, port_labels, selected_map
from complete_angular_modes import mouth_channels, reflection_parities
from exterior_boundary_map import _complex_json
import threshold_return_bound as previous


N=4
J=9
RHO=32.
ELL=.2*math.hypot(1.,RHO)
KMAX=1e-6/ELL
DESTINATION=project_path('complete_threshold_bound_results.json')
INPUTS=tuple(dict.fromkeys(previous.INPUTS+('threshold_return_bound.py',
    'threshold_return_bound_results.json','complete_angular_modes.py')))
PROOF_FIELDS=('parameters','open_gram','open_right','open_left','tail','core',
    'core_controls','tail_green_control','probes','gram_conversion','phases',
    'matrices','scope','seconds','tests','passed','certificate_gate_passed')


def hashes():
    return {name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in INPUTS}


def exact_sector_core(mu,raw,pressure,core):
    """Independent port-label construction checks the selected full map."""
    ns=mouth_channels(mu,N);incoming,outgoing=port_labels(mu,N)
    zero=mu==.5;c=len(incoming);rows=[]
    maximum_map=maximum_pressure=maximum_trace=0.
    analytic_norms=np.array([4*sum(core['records'][int(n)-1]['probes'][j]['zero_derivative_norm'] for n in ns)
                             for j in range(2)])
    slope=np.array([4*sum(core['records'][int(n)-1]['probes'][j]['mixed_frequency_metric_norm_upper'] for n in ns)
                    for j in range(2)])
    for angle in (0.,math.pi):
        direct=np.zeros((c,c),complex);dp=np.zeros((2,c,c),complex)
        if zero:direct[incoming.index(('upper',0)),outgoing.index(('lower',0))]=np.exp(1j*angle)
        for n,parity in zip(ns,reflection_parities(mu,N)):
            n=int(n);item=core['records'][n-1]
            action=item['mass_action'];t=item['transmission'];reflection=-1j*math.tanh(action)
            iu,il=incoming.index(('upper',n)),incoming.index(('lower',n))
            ou,ol=outgoing.index(('upper',n)),outgoing.index(('lower',n))
            forward=parity*np.exp(1j*angle);reverse=parity*np.exp(-1j*angle)
            direct[iu,ou]=direct[il,ol]=reflection
            direct[iu,ol]=forward*t;direct[il,ou]=reverse*t
            for j in range(2):
                mass_weight=item['probes'][j]['J_f']
                dp[j,iu,ou]=dp[j,il,ol]=1j*mass_weight*t*t
                dp[j,iu,ol]=forward*mass_weight*t*math.tanh(action)
                dp[j,il,ou]=reverse*mass_weight*t*math.tanh(action)
        selected=selected_map(raw,mu,angle,N)
        selected_pressure=np.array([selected_map(v,mu,angle,N) for v in pressure])
        maximum_map=max(maximum_map,float(np.max(abs(direct-selected))))
        maximum_pressure=max(maximum_pressure,float(np.max(abs(dp-selected_pressure))))
        numerical_norms=[float(np.linalg.norm(block_diag(v,v.conj().T),ord='nuc')) for v in selected_pressure]
        maximum_trace=max(maximum_trace,float(np.max(abs(np.array(numerical_norms)-analytic_norms))))

        # Canonical cap map has the opposite input/output ordering from B.
        cap=np.zeros((c,c),complex)
        for n in ns:
            for port in ('upper','lower'):
                cap[outgoing.index((port,int(n))),incoming.index((port,int(n)))]=-1j
        reference=block_diag(cap,cap.conj().T)
        if zero:
            up=incoming.index(('upper',0));lo=outgoing.index(('lower',0))
            reference[c+up,up]=-1j;reference[lo,c+lo]=-1j
        B=block_diag(selected,selected.conj().T)
        gap=float(np.linalg.svd(np.eye(2*c)-reference@B,compute_uv=False)[-1])
        expected=math.sqrt(2) if zero else math.sqrt(2*(1+math.tanh(core['records'][int(ns[0])-1]['mass_action'])))
        rows.append(dict(phase=angle,reference_zero_gap=gap,exact_reference_gap=expected,
            reference_gap_error=abs(gap-expected),B_unitarity_error=float(np.max(abs(B.conj().T@B-np.eye(2*c))))))
    return dict(mouth_channels=ns.tolist(),incoming_labels=[list(v) for v in incoming],outgoing_labels=[list(v) for v in outgoing],
        mouth_dimension_per_weyl=c,has_zero_mode=zero,
        pressure_trace_norm_at_zero=analytic_norms.tolist(),paired_pressure_trace_norm_slope=slope.tolist(),
        direct_map_error=maximum_map,direct_pressure_error=maximum_pressure,
        pressure_trace_norm_error=maximum_trace,phase_controls=rows,
        passed=bool(maximum_map<1e-13 and maximum_pressure<1e-13 and maximum_trace<1e-13
            and max(v['reference_gap_error'] for v in rows)<1e-13
            and max(v['B_unitarity_error'] for v in rows)<1e-13))


def sector_case(m,raw,pressure,core):
    started=time.perf_counter();mu=m+.5
    sector=CompleteSector(mu,h=.05,nmax=N,outer_nmax=J,rho_c=RHO)
    selected=exact_sector_core(mu,raw,pressure,core)
    c=selected['mouth_dimension_per_weyl'];p=len(sector.config['outer_channels'])
    gram=previous.open_gram_check(sector)
    SR,right=previous.open_response(sector)
    SL,left=previous.open_response(sector,True)
    Z,tail=previous.zero_tail(mu=mu,J=J);enclosure=tail['fine']
    if not enclosure['passed']:
        return dict(m=m,mu=mu,config=sector.config,passed=False,certificate_gate_passed=False,
            failure='Inward tail enclosure failed; no sewn bound asserted.',tail=tail,
            sector_core=selected,open_gram=gram,open_right=right,open_left=left,
            tests=dict(sector_core=selected['passed'],inward_comparison=False),
            phases=[],seconds=time.perf_counter()-started)
    rf=sector.right_factors['outer'];lf=sector.left_factors['outer']
    ZR=previous.convert(rf[1],Z,rf[0]);ZL=previous.convert(lf[1],Z.conj().T,lf[0])
    gamma_core=max(float(np.linalg.norm(sector.Fin,2)*np.linalg.norm(np.linalg.inv(sector.Fout),2)),
                   float(np.linalg.norm(sector.Lin,2)*np.linalg.norm(np.linalg.inv(sector.Lout),2)))
    gamma_outer=max(float(np.linalg.norm(rf[1],2)*np.linalg.norm(np.linalg.inv(rf[0]),2)),
                    float(np.linalg.norm(lf[1],2)*np.linalg.norm(np.linalg.inv(lf[0]),2)))
    S=block_diag(SR,SL);eye=np.eye(len(S))
    lopen=max(right['response_lipschitz_upper'],left['response_lipschitz_upper'])
    mouth_indices=np.r_[np.arange(c),c+p+np.arange(c)]
    outer_indices=np.r_[c+np.arange(p),2*c+p+np.arange(p)]
    phases=[]
    for angle in (0.,math.pi):
        b=selected_map(raw,mu,angle,N)
        br=previous.convert(sector.Fin,b,sector.Fout)
        bl=previous.convert(sector.Lin,b.conj().T,sector.Lout)
        B=block_diag(br,ZR,bl,ZL);M=eye-S@B
        try:
            inverse=np.linalg.solve(M,eye)
        except np.linalg.LinAlgError as error:
            phases.append(dict(phase=angle,neumann_passed=False,neumann_epsilon=None,
                minimum_nominal_sewing_singular_value=float(np.linalg.svd(M,compute_uv=False)[-1]),
                failure='Nominal sewing inverse failed: '+str(error),B_nominal=_complex_json(B)))
            continue
        inverse_residual=previous.upper_norm(eye-inverse@M)
        rinv=previous.upper_norm(inverse);RS=inverse@S
        tail_error=gamma_outer*(enclosure['Z_error_upper']+KMAX*enclosure['energy_upper'])
        core_error=KMAX*gamma_core*core['frequency_lipschitz']
        Bsup=float(np.linalg.norm(B,2))+max(tail_error,core_error)
        bulk_term=rinv*KMAX*lopen*Bsup
        core_term=previous.upper_norm(RS[:,mouth_indices])*core_error
        tail_term=previous.upper_norm(RS[:,outer_indices])*tail_error
        neumann=inverse_residual+bulk_term+core_term+tail_term
        inverse_upper=rinv/(1-neumann) if neumann<1 else None
        Dfull=[];pressures=[];numerical_norms=[]
        for v in pressure:
            vv=selected_map(v,mu,angle,N)
            d=block_diag(previous.convert(sector.Fin,vv,sector.Fout),np.zeros((p,p)),
                previous.convert(sector.Lin,vv.conj().T,sector.Lout),np.zeros((p,p)))
            Dfull.append(d);pressures.append(float(np.trace(RS@d).real/math.pi))
            numerical_norms.append(float(np.linalg.norm(block_diag(vv,vv.conj().T),ord='nuc')))
        # The same conservative sqrt(2) reference gap is valid in every sector;
        # nonzero-m sectors have the larger exact massive gap checked above.
        ref_gap=math.sqrt(2)-KMAX*(core['frequency_lipschitz']+2*ELL*(math.pi/2-math.atan(RHO)))
        bounds=[]
        for j,probe in enumerate(sector.probes):
            zero_slope=4*probe['I_f'] if selected['has_zero_mode'] else 0.
            pnorm=numerical_norms[j]+KMAX*(selected['paired_pressure_trace_norm_slope'][j]+zero_slope)
            snorm=KMAX*probe['I_f']*2*c
            factor=((inverse_upper*(1+KMAX*lopen)*gamma_core+1/ref_gap)/math.pi
                    if inverse_upper is not None else None)
            bounds.append(dict(pressure_derivative_trace_norm_upper=pnorm,
                shift_derivative_trace_norm_upper=snorm,
                pressure_integral_absolute_upper=KMAX*factor*pnorm if factor is not None else None,
                shift_integral_absolute_upper=KMAX*factor*snorm if factor is not None else None,
                either_null_integral_absolute_upper=KMAX*factor*(pnorm+2*snorm) if factor is not None else None))
        phases.append(dict(phase=angle,minimum_nominal_sewing_singular_value=float(np.linalg.svd(M,compute_uv=False)[-1]),
            nominal_inverse_frobenius_upper=rinv,inverse_residual=inverse_residual,
            core_map_error_upper=core_error,tail_map_error_upper=tail_error,
            bulk_neumann_term=bulk_term,core_neumann_term=core_term,tail_neumann_term=tail_term,
            neumann_epsilon=neumann,neumann_passed=bool(neumann<1),inverse_upper=inverse_upper,
            reference_inverse_gap_lower=ref_gap,nominal_connected_pressure_at_zero=pressures,
            conditional_source_remainder=bounds,B_nominal=_complex_json(B),pressure_derivatives=_complex_json(Dfull)))
    tests=dict(sector_core=selected['passed'],far_volterra=tail['far']['passed'],
        inward_comparison=enclosure['passed'],open_dg_gram_identity=gram['relative_residual']<1e-12,
        positive_mass=gram['minimum_cell_mass_eigenvalue']>0,
        open_solves=max(right['primal_relative_residual'],right['dual_relative_residual'],left['primal_relative_residual'],left['dual_relative_residual'])<1e-10,
        opposite_open_adjoint=float(np.max(abs(SL-SR.conj().T)))<1e-10,
        nominal_tail_unitarity=enclosure['nominal_unitarity_error']<1e-9,
        actual_outer_labels=sector.config['outer_channels']==list(range(int(abs(mu)+.5),J+1)),
        actual_mouth_dimension=c==2*len(selected['mouth_channels'])+int(m==0)==sector.Gm.shape[1]==sector.Cm.shape[0],
        actual_tail_dimension=Z.shape==(p,p))
    return dict(m=m,mu=mu,config=sector.config,passed=all(tests.values()),tests=tests,
        certificate_gate_passed=all(v['neumann_passed'] for v in phases),
        origin='Fresh zero-frequency DG response in this sector; no frequency integration.',
        sector_core=selected,open_gram=gram,open_right=right,open_left=left,tail=tail,
        gram_conversion=dict(core_norm_factor=gamma_core,outer_norm_factor=gamma_outer),phases=phases,
        matrices=dict(S_open_right=_complex_json(SR),S_open_left=_complex_json(SL),
            Z_zero_nominal=_complex_json(Z),F_mouth_right_in=_complex_json(sector.Fin),
            F_mouth_right_out=_complex_json(sector.Fout),F_mouth_left_in=_complex_json(sector.Lin),
            F_mouth_left_out=_complex_json(sector.Lout),F_outer_right=[_complex_json(v) for v in rf],
            F_outer_left=[_complex_json(v) for v in lf]),seconds=time.perf_counter()-started)


def aggregate(sectors):
    phases=[]
    complete=[v['m'] for v in sectors]==list(range(-N,N+1))
    for index,angle in enumerate((0.,math.pi)):
        valid=[v for v in sectors if len(v.get('phases',[]))==2 and v['phases'][index]['neumann_passed']]
        failed=[v['m'] for v in sectors if v not in valid]
        fields=('pressure_integral_absolute_upper','shift_integral_absolute_upper','either_null_integral_absolute_upper')
        sums={key:[sum(v['phases'][index]['conditional_source_remainder'][j][key] for v in valid)
                    for j in range(2)] for key in fields}
        phases.append(dict(phase=angle,complete=complete,valid_sectors=[v['m'] for v in valid],failed_sectors=failed,
            all_sector_gate_passed=bool(complete and not failed),
            conditional_aggregate_bounds=sums if complete and not failed else None,
            passing_sector_partial_bounds=sums,
            maximum_neumann_epsilon=max((v['phases'][index]['neumann_epsilon'] for v in sectors
                if v.get('phases') and v['phases'][index]['neumann_epsilon'] is not None),default=None)))
    return phases


def save(data):
    DESTINATION.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(DESTINATION.read_text(encoding='utf-8'))==data


def run_checks():
    started=time.perf_counter();inputs=hashes()
    baseline=json.loads(project_path('threshold_return_bound_results.json').read_text(encoding='utf-8'))
    if not source_matches('threshold_return_bound.py', baseline['source_sha256'], inputs['threshold_return_bound.py']):
        raise ValueError('Frozen104 execution source does not match its recorded hash')
    if not baseline['passed'] or not baseline['certificate_gate_passed']:
        raise ValueError('Frozen104 m=0 evidence does not pass its own gates')
    raw,pressure,core=previous.zero_core(N,RHO)
    zero=copy.deepcopy({name:baseline[name] for name in PROOF_FIELDS})
    zero.update(m=0,mu=.5,origin='Exact saved104 numerical/proof payload; no new m=0 exterior solve.',
        sector_core=exact_sector_core(.5,raw,pressure,core),
        reused_result_sha256=inputs['threshold_return_bound_results.json'],
        reused_execution_source_sha256=baseline['source_sha256'])
    data=dict(stage='running',passed=False,certificate_gate_passed=False,
        parameters=dict(q=1,b=.2,d=1,rho_c=RHO,ell=ELL,h=.05,R=2.,L=40.,nmax=N,outer_nmax=J,angular_order=64,
            m_values=list(range(-N,N+1)),kappa_interval=[0.,KMAX],scaled_x_interval=[0.,1e-6],
            fresh_sector_count=8,reused_sector_count=1),
        core=core,probes=baseline['probes'],sectors=[zero],aggregate=[],
        input_sha256=inputs,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope=[
            'Complete m=-N..N only at the fixed mouth cutoff N4; not the infinite angular source.',
            'Both actual Weyls and the declared reflected angular/spin connection are used in every sector. No degeneracy or silent-chirality multiplier.',
            'Actual pair tail to infinity in the fixed J9, GL64 implemented radial equation, with the far/inward enclosures from104.',
            'Conditional analytic finite-DG bounds evaluated in floating arithmetic; LU residuals and analytic allowances are not verified interval roundoff.',
            'No bound on the remaining continuum, finite-J/GL64 errors, smooth junction or self-consistent Einstein solution is asserted here.',
            'This infrared bound does not replace or repair the finite-L40 truncation of the saved nonzero-frequency integration.',
            'Every failed Neumann gate is retained; an aggregate bound is asserted only if all nine sectors pass.'])
    save(data)
    for m in range(-N,N+1):
        if m==0:continue
        print(f'105: begin m={m:+d}, mu={m+.5:+g}',flush=True)
        row=sector_case(m,raw,pressure,core)
        data['sectors'].append(row);data['sectors'].sort(key=lambda v:v['m'])
        data['aggregate']=aggregate(data['sectors']);data['seconds']=time.perf_counter()-started
        save(data)
        print(f"105: m={m:+d}, controls={row['passed']}, Neumann={row['certificate_gate_passed']}, eps={[v['neumann_epsilon'] for v in row['phases']]}, seconds={row['seconds']:.3f}",flush=True)
    counts={str(n):sum(n in v['sector_core']['mouth_channels'] for v in data['sectors']) for n in range(1,N+1)}
    paired=sum(len(v['sector_core']['mouth_channels']) for v in data['sectors'])
    mouth=sum(v['sector_core']['mouth_dimension_per_weyl'] for v in data['sectors'])
    original=next(v for v in data['sectors'] if v['m']==0)
    unchanged=all(original[name]==baseline[name] for name in PROOF_FIELDS)
    tests=dict(all_nine_sectors=[v['m'] for v in data['sectors']]==list(range(-N,N+1)),
        sector_controls=all(v['passed'] for v in data['sectors']),
        selected_core_controls=all(v['sector_core']['passed'] for v in data['sectors']),
        full_multiplet_counts=all(counts[str(n)]==2*n+1 for n in range(1,N+1)),
        paired_count=paired==N*(N+2),mouth_dimension=mouth==2*N*(N+2)+1,
        unique_zero_sector=[v['m'] for v in data['sectors'] if v['sector_core']['has_zero_mode']]==[0],
        m0_saved_evidence_exact=unchanged,inputs_preserved=hashes()==inputs,
        all_summands_nonnegative=all(x>=0 for a in data['aggregate'] for values in a['passing_sector_partial_bounds'].values() for x in values))
    data.update(stage='completed',tests=tests,passed=all(tests.values()),
        certificate_gate_passed=all(a['all_sector_gate_passed'] for a in data['aggregate']),
        completeness=dict(per_n_sector_counts=counts,paired_modes=paired,mouth_dimension_per_weyl_sum=mouth,
            both_weyl_mouth_dimension_sum=2*mouth),seconds=time.perf_counter()-started)
    save(data)
    return data


if __name__=='__main__':
    result=run_checks()
    print(json.dumps(dict(passed=result['passed'],certificate_gate_passed=result['certificate_gate_passed'],
        seconds=result['seconds'],aggregate=result['aggregate']),indent=2))
    if not result['passed']:raise AssertionError('Complete threshold control failed; diagnostics saved')
