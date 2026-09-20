"""Finite-window complete N2/N4 source integrals for the rho_c=32 throat.

This reuses the frozen exterior, paired Weyl, radial variation and AdS-cap
implementations. Integration is in x=kappa*ell: saved integrands are the
physical kappa kernels divided by ell, and saved weights are x weights.
The rho_c=32 geometry/probes and three-node witness come from increment102.
No old source/data is edited, and no angular degeneracy multiplier is used.

Observed endpoint, panel, angular-cut and mesh changes are separate controls;
none is a certified omitted-tail or continuum bound. The sharp junction and
the self-consistent Einstein source remain unsolved.

Run: .venv\\Scripts\\python.exe -B -W error long_throat_integral.py
"""
from project_paths import project_path, source_matches
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

from ads_reference_source import reference_relative
from angular_reference_sum import comparisons, sum_records
from complete_angular_modes import mouth_channels, reflection_parities
from complete_angular_source import CompleteSector, FIELDS, PHASES, _decode


PILOT = 'long_throat_source_results.json'
OUTPUT = project_path('long_throat_integral_results.json')
X_EDGES = (1e-6,.02,.2,1.,4.,12.,24.,48.)
ENDPOINT_X = (1e-7,1e-6,48.,96.)
RECIPES = ('coarse','fine','finer_mesh')
ELL = .2*math.hypot(1.,32.)
SECTOR_TARGET = 2e-6


def hashes(names):
    return {name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names}


def write_result(result):
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if json.loads(OUTPUT.read_text(encoding='utf-8'))!=result:
        raise ValueError('Saved result did not round-trip exactly')


def reduced_config(config, nmax):
    mu=config['mu']
    return dict(config,nmax=nmax,mouth_channels=mouth_channels(mu,nmax).tolist(),
        reflection_cross_signs=reflection_parities(mu,nmax).tolist())


def integrate_cuts(case, order):
    """Share each physical exterior response; independently sew lower N."""
    started=time.perf_counter()
    nodes,weights=leggauss(order)
    cuts=(2,4) if abs(case.config['m'])<=2 else (4,)
    records={str(n):dict(config=reduced_config(case.config,n),order=order,
        edges=list(X_EDGES),frequency_coordinate='x=kappa*ell',ell=ELL,panels=[]) for n in cuts}
    minima={str(n):np.full((2,2),np.inf) for n in cuts}
    conditions={str(n):np.zeros((2,2)) for n in cuts}
    errors={str(n):dict(equation=0.,identity=0.,reference_shift=0.,phase=0.,jacobian=0.) for n in cuts}
    for lo,hi in zip(X_EDGES[:-1],X_EDGES[1:]):
        xs=(lo+hi)/2+(hi-lo)*nodes/2
        wx=(hi-lo)*weights/2
        raw={str(n):{p:{f:[] for f in FIELDS} for p in PHASES} for n in cuts}
        singular={str(n):[] for n in cuts}
        condition_samples={str(n):[] for n in cuts}
        for x in xs:
            high=case.evaluate(float(x/ELL),lower_nmax=2)
            for n in cuts:
                key=str(n)
                row=high if n==4 else high['lower_cut_response']
                errors[key]['equation']=max(errors[key]['equation'],*row['exterior_equation_relative_residuals'])
                sv=[p['minimum_singular_values'] for p in row['phases']]
                cc=[p['condition_numbers'] for p in row['phases']]
                minima[key]=np.minimum(minima[key],sv)
                conditions[key]=np.maximum(conditions[key],cc)
                singular[key].append(sv)
                condition_samples[key].append(cc)
                for ip,p in enumerate(PHASES):
                    sample=row['phases'][ip]
                    for f in FIELDS:
                        raw[key][p][f].append(sample[f+'_integrands'])
                    errors[key]['identity']=max(errors[key]['identity'],sample['resolvent_difference_identity_error'])
                    errors[key]['reference_shift']=max(errors[key]['reference_shift'],*abs(np.asarray(sample['reference_shift_trace'])))
                errors[key]['phase']=max(errors[key]['phase'],row['reference_phase_pressure_error'],row['reference_phase_shift_error'])
        for n in cuts:
            key=str(n)
            scaled={p:{f:(np.asarray(raw[key][p][f])/ELL).tolist() for f in FIELDS} for p in PHASES}
            responses={p:{f:np.einsum('i,ij->j',wx,scaled[p][f]).tolist() for f in FIELDS} for p in PHASES}
            for p in PHASES:
                for f in FIELDS:
                    alternate=np.einsum('i,ij->j',wx/ELL,raw[key][p][f])
                    errors[key]['jacobian']=max(errors[key]['jacobian'],float(np.max(abs(alternate-responses[p][f]))))
            records[key]['panels'].append(dict(bounds=[lo,hi],x=xs.tolist(),
                kappa=(xs/ELL).tolist(),weights=wx.tolist(),kappa_weights=(wx/ELL).tolist(),
                kappa_bounds=[lo/ELL,hi/ELL],integrands=scaled,raw_kappa_integrands=raw[key],
                responses=responses,minimum_singular_values=singular[key],condition_numbers=condition_samples[key]))
        if hi in (1.,12.,48.):
            print(json.dumps(dict(stage='integration_panel',m=case.config['m'],h=case.config['h'],
                order=order,through_x=hi,seconds=time.perf_counter()-started)),flush=True)
    for n in cuts:
        key=str(n)
        record=records[key]
        record['totals']={p:{f:np.sum([v['responses'][p][f] for v in record['panels']],axis=0).tolist()
            for f in FIELDS} for p in PHASES}
        record['diagnostics']=dict(minimum_singular_values=minima[key].tolist(),
            maximum_condition_numbers=conditions[key].tolist(),
            maximum_equation_relative_residual=errors[key]['equation'],
            resolvent_identity_error=errors[key]['identity'],
            reference_shift_trace_error=errors[key]['reference_shift'],
            reference_phase_error=errors[key]['phase'],jacobian_reintegration_error=errors[key]['jacobian'])
        record['seconds_shared_exterior']=time.perf_counter()-started
        record['exterior_reuse']=('Fresh N4 mouth-response columns' if n==4 else
            'Exact output-Gram restriction of the same N4 environment, followed by independent N2 core/reference sewing')
    return records


def record_controls(record):
    null_error=aggregate_error=jacobian_error=0.
    for recipe in RECIPES:
        row=record[recipe]
        if row['edges']!=list(X_EDGES) or row['ell']!=ELL:
            raise ValueError('Mismatched scaled-frequency integration window')
        for p in PHASES:
            pressure=np.asarray(row['totals'][p]['pressure'])
            shift=np.asarray(row['totals'][p]['shift'])
            null_error=max(null_error,float(np.max(abs(np.asarray(row['totals'][p]['null_plus'])-pressure+2*shift))),
                float(np.max(abs(np.asarray(row['totals'][p]['null_minus'])-pressure-2*shift))))
            for f in FIELDS:
                value=sum((np.einsum('i,ij->j',v['weights'],v['integrands'][p][f]) for v in row['panels']),np.zeros(2))
                alternate=sum((np.einsum('i,ij->j',v['kappa_weights'],v['raw_kappa_integrands'][p][f])
                    for v in row['panels']),np.zeros(2))
                aggregate_error=max(aggregate_error,float(np.max(abs(value-row['totals'][p][f]))))
                jacobian_error=max(jacobian_error,float(np.max(abs(value-alternate))))
    records=[record[r] for r in RECIPES]
    envelope=max(v for phase in record['comparisons']['quadrature_panel_envelope'].values()
        for values in phase.values() for v in values)
    tests=dict(exterior_equations=max(r['diagnostics']['maximum_equation_relative_residual'] for r in records)<1e-8,
        sewing_nonsingular=min(v for r in records for pair in r['diagnostics']['minimum_singular_values'] for v in pair)>1e-8,
        resolvent_identity=max(r['diagnostics']['resolvent_identity_error'] for r in records)<2e-9,
        reference_shift_zero=max(r['diagnostics']['reference_shift_trace_error'] for r in records)<2e-9,
        observed_quadrature_panel_envelope=envelope<SECTOR_TARGET,
        null_identities=null_error<1e-12,independent_reaggregation=aggregate_error<1e-12,
        frequency_jacobian=jacobian_error<1e-12)
    return dict(passed=all(tests.values()),tests={k:bool(v) for k,v in tests.items()},
        quadrature_target=SECTOR_TARGET,maximum_panel_envelope=envelope,
        null_identity_error=null_error,reaggregation_error=aggregate_error,jacobian_error=jacobian_error)


def witness_checks(pilot):
    """Reconstruct the frozen102 middle matrices, without an exterior solve."""
    group=next(g for g in pilot['geometries'] if g['geometry']['rho_c']==32.)
    rows=[]
    for sector in group['sectors']:
        old=sector['samples'][1]
        mats=old['matrices']
        connected=_decode(mats['S_connected'])
        reference=_decode(mats['S_reference'])
        error=0.
        for expected,saved in zip(old['phases'],mats['phases']):
            calculated=reference_relative(connected,reference,_decode(saved['B']),
                _decode(saved['pressure']),_decode(saved['shift']))
            error=max(error,*[float(np.max(abs(np.asarray(calculated[f+'_integrands'])-expected[f+'_integrands']))) for f in FIELDS])
        rows.append(dict(m=sector['m'],x=old['x'],kappa=old['kappa'],kernel_reconstruction_error=error))
    return dict(rows=rows,passed=max(r['kernel_reconstruction_error'] for r in rows)<1e-12,
        frozen_source_sha256=pilot['source_sha256'],frozen_tests=pilot['tests'],
        geometry=group['geometry'],samples=group['sums'],
        scope='All three frozen source samples are retained as a witness; only their already saved middle-node matrices are algebraically reconstructed. No pilot solve is repeated.')


def endpoint_samples(case):
    samples=[]
    for x in ENDPOINT_X:
        row=case.evaluate(x/ELL,lower_nmax=2)
        row['x']=x
        samples.append(row)
    return dict(config=case.config,samples=samples)


def sum_endpoints(sectors):
    result={}
    for nmax in (2,4):
        selected=[s for s in sectors if abs(s['config']['m'])<=nmax]
        rows=[]
        for index,x in enumerate(ENDPOINT_X):
            data=[s['samples'][index] if nmax==4 else s['samples'][index]['lower_cut_response'] for s in selected]
            raw={p:{f:np.sum([r['phases'][ip][f+'_integrands'] for r in data],axis=0).tolist()
                for f in FIELDS} for ip,p in enumerate(PHASES)}
            rows.append(dict(x=x,kappa=x/ELL,raw_kappa_kernels=raw,
                scaled_kernels={p:{f:(np.asarray(raw[p][f])/ELL).tolist() for f in FIELDS} for p in PHASES},
                minimum_sewing_singular_value=min(v for r in data for ph in r['phases'] for v in ph['minimum_singular_values']),
                maximum_equation_residual=max(v for r in data for v in r['exterior_equation_relative_residuals'])))
        result[str(nmax)]=rows
    return result


def outer_samples(base):
    """Same actual x and h, separate J; no unmatched old-node substitution."""
    other=CompleteSector(base.config['mu'],nmax=4,outer_nmax=11,rho_c=32.)
    samples=[]
    for x in (.4,4.):
        a=base.evaluate(x/ELL,lower_nmax=2)
        b=other.evaluate(x/ELL,lower_nmax=2)
        changes={p:{f:((np.asarray(b['phases'][ip][f+'_integrands'])-
            a['phases'][ip][f+'_integrands'])/ELL).tolist() for f in FIELDS} for ip,p in enumerate(PHASES)}
        samples.append(dict(x=x,kappa=x/ELL,J9=a,J11=b,scaled_J11_minus_J9=changes))
    return dict(m=base.config['m'],h=.05,nmax=4,samples=samples)


def aggregate_records(sectors,recipe):
    if any(s[recipe]['config']['outer_nmax']!=9 for s in sectors):
        raise ValueError('Cannot aggregate unequal outer angular cutoffs')
    summed=sum_records(sectors,recipe)
    for key in ('m','mouth_channels','outer_channels','has_zero_mode'):
        summed['config'].pop(key,None)
    summed.update(frequency_coordinate='x=kappa*ell',ell=ELL,m_values=[s['m'] for s in sectors])
    for ip,panel in enumerate(summed['panels']):
        first=sectors[0][recipe]['panels'][ip]
        for key in ('x','kappa_weights','kappa_bounds'):
            panel[key]=first[key]
        panel['raw_kappa_integrands']={p:{f:np.sum([s[recipe]['panels'][ip]['raw_kappa_integrands'][p][f]
            for s in sectors],axis=0).tolist() for f in FIELDS} for p in PHASES}
    return summed


def finalize(result):
    integration=result['integration']
    for key,cut in integration['cutoffs'].items():
        nmax=int(key)
        if [s['m'] for s in cut['sectors']]!=list(range(-nmax,nmax+1)):
            raise ValueError('Incomplete finite angular multiplet')
        sums={recipe:aggregate_records(cut['sectors'],recipe) for recipe in RECIPES}
        cut['partial_sums']=sums
        cut['comparisons']=comparisons(sums)
        cut['sum_of_sector_quadrature_envelopes']={p:{f:np.sum([
            s['comparisons']['quadrature_panel_envelope'][p][f] for s in cut['sectors']],axis=0).tolist()
            for f in FIELDS} for p in PHASES}
    cuts=integration['cutoffs']
    integration['N4_minus_N2']={r:{p:{f:(np.asarray(cuts['4']['partial_sums'][r]['totals'][p][f])-
        cuts['2']['partial_sums'][r]['totals'][p][f]).tolist() for f in FIELDS} for p in PHASES} for r in RECIPES}
    endpoints=result['frequency_controls']
    endpoints['sums']=sum_endpoints(endpoints['sectors'])
    endpoints['sampled_low_rectangle_scale']={key:{p:{f:(X_EDGES[0]*np.maximum(
        abs(np.asarray(rows[0]['scaled_kernels'][p][f])),abs(np.asarray(rows[1]['scaled_kernels'][p][f])))).tolist()
        for f in FIELDS} for p in PHASES} for key,rows in endpoints['sums'].items()}
    endpoints['low_scale_scope']='Window width times the larger of two sampled low-frequency magnitudes; NOT a supremum or rigorous omitted-IR bound.'
    endpoints['high_scope']='x96 lies beyond the x48 integration cutoff. Its measured decay or sign and the final panel do NOT bound the omitted UV integral.'
    aggregate_error=jacobian_error=null_error=0.
    for cut in cuts.values():
        for record in cut['partial_sums'].values():
            for p in PHASES:
                pressure=np.asarray(record['totals'][p]['pressure']);shift=np.asarray(record['totals'][p]['shift'])
                null_error=max(null_error,float(np.max(abs(np.asarray(record['totals'][p]['null_plus'])-pressure+2*shift))),
                    float(np.max(abs(np.asarray(record['totals'][p]['null_minus'])-pressure-2*shift))))
                for f in FIELDS:
                    a=sum((np.einsum('i,ij->j',v['weights'],v['integrands'][p][f]) for v in record['panels']),np.zeros(2))
                    b=sum((np.einsum('i,ij->j',v['kappa_weights'],v['raw_kappa_integrands'][p][f]) for v in record['panels']),np.zeros(2))
                    aggregate_error=max(aggregate_error,float(np.max(abs(a-record['totals'][p][f]))))
                    jacobian_error=max(jacobian_error,float(np.max(abs(a-b))))
    extra=[r for s in endpoints['sectors'] for r in s['samples']]
    extra += [r[name] for s in result['outer_cutoff_controls'] for r in s['samples'] for name in ('J9','J11')]
    extra_equation=max(v for r in extra for v in r['exterior_equation_relative_residuals'])
    extra_singular=min(v for r in extra for p in r['phases'] for v in p['minimum_singular_values'])
    envelope={key:max(v for ph in cut['sum_of_sector_quadrature_envelopes'].values() for values in ph.values() for v in values)
        for key,cut in cuts.items()}
    tests=dict(frozen102_passed=result['prior_controls']['passed'],
        saved_matrix_witness=result['witness']['passed'],
        geometry_and_probes_exact=result['geometry']==result['witness']['geometry'],
        all_complete_sector_controls=all(s['controls']['passed'] for c in cuts.values() for s in c['sectors']),
        N2_observed_quadrature=envelope['2']<5*SECTOR_TARGET,
        N4_observed_quadrature=envelope['4']<9*SECTOR_TARGET,
        extra_equation_controls=extra_equation<1e-8,extra_sewing_nonsingular=extra_singular>1e-8,
        summed_null_identities=null_error<1e-12,independent_reaggregation=aggregate_error<1e-12,
        correct_frequency_jacobian=jacobian_error<1e-12,
        inputs_preserved=result['input_sha256']==hashes(result['input_sha256']),
        execution_source_preserved=result['source_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result.update(passed=all(tests.values()),tests={k:bool(v) for k,v in tests.items()})
    result['diagnostics']=dict(sum_null_error=null_error,independent_reaggregation_error=aggregate_error,
        frequency_jacobian_error=jacobian_error,summed_quadrature_envelope=envelope,
        extra_maximum_equation_residual=extra_equation,extra_minimum_sewing_singular_value=extra_singular,
        main_exterior_frequency_responses=9*7*(8+12+8),extra_endpoint_responses=9*4,
        extra_outer_cutoff_responses=3*2*2,N2_extra_exterior_solves=0,
        omitted_low_frequency_bound=None,omitted_high_frequency_bound=None,
        omitted_angular_tail_bound=None,continuum_error_bound=None)


def run_checks():
    started=time.perf_counter()
    prior=json.loads(project_path(PILOT).read_text(encoding='utf-8'))
    names=sorted(set(prior['input_sha256'])|{'long_throat_source.py',PILOT})
    inputs=hashes(names)
    result=dict(passed=False,tests={},stage='Preparing finite-window long-throat integral',seconds=0.,
        input_sha256=inputs,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parameters=dict(q=1,b=.2,d=1.,R=2.,L=40.,rho_c=32.,ell=ELL,nmax=4,lower_nmax=2,outer_nmax=9,
            x_edges=list(X_EDGES),kappa_edges=[x/ELL for x in X_EDGES],
            recipes=dict(coarse=dict(h=.05,order=8),fine=dict(h=.05,order=12),finer_mesh=dict(h=.035,order=8)),
            primary_recipe='fine',primary_phase=math.pi,phase_settings=[0.,math.pi],
            individual_quadrature_target=SECTOR_TARGET,summed_N4_target=9*SECTOR_TARGET,
            sign_is_not_an_acceptance_gate=True,units='hbar=c=d=1'),
        prior_controls=dict(passed=prior['passed'],tests=prior['tests'],source_sha256=prior['source_sha256']),
        integration=dict(cutoffs={'2':dict(sectors=[]),'4':dict(sectors=[])}),
        frequency_controls=dict(x=list(ENDPOINT_X),sectors=[]),outer_cutoff_controls=[],
        observable_definitions=dict(integrands='Q(kappa)/ell, with x=kappa*ell and dx weights',
            raw_kappa_integrands='Q(kappa), recorded independently with weights dx/ell',
            pressure='Smeared connected-minus-invariant-AdS-reference rho+p_r',
            shift='Negative smeared lower orthonormal radial flux',null_plus='pressure-2*shift',null_minus='pressure+2*shift'),
        sources=prior.get('sources',[]),
        scope=prior['scope']+[
            'This increment integrates only the finite x window [1e-6,48], equivalently the recorded kappa window. It does not integrate the saved three-node pilot.',
            'N2/N4 sums each complete their own finite angular multiplets, with both actual Weyl responses; no full angular completion or arbitrary degeneracy scaling.',
            'The observed quadrature envelope, matched-order mesh difference, N4-minus-N2 change and selected outer J11-minus-J9 change are distinct diagnostics.',
            'The two low-frequency samples and high-frequency x96 sample are actual source evaluations, not assumed extrapolations or certified omitted-integral bounds.',
            'No tail-cache hook is added to the frozen solver. Existing radial caches are reused; no global monkeypatching or historical source changes occur.'])
    # Replace the inherited pilot-only statement while preserving its scope in the witness.
    result['scope']=[s for s in result['scope'] if not s.startswith('Three matched x=')]
    try:
        if not prior['passed'] or not all(prior['tests'].values()):
            raise ValueError('Frozen102 must pass before this integration')
        if not source_matches('long_throat_source.py', prior['source_sha256'], inputs['long_throat_source.py']):
            raise ValueError('Frozen102 execution source changed')
        result['witness']=witness_checks(prior)
        result['geometry']=result['witness']['geometry']
        if result['geometry']['ell']!=ELL or not result['witness']['passed']:
            raise ValueError('Frozen geometry or matrix witness failed')
        write_result(result)
        for m in range(-4,5):
            tick=time.perf_counter()
            case=CompleteSector(m+.5,nmax=4,outer_nmax=9,rho_c=32.)
            coarse=integrate_cuts(case,8)
            fine=integrate_cuts(case,12)
            result['frequency_controls']['sectors'].append(endpoint_samples(case))
            if m in (-4,0,4):
                result['outer_cutoff_controls'].append(outer_samples(case))
            small=CompleteSector(m+.5,h=.035,nmax=4,outer_nmax=9,rho_c=32.)
            finer=integrate_cuts(small,8)
            for key in coarse:
                record=dict(mu=m+.5,m=m,coarse=coarse[key],fine=fine[key],finer_mesh=finer[key])
                record['comparisons']=comparisons(record)
                record['controls']=record_controls(record)
                result['integration']['cutoffs'][key]['sectors'].append(record)
            result.update(stage=f'Completed finite-window sector m={m}',seconds=time.perf_counter()-started)
            write_result(result)
            print(json.dumps(dict(stage='sector_integrated',m=m,seconds=time.perf_counter()-tick,
                primary_pi=fine['4']['totals']['pi'],finer_mesh_pi=finer['4']['totals']['pi'],
                controls=result['integration']['cutoffs']['4']['sectors'][-1]['controls'])),flush=True)
            del case,small
        finalize(result)
        result.update(stage='Completed rho32 finite-window complete N2/N4 source integrals',seconds=time.perf_counter()-started)
        write_result(result)
    except Exception as exc:
        result.update(passed=False,stage='Stopped with retained diagnostics',failure=dict(type=type(exc).__name__,message=str(exc)),
            seconds=time.perf_counter()-started)
        write_result(result)
        raise
    return result


if __name__=='__main__':
    result=run_checks()
    print(json.dumps(dict(passed=result['passed'],tests=result['tests'],seconds=result['seconds'],
        primary=result['integration']['cutoffs'].get('4',{}).get('partial_sums',{}).get('fine',{}).get('totals'))),flush=True)
    if not result['passed']:
        raise SystemExit('A numerical gate failed; calculated diagnostics were saved.')
