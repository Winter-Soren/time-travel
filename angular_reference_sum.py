"""Explicit three-sector sum of the compact AdS-reference source diagnostic.

The mu=1/2 quadratures are reused byte-for-value from frozen increment99.
Only mu=-1/2 and3/2 are newly solved, each with both actual Weyl responses.
All use the same reflected endpoint connection, frequency panels and mesh
recipes.  This completes the n=1 multiplet but not higher angular multiplets.
It is not a full absolute stress, an NEC proof or a self-consistent medium.

Run: .venv\\Scripts\\python.exe -B -W error angular_reference_sum.py
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from ads_reference_source import integrate_reference,saved_case,FIELDS,PHASES
from azimuthal_source import AzimuthalSector


EDGES=(1e-6,.02,.2,1.,4.,12.,24.,48.,96.)
RECIPES=('coarse','fine','finer_mesh')
QUADRATURE_TARGET=2e-6


def difference(a,b):
    return {phase:{field:(np.asarray(b['totals'][phase][field])-a['totals'][phase][field]).tolist()
                  for field in FIELDS} for phase in PHASES}


def comparisons(record):
    coarse,fine,finer=(record[name] for name in RECIPES)
    envelope={phase:{field:np.sum([abs(np.asarray(b['responses'][phase][field])-a['responses'][phase][field])
                       for a,b in zip(coarse['panels'],fine['panels'])],axis=0).tolist()
                    for field in FIELDS} for phase in PHASES}
    return dict(quadrature_difference=difference(coarse,fine),quadrature_panel_envelope=envelope,
        matched_order_mesh_difference=difference(coarse,finer),
        final_panel_contribution=fine['panels'][-1]['responses'],
        final_two_panel_contribution={phase:{field:(np.asarray(fine['panels'][-2]['responses'][phase][field])+
            fine['panels'][-1]['responses'][phase][field]).tolist() for field in FIELDS} for phase in PHASES})


def record_controls(record):
    null_error=reaggregation_error=0.
    for name in RECIPES:
        r=record[name]
        if r['edges']!=list(EDGES):raise ValueError('Mismatched frequency window')
        for phase in PHASES:
            p=np.asarray(r['totals'][phase]['pressure']);s=np.asarray(r['totals'][phase]['shift'])
            null_error=max(null_error,float(np.max(abs(np.asarray(r['totals'][phase]['null_plus'])-p+2*s))),
                           float(np.max(abs(np.asarray(r['totals'][phase]['null_minus'])-p-2*s))))
            for field in FIELDS:
                weighted=sum((np.einsum('i,ij->j',panel['weights'],panel['integrands'][phase][field])
                              for panel in r['panels']),np.zeros(2))
                reaggregation_error=max(reaggregation_error,float(np.max(abs(weighted-r['totals'][phase][field]))))
    records=[record[name] for name in RECIPES]
    envelope=max(v for phase in record['comparisons']['quadrature_panel_envelope'].values() for values in phase.values() for v in values)
    checks=dict(exterior_equation=max(r['diagnostics']['maximum_equation_relative_residual'] for r in records)<1e-8,
        sewing_nonsingular=min(v for r in records for pair in r['diagnostics']['minimum_singular_values'] for v in pair)>1e-8,
        resolvent_identity=max(r['diagnostics']['resolvent_identity_error'] for r in records)<2e-9,
        reference_shift_zero=max(r['diagnostics']['reference_shift_trace_error'] for r in records)<2e-9,
        quadrature_panel_envelope=envelope<QUADRATURE_TARGET,
        null_identities=null_error<1e-12,independent_reaggregation=reaggregation_error<1e-12)
    return dict(passed=all(checks.values()),tests={k:bool(v) for k,v in checks.items()},
        quadrature_target=QUADRATURE_TARGET,maximum_panel_envelope=envelope,
        null_identity_error=null_error,reaggregation_error=reaggregation_error)


def sum_records(sectors,key):
    """Sum actual sectors only after matching geometry, quadrature and nodes."""
    records=[s[key] for s in sectors];base=records[0]
    for r in records[1:]:
        for name in ('h','R','nmax','L','rho_c'):
            if r['config'][name]!=base['config'][name]:raise ValueError('Unmatched source geometry: '+name)
        if r['order']!=base['order'] or r['edges']!=base['edges']:raise ValueError('Unmatched quadrature recipe')
    panels=[]
    for ip,p in enumerate(base['panels']):
        group=[r['panels'][ip] for r in records]
        if any(g['kappa']!=p['kappa'] or g['weights']!=p['weights'] or g['bounds']!=p['bounds'] for g in group):
            raise ValueError('Sector sums require exactly matching nodes')
        integrands={phase:{field:np.sum([g['integrands'][phase][field] for g in group],axis=0).tolist()
                           for field in FIELDS} for phase in PHASES}
        responses={phase:{field:np.sum([g['responses'][phase][field] for g in group],axis=0).tolist()
                          for field in FIELDS} for phase in PHASES}
        panels.append(dict(bounds=p['bounds'],kappa=p['kappa'],weights=p['weights'],integrands=integrands,responses=responses))
    totals={phase:{field:np.sum([r['totals'][phase][field] for r in records],axis=0).tolist()
                   for field in FIELDS} for phase in PHASES}
    return dict(config={k:v for k,v in base['config'].items() if k not in ('mu','reflection_cross_signs')},
        sectors=[s['mu'] for s in sectors],order=base['order'],edges=base['edges'],panels=panels,totals=totals)


def write_snapshot(result):
    path=Path(__file__).with_name('angular_reference_sum_results.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(path.read_text(encoding='utf-8'))==result


def checks():
    started=time.perf_counter();prior_path=Path('ads_reference_source_results.json')
    prior=json.loads(prior_path.read_text(encoding='utf-8'))
    if not prior['passed'] or not all(prior['tests'].values()):raise ValueError('Frozen99 acceptance record is not passing')
    input_names=set(prior['input_sha256'])|{'ads_reference_source.py','ads_reference_source_results.json'}
    hashes={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in sorted(input_names)}
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    reused=dict(mu=.5,reused_from='ads_reference_source_results.json',computed_here=False,
        **{key:prior['integration'][key] for key in RECIPES})
    reused['comparisons']=comparisons(reused);reused['controls']=record_controls(reused)
    result=dict(passed=False,tests={},stage='Running: frozen mu=1/2 records loaded, new sectors pending',
        parameters=dict(q=1,sectors=[-.5,.5,1.5],b=.2,d=1.,R=2.,nmax=8,L=40.,rho_c=2.,
            units='hbar=c=d=1',primary_phase=math.pi,primary_recipe='fine',connection='explicit reflected spin/gauge endpoint branch'),
        sectors=[reused],probes=prior['probes'],reused_reference_controls=prior['tests'],
        sources=prior['sources'],input_sha256=hashes,source_sha256=source_hash,seconds=0.,
        observable_definitions=prior['observable_definitions'],
        scope=['Every newly added sector is explicitly computed with both actual physical Weyl responses; no angular degeneracy or chirality multiplier is used.',
            'The mu=1/2 integration is reused exactly from frozen99 at identical nodes, mesh, angular cutoff, connection and phase; it is not recomputed.',
            'The three sectors complete only the n=1 multiplet. Each retains coupled n1..8 and exterior mixing; higher angular multiplets are incomplete.',
            'The source target is connected-minus-standard-invariant-AdS-reference at each separately declared phase, not phase pi-minus-zero.',
            'The invariant reference has zero radial null contraction by symmetry. A full absolute local source still requires the complete angular/frequency limits, correctly matched states and continuum sewing, and common covariant renormalization.',
            'The frequency window [1e-6,96] and observed mesh/quadrature changes are finite numerical diagnostics, not certified missing-tail or continuum bounds.',
            'The prescribed product throat has G_kk=0. A nonzero completed source would require geometry response; this calculation supplies no matched collar, formation mechanism, stable medium or time machine.'])
    write_snapshot(result)
    try:
        for mu in (-.5,1.5):
            tick=time.perf_counter();case=AzimuthalSector(mu,h=.05,R=2.,nmax=8,L=40.,rho_c=2.)
            witness=saved_case(case.config,case.evaluate(.4,save_matrices=True,adjoint_control=False))
            baseline=next(row for group in prior['mesh_refinement'] if group['config']['mu']==mu and group['config']['h']==.05
                          for row in group['samples'] if row['kappa']==.4)
            witness_error=max(abs(np.asarray(a[field+'_integrands'])-b[field+'_integrands']).max()
                for a,b in zip(witness['phases'],baseline['phases']) for field in FIELDS)
            coarse=integrate_reference(case,8,EDGES);fine=integrate_reference(case,12,EDGES)
            small=AzimuthalSector(mu,h=.035,R=2.,nmax=8,L=40.,rho_c=2.)
            finer=integrate_reference(small,8,EDGES)
            record=dict(mu=mu,computed_here=True,coarse=coarse,fine=fine,finer_mesh=finer,
                saved_pilot_reproduction_error=float(witness_error),seconds=time.perf_counter()-tick)
            record['comparisons']=comparisons(record);record['controls']=record_controls(record)
            result['sectors'].append(record);result['sectors'].sort(key=lambda r:r['mu'])
            result['stage']='Completed mu='+str(mu)+'; '+('one new sector remains' if mu==-.5 else 'final sum/checks pending')
            result['seconds']=time.perf_counter()-started;write_snapshot(result)
            print(json.dumps(dict(completed_mu=mu,primary=record['fine']['totals'],
                finer_mesh=record['finer_mesh']['totals'],controls=record['controls'],seconds=record['seconds'])),flush=True)
            del case,small
        sums={key:sum_records(result['sectors'],key) for key in RECIPES}
        total_comparisons=comparisons(sums)
        conservative={phase:{field:np.sum([r['comparisons']['quadrature_panel_envelope'][phase][field] for r in result['sectors']],axis=0).tolist()
                            for field in FIELDS} for phase in PHASES}
        null_error=aggregation_error=0.
        for record in sums.values():
            for phase in PHASES:
                p=np.asarray(record['totals'][phase]['pressure']);s=np.asarray(record['totals'][phase]['shift'])
                null_error=max(null_error,float(np.max(abs(np.asarray(record['totals'][phase]['null_plus'])-p+2*s))),
                               float(np.max(abs(np.asarray(record['totals'][phase]['null_minus'])-p-2*s))))
                for field in FIELDS:
                    computed=sum((np.einsum('i,ij->j',panel['weights'],panel['integrands'][phase][field]) for panel in record['panels']),np.zeros(2))
                    aggregation_error=max(aggregation_error,float(np.max(abs(computed-record['totals'][phase][field]))))
        tests=dict(frozen99_prior_passed=prior['passed'] and all(prior['tests'].values()),
            frozen_mu05_records_exact=all(reused[key]==prior['integration'][key] for key in RECIPES),
            both_added_sectors_completed=len(result['sectors'])==3,
            selected_saved_pilot_reproduction=max(r.get('saved_pilot_reproduction_error',0.) for r in result['sectors'])<2e-11,
            all_individual_source_controls=all(r['controls']['passed'] for r in result['sectors']),
            combined_observed_quadrature=max(v for phase in conservative.values() for values in phase.values() for v in values)<3*QUADRATURE_TARGET,
            summed_null_identities=null_error<1e-12,summed_nodes_and_integrals=aggregation_error<1e-12,
            all_inputs_preserved=hashes=={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in hashes},
            executed_source_preserved=source_hash==hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        result.update(passed=all(tests.values()),tests={k:bool(v) for k,v in tests.items()},
            stage='Completed explicit three-sector AdS-reference source sum with observed controls',
            partial_sums=sums,comparisons=total_comparisons,
            diagnostics=dict(sum_of_individual_quadrature_envelopes=conservative,
                individual_quadrature_target=QUADRATURE_TARGET,combined_quadrature_target=3*QUADRATURE_TARGET,
                summed_null_identity_error=null_error,independent_aggregate_error=aggregation_error,
                omitted_frequency_tail_bound=None,full_angular_remainder_bound=None,continuum_error_bound=None),
            seconds=time.perf_counter()-started)
        write_snapshot(result)
        return result
    except Exception as exc:
        result.update(passed=False,stage='Stopped with recorded failure',failure=dict(type=type(exc).__name__,message=str(exc)),seconds=time.perf_counter()-started)
        write_snapshot(result)
        raise


if __name__=='__main__':
    result=checks()
    print(json.dumps(dict(passed=result['passed'],tests=result['tests'],seconds=result['seconds'],
        primary=result['partial_sums']['fine']['totals'],source_sha256=result['source_sha256'])),flush=True)
    if not result['passed']:raise SystemExit('An angular source-sum gate failed; numerical evidence is retained.')
