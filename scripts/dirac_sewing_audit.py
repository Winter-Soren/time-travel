"""Opposite-Weyl assembly and physical versus adjoint-paired facet sewing.

The second Weyl transport is assembled independently on the SAME mesh and
regular-axis domain: reverse all spatial Weyl terms and all upwind senses,
but keep the positive imaginary-frequency mass. Spherical lead bases are
swapped and reprojected; their physical Gram-converted maps need not equal
the adjoints of the first finite-facet maps. This distinction is measured,
never removed by clipping or by silently doubling a one-chirality result.

Run: .venv\\Scripts\\python.exe -B -W error dirac_sewing_audit.py
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from exterior_boundary_map import ExteriorDG,ID,SX,SY,SZ,_complex_json,_quadrature,pair_potential
from exterior_tail_transport import facet_map
from imaginary_weyl_leads import tail_map
from phase_vacuum_contrast import PhaseContrast
from spherical_weyl_leads import mouth_characteristics,outer_characteristics,zero_mode
from throat_null_response import raw_variation,shift_derivatives,phase_map,probe_geometry


def opposite_assembly(model):
    """Independent E=0 left-Weyl volume and upwind faces; no matrix transpose."""
    rows=[];cols=[];values=[]
    bary,qw,_,_=_quadrature()
    def block(tl,tr,value):
        ii=6*tl+np.arange(6);jj=6*tr+np.arange(6)
        rows.extend(np.repeat(ii,6));cols.extend(np.tile(jj,6));values.extend(value.ravel())
    for it,tri in enumerate(model.triangles):
        xy=model.vertices[tri]
        jac=np.linalg.det(np.stack((xy[1]-xy[0],xy[2]-xy[0])))
        grad=np.linalg.inv(np.column_stack((np.ones(3),xy)))[1:].T
        points=bary@xy;rho,z=points.T;weights=qw*jac
        moments=np.einsum('p,pa->a',weights*rho,bary)
        plain=np.einsum('p,pa,pb->ab',weights,bary,bary)
        azimuth=model.mu-pair_potential(rho,z,model.q,model.d)
        angular=np.einsum('p,pa,pb->ab',weights*azimuth,bary,bary)
        volume=np.empty((3,2,3,2),complex)
        for a in range(3):
            for b in range(3):
                volume[a,:,b,:]=(-moments[a]*(grad[b,0]*SX+grad[b,1]*SZ)
                                -plain[a,b]*SX/2-1j*angular[a,b]*SY)
        block(it,it,volume.reshape(6,6))
    for edge in model.interior:
        it,jt=edge['triangle'],edge['other'];w=edge['weights']
        left,right=edge['basis'],edge['other_basis']
        an=edge['normal'][0]*SX+edge['normal'][1]*SZ
        incoming_left=(ID+an)/2;incoming_right=(ID-an)/2
        ll=np.einsum('p,pa,pb->ab',w,left,left)
        lr=np.einsum('p,pa,pb->ab',w,left,right)
        rr=np.einsum('p,pa,pb->ab',w,right,right)
        block(it,it,np.kron(ll,incoming_left));block(it,jt,-np.kron(lr,incoming_left))
        block(jt,jt,np.kron(rr,incoming_right));block(jt,it,-np.kron(lr.T,incoming_right))
    boundary=[]
    for edge in model.boundary:
        an=edge['normal'][0]*SX+edge['normal'][1]*SZ
        incoming=(ID+an)/2;outgoing=(ID-an)/2
        fm=np.einsum('p,pa,pb->ab',edge['weights'],edge['basis'],edge['basis'])
        block(edge['triangle'],edge['triangle'],np.kron(fm,incoming))
        boundary.append(dict(edge,pm=incoming,pp=outgoing))
    matrix=coo_matrix((values,(rows,cols)),shape=(model.ndof,model.ndof)).tocsc()
    proxy=SimpleNamespace(ndof=model.ndof,boundary=boundary,active=model.active.copy())
    return matrix,proxy


def raw_basis(kind,points,R,nmax):
    if kind=='outer':
        return outer_characteristics(.5,np.arctan2(points[:,0],points[:,1]),R,nmax)
    theta=np.arctan2(points[:,0],points[:,1]-(.5 if kind=='upper' else -.5))
    plus,minus=mouth_characteristics(theta,.2,kind,nmax)
    zero=zero_mode(theta,.2,kind)[:,:,None]
    if kind=='upper':minus=np.concatenate((zero,minus),axis=2)
    else:plus=np.concatenate((zero,plus),axis=2)
    return plus,minus


def port_factors(port,basis_function):
    """Rebuild physical spherical-to-facet factors for a possibly unpaired port."""
    factors=[];error=0.
    for index,sense,projector in ((0,'plus','pp'),(1,'minus','pm')):
        raw=np.concatenate([np.einsum('ab,pbc->pac',e[projector],
            np.moveaxis(basis_function(e['points'])[index],0,1)) for e in port['edges']])
        gram=np.einsum('p,pai,paj->ij',port['weights'],raw.conj(),raw)
        F=np.linalg.cholesky((gram+gram.conj().T)/2).conj().T
        normalized=np.linalg.solve(F.T,raw.reshape(-1,len(F)).T).T.reshape(raw.shape)
        error=max(error,float(np.max(abs(normalized-port[sense]))));factors.append(F)
    return factors,error


def convert_map(raw,Fin,Fout):
    return np.linalg.solve(Fout.T,(Fin@raw).T).T


def exterior_response(matrix,active,Gm,Cm,Go,Co,Z):
    """Eliminate only the outer channels, leaving every mouth input independent."""
    A=matrix[active][:,active].tocsc();lu=splu(A)
    gm,go=Gm[active],Go[active];cm,co=Cm[:,active],Co[:,active]
    count=gm.shape[1];lifted=lu.solve(np.column_stack((gm,go)))
    Xm,Xo=lifted[:,:count],lifted[:,count:]
    M=np.eye(len(Z))-co@Xo@Z
    X=Xm+Xo@Z@np.linalg.solve(M,co@Xm)
    residual=A@X-gm-go@Z@(co@X)
    return cm@X,float(np.linalg.norm(residual)/np.linalg.norm(gm))


def observables(S,Bs,pressures,shifts):
    logabs=[];phase=[];minimum=[];condition=[];p=[];s=[]
    for B,db,ds in zip(Bs,pressures,shifts):
        M=np.eye(len(B))-S@B;sign,logdet=np.linalg.slogdet(M)
        sv=np.linalg.svd(M,compute_uv=False)
        if sign==0 or not math.isfinite(logdet):raise ArithmeticError('Singular sampled sewing')
        logabs.append(float(logdet));phase.append(_complex_json(sign))
        minimum.append(float(sv[-1]));condition.append(float(sv[0]/sv[-1]))
        p.append(np.array([np.trace(np.linalg.solve(M,S@d)) for d in db]))
        s.append(np.array([np.trace(np.linalg.solve(M,S@d)) for d in ds]))
    pressure=(p[1]-p[0]).real/math.pi;shift=(s[1]-s[0]).real/math.pi
    contrast=logabs[1]-logabs[0]
    return dict(log_determinants=logabs,determinant_phases=phase,delta_logabs=contrast,
        energy_integrand=-contrast/math.pi,pressure_integrands=pressure.tolist(),
        shift_integrands=shift.tolist(),null_plus_integrands=(pressure-2*shift).tolist(),
        null_minus_integrands=(pressure+2*shift).tolist(),minimum_singular_values=minimum,
        condition_numbers=condition)


class DiracAudit:
    def __init__(self,h=.08,R=2.,nmax=8,L=40.,rho_c=2.):
        self.right=PhaseContrast(h,R,nmax,L,rho_c);r=self.right
        self.config=r.config.copy();self.active=r.model.active
        self.left_zero,self.proxy=opposite_assembly(r.model)
        self.left_ports=[];self.left_factors={};self.port_checks=[]
        for kind,right_port in zip(('upper','lower','outer'),(r.upper,r.lower,r.outer)):
            right_basis=lambda points,kind=kind:raw_basis(kind,points,R,nmax)
            left_basis=lambda points,base=right_basis:base(points)[::-1]
            port=ExteriorDG.modal_port(self.proxy,kind,left_basis)
            factors,error=port_factors(port,left_basis)
            right_factors,right_error=port_factors(right_port,right_basis)
            self.left_ports.append(port);self.left_factors[kind]=factors
            self.port_checks.append(dict(kind=kind,reproduction_error=max(error,right_error),
                incoming_vs_right_outgoing_adjoint=float(np.max(abs(port['G']-right_port['C'].conj().T))),
                outgoing_vs_right_incoming_adjoint=float(np.max(abs(port['C']-right_port['G'].conj().T))),
                swapped_factor_error=max(float(np.max(abs(factors[j]-right_factors[1-j]))) for j in (0,1))))
        upper,lower,outer=self.left_ports
        self.Gm=np.column_stack((upper['G'],lower['G']));self.Cm=np.vstack((upper['C'],lower['C']))
        self.Fin=block_diag(self.left_factors['upper'][1],self.left_factors['lower'][1])
        self.Fout=block_diag(self.left_factors['upper'][0],self.left_factors['lower'][0])
        difference=(self.left_zero-r.unshifted.conj().T).tocoo()
        scale=np.max(abs(r.unshifted.data))
        self.assembly=dict(triangles=len(r.model.triangles),active_dof=len(self.active),
            opposite_adjoint_maximum_error=float(np.max(abs(difference.data))) if difference.nnz else 0.,
            opposite_adjoint_relative_error=(float(np.max(abs(difference.data)))/scale if difference.nnz else 0.),
            axis_domain_equal=bool(np.array_equal(self.proxy.active,r.model.active)),
            description='Opposite volume and every upwind face rebuilt independently; matrix conjugate transpose used only afterwards as a test.')
        self.probes=probe_geometry(.2*math.hypot(1.,rho_c),rho_c);self.cache={}

    def evaluate(self,kappa,save_matrices=False,adjoint_control=True):
        kappa=float(kappa)
        if not math.isfinite(kappa) or kappa<=0:raise ValueError('Positive finite kappa required')
        if save_matrices and not adjoint_control:raise ValueError('Saved matrix audit includes the adjoint control')
        key=(kappa,bool(save_matrices),bool(adjoint_control))
        if key in self.cache:return self.cache[key]
        started=time.perf_counter();r=self.right;N=self.config['nmax']
        rawZ,tail=tail_map(kappa,self.config['R'],self.config['L'],nmax=N)
        ZR,FR,_=facet_map(r.outer,r.outer_basis,rawZ)
        left_outer_basis=lambda points:r.outer_basis(points)[::-1]
        ZL,FL,_=facet_map(self.left_ports[2],left_outer_basis,rawZ.conj().T)
        AR=r.unshifted+kappa*r.mass;AL=self.left_zero+kappa*r.mass
        SR,resR=exterior_response(AR,self.active,r.Gm,r.Cm,r.outer['G'],r.outer['C'],ZR)
        SL,resL=exterior_response(AL,self.active,self.Gm,self.Cm,self.left_ports[2]['G'],self.left_ports[2]['C'],ZL)
        if adjoint_control:
            SA,resA=exterior_response(AL,self.active,self.Gm,self.Cm,self.left_ports[2]['G'],self.left_ports[2]['C'],ZR.conj().T)
        else:SA,resA=None,None
        raw,dp,modes=raw_variation(kappa,N,self.config['rho_c'])
        ds=shift_derivatives(raw,kappa,N,self.probes)
        BR=[];PR=[];QR=[];BL=[];PL=[];QL=[]
        for angle in (0.,math.pi):
            b=phase_map(raw,angle,N)
            pb=[phase_map(v,angle,N) for v in dp];qb=[phase_map(v,angle,N) for v in ds]
            BR.append(convert_map(b,r.Fin,r.Fout));PR.append([convert_map(v,r.Fin,r.Fout) for v in pb])
            QR.append([convert_map(v,r.Fin,r.Fout) for v in qb])
            BL.append(convert_map(b.conj().T,self.Fin,self.Fout))
            PL.append([convert_map(v.conj().T,self.Fin,self.Fout) for v in pb])
            QL.append([-convert_map(v.conj().T,self.Fin,self.Fout) for v in qb])
        BA=[b.conj().T for b in BR];PA=[[v.conj().T for v in p] for p in PR]
        QA=[[-v.conj().T for v in q] for q in QR]
        right=observables(SR,BR,PR,QR);physical=observables(SL,BL,PL,QL)
        adjoint=observables(SA,BA,PA,QA) if adjoint_control else None
        compare=lambda a,b,sign=1:[float(x-sign*y) for x,y in zip(a,b)]
        paired=dict(energy_integrand=right['energy_integrand']+physical['energy_integrand'],
            pressure_integrands=compare(right['pressure_integrands'],physical['pressure_integrands'],-1),
            shift_integrands=compare(right['shift_integrands'],physical['shift_integrands'],-1),
            null_plus_integrands=compare(right['null_plus_integrands'],physical['null_plus_integrands'],-1),
            null_minus_integrands=compare(right['null_minus_integrands'],physical['null_minus_integrands'],-1))
        discrepancies=dict(energy_physical_minus_right=physical['energy_integrand']-right['energy_integrand'],
            pressure_physical_minus_right=compare(physical['pressure_integrands'],right['pressure_integrands']),
            shift_physical_plus_right=compare(physical['shift_integrands'],right['shift_integrands'],-1),
            exterior_physical_minus_right_adjoint=float(np.linalg.norm(SL-SR.conj().T)),
            exterior_adjoint_prescription_error=float(np.linalg.norm(SA-SR.conj().T)) if adjoint_control else None,
            outer_physical_map_adjoint_defect=float(np.linalg.norm(ZL-ZR.conj().T)),
            mouth_physical_map_adjoint_defect=max(float(np.linalg.norm(a-b)) for a,b in zip(BL,BA)),
            adjoint_energy_error=abs(adjoint['energy_integrand']-right['energy_integrand']) if adjoint_control else None,
            adjoint_pressure_error=max(abs(np.array(adjoint['pressure_integrands'])-right['pressure_integrands'])) if adjoint_control else None,
            adjoint_shift_error=max(abs(np.array(adjoint['shift_integrands'])+right['shift_integrands'])) if adjoint_control else None)
        discrepancies={k:(float(v) if np.isscalar(v) else v) for k,v in discrepancies.items()}
        row=dict(kappa=kappa,right=right,left_physical=physical,left_adjoint_prescription=adjoint,
            actual_two_chirality_sum=paired,discrepancies=discrepancies,
            exterior_equation_relative_residuals=([resR,resL,resA] if adjoint_control else [resR,resL]),
            tail_accepted_steps=tail['accepted_steps'],seconds=time.perf_counter()-started)
        if save_matrices:
            row['matrices']={key:_complex_json(value) for key,value in dict(
                S_right=SR,S_left_physical=SL,S_left_adjoint_prescription=SA,
                Z_raw_right=rawZ,Z_facet_right=ZR,Z_facet_left_physical=ZL,
                F_outer_right_plus=FR[0],F_outer_right_minus=FR[1],F_outer_left_plus=FL[0],F_outer_left_minus=FL[1],
                F_mouth_right_in=r.Fin,F_mouth_right_out=r.Fout,F_mouth_left_in=self.Fin,F_mouth_left_out=self.Fout,
                B_right=BR,B_left_physical=BL,B_left_adjoint_prescription=BA,
                pressure_right=PR,pressure_left_physical=PL,shift_right=QR,shift_left_physical=QL,
                B_raw_right=raw,pressure_raw_right=dp,shift_raw_right=ds).items()}
        self.cache[key]=row;return row


def integrate(evaluator,order=8,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
    """Integrate two separately solved physical chiralities, never twice a saved value."""
    nodes,weights=np.polynomial.legendre.leggauss(order);edges=np.asarray(edges,float)
    fields=dict(energy='energy_integrand',rho_plus_pr='pressure_integrands',shift='shift_integrands',
                null_plus='null_plus_integrands',null_minus='null_minus_integrands')
    cases=dict(right='right',left_physical='left_physical',paired='actual_two_chirality_sum')
    panels=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        frequencies=(lo+hi)/2+(hi-lo)/2*nodes
        rows=[evaluator.evaluate(float(k),adjoint_control=False) for k in frequencies]
        values={name:{field:np.array([r[key][source] for r in rows]) for field,source in fields.items()}
                for name,key in cases.items()}
        integrated={name:{field:((hi-lo)/2*(weights@array)).tolist() for field,array in data.items()} for name,data in values.items()}
        panels.append(dict(lower=float(lo),upper=float(hi),kappa=frequencies.tolist(),
            integrands={name:{field:a.tolist() for field,a in data.items()} for name,data in values.items()},
            responses=integrated,minimum_sewing_singular_value=min(min(r[c]['minimum_singular_values'])
                for r in rows for c in ('right','left_physical'))))
    totals={name:{field:np.sum([p['responses'][name][field] for p in panels],axis=0).tolist() for field in fields} for name in cases}
    return dict(config=evaluator.config.copy(),order=int(order),edges=edges.tolist(),panels=panels,totals=totals,
        scope='Actual retained two-chirality response over this finite frequency window; finite dissipative approximation and physical raw-sphere maps.')


def resolve_quadrature(evaluator,coarse,fine,target=2e-7):
    """If needed, refine ONLY the worst panel; retain the failed coarse evidence.

    The sum of observed absolute panel changes is an operational convergence
    criterion, not a certified integration or continuum error bound.
    """
    fields=tuple(fine['totals']['paired'])
    changes=[{field:abs(np.asarray(b['responses']['paired'][field])-a['responses']['paired'][field])
              for field in fields} for a,b in zip(coarse['panels'],fine['panels'])]
    envelope=lambda rows:{field:np.sum([r[field] for r in rows],axis=0) for field in fields}
    maximum=lambda values:max(float(np.max(v)) for v in values.values())
    initial=envelope(changes);primary=fine;refinement=None
    orders=[fine['order']]*len(fine['panels'])
    if maximum(initial)>target:
        index=int(np.argmax([maximum(r) for r in changes]));panel=fine['panels'][index]
        stamp=time.perf_counter()
        refined=integrate(evaluator,24,(panel['lower'],panel['upper']))
        refined['seconds']=time.perf_counter()-stamp;new=refined['panels'][0]
        changes[index]={field:abs(np.asarray(new['responses']['paired'][field])-panel['responses']['paired'][field]) for field in fields}
        panels=list(fine['panels']);panels[index]=new;orders[index]=24
        totals={name:{field:np.sum([p['responses'][name][field] for p in panels],axis=0).tolist()
                      for field in fields} for name in fine['totals']}
        primary=dict(fine,order=None,panel_orders=orders,panels=panels,totals=totals,
            scope='Mixed quadrature: order16 except the explicitly identified order24 panel; actual separately solved chiralities.')
        refinement=dict(panel_index=index,lower=panel['lower'],upper=panel['upper'],before_order=16,
                        after_order=24,record=refined)
        print(json.dumps(dict(quadrature_refined_panel=index,edges=[panel['lower'],panel['upper']],
            order=24,totals=totals,seconds=refined['seconds'])),flush=True)
    final=envelope(changes)
    return primary,dict(target=target,initial_coarse_gate_passed=maximum(initial)<=target,
        initial_observed_panel_change_sum={k:v.tolist() for k,v in initial.items()},
        targeted_refinement=refinement,final_panel_orders=orders,
        final_observed_panel_change_sum={k:v.tolist() for k,v in final.items()},
        final_gate_passed=maximum(final)<=target,
        scope='Observed finite-window panel convergence only, without interval arithmetic or a bound on omitted frequency tails.')


def checks():
    start=time.perf_counter();root=Path(__file__).resolve().parent
    inputs=('throat_null_response.py','throat_null_response_results.json','phase_vacuum_contrast.py',
        'curved_throat_modes.py','curved_exterior_closure.py','imaginary_exterior.py','imaginary_weyl_leads.py',
        'exterior_boundary_map.py','exterior_tail_transport.py','magnetic_tail_map.py','spherical_weyl_leads.py','two_mouth_mesh.py')
    hashes={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs}
    records=[];evaluators={}
    for h in (.08,.05,.035):
        case=DiracAudit(h=h);evaluators[h]=case
        rows=[case.evaluate(k,save_matrices=True) for k in (.01,.4,2.)]
        records.append(dict(config=case.config,assembly=case.assembly,ports=case.port_checks,samples=rows))
        print(json.dumps(dict(h=h,assembly=case.assembly,discrepancies=[r['discrepancies'] for r in rows])),flush=True)
    all_rows=[r for c in records for r in c['samples']]
    previous=json.loads((project_path('throat_null_response_results.json')).read_text(encoding='utf-8'))
    preservation=[]
    for now,old in zip(records[1]['samples'],previous['finite_difference_checks']):
        if now['kappa']!=old['kappa']:raise AssertionError('Frozen96 comparison frequencies changed')
        oldS=old['audit']['matrices']['S_exterior']
        oldS=np.asarray(oldS['real'])+1j*np.asarray(oldS['imag'])
        newS=now['matrices']['S_right'];newS=np.asarray(newS['real'])+1j*np.asarray(newS['imag'])
        preservation.append(dict(kappa=now['kappa'],exterior_matrix_error=float(np.max(abs(newS-oldS))),
            rho_plus_pr_error=max(abs(np.array(now['right']['pressure_integrands'])-[p['response_integrand'] for p in old['probes']])),
            shift_error=max(abs(np.array(now['right']['shift_integrands'])-[p['shift_response_integrand'] for p in old['probes']]))))
    preservation=[{k:float(v) for k,v in row.items()} for row in preservation]
    integrations=[]
    for h,order in ((.05,8),(.05,16),(.035,8)):
        stamp=time.perf_counter();integrated=integrate(evaluators[h],order)
        integrated['seconds']=time.perf_counter()-stamp;integrations.append(integrated)
        print(json.dumps(dict(h=h,order=order,totals=integrated['totals'],seconds=integrated['seconds'])),flush=True)
    coarse,fine,finemesh=integrations
    primary,resolved=resolve_quadrature(evaluators[.05],coarse,fine)
    integral_fields=('energy','rho_plus_pr','shift','null_plus','null_minus')
    quadrature={field:(np.asarray(fine['totals']['paired'][field])-coarse['totals']['paired'][field]).tolist() for field in integral_fields}
    mesh_difference={field:(np.asarray(finemesh['totals']['paired'][field])-coarse['totals']['paired'][field]).tolist() for field in integral_fields}
    every_row=[r for case in evaluators.values() for r in case.cache.values()]
    trend=[]
    for record in records:
        samples=record['samples']
        trend.append(dict(h=record['config']['h'],
            maximum_energy_mismatch=max(abs(r['discrepancies']['energy_physical_minus_right']) for r in samples),
            maximum_rho_plus_pr_mismatch=max(abs(v) for r in samples for v in r['discrepancies']['pressure_physical_minus_right']),
            maximum_shift_sum=max(abs(v) for r in samples for v in r['discrepancies']['shift_physical_plus_right']),
            maximum_outer_Gram_defect=max(r['discrepancies']['outer_physical_map_adjoint_defect'] for r in samples),
            maximum_mouth_Gram_defect=max(r['discrepancies']['mouth_physical_map_adjoint_defect'] for r in samples)))
    tests=dict(independent_opposite_assembly=max(c['assembly']['opposite_adjoint_relative_error'] for c in records)<1e-12,
        identical_regular_axis_domain=all(c['assembly']['axis_domain_equal'] for c in records),
        independently_rebuilt_swapped_ports=max(max(p['incoming_vs_right_outgoing_adjoint'],p['outgoing_vs_right_incoming_adjoint'],p['swapped_factor_error'],p['reproduction_error'])
            for c in records for p in c['ports'])<1e-12,
        direct_sparse_response_equations=max(max(r['exterior_equation_relative_residuals']) for r in every_row)<1e-11,
        deliberately_adjoint_paired_finite_identity=max(max(r['discrepancies'][k] for k in (
            'exterior_adjoint_prescription_error','adjoint_energy_error','adjoint_pressure_error','adjoint_shift_error')) for r in all_rows)<2e-9,
        sampled_physical_sewings_nonsingular=min(min(r['left_physical']['minimum_singular_values']) for r in every_row)>1e-7,
        frozen96_right_response_reproduced=max(max(r[k] for k in ('exterior_matrix_error','rho_plus_pr_error','shift_error')) for r in preservation)<1e-12,
        finite_window_quadrature=resolved['final_gate_passed'],
        physical_sewing_mismatch_envelopes_shrink=all(trend[j+1][field]<trend[j][field] for j in (0,1)
            for field in ('maximum_energy_mismatch','maximum_rho_plus_pr_mismatch','maximum_shift_sum','maximum_outer_Gram_defect','maximum_mouth_Gram_defect')),
        frozen_dependencies_preserved=hashes=={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs})
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,stage='Selected-frequency physical-sewing audit plus an actual finite-window two-chirality integral.',
        parameters=dict(q=1,mu=.5,b=.2,d=1.,R=2.,nmax=8,L=40.,rho_c=2.,frequencies=[.01,.4,2.],
            phase_difference='pi minus0',units='hbar=c=d=1'),mesh_refinement=records,
        representative=records[1]['samples'][1],frozen96_reproduction=preservation,
        quadrature_refinement=integrations[:2],primary_quadrature=primary,quadrature_resolution=resolved,
        integrated_mesh_refinement=[integrations[0],integrations[2]],
        pointwise_mesh_trend=trend,
        diagnostics=dict(primary_two_chirality_integrals=primary['totals']['paired'],
            finer_mesh_two_chirality_integrals=finemesh['totals']['paired'],
            frequency_quadrature_differences=quadrature,integrated_mesh_differences=mesh_difference,
            minimum_physical_sewing_singular_value=min(min(r['left_physical']['minimum_singular_values']) for r in every_row),
            maximum_sparse_equation_residual=max(max(r['exterior_equation_relative_residuals']) for r in every_row),
            frequency_geometry_evaluations=len(every_row),combined_continuum_bound=None,omitted_IR_bound=None,omitted_UV_bound=None),
        observable_definitions=dict(pressure_integrands='Smeared energy-plus-radial-pressure PHASE DIFFERENCE, not radial pressure alone.',
            shift_integrands='Minus the smeared lower orthonormal radial energy-flux PHASE DIFFERENCE.',
            null_plus_integrands='pressure_integrands - 2*shift_integrands, k_hat=(1,+1).',
            null_minus_integrands='pressure_integrands + 2*shift_integrands, k_hat=(1,-1).'),
        equations=dict(bare='A_L(kappa) is independently assembled; compare with A_R(kappa)^dagger only as a check.',
            physical_outer='Z_L,facet=F_R,plus Z_R,raw^dagger F_R,minus^(-1).',
            physical_mouth='B_L,facet=F_R,out B_R,raw^dagger F_R,in^(-1).',
            derivatives='dB_L,pressure=(dB_R,pressure)^dagger in RAW coordinates; dB_L,shift=-(dB_R,shift)^dagger, then map factors.',
            adjoint_prescription='Separately: Z_L=Z_R,facet^dagger and B_L=B_R,facet^dagger. Exact finite algebra does not make this the physical projected-sphere map.'),
        scope=['Opposite transport is actually assembled and solved; its response is not obtained by doubling the saved right-Weyl values.',
            'Physical raw spherical maps are adjointed before independent facet Gram conversion. The alternative finite adjoint prescription is clearly separated.',
            'The same-sign facet projection neglects cross-sign chord geometry effects; their continuum removal needs mesh control.',
            'Three frequencies have full saved matrix audits; a separate finite-window integral solves BOTH physical chiralities. It is not obtained by doubling increment96.',
            'The failed order8 energy-resolution gate is retained. Only its dominant panel is refined to order24; the exact mixed primary quadrature and remaining observed panel changes are explicit.',
            'The integral keeps one q=1, mu=.5 sector and a finite angular cutoff. Other azimuthal sectors, the absolute source and the complete Dirac determinant remain missing.',
            'Shrinking signed physical-map discrepancies have cancellations; no convergence order or rigorous continuum/tail error bound is inferred. The adjoint-paired finite prescription alone does not prove self-adjointness.',
            'The synthetic throat, endpoint identification, magnetic exterior, state and source limitations of increments93--96 remain; no physical junction, absolute stress, stable medium or time travel is established.'],
        input_source_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.perf_counter()-start)
    return result


if __name__=='__main__':
    result=checks();path=project_path('dirac_sewing_audit_results.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(path.read_text(encoding='utf-8'))==result
    print(json.dumps({k:result[k] for k in ('passed','tests','seconds','source_sha256')},indent=2))
    if not result['passed']:raise SystemExit('An opposite-Weyl sewing control failed; diagnostics retained.')
