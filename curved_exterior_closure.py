"""Finite curved-throat mode sewing on the computed imaginary exterior.

The synthetic throat has sphere radius b and lapse one at its two endpoints.
This fixes its induced endpoint metric only, not a physical junction or the
Einstein/matter equations.  A unit boundary forcing remains after sewing, so
the correct attenuation ledger uses source_work=1+2 Re(looped_upper_input).
"""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import block_diag

from exterior_boundary_map import _complex_json
from imaginary_exterior import build_problem, attenuation_ledger
from spherical_weyl_leads import mouth_characteristics, zero_mode


def mouth_factors(problem):
    """Reproduce frozen endpoint Gram coordinates without changing them."""
    b, d, nmax = .2, 1., problem['nmax']
    factors = {}
    reproduction_error = 0.
    for port in problem['ports'][:2]:
        raw = [[], []]
        for edge in port['edges']:
            points = edge['points']
            theta = np.arctan2(points[:,0], points[:,1]-(d/2 if port['kind']=='upper' else -d/2))
            plus,minus = mouth_characteristics(theta,b,port['kind'],nmax)
            zero = zero_mode(theta,b,port['kind'])[:,:,None]
            if port['kind']=='upper':
                minus = np.concatenate((zero,minus),axis=2)
            else:
                plus = np.concatenate((zero,plus),axis=2)
            for index,(basis,projector) in enumerate(((plus,'pp'),(minus,'pm'))):
                raw[index].append(np.einsum('ab,pbc->pac',edge[projector],np.moveaxis(basis,0,1)))
        for index,sense in enumerate(('plus','minus')):
            basis = np.concatenate(raw[index])
            gram = np.einsum('p,pai,paj->ij',port['weights'],basis.conj(),basis)
            factor = np.linalg.cholesky((gram+gram.conj().T)/2).conj().T
            transformed = basis@np.linalg.solve(factor,np.eye(len(factor)))
            reproduction_error = max(reproduction_error,float(np.max(abs(transformed-port[sense]))))
            factors[port['kind']+'_'+sense] = factor
    incoming = block_diag(factors['upper_minus'],factors['lower_minus'])
    outgoing = block_diag(factors['upper_plus'],factors['lower_plus'])
    return incoming,outgoing,reproduction_error


@lru_cache(maxsize=32)
def throat_blocks(kappa,nmax,rho_c=2.,profile='secant'):
    from curved_throat_modes import throat_scattering
    if not math.isfinite(rho_c) or rho_c<=0:
        raise ValueError('rho_c must be finite and positive')
    b=.2
    ell=b*math.sqrt(1+rho_c*rho_c)
    length=2*ell*math.atan(rho_c)
    rows=[]
    for n in range(1,nmax+1):
        S,diag=throat_scattering(1j*kappa,n,ell=ell,rho_c=rho_c,q=1,profile=profile,
                               b=b if profile=='constant' else None)
        if not np.all(np.isfinite(S)) or np.shape(S)!=(2,2):
            raise ArithmeticError('Invalid paired throat scattering matrix')
        rows.append((S,diag))
    return rows,dict(b=b,rho_c=rho_c,ell=ell,optical_length=length,
                     endpoint_lapse=b*math.sqrt(1+rho_c*rho_c)/ell,profile=profile)


def raw_mouth_map(problem,phase=0.,rho_c=2.,kind='curved'):
    """Canonical map: outgoing [upper_n,lower_0,lower_n] -> incoming
    [upper_0,upper_n,lower_n].  Endpoint spin/gauge identification is assumed.
    """
    if kind not in ('curved','zero_only','constant') or not math.isfinite(phase):
        raise ValueError('Choose curved, zero_only or constant and finite phase')
    kappa,nmax=problem['kappa'],problem['nmax']
    b=.2
    ell=b*math.sqrt(1+rho_c*rho_c)
    length=2*ell*math.atan(rho_c)
    B=np.zeros((2*nmax+1,2*nmax+1),complex)
    forward=np.exp(1j*phase)
    B[0,nmax]=forward*np.exp(-kappa*length)
    scattering=[]
    if kind=='zero_only':
        for n,beta in enumerate(problem['beta']):
            B[n+1,n]=beta
            B[nmax+1+n,nmax+1+n]=beta
        geometry=dict(b=b,rho_c=rho_c,ell=ell,optical_length=length,endpoint_lapse=1.,
                      profile='zero mode only; old semi-infinite product relation for paired modes')
    else:
        rows,geometry=throat_blocks(kappa,nmax,rho_c,'secant' if kind=='curved' else 'constant')
        for n,(S,diag) in enumerate(rows):
            # S maps [lower outgoing,upper outgoing] to
            # [lower incoming,upper incoming].  Reflections stay unphased.
            B[n+1,n]=S[1,1]
            B[n+1,nmax+1+n]=forward*S[1,0]
            B[nmax+1+n,n]=forward.conjugate()*S[0,1]
            B[nmax+1+n,nmax+1+n]=S[0,0]
            scattering.append(dict(n=n+1,S=_complex_json(S),diagnostics=diag))
    return B,geometry,scattering


def solve_closure(problem,phase=0.,rho_c=2.,kind='curved',save_matrices=False):
    """Solve a driven sewn resolvent, never a unit transmission experiment."""
    started=time.perf_counter()
    model=problem['model']
    upper,lower,outer=problem['ports']
    nmax=problem['nmax']
    Braw,geometry,paired=raw_mouth_map(problem,phase,rho_c,kind)
    Fin,Fout,reproduction=mouth_factors(problem)
    B=np.linalg.solve(Fout.T,(Fin@Braw).T).T
    Um,sm,Vhm=np.linalg.svd(B)
    Ut,st,Vht=np.linalg.svd(problem['facet_Z'])
    Gmouth=np.column_stack((upper['G'],lower['G']))
    Cmouth=np.vstack((upper['C'],lower['C']))
    G=np.column_stack((Gmouth@Um,outer['G']@Ut))
    C=np.vstack((Vhm@Cmouth,Vht@outer['C']))
    singular=np.r_[sm,st]
    solution=model.solve(problem['source'],(G,singular,C))
    outgoing=[p['C']@solution for p in problem['ports']]
    c=np.r_[outgoing[0],outgoing[1]]
    looped=B@c
    incoming_mouth=looped.copy()
    incoming_mouth[0]+=1.
    incoming=[incoming_mouth[:nmax+1],incoming_mouth[nmax+1:],problem['facet_Z']@outgoing[2]]
    data=model.modal_incoming(problem['ports'],incoming)
    rhs=sum((p['G']@v for p,v in zip(problem['ports'],incoming)),np.zeros(model.ndof,complex))
    ledger=attenuation_ledger(model,problem['mass'],problem['kappa'],solution,data,rhs)
    omitted={p['kind']:float(ledger['outgoing'][p['kind']]-np.vdot(v,v).real)
             for p,v in zip(problem['ports'],outgoing)}
    source_work=float(1+2*looped[0].real)
    mouth_net=float(np.vdot(c,c).real-np.vdot(looped,looped).real)
    raw_out=np.linalg.solve(Fout,c)
    raw_in=Braw@raw_out
    raw_net=float(np.vdot(raw_out,raw_out).real-np.vdot(raw_in,raw_in).real)
    mouth_geometry=float((np.vdot(raw_out,(Fout.conj().T@Fout-np.eye(len(B)))@raw_out)
                           -np.vdot(raw_in,(Fin.conj().T@Fin-np.eye(len(B)))@raw_in)).real)
    outer_net=float(np.vdot(outgoing[2],outgoing[2]).real-np.vdot(incoming[2],incoming[2]).real)
    out_sphere=np.linalg.solve(problem['gram_factors'][0],outgoing[2])
    in_sphere=problem['spherical_Z']@out_sphere
    raw_outer_net=float(np.vdot(out_sphere,out_sphere).real-np.vdot(in_sphere,in_sphere).real)
    outer_geometry=outer_net-raw_outer_net
    reduced=mouth_net+outer_net+sum(omitted.values())+ledger['numerical_loss']+ledger['bulk_attenuation']-source_work
    minraw=float(np.linalg.eigvalsh(np.eye(len(Braw))-Braw.conj().T@Braw)[0])
    minfacet=float(np.linalg.eigvalsh(np.eye(len(B))-B.conj().T@B)[0])
    response=problem['receiver']@solution
    result=dict(kappa=problem['kappa'],h=problem['h'],R=problem['R'],nmax=nmax,L=problem['L'],
        phase=phase,kind=kind,geometry=geometry,driven_lower_response=_complex_json(response),
        driven_response_magnitude_squared=float(abs(response)**2),
        looped_upper_input=_complex_json(looped[0]),source_work=source_work,
        facet_mouth_net=mouth_net,canonical_throat_net=raw_net,
        signed_mouth_geometry_difference=mouth_geometry,mouth_geometry_identity_error=abs(mouth_net-raw_net-mouth_geometry),
        facet_outer_net=outer_net,spherical_outer_net=raw_outer_net,signed_outer_geometry_difference=outer_geometry,
        omitted_outgoing_norm=omitted,ledger=ledger,reduced_source_work_balance_residual=float(reduced),
        raw_mouth_minimum_attenuation_eigenvalue=minraw,facet_mouth_minimum_attenuation_eigenvalue=minfacet,
        mouth_basis_reproduction_error=reproduction,
        mouth_coordinate_identity_error=float(np.linalg.norm(B@Fout-Fin@Braw)),
        full_feedback_svd_error=float(np.linalg.norm((Um*sm)@Vhm-B)+np.linalg.norm((Ut*st)@Vht-problem['facet_Z'])),
        outgoing_coefficients={p['kind']:_complex_json(v) for p,v in zip(problem['ports'],outgoing)},
        incoming_coefficients={p['kind']:_complex_json(v) for p,v in zip(problem['ports'],incoming)},
        seconds=time.perf_counter()-started,
        interpretation='Driven imaginary-energy sewn response; squared magnitudes are NOT probabilities.')
    if save_matrices:
        result['boundary_matrices']=dict(B_raw=_complex_json(Braw),B_facet=_complex_json(B),
            F_in=_complex_json(Fin),F_out=_complex_json(Fout),paired_scattering=paired,
            outgoing_order='upper paired n1..N; lower zero; lower paired n1..N',
            incoming_order='upper zero; upper paired n1..N; lower paired n1..N',
            phase_forward='lower to upper exp(+i*phase)',phase_reverse='upper to lower exp(-i*phase)')
    return result


def response(row):
    v=row['driven_lower_response']
    return complex(v['real'],v['imag'])


def phase_placement_check(problem,rho_c=2.):
    B0,_,_=raw_mouth_map(problem,0.,rho_c)
    phase=.371
    B,_,_=raw_mouth_map(problem,phase,rho_c)
    N=problem['nmax']
    inc=np.r_[np.full(N+1,np.exp(1j*phase)),np.ones(N)]
    out=np.r_[np.full(N,np.exp(1j*phase)),np.ones(N+1)]
    target=inc[:,None]*B0*out.conj()[None,:]
    return float(np.max(abs(B-target)))


def checks():
    started=time.perf_counter()
    problems={}
    records={}

    def problem(kappa=.4,h=.05,R=2.,nmax=12):
        key=(kappa,h,R,nmax)
        if key not in problems:
            problems[key]=build_problem(kappa,h,R,nmax,40.)
        return problems[key]

    def case(kappa=.4,h=.05,R=2.,nmax=12,phase=np.pi,kind='curved',save=False):
        key=(kappa,h,R,nmax,phase,kind,save)
        if key not in records:
            records[key]=solve_closure(problem(kappa,h,R,nmax),phase,2.,kind,save)
            r=records[key]
            print(f"kappa={kappa:g} h={h:g} R={R:g} N={nmax} phase={phase:g} {kind}: response={response(r)}",flush=True)
        return records[key]

    comparisons=[]
    for kappa,phase in ((.1,np.pi),(.4,np.pi),(1.,np.pi),(.4,0.)):
        curved=case(kappa=kappa,phase=phase)
        zero=case(kappa=kappa,phase=phase,kind='zero_only')
        comparisons.append(dict(kappa=kappa,phase=phase,curved=curved,zero_only=zero,
            complex_response_difference=abs(response(curved)-response(zero)),
            relative_response_difference=abs(response(curved)-response(zero))/abs(response(zero))))
    mesh_rows=[case(h=h) for h in (.08,.05,.035)]
    paired_mesh=[]
    for curved in mesh_rows:
        zero=case(h=curved['h'],kind='zero_only')
        paired_mesh.append(dict(h=curved['h'],curved=curved,zero_only=zero,
            complex_response_difference=abs(response(curved)-response(zero)),
            signed_real_response_difference=float((response(curved)-response(zero)).real)))
    radius_rows=[case(R=R) for R in (2.,4.)]
    angular_rows=[case(nmax=N) for N in (8,12)]
    representative=case(h=.035,R=4.,save=True)
    representative_zero=case(h=.035,R=4.,kind='zero_only')
    constant=case(kind='constant')
    phase_error=phase_placement_check(problem())
    rows=list(records.values())
    tests=dict(
        actual_input_attenuation_identity=max(r['ledger']['attenuation_identity_error'] for r in rows)<1e-10,
        driven_source_work_identity=max(abs(r['reduced_source_work_balance_residual']) for r in rows)<1e-10,
        sparse_equation_residual=max(r['ledger']['linear_residual_relative'] for r in rows)<1e-10,
        canonical_throat_contraction=min(r['raw_mouth_minimum_attenuation_eigenvalue'] for r in rows)>-1e-10,
        signed_mouth_geometry_identity=max(r['mouth_geometry_identity_error'] for r in rows)<1e-11,
        reproduced_frozen_mouth_bases=max(r['mouth_basis_reproduction_error'] for r in rows)<1e-11,
        full_coordinate_map=max(r['mouth_coordinate_identity_error'] for r in rows)<1e-11,
        full_matrix_svd=max(r['full_feedback_svd_error'] for r in rows)<1e-11,
        opposite_cross_transmission_phases=phase_error<1e-12,
        positive_source_work=all(r['source_work']>0 for r in rows),
        modal_remainders_nonnegative=min(min(r['omitted_outgoing_norm'].values()) for r in rows)>-1e-10,
        mesh_loss_refines=all(a['ledger']['numerical_loss']>b['ledger']['numerical_loss'] for a,b in zip(mesh_rows[:-1],mesh_rows[1:])),
        actual_paired_curved_response_differs=max(r['complex_response_difference'] for r in comparisons)>1e-7,
    )
    assert all(tests.values()),tests
    root=Path(__file__).resolve().parent
    inputs=('imaginary_exterior.py','curved_throat_modes.py','imaginary_weyl_leads.py',
            'exterior_boundary_map.py','exterior_tail_transport.py','spherical_weyl_leads.py','two_mouth_mesh.py')
    return dict(passed=True,tests=tests,parameters=dict(q=1,mu=.5,b=.2,d=1.,rho_c=2.,energy='i*kappa'),
        representative=representative,representative_zero_only=representative_zero,
        curved_vs_zero_only=comparisons,mesh_refinement=mesh_rows,paired_mesh_refinement=paired_mesh,
        radius_refinement=radius_rows,angular_refinement=angular_rows,finite_constant_profile=constant,
        phase_placement_error=phase_error,
        sensitivities=dict(fine_mesh_change=abs(response(mesh_rows[-1])-response(mesh_rows[-2])),
            matching_radius_change=abs(response(radius_rows[1])-response(radius_rows[0])),
            angular_change=abs(response(angular_rows[1])-response(angular_rows[0])),
            curved_vs_finite_constant_change=abs(response(constant)-response(case())),
            certified_continuum_or_infinite_tail_bound=False),
        equations=dict(source='incoming_mouth=e_upper_zero+B_facet*outgoing_mouth',
            map='B_facet=F_in*B_raw*F_out^(-1)',
            source_work='1+2*Re((B_facet*outgoing_mouth)[upper_zero])',
            balance='source_work=mouth_net+outer_net+omitted_modes+DG_loss+2*kappa*bulk_norm',
            optical='ell=b*sqrt(1+rho_c^2); Lopt=2*ell*atan(rho_c); Ncut=1',
            phases='lower->upper exp(+i*phase), upper->lower exp(-i*phase), reflections unchanged'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs},
        seconds=time.perf_counter()-started,
        scope=['Finite curved throat paired modes coupled to the actual connected flat q=1 exterior in the specified approximation.',
          'Synthetic b=.2 throat matches endpoint lapse and sphere radius only; its junction stress and bulk Einstein source are NOT supplied.',
          'Endpoint spin/gauge/angular identification and phases 0/pi are declared reference choices, not derived physical collar data.',
          'Both paired-mode and zero-only comparisons use the same optical length, phase, forcing, tail and Gram conversion.',
          'The field is driven after sewing: source work, not a unit incoming-flux normalization, closes the attenuation ledger.',
          'Positive imaginary-energy squared response can exceed one and is not a probability.',
          'Tail truncation, same-sign polygon projection, omitted angular modes and mesh errors remain recorded approximations.',
          'No renormalized quantum determinant/stress, physical MMP junction, stable material, finite payload traversal or backward time travel is proved.'])


if __name__=='__main__':
    result=checks()
    target=Path(__file__).with_name('curved_exterior_closure_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({key:result[key] for key in ('passed','tests','seconds','source_sha256')},indent=2))
