"""Compact lapse variation of the retained phase-contrast determinant.

Spatial geometry, magnetic potential, cuts, throat and endpoint phases are
fixed. N=1+eta*f equals one near every cut. After the Weyl congruence the
weak operator is D+kappa*M_(1/N), so its derivative is -kappa*M_f.
This computes a smeared PHASE DIFFERENCE, not an absolute stress, pressure,
NEC test, physical tracer calculation or self-consistent supporting source.
The finite DG sector and its continuum/renormalization limitations remain.

Run: .venv\\Scripts\\python.exe -B -W error local_lapse_response.py
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.sparse import block_diag
from scipy.sparse.linalg import splu

from curved_exterior_closure import raw_mouth_map
from exterior_boundary_map import ID, _complex_json, _quadrature
from exterior_tail_transport import facet_map
from imaginary_weyl_leads import tail_map
from phase_vacuum_contrast import PhaseContrast


REGIONS=(dict(name='between_mouths',center=(.35,0.),radius=.20),
         dict(name='beside_upper_mouth',center=(.50,.50),radius=.18),
         dict(name='outer_midplane',center=(.80,0.),radius=.25))


def bump(points,region):
    """Peak-one C-infinity bump exp(1-1/(1-u)), with zero outside u<1."""
    points=np.asarray(points,float)
    u=np.sum((points-np.asarray(region['center']))**2,axis=-1)/region['radius']**2
    result=np.zeros_like(u)
    inside=u<1.
    result[inside]=np.exp(1-1/(1-u[inside]))
    return result


def weighted_mass(model,weight,order=12):
    """Assemble the actual P1 spinor mass integral rho*weight, not nodal lumping."""
    bary,qw,_,_=_quadrature(order)
    blocks=[]
    integrated=0.
    for tri in model.triangles:
        xy=model.vertices[tri]
        jac=np.linalg.det(np.stack((xy[1]-xy[0],xy[2]-xy[0])))
        points=bary@xy
        values=np.asarray(weight(points),float)
        if values.shape!=(len(bary),) or not np.all(np.isfinite(values)):
            raise ValueError('Mass weight must be finite at every quadrature point')
        weights=qw*jac*points[:,0]*values
        local=np.einsum('p,pa,pb->ab',weights,bary,bary)
        blocks.append(np.kron(local,ID))
        integrated+=float(np.sum(weights))
    matrix=block_diag(blocks,format='csc')
    matrix.eliminate_zeros()
    return matrix,integrated


def _segment_distance(point,ends):
    direction=ends[1]-ends[0]
    t=np.dot(point-ends[0],direction)/np.dot(direction,direction)
    return float(np.linalg.norm(point-(ends[0]+np.clip(t,0.,1.)*direction)))


def region_geometry(model,region):
    center=np.asarray(region['center'])
    clearances={kind:math.inf for kind in ('axis','upper','lower','outer')}
    for edge in model.boundary:
        distance=_segment_distance(center,model.vertices[list(edge['vertices'])])-region['radius']
        clearances[edge['kind']]=min(clearances[edge['kind']],distance)
    if min(clearances.values())<=0:
        raise ValueError('The entire bump support must lie strictly inside the numerical domain')
    radial=quad(lambda u:math.exp(1-1/(1-u)) if u<1 else 0.,0.,1.,epsabs=2e-14,epsrel=2e-14)[0]
    reduced=math.pi*center[0]*region['radius']**2*radial
    return dict(name=region['name'],center=list(region['center']),radius=region['radius'],
        formula='u=((rho-rho0)^2+(z-z0)^2)/a^2; f=exp(1-1/(1-u)) for u<1, otherwise0',
        peak_weight=1.,polygon_clearances=clearances,
        exact_meridional_weight=reduced,exact_full_three_dimensional_weight=2*math.pi*reduced,
        normalization='The fixed-mu spinor already has normalized azimuthal dependence; do not multiply its trace by2pi or divide by this volume to claim stress.')


class LocalLapse:
    """Reuse one frozen exterior geometry and both Green factors for all bumps."""

    def __init__(self,h=.05,R=2.,nmax=8,L=40.,mass_order=16):
        if isinstance(mass_order,bool) or int(mass_order)!=mass_order or mass_order<4:
            raise ValueError('mass_order must be an integer >=4')
        self.base=PhaseContrast(h,R,nmax,L)
        self.config=dict(self.base.config,mass_order=int(mass_order))
        self.active=self.base.model.active
        self.region_info=[region_geometry(self.base.model,r) for r in REGIONS]
        self.masses=[]
        for region,info in zip(REGIONS,self.region_info):
            matrix,weight=weighted_mass(self.base.model,lambda p:bump(p,region),mass_order)
            self.masses.append(matrix[self.active][:,self.active])
            info['quadrature_meridional_weight']=weight
            info['weight_relative_error']=abs(weight/info['exact_meridional_weight']-1)
        self.cache={}

    def _maps(self,kappa):
        base=self.base
        rawZ,diag=tail_map(kappa,base.config['R'],base.config['L'],nmax=base.config['nmax'])
        Z,_,geometry=facet_map(base.outer,base.outer_basis,rawZ)
        base.problem['kappa']=kappa
        Bs=[]
        for phase in (0.,math.pi):
            rawB,_,_=raw_mouth_map(base.problem,phase,base.config['rho_c'],'curved')
            Bs.append(np.linalg.solve(base.Fout.T,(base.Fin@rawB).T).T)
        return Z,Bs,diag,geometry

    def _green(self,matrix,Z,transpose=False):
        base=self.base
        lu=splu(matrix)
        Gm,Go=base.Gm[self.active],base.outer['G'][self.active]
        Cm,Co=base.Cm[:,self.active],base.outer['C'][:,self.active]
        count=Gm.shape[1]
        lifted=lu.solve(np.column_stack((Gm,Go)))
        Xm,Xo=lifted[:,:count],lifted[:,count:]
        Qmo,Qoo=Cm@Xo,Co@Xo
        M=np.eye(len(Z))-Qoo@Z
        X=Xm+Xo@Z@np.linalg.solve(M,Co@Xm)
        S=Cm@X
        row=dict(lu=lu,X=X,S=S)
        if transpose:
            # Ordinary transpose is required for the row of an inverse.
            # A fully conjugated adjoint solve is equivalent; Y=X^dagger is not.
            rows=lu.solve(np.vstack((Cm,Co)).T,trans='T').T
            Ym,Yo=rows[:count],rows[count:]
            Y=Ym+Qmo@Z@np.linalg.solve(M,Yo)
            column_error=matrix@X-Gm-Go@Z@(Co@X)
            row_error=(matrix.T@Y.T).T-Cm-(Y@Go)@Z@Co
            row.update(Y=Y,column_residual=float(np.linalg.norm(column_error)/np.linalg.norm(Gm)),
                transpose_residual=float(np.linalg.norm(row_error)/np.linalg.norm(Cm)))
        return row

    def _determinants(self,S,Bs):
        output=[]
        for B in Bs:
            M=np.eye(len(B))-S@B
            sign,logabs=np.linalg.slogdet(M)
            singular=np.linalg.svd(M,compute_uv=False)
            if sign==0 or not math.isfinite(logabs):
                raise ArithmeticError('Singular or nonfinite sewing determinant')
            output.append(dict(logabs=float(logabs),phase=_complex_json(sign),
                minimum_singular_value=float(singular[-1]),condition_number=float(singular[0]/singular[-1])))
        return output

    def evaluate(self,kappa,audit=False):
        kappa=float(kappa)
        if not math.isfinite(kappa) or kappa<=0:
            raise ValueError('kappa must be finite and positive')
        key=(kappa,bool(audit))
        if key in self.cache:return self.cache[key]
        started=time.perf_counter()
        base=self.base
        matrix=(base.unshifted+kappa*base.mass)[self.active][:,self.active].tocsc()
        Z,Bs,tail,geometry=self._maps(kappa)
        green=self._green(matrix,Z,transpose=True)
        X,Y,S=(green[k] for k in ('X','Y','S'))
        Ks=[np.linalg.solve((np.eye(len(B))-S@B).T,B.T).T for B in Bs]
        K=Ks[1]-Ks[0]
        reduced=[Y@(mass@X) for mass in self.masses]
        traces=[np.trace(K@small) for small in reduced]
        values=[float(kappa/math.pi*t.real) for t in traces]
        row=dict(kappa=kappa,regions=[dict(name=r['name'],green_trace=_complex_json(t),response_integrand=value)
                      for r,t,value in zip(REGIONS,traces,values)],
            column_equation_relative_residual=green['column_residual'],
            transpose_equation_relative_residual=green['transpose_residual'],
            determinants=self._determinants(S,Bs),outer_geometry=geometry,
            tail_accepted_steps=tail['accepted_steps'],seconds=time.perf_counter()-started)
        if audit:
            fd=[]
            for region,exact,reduced_matrix in zip(REGIONS,values,reduced):
                estimates=[]
                for eta in (.02,.01,.005):
                    logs=[];Smaps=[]
                    for signed in (-eta,eta):
                        # Stable correction to the exact baseline mass. Use
                        # the same non-lumped quadrature as the derivative.
                        delta,_=weighted_mass(base.model,
                            lambda p:-signed*bump(p,region)/(1+signed*bump(p,region)),self.config['mass_order'])
                        perturbed=matrix+kappa*delta[self.active][:,self.active]
                        trial=self._green(perturbed,Z)
                        det=self._determinants(trial['S'],Bs)
                        logs.append(det[1]['logabs']-det[0]['logabs'])
                        Smaps.append(trial['S'])
                    derivative=-(logs[1]-logs[0])/(2*eta*math.pi)
                    map_fd=(Smaps[1]-Smaps[0])/(2*eta)
                    estimates.append(dict(eta=eta,minimum_lapse=1-eta,
                        response_integrand_derivative=float(derivative),
                        absolute_error=abs(float(derivative)-exact),
                        relative_error=abs(float(derivative)-exact)/max(abs(exact),1e-15),
                        exterior_map_derivative_relative_error=float(np.linalg.norm(map_fd-kappa*reduced_matrix)
                            /max(np.linalg.norm(kappa*reduced_matrix),1e-15))))
                fd.append(dict(name=region['name'],analytic_response=exact,finite_differences=estimates))
                richardson=(4*estimates[-1]['response_integrand_derivative']
                            -estimates[-2]['response_integrand_derivative'])/3
                fd[-1].update(richardson_response=float(richardson),richardson_absolute_error=abs(richardson-exact),
                    richardson_relative_error=abs(richardson-exact)/max(abs(exact),1e-15),
                    successive_absolute_error_ratios=[estimates[i]['absolute_error']/max(estimates[i+1]['absolute_error'],1e-30)
                                                       for i in range(len(estimates)-1)])
            # Independent direct trace contraction in the DG space, without
            # materializing the full difference of Green matrices.
            contractions=[]
            for mass,small,trace in zip(self.masses,reduced,traces):
                alternative=np.einsum('ij,ji->',X@K,Y@mass)
                contractions.append(float(abs(alternative-trace)))
            Cm=base.Cm[:,self.active]
            ordinary=green['lu'].solve(Cm.T,trans='T').T
            adjoint=green['lu'].solve(Cm.conj().T,trans='H').conj().T
            row['audit']=dict(finite_differences=fd,
                direct_DG_trace_errors=contractions,
                transpose_vs_fully_conjugated_adjoint=float(np.linalg.norm(ordinary-adjoint)/np.linalg.norm(ordinary)),
                matrices=dict(S_exterior=_complex_json(S),K_zero=_complex_json(Ks[0]),K_pi=_complex_json(Ks[1]),
                    K_difference=_complex_json(K),B_zero=_complex_json(Bs[0]),B_pi=_complex_json(Bs[1]),
                    regional_Y_Mf_X={r['name']:_complex_json(v) for r,v in zip(REGIONS,reduced)}))
        row['seconds']=time.perf_counter()-started
        self.cache[key]=row
        return row


def mass_controls(evaluator,kappa=.4):
    """Constant polynomial, circular-support and nonuniform resolvent checks."""
    base=evaluator.base;model=base.model
    constant,_=weighted_mass(model,lambda p:np.ones(len(p)),4)
    constant_error=np.linalg.norm((constant-base.mass).data)/np.linalg.norm(base.mass.data)
    analytic_error=0.
    for tri,cell in zip(model.triangles,model.cells):
        xy=model.vertices[tri]
        area=np.linalg.det(np.stack((xy[1]-xy[0],xy[2]-xy[0])))/2
        rho=xy[:,0]
        exact=np.empty((3,3))
        for i in range(3):
            for j in range(3):
                if i==j:exact[i,j]=area*(2*rho[i]+sum(rho))/30
                else:exact[i,j]=area*(rho[i]+rho[j]+sum(rho))/60
        analytic_error=max(analytic_error,float(np.max(abs(exact-cell['mass']))))
    Z,Bs,_,_=evaluator._maps(kappa)
    matrix=(base.unshifted+kappa*base.mass)[evaluator.active][:,evaluator.active].tocsc()
    green=evaluator._green(matrix,Z,True)
    X,Y,S=(green[k] for k in ('X','Y','S'))
    Ks=[np.linalg.solve((np.eye(len(B))-S@B).T,B.T).T for B in Bs]
    K=Ks[1]-Ks[0]

    def affine_spinor(points):
        rho,z=np.asarray(points).T
        return np.column_stack((1+.2j+.3*rho-.2j*z,.4-.1j+.1j*rho+.25*z))
    vector=affine_spinor(model.vertices[model.triangles].reshape(-1,2)).reshape(-1)
    rows=[]
    for region,info in zip(REGIONS,evaluator.region_info):
        reference,reference_weight=weighted_mass(model,lambda p:bump(p,region),24)
        restricted=reference[evaluator.active][:,evaluator.active]
        reference_response=float(kappa/math.pi*np.trace(K@(Y@(restricted@X))).real)
        angles=np.arange(16)*(2*math.pi/16)
        directions=np.column_stack((np.cos(angles),np.sin(angles)))
        center=np.asarray(region['center']);radius=region['radius']
        # An independent polar integral over the exact support disk. The
        # angular trapezoid is exact for this degree-three affine norm weight.
        def radial_integrand(r):
            points=center+r*directions
            weight=float(bump(np.array([[center[0]+r,center[1]]]),region)[0])
            density=points[:,0]*np.sum(abs(affine_spinor(points))**2,axis=1)
            return float(2*math.pi*r*weight*np.mean(density))
        exact_form=quad(radial_integrand,0,radius,epsabs=2e-14,epsrel=2e-13)[0]
        orders=[]
        for order in (4,8,12,16,24):
            mass,volume=(reference,reference_weight) if order==24 else weighted_mass(model,lambda p:bump(p,region),order)
            selected=mass[evaluator.active][:,evaluator.active]
            response=float(kappa/math.pi*np.trace(K@(Y@(selected@X))).real)
            form=float(np.vdot(vector,mass@vector).real)
            difference=mass-reference
            orders.append(dict(order=order,meridional_weight=volume,
                analytic_weight_relative_error=abs(volume/info['exact_meridional_weight']-1),
                matrix_relative_difference_from24=float(np.linalg.norm(difference.data)/np.linalg.norm(reference.data)),
                response_integrand=response,response_relative_difference_from24=abs(response/reference_response-1),
                affine_spinor_quadratic_form=form,independent_polar_quadratic_form=exact_form,
                affine_quadratic_form_relative_error=abs(form/exact_form-1)))
        rows.append(dict(name=region['name'],kappa=kappa,orders=orders))
    return dict(constant_mass_relative_error=float(constant_error),
        exact_triangle_scalar_mass_maximum_error=analytic_error,regions=rows,
        scope='Order refinement on one mesh; the order24 response is a numerical reference, not an interval error bound.')


def integrate_response(evaluator,order=8,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
    if isinstance(order,bool) or int(order)!=order or order<2:
        raise ValueError('Require integer order >=2')
    edges=np.asarray(edges,float)
    if len(edges)<2 or not np.all(np.isfinite(edges)) or edges[0]<=0 or np.any(np.diff(edges)<=0):
        raise ValueError('Panel edges must be finite, increasing and positive')
    nodes,weights=np.polynomial.legendre.leggauss(int(order))
    panels=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        kappa=lo+(nodes+1)*(hi-lo)/2
        rows=[evaluator.evaluate(k) for k in kappa]
        values=np.array([[r['response_integrand'] for r in row['regions']] for row in rows])
        responses=(hi-lo)/2*(weights@values)
        panels.append(dict(lower=float(lo),upper=float(hi),kappa=kappa.tolist(),
            response_integrands=values.tolist(),responses=responses.tolist(),
            minimum_sewing_singular_value=min(d['minimum_singular_value'] for row in rows for d in row['determinants'])))
    return dict(config=evaluator.config.copy(),order=int(order),edges=edges.tolist(),
        region_order=[r['name'] for r in REGIONS],panels=panels,
        total_responses=np.sum([p['responses'] for p in panels],axis=0).tolist(),
        scope='Smeared phase-difference response over only the displayed frequency interval.')


def independent_lapse_matrix_audit():
    """Independent non-Hermitian gate supplied and checked by formation_updates."""
    rng=np.random.default_rng(950731)
    def rnd(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
    n,m,kappa=7,3,.73
    h=rnd((n,n));h=(h+h.conj().T)/4
    A=kappa*np.eye(n)+1j*h
    G=rnd((n,m))/3;C=G.conj().T
    B0=rnd((m,m))/12;Bpi=np.diag([-1,1,-1])@B0@np.diag([1,-1,1])
    f=np.array([0.,.2,.8,1.,.35,0.,0.]);Mf=np.diag(f)
    X=np.linalg.solve(A,G);Y=np.linalg.solve(A.T,C.T).T
    wrongY=np.linalg.solve(A.conj().T,C.T).T
    adjointY=np.linalg.solve(A.conj().T,C.conj().T).conj().T
    S=C@X
    def part(B):return np.linalg.solve((np.eye(m)-S@B).T,B.T).T
    K=part(Bpi)-part(B0)
    difference=np.linalg.solve(A-G@Bpi@C,np.eye(n))-np.linalg.solve(A-G@B0@C,np.eye(n))
    small=np.trace(K@Y@Mf@X);direct=np.trace(difference@Mf)
    wrongK=np.linalg.solve(np.eye(m)-S@Bpi,Bpi)-np.linalg.solve(np.eye(m)-S@B0,B0)
    def logcontrast(eta):
        Ae=A+kappa*np.diag(1/(1+eta*f)-1)
        return np.linalg.slogdet(Ae-G@Bpi@C)[1]-np.linalg.slogdet(Ae-G@B0@C)[1]
    analytic=-kappa*small.real
    estimates=[float((logcontrast(e)-logcontrast(-e))/(2*e)) for e in (1e-3,5e-4)]
    result=dict(seed=950731,woodbury_error=float(np.linalg.norm(difference-X@K@Y)),
        trace_error=float(abs(small-direct)),correct_adjoint_error=float(np.linalg.norm(Y-adjointY)),
        mixed_adjoint_error=float(abs(np.trace(K@wrongY@Mf@X)-direct)),
        wrong_order_error=float(abs(np.trace(wrongK@Y@Mf@X)-direct)),
        analytic_derivative=float(analytic),finite_difference_estimates=estimates,
        richardson_error=float(abs(analytic-(4*estimates[1]-estimates[0])/3)))
    result['passed']=(result['woodbury_error']<1e-13 and result['trace_error']<1e-13
        and result['correct_adjoint_error']<1e-13 and result['richardson_error']<1e-10
        and result['mixed_adjoint_error']>1e-3 and result['wrong_order_error']>1e-3)
    return result


def checks():
    started=time.perf_counter();root=Path(__file__).resolve().parent
    inputs=('phase_vacuum_contrast.py','phase_vacuum_contrast_results.json',
        'curved_exterior_closure.py','curved_throat_modes.py','imaginary_exterior.py',
        'imaginary_weyl_leads.py','exterior_boundary_map.py','exterior_tail_transport.py',
        'magnetic_tail_map.py','spherical_weyl_leads.py','two_mouth_mesh.py')
    hashes={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs}
    evaluators={};integrations={}

    def get(h=.05,R=2.,N=8,L=40.):
        key=(h,R,N,L)
        if key not in evaluators:evaluators[key]=LocalLapse(h,R,N,L)
        return evaluators[key]

    def integral(h=.05,order=8,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
        key=(h,order,tuple(edges))
        if key not in integrations:
            stamp=time.perf_counter()
            record=integrate_response(get(h=h),order,edges)
            record['seconds']=time.perf_counter()-stamp;integrations[key]=record
            print(json.dumps(dict(h=h,order=order,edges=edges,responses=record['total_responses'],seconds=record['seconds'])),flush=True)
        return integrations[key]

    primary=get()
    matrix_audit=independent_lapse_matrix_audit()
    mass=mass_controls(primary)
    coarse=integral();fine=integral(order=16)
    mesh=[integral(h=h) for h in (.08,.05,.035)]
    fd_cases=[primary.evaluate(k,audit=True) for k in (.01,.4,2.)]
    frequencies=(1e-4,.01,.4,2.,8.)
    reference=[primary.evaluate(k) for k in frequencies]
    point_controls=[]
    for name,config in (('angular',dict(N=12)),('matching_radius',dict(R=4.)),('far_tail',dict(L=80.))):
        evaluator=get(**config)
        rows=[evaluator.evaluate(k) for k in frequencies]
        differences=[[trial['response_integrand']-base['response_integrand']
                      for trial,base in zip(row['regions'],ref['regions'])] for row,ref in zip(rows,reference)]
        point_controls.append(dict(name=name,config=evaluator.config,frequencies=list(frequencies),
            samples=rows,response_differences=differences,
            maximum_absolute_difference=float(np.max(abs(np.array(differences)))),
            scope='Selected-frequency sensitivity only; this is not an integrated error bound.'))
    ir=integral(edges=(1e-8,1e-7,1e-6))
    uv=integral(edges=(12.,16.))
    ir_mesh=[dict(h=h,samples=[get(h=h).evaluate(k) for k in (1e-8,1e-6,1e-4,.001,.01)]) for h in (.08,.05,.035)]
    all_rows=[row for evaluator in evaluators.values() for row in evaluator.cache.values()]
    fine_mass=[next(v for v in region['orders'] if v['order']==16) for region in mass['regions']]
    fd=[r for row in fd_cases for r in row['audit']['finite_differences']]
    quadrature=np.array(fine['total_responses'])-np.array(coarse['total_responses'])
    mesh_changes=[(np.array(b['total_responses'])-np.array(a['total_responses'])).tolist() for a,b in zip(mesh[:-1],mesh[1:])]
    tests=dict(
        independent_nonhermitian_order_and_sign=matrix_audit['passed'],
        positive_compact_lapse=min(v['minimum_lapse'] for r in fd for v in r['finite_differences'])>0,
        compact_supports_avoid_every_polygon_cut=min(min(info['polygon_clearances'].values())
                      for evaluator in evaluators.values() for info in evaluator.region_info)>0,
        exact_constant_mass=mass['constant_mass_relative_error']<1e-12,
        exact_triangle_mass=mass['exact_triangle_scalar_mass_maximum_error']<1e-14,
        independent_affine_polar_mass=max(r['affine_quadratic_form_relative_error'] for r in fine_mass)<2e-6,
        weighted_resolvent_mass_refinement=max(r['response_relative_difference_from24'] for r in fine_mass)<1e-6,
        column_equations=max(r['column_equation_relative_residual'] for r in all_rows)<1e-11,
        ordinary_transpose_equations=max(r['transpose_equation_relative_residual'] for r in all_rows)<1e-11,
        direct_DG_trace=max(max(r['audit']['direct_DG_trace_errors']) for r in fd_cases)<1e-10,
        fully_conjugated_adjoint_equivalence=max(r['audit']['transpose_vs_fully_conjugated_adjoint'] for r in fd_cases)<1e-12,
        direct_lapse_determinant_Richardson=max(r['richardson_relative_error'] for r in fd)<1e-6,
        direct_lapse_second_order=min(min(r['successive_absolute_error_ratios']) for r in fd)>3.7
                and max(max(r['successive_absolute_error_ratios']) for r in fd)<4.3,
        finite_frequency_quadrature=max(abs(quadrature))<1e-7,
        sampled_sewings_nonsingular=min(d['minimum_singular_value'] for r in all_rows for d in r['determinants'])>1e-7,
        frozen_dependencies_preserved=hashes=={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs},
    )
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,regions=primary.region_info,
        parameters=dict(q=1,mu=.5,b=.2,d=1.,rho_c=2.,phase_difference='pi minus0',
            units='hbar=c=d=1; eta and f are dimensionless, responses restore units hbar*c/d',
            source='Compact change of the lapse only; fixed flat spatial metric and magnetic A_i, A_t=0.'),
        mass_controls=mass,independent_matrix_audit=matrix_audit,
        quadrature_refinement=[coarse,fine],mesh_refinement=mesh,
        representative=fd_cases[1],finite_difference_checks=fd_cases,
        selected_frequency_reference=reference,selected_frequency_sensitivities=point_controls,
        infrared_extension=ir,ultraviolet_extension=uv,infrared_mesh_samples=ir_mesh,
        diagnostics=dict(finite_window_responses=fine['total_responses'],
            frequency_quadrature_differences=quadrature.tolist(),mesh_successive_changes=mesh_changes,
            maximum_Richardson_relative_error=max(r['richardson_relative_error'] for r in fd),
            maximum_column_residual=max(r['column_equation_relative_residual'] for r in all_rows),
            maximum_transpose_residual=max(r['transpose_equation_relative_residual'] for r in all_rows),
            minimum_sampled_sewing_singular_value=min(d['minimum_singular_value'] for r in all_rows for d in r['determinants']),
            maximum_sampled_sewing_condition_number=max(d['condition_number'] for r in all_rows for d in r['determinants']),
            frequency_geometry_evaluations=len(all_rows),
            combined_continuum_bound=None,omitted_IR_bound=None,omitted_UV_bound=None,
            integrated_N_R_L_sensitivity_bound=None),
        equations=dict(canonical='H[N]=sqrt(N) H0 sqrt(N), delta H={f,H0}/2.',
            weak='A_kappa[N]=D+kappa*M_(1/N); delta A=-kappa*M_f.',
            green='X=A_out^(-1)G_m; Y=C_m A_out^(-1); S=C_m X.',
            sewing='K_phi=B_phi(I-S B_phi)^(-1); G_pi-G_0=X(K_pi-K_0)Y.',
            response='d Delta E/d eta=(1/pi) integral kappa ReTr[(K_pi-K_0) Y M_f X] d kappa.',
            continuum='Conditionally, at N=1 this is integral f(x) Delta<T_hat00(x)> d^3x.',
            volume='Meridional bump weight=pi*rho0*a^2*(1-e*E1(1)); full3D weight=2pi times this. The normalized fixed-mu trace has NO extra2pi.'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_source_sha256=hashes,
        seconds=time.perf_counter()-started,
        scope=['Smeared phase-difference energy-density diagnostic within the retained dissipative q=1, mu=.5 sector, on a prescribed synthetic geometry.',
            'All local N-dependent determinant factors and geometries are common to the two phases. The compact support leaves every mouth, tail and throat map unchanged.',
            'The derivative is tested by independently rebuilding the finite-lapse operator and computing its phase determinants, not by perturbing saved probabilities.',
            'Bumps are nonnegative peak-one weights, not point probes. No response is divided by a volume to claim a pointwise physical stress.',
            'The continuum energy-density interpretation requires self-adjoint sewing, common renormalization and continuum convergence, which this DG calculation does not prove.',
            'Only phase differences in one retained sector are computed; absolute energy, other azimuthal sectors/chirality and radial pressure/NEC are absent.',
            'Mesh and mass/frequency quadratures are separately checked. N/R/L controls are sampled frequencies only, not bounds on their integrated effects.',
            'The IR mesh plateau is numerical and not evidence of a physical gap. Uncomputed frequency tails and a combined continuum bound remain unresolved.',
            'This is a mathematical source diagnostic, not a modeled abiotic object, prepared state, corrected Einstein solution, stable medium, or time-travel demonstration.'],
        sources=['https://arxiv.org/pdf/1106.0173','https://arxiv.org/pdf/1606.07834',
                 'https://arxiv.org/pdf/0711.2395','https://arxiv.org/html/1807.04726v3'])
    return result


if __name__=='__main__':
    result=checks()
    target=Path(__file__).with_name('local_lapse_response_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({k:result[k] for k in ('passed','tests','diagnostics','seconds','source_sha256')},indent=2))
    if not result['passed']:raise SystemExit('A local-lapse control failed; diagnostics retained.')
