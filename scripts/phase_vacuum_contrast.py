"""Phase contrast of the retained connected exterior/throat determinant.

This evaluates a finite-sector, mesh-dependent imaginary-frequency functional
on the synthetic q=1 geometry of increments 90--93.  It is not a supporting
stress tensor or a renormalized four-dimensional determinant.  For a continuum
self-adjoint sewing with a common local subtraction, -1/pi times this integral
has the symmetric-sea (-Tr|H|/2) normalization for one complex Weyl sector.
The dissipative DG approximation does not itself prove those hypotheses.

The two sewings differ only by endpoint holonomy phi=pi versus phi=0.  No
disconnected reference, extra chirality, length variation, or quantum state
preparation is inferred.  Primary spectral-determinant foundations:
https://arxiv.org/pdf/1606.07834 ; https://arxiv.org/pdf/0711.2395 .

Run: .venv\\Scripts\\python.exe -B -W error phase_vacuum_contrast.py
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.sparse.linalg import splu

from curved_exterior_closure import mouth_factors, raw_mouth_map
from exterior_boundary_map import _complex_json
from exterior_tail_transport import facet_map
from imaginary_exterior import build_problem
from imaginary_weyl_leads import tail_map
from spherical_weyl_leads import outer_characteristics


class PhaseContrast:
    """One fixed geometry; only the positive imaginary frequency is changed."""

    def __init__(self, h=.08, R=2., nmax=8, L=40., rho_c=2.):
        if (not all(math.isfinite(float(v)) for v in (h,R,L,rho_c))
                or h<=0 or not L>R>.7 or rho_c<=0
                or isinstance(nmax,bool) or int(nmax)!=nmax or nmax<1):
            raise ValueError('Require positive h/rho_c, L>R>.7 and integer nmax>=1')
        self.config=dict(h=float(h),R=float(R),nmax=int(nmax),L=float(L),rho_c=float(rho_c))
        self.problem=build_problem(.4,h,R,nmax,L)
        self.model=self.problem['model']
        self.mass=self.problem['mass']
        self.unshifted=self.model.matrix_full-.4*self.mass
        self.upper,self.lower,self.outer=self.problem['ports']
        self.Gm=np.column_stack((self.upper['G'],self.lower['G']))
        self.Cm=np.vstack((self.upper['C'],self.lower['C']))
        self.Fin,self.Fout,self.basis_error=mouth_factors(self.problem)
        self.outer_basis=lambda points:outer_characteristics(.5,
            np.arctan2(points[:,0],points[:,1]),R,nmax)
        self.cache={}

    def evaluate(self,kappa,audit=False):
        kappa=float(kappa)
        if not math.isfinite(kappa) or kappa<=0:
            raise ValueError('kappa must be finite and positive')
        key=(kappa,bool(audit))
        if key in self.cache:
            return self.cache[key]
        started=time.perf_counter()
        model=self.model
        model.matrix_full=self.unshifted+kappa*self.mass
        model.matrix=model.matrix_full[model.active][:,model.active].tocsc()
        model.lu=splu(model.matrix)
        self.problem['kappa']=kappa
        rawZ,tail=tail_map(kappa,self.config['R'],self.config['L'],nmax=self.config['nmax'])
        Z,_,geometry=facet_map(self.outer,self.outer_basis,rawZ)
        Go,Co=self.outer['G'],self.outer['C']
        # Eliminate the outer modes only. Mouth inputs remain independent.
        bare=model.solve(np.column_stack((self.Gm,Go)))
        count=self.Gm.shape[1]
        Xm,Xo=bare[:,:count],bare[:,count:]
        Qmm,Qmo=self.Cm@Xm,self.Cm@Xo
        Qom,Qoo=Co@Xm,Co@Xo
        Mouter=np.eye(len(Z))-Qoo@Z
        coefficients=np.linalg.solve(Mouter,Qom)
        solutions=Xm+Xo@Z@coefficients
        S=Qmm+Qmo@Z@coefficients
        sparse_residual=(model.matrix@solutions[model.active]-self.Gm[model.active]
                         -(Go@Z@(Co@solutions))[model.active])
        residual=float(np.linalg.norm(sparse_residual)/np.linalg.norm(self.Gm[model.active]))
        Sraw=np.linalg.solve(self.Fout,S@self.Fin)
        records=[]
        maps=[]
        for phase in (0.,math.pi):
            rawB,throat,_=raw_mouth_map(self.problem,phase,self.config['rho_c'],'curved')
            B=np.linalg.solve(self.Fout.T,(self.Fin@rawB).T).T
            matrix=np.eye(count)-S@B
            sign,logabs=np.linalg.slogdet(matrix)
            if sign==0 or not math.isfinite(logabs):
                raise ArithmeticError('A sampled sewing determinant is singular')
            raw_sign,raw_log=np.linalg.slogdet(np.eye(count)-Sraw@rawB)
            sewing_singular=np.linalg.svd(matrix,compute_uv=False)
            records.append(dict(phase=phase,log_absolute_determinant=float(logabs),
                determinant_phase=_complex_json(sign),
                minimum_sewing_singular_value=float(sewing_singular[-1]),
                sewing_condition_number=float(sewing_singular[0]/sewing_singular[-1]),
                canonical_logdet_difference=float(abs(raw_log-logabs)),
                canonical_phase_difference=float(abs(raw_sign-sign))))
            maps.append((rawB,B,matrix))
        row=dict(kappa=kappa,delta_log_absolute_determinant=records[1]['log_absolute_determinant']
                 -records[0]['log_absolute_determinant'],phases=records,
            sparse_multi_input_relative_residual=residual,
            exterior_map_coordinate_error=float(np.linalg.norm(self.Fout@Sraw-S@self.Fin)),
            exterior_map_maximum_singular_value=float(np.linalg.svd(S,compute_uv=False)[0]),
            outer_coordinate_diagnostics=geometry,tail_accepted_steps=tail['accepted_steps'],
            throat_optical_length=throat['optical_length'],seconds=time.perf_counter()-started)
        if audit:
            outer_sign,outer_log=np.linalg.slogdet(Mouter)
            direct_errors=[]
            for record,(_,B,matrix) in zip(records,maps):
                full=np.block([[np.eye(count)-Qmm@B,-Qmo@Z],
                               [-Qom@B,Mouter]])
                full_sign,full_log=np.linalg.slogdet(full)
                determinant_phase=complex(record['determinant_phase']['real'],record['determinant_phase']['imag'])
                direct_errors.append(dict(logabs_error=float(abs(full_log-outer_log-record['log_absolute_determinant'])),
                    phase_error=float(abs(full_sign/outer_sign-determinant_phase))))
            left,singular,right=np.linalg.svd(Z)
            independent=model.solve(self.Gm,(Go@left,singular,right@Co))
            phi=.371
            rawB,_,_=raw_mouth_map(self.problem,phi,self.config['rho_c'],'curved')
            n=self.config['nmax']
            incoming_generator=np.r_[np.ones(n+1),np.zeros(n)]
            outgoing_generator=np.r_[np.ones(n),np.zeros(n+1)]
            derivative=1j*(incoming_generator[:,None]-outgoing_generator[None,:])*rawB
            B=np.linalg.solve(self.Fout.T,(self.Fin@rawB).T).T
            dB=np.linalg.solve(self.Fout.T,(self.Fin@derivative).T).T
            analytic=-np.trace(np.linalg.solve(np.eye(count)-S@B,S@dB)).real
            estimates=[]
            for step in (1e-4,5e-5):
                values=[]
                for angle in (phi-step,phi+step):
                    rb,_,_=raw_mouth_map(self.problem,angle,self.config['rho_c'],'curved')
                    bb=np.linalg.solve(self.Fout.T,(self.Fin@rb).T).T
                    values.append(np.linalg.slogdet(np.eye(count)-S@bb)[1])
                estimates.append(float((values[1]-values[0])/(2*step)))
            row['audit']=dict(full_boundary_schur=direct_errors,
                independent_woodbury_relative_difference=float(np.linalg.norm(independent-solutions)/np.linalg.norm(solutions)),
                derivative_phase=phi,analytic_phase_derivative=float(analytic),
                centered_phase_derivatives=estimates,
                phase_derivative_error=abs(float(analytic)-estimates[-1]),
                matrices=dict(S_exterior_facet=_complex_json(S),S_exterior_canonical=_complex_json(Sraw),
                    F_in=_complex_json(self.Fin),F_out=_complex_json(self.Fout),
                    B_zero_raw=_complex_json(maps[0][0]),B_pi_raw=_complex_json(maps[1][0]),
                    B_zero_facet=_complex_json(maps[0][1]),B_pi_facet=_complex_json(maps[1][1]),
                    Z_outer_facet=_complex_json(Z)),
                incoming_order='upper zero; upper paired n1..N; lower paired n1..N',
                outgoing_order='upper paired n1..N; lower zero; lower paired n1..N')
        self.cache[key]=row
        return row


def integrate(evaluator,order=12,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
    """Gauss integration over explicit positive-frequency panels; no tail fit."""
    if isinstance(order,bool) or int(order)!=order or order<2:
        raise ValueError('order must be an integer >=2')
    edges=np.asarray(edges,float)
    if not np.all(np.isfinite(edges)) or edges[0]<=0 or np.any(np.diff(edges)<=0):
        raise ValueError('Positive, finite, increasing panel edges required')
    nodes,weights=np.polynomial.legendre.leggauss(int(order))
    panels=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        kappa=lo+(nodes+1)*(hi-lo)/2
        rows=[evaluator.evaluate(k) for k in kappa]
        values=np.array([r['delta_log_absolute_determinant'] for r in rows])
        integral=float((hi-lo)/2*(weights@values))
        panels.append(dict(lower=float(lo),upper=float(hi),integral=integral,
            phase_contrast=-integral/math.pi,kappa=kappa.tolist(),integrand=values.tolist(),
            minimum_sewing_singular_value=min(p['minimum_sewing_singular_value'] for r in rows for p in r['phases'])))
    total=sum(panel['integral'] for panel in panels)
    return dict(config=evaluator.config.copy(),order=int(order),edges=edges.tolist(),panels=panels,
                integral=total,phase_contrast=-total/math.pi,
                frequency_scope='Only the displayed positive-frequency interval; omitted tails are separate diagnostics.')


def normalization_check():
    """Circle NS minus periodic: one complex chiral symmetric sea, not Dirac."""
    from scattering_vacuum import LoopReservoir
    with mp.workdps(50):
        length=mp.mpf('1.3')
        # D_pi/D_0=(1+exp(-k L))/(1-exp(-k L)); its endpoint log is integrable.
        direct=-mp.quad(lambda y:mp.log1p(mp.exp(-y))-mp.log(-mp.expm1(-y)),
                        [0,1,mp.inf])/(mp.pi*length)
        expected=-mp.pi/(4*length)
        ns_zeta=-2*mp.pi*mp.zeta(-1,mp.mpf('.5'))/length
        periodic_zeta=-2*mp.pi*mp.zeta(-1)/length
        previous=(LoopReservoir(1.,float(length),0.).vacuum_energy('complex_chiral',False)
                  -LoopReservoir(1.,float(length),math.pi).vacuum_energy('complex_chiral',False))
        return dict(length=float(length),circle_ns_minus_periodic=float(direct),expected=float(expected),
            analytic_error=float(abs(direct-expected)),previous_loop_reservoir_difference=float(previous),
            previous_model_difference=float(abs(previous-direct)),
            independent_zeta_contrast=float(ns_zeta-periodic_zeta),
            independent_zeta_error=float(abs(ns_zeta-periodic_zeta-direct)),
            ns_chiral_energy=float(-mp.pi/(12*length)),periodic_chiral_energy=float(mp.pi/(6*length)),
            interpretation='One complex chiral symmetric-sea normalization; no second chirality is silently added.')


def checks():
    started=time.perf_counter()
    root=Path(__file__).resolve().parent
    inputs=('curved_exterior_closure.py','curved_throat_modes.py','imaginary_exterior.py',
            'imaginary_weyl_leads.py','exterior_boundary_map.py','exterior_tail_transport.py',
            'magnetic_tail_map.py','spherical_weyl_leads.py','two_mouth_mesh.py','scattering_vacuum.py')
    input_hashes={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs}
    evaluators={}
    integrations={}

    def get(h=.05,R=2.,N=8,L=40.):
        key=(h,R,N,L)
        if key not in evaluators:
            evaluators[key]=PhaseContrast(h,R,N,L)
        return evaluators[key]

    def integral(h=.05,R=2.,N=8,L=40.,order=12,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
        key=(h,R,N,L,order,tuple(edges))
        if key not in integrations:
            stamp=time.perf_counter()
            result=integrate(get(h,R,N,L),order,edges)
            result['seconds']=time.perf_counter()-stamp
            integrations[key]=result
            print(json.dumps(dict(h=h,R=R,N=N,L=L,order=order,edges=edges,
                phase_contrast=result['phase_contrast'],seconds=result['seconds'])),flush=True)
        return integrations[key]

    # Same panels and parameters for quadrature and each geometry comparison.
    coarse=integral()
    fine=integral(order=24)
    mesh=[integral(h=h) for h in (.08,.05,.035)]
    angular=[coarse,integral(N=12)]
    radius=[coarse,integral(R=4.)]
    radial_tail=[coarse,integral(L=80.)]
    representative=get().evaluate(.4,audit=True)
    infrared=integral(order=24,edges=(1e-8,1e-7,1e-6))
    ultraviolet=integral(order=12,edges=(12.,16.))
    ir_samples=[]
    for h in (.08,.05,.035):
        evaluator=get(h=h)
        ir_samples.append(dict(h=h,samples=[evaluator.evaluate(k) for k in (1e-8,1e-6,1e-4,.001,.01)]))
    normalization=normalization_check()
    all_rows=[row for evaluator in evaluators.values() for row in evaluator.cache.values()]
    audit=representative['audit']
    quadrature_error=abs(fine['phase_contrast']-coarse['phase_contrast'])
    mesh_changes=[abs(b['phase_contrast']-a['phase_contrast']) for a,b in zip(mesh[:-1],mesh[1:])]
    tests=dict(
        all_samples_finite=all(math.isfinite(row['delta_log_absolute_determinant']) for row in all_rows),
        independent_circle_integral=normalization['analytic_error']<1e-40,
        independent_circle_zeta=normalization['independent_zeta_error']<1e-40,
        previous_loop_normalization=normalization['previous_model_difference']<1e-14,
        multi_input_sparse_equations=max(row['sparse_multi_input_relative_residual'] for row in all_rows)<1e-11,
        mouth_coordinate_invariance=max(p['canonical_logdet_difference'] for row in all_rows for p in row['phases'])<1e-10,
        full_boundary_schur=max(r['logabs_error'] for r in audit['full_boundary_schur'])<1e-11,
        full_boundary_phase=max(r['phase_error'] for r in audit['full_boundary_schur'])<1e-11,
        independent_woodbury=audit['independent_woodbury_relative_difference']<1e-11,
        analytic_phase_derivative=audit['phase_derivative_error']<2e-8,
        frequency_quadrature_refinement=quadrature_error<1e-6,
        sampled_positive_kappa_sewings_nonsingular=min(p['minimum_sewing_singular_value'] for row in all_rows for p in row['phases'])>1e-7,
        frozen_inputs_preserved=input_hashes=={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs},
    )
    result=dict(passed=all(tests.values()),tests=tests,
        parameters=dict(q=1,mu=.5,b=.2,d=1.,rho_c=2.,energy='i*kappa',phase_difference='pi minus 0',
            units='hbar=c=d=1; this functional has units hbar*c/d after restoring separation'),
        normalization=normalization,representative=representative,
        quadrature_refinement=[coarse,fine],mesh_refinement=mesh,angular_refinement=angular,
        matching_radius_refinement=radius,far_tail_refinement=radial_tail,
        infrared_extension=infrared,ultraviolet_extension=ultraviolet,infrared_mesh_samples=ir_samples,
        diagnostics=dict(finite_window_phase_contrast=fine['phase_contrast'],
            frequency_quadrature_difference=quadrature_error,mesh_successive_changes=mesh_changes,
            angular_change=abs(angular[1]['phase_contrast']-angular[0]['phase_contrast']),
            matching_radius_change=abs(radius[1]['phase_contrast']-radius[0]['phase_contrast']),
            far_tail_change=abs(radial_tail[1]['phase_contrast']-radial_tail[0]['phase_contrast']),
            minimum_sampled_sewing_singular_value=min(p['minimum_sewing_singular_value'] for row in all_rows for p in row['phases']),
            maximum_sampled_sewing_condition_number=max(p['sewing_condition_number'] for row in all_rows for p in row['phases']),
            maximum_sparse_equation_residual=max(row['sparse_multi_input_relative_residual'] for row in all_rows),
            evaluated_frequency_geometry_pairs=len(all_rows),
            combined_continuum_error_bound=None,omitted_zero_to_1e_minus8_bound=None,omitted_above_16_bound=None),
        equations=dict(exterior='S_ext=C_m*(A_kappa-G_o*Z_facet*C_o)^(-1)*G_m; mouth feedback excluded.',
            canonical='S_raw=F_out^(-1)*S_ext*F_in; B_facet=F_in*B_raw*F_out^(-1).',
            determinant='D_phi=det(I-S_ext*B_phi).',
            integrand='delta_log=log|D_pi|-log|D_0|, at exactly the same local geometry and reference.',
            functional='Delta E_sector=-1/pi integral(delta_log d kappa), symmetric-sea one-complex-Weyl normalization.',
            derivative='d_phi log|D_phi|=-Re Tr[(I-S_ext B_phi)^(-1) S_ext d_phi B_phi].',
            cross_phases='lower->upper exp(+i*phi); upper->lower exp(-i*phi); same-end reflections unchanged.'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_source_sha256=input_hashes,
        seconds=time.perf_counter()-started,
        scope=['Retained-sector phase-contrast determinant functional of the actual connected finite DG exterior and synthetic finite curved throat.',
            'The finite window [1e-6,12] is reported separately from computed IR/UV extensions; no uncomputed tail is set to zero or certified by a fitted decay.',
            'DG upwinding is dissipative. Continuum self-adjoint sewing and full spectral convergence are not established by this finite calculation.',
            'For a self-adjoint continuum operator with common local subtraction, the stated prefactor follows from -Tr|H|/2. Lower-axis adjointness does not presume a second physical chirality or +/-E symmetry at fixed mu.',
            'Only one complex Weyl sector q=1, mu=.5 and finitely many paired modes is retained. Other azimuthal sectors and Dirac chirality are not silently included.',
            'The same-end throat factors and all local geometry are common to both phases; this is a holonomy comparison, not an invented disconnected subtraction.',
            'Endpoint spin/gauge matching, NS versus periodic phase and the synthetic metric are assumptions. Neither the junction stress nor its semiclassical Einstein source is supplied.',
            'The finite mesh regulates the near-zero periodic behavior; the separate circle control retains its integrable log(kappa) singularity. No uniform continuum IR bound is claimed.',
            'A phase derivative of this modulus functional is a consistency diagnostic, not an eta invariant, physical charge preparation, local stress, radial NEC, or supporting-energy density.',
            'No length derivative is taken. This cannot repair the previously established phase-independent geometry/source mismatch, prove stability, or establish time travel.'],
        sources=['https://arxiv.org/pdf/1606.07834','https://arxiv.org/pdf/0711.2395',
                 'https://arxiv.org/html/1807.04726v3'])
    return result


if __name__=='__main__':
    result=checks()
    # Retain diagnostic output even if a substantive control fails.
    target=project_path('phase_vacuum_contrast_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({key:result[key] for key in ('passed','tests','diagnostics','seconds','source_sha256')},indent=2))
    if not result['passed']:
        raise SystemExit('A phase-contrast control failed; diagnostics retained.')
