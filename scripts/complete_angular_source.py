"""Complete finite q=1 mouth multiplets, with a separate outer cutoff.

This finite-cut calculation reuses the frozen DG bulk, opposite-Weyl
assembly, radial throat variations, magnetic tail, and AdS reference kernels.
Only the angular/port adapter is new.  The older adapter restricts mu to
three values and ties the two cutoffs together, so its short port assembly
is composed here without monkeypatching or changing historical sources.

The observables are compact connected-minus-AdS-reference source kernels
in a prescribed product throat.  Completing n<=N does not complete the
angular sum, frequency integral, continuum limit, or Einstein--matter join.

Run: .venv\\Scripts\\python.exe -B -W error complete_angular_source.py
"""
from project_paths import project_path, source_matches
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import block_diag as sparse_block_diag
from scipy.special import jv, jvp
from numpy.polynomial.legendre import leggauss

from ads_reference_source import cap_reflection, reference_map, reference_relative
from angular_reference_sum import comparisons, record_controls, sum_records
from complete_angular_modes import (basis_checks, mouth_channels,
    mouth_characteristics, reflection_parities, zero_mode)
from dirac_sewing_audit import opposite_assembly, port_factors, exterior_response
from exterior_boundary_map import ExteriorDG, ID, SX, SY, SZ, _complex_json
from exterior_tail_transport import facet_map
from imaginary_weyl_leads import tail_map
from spherical_weyl_leads import outer_characteristics
from throat_null_response import raw_variation, phase_map, shift_derivatives, probe_geometry
from two_mouth_mesh import mesh


FIELDS=('pressure','shift','null_plus','null_minus')
PHASES=('0','pi')
INPUTS=('complete_angular_modes.py','ads_reference_source.py',
    'ads_reference_source_results.json','azimuthal_modes.py','azimuthal_source.py',
    'dirac_sewing_audit.py','exterior_boundary_map.py','exterior_tail_transport.py',
    'imaginary_weyl_leads.py','magnetic_tail_map.py','spherical_weyl_leads.py',
    'throat_null_response.py','curved_throat_modes.py','two_mouth_mesh.py')
EDGES=(1e-6,.02,.2,1.,4.,12.,24.,48.,96.)
RECIPES=('coarse','fine','finer_mesh')


def _decode(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def _hashes():
    return {name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in INPUTS}


def configuration(mu,h=.05,R=2.,nmax=2,outer_nmax=9,L=40.,rho_c=2.):
    """Validate before allocating a mesh; do not silently change either cutoff."""
    if not math.isfinite(float(mu)) or abs(2*mu-round(2*mu))>1e-12 or round(2*mu)%2!=1:
        raise ValueError('mu must be a finite half-odd-integer')
    m=int(round(mu-.5));mu=m+.5
    for value in (nmax,outer_nmax):
        if isinstance(value,bool) or not isinstance(value,(int,np.integer)):
            raise ValueError('Both angular cutoffs must be integers')
    nmax,outer_nmax=int(nmax),int(outer_nmax)
    ns=mouth_channels(mu,nmax)
    outer_first=int(round(abs(mu)+.5))
    if not outer_first<=outer_nmax<64:
        raise ValueError('Nonempty outer channels and J<64 are required by this default tail quadrature')
    if (not all(math.isfinite(float(v)) for v in (h,R,L,rho_c)) or h<=0
            or not L>R>.7 or rho_c<=0):
        raise ValueError('Require h/rho_c>0 and L>R>.7')
    return dict(mu=float(mu),m=m,h=float(h),R=float(R),nmax=nmax,
        outer_nmax=outer_nmax,L=float(L),rho_c=float(rho_c),
        mouth_channels=ns.tolist(),outer_channels=list(range(outer_first,outer_nmax+1)),
        has_zero_mode=m==0,reflection_cross_signs=reflection_parities(mu,nmax).tolist())


def selectors(mu,nmax):
    """B maps core-entering columns to core-exiting rows (exterior in)."""
    ns=mouth_channels(mu,nmax)
    if mu==.5:
        indices=np.arange(2*nmax+1)
        return indices,indices
    return np.r_[ns,nmax+ns],np.r_[ns-1,nmax+ns]


def selected_map(value,mu,angle,nmax):
    """Apply physical endpoint phases/parities before selecting absolute n."""
    ns=mouth_channels(mu,nmax);parity=reflection_parities(mu,nmax)
    result=phase_map(value,angle,nmax).copy()
    result[ns,nmax+ns]*=parity
    result[nmax+ns,ns-1]*=parity
    rows,cols=selectors(mu,nmax)
    return result[np.ix_(rows,cols)]


def complete_reference(kappa,mu,nmax,rho_c=2.):
    """AdS cap map: its input/output selectors are the reverse of B's."""
    if mu==.5:
        return reference_map(kappa,nmax,mu,rho_c)
    ell=.2*math.hypot(1.,rho_c);caps=[];values=[]
    for n in mouth_channels(mu,nmax):
        value,diag=cap_reflection(1j*kappa,math.sqrt(n*(n+1)),ell,rho_c)
        values.append(value);caps.append(dict(n=int(n),reflection=_complex_json(value),**diag))
    right=np.diag(np.r_[values,values])
    return block_diag(right,right.conj().T),dict(ell=ell,
        half_cap_optical_length=ell*(math.pi/2-math.atan(rho_c)),
        global_ads_optical_length=math.pi*ell,has_cross_chirality_zero=False,caps=caps)


def port_labels(mu,nmax):
    ns=mouth_channels(mu,nmax).tolist()
    incoming=[('upper',n) for n in ([0]+ns if mu==.5 else ns)]+[('lower',n) for n in ns]
    outgoing=[('upper',n) for n in ns]+[('lower',n) for n in ([0]+ns if mu==.5 else ns)]
    return incoming,outgoing


def reproject_environment(connected,mu,high_nmax,low_nmax,right_out,left_out):
    """Exactly select a smaller mouth space on the SAME bulk/outer system.

    S_raw=W_out^-1 T_raw. Thus S_low=W_low^-1 Pout W_high S_high Pin^T,
    W_low=Pout W_high Pout^T. A naive S_raw submatrix omits this Gram
    correction. Left-Weyl input/output selectors swap the right selectors.
    The smaller nonlinear core/reference sewing must then be solved afresh.
    """
    if low_nmax>high_nmax:raise ValueError('Cannot expand an exterior response')
    hi_in,hi_out=port_labels(mu,high_nmax);lo_in,lo_out=port_labels(mu,low_nmax)
    pin=np.array([hi_in.index(v) for v in lo_in]);pout=np.array([hi_out.index(v) for v in lo_out])
    count=len(hi_in);ins=np.r_[pin,count+pout];outs=np.r_[pout,count+pin]
    Fout=block_diag(right_out,left_out);weight=Fout.conj().T@Fout
    low_weight=weight[np.ix_(outs,outs)]
    response=np.linalg.solve(low_weight,(weight@connected)[np.ix_(outs,ins)])
    return response,dict(incoming_indices=ins.tolist(),outgoing_indices=outs.tolist(),
        output_gram_condition=float(np.linalg.cond(low_weight)),
        naive_submatrix_difference=float(np.max(abs(response-connected[np.ix_(outs,ins)]))))


def source_response(connected,kappa,mu,nmax,rho_c,probes):
    """Unchanged radial kernels, independently sewn at the selected cutoff."""
    reference,capdiag=complete_reference(kappa,mu,nmax,rho_c)
    raw,pressure,_=raw_variation(kappa,nmax,rho_c)
    shift=shift_derivatives(raw,kappa,nmax,probes);phases=[]
    for angle in (0.,math.pi):
        b=selected_map(raw,mu,angle,nmax)
        dp=[selected_map(v,mu,angle,nmax) for v in pressure]
        ds=[selected_map(v,mu,angle,nmax) for v in shift]
        B=block_diag(b,b.conj().T)
        P=[block_diag(v,v.conj().T) for v in dp]
        Q=[block_diag(v,-v.conj().T) for v in ds]
        phases.append(dict(phase=angle,**reference_relative(connected,reference,B,P,Q)))
    return dict(kappa=kappa,phases=phases,cap_diagnostics=capdiag,
        reference_phase_pressure_error=float(np.max(abs(np.asarray(phases[0]['reference_pressure_trace'])-phases[1]['reference_pressure_trace']))),
        reference_phase_shift_error=float(np.max(abs(np.asarray(phases[0]['reference_shift_trace'])-phases[1]['reference_shift_trace']))))


class CompleteSector:
    """One physical mu sector, with independent mouth N and outer J."""

    def __init__(self,mu,h=.05,R=2.,nmax=2,outer_nmax=9,L=40.,rho_c=2.):
        self.config=configuration(mu,h,R,nmax,outer_nmax,L,rho_c)
        mu=self.config['mu'];self.ns=np.array(self.config['mouth_channels'])
        self.model=ExteriorDG(mesh(b=.2,d=1.,R=R,h=h),q=1,mu=mu,energy=0.)
        self.active=self.model.active;self.unshifted=self.model.matrix_full
        self.mass=sparse_block_diag([np.kron(c['mass'],ID) for c in self.model.cells],format='csc')
        self.left_zero,self.proxy=opposite_assembly(self.model)
        self.right_ports=[];self.left_ports=[];self.right_factors={};self.left_factors={}
        self.port_checks=[];self.bases={}
        for kind in ('upper','lower','outer'):
            def basis(points,kind=kind):
                if kind=='outer':
                    return outer_characteristics(mu,np.arctan2(points[:,0],points[:,1]),R,outer_nmax)
                theta=np.arctan2(points[:,0],points[:,1]-(.5 if kind=='upper' else -.5))
                plus,minus=mouth_characteristics(mu,theta,.2,kind,nmax)
                zero=zero_mode(mu,theta,.2,kind)
                if (zero is not None)!=(mu==.5):
                    raise AssertionError('Incorrect zero-channel count')
                if zero is not None:
                    if kind=='upper':minus=np.concatenate((zero[:,:,None],minus),axis=2)
                    else:plus=np.concatenate((zero[:,:,None],plus),axis=2)
                return plus,minus
            swapped=lambda points,base=basis:base(points)[::-1]
            right=ExteriorDG.modal_port(self.model,kind,basis)
            left=ExteriorDG.modal_port(self.proxy,kind,swapped)
            rf,re=port_factors(right,basis);lf,le=port_factors(left,swapped)
            self.right_ports.append(right);self.left_ports.append(left)
            self.right_factors[kind]=rf;self.left_factors[kind]=lf;self.bases[kind]=basis
            self.port_checks.append(dict(kind=kind,incoming_count=right['G'].shape[1],
                outgoing_count=right['C'].shape[0],basis_reproduction_error=max(re,le),
                incoming_adjoint_error=float(np.max(abs(left['G']-right['C'].conj().T))),
                outgoing_adjoint_error=float(np.max(abs(left['C']-right['G'].conj().T))),
                factor_swap_error=max(float(np.max(abs(lf[j]-rf[1-j]))) for j in (0,1)),
                gram_condition_numbers=[float(np.linalg.cond(f.conj().T@f)) for f in rf]))
        ur,lr,self.outer=self.right_ports;ul,ll,self.left_outer=self.left_ports
        self.Gm=np.column_stack((ur['G'],lr['G']));self.Cm=np.vstack((ur['C'],lr['C']))
        self.GL=np.column_stack((ul['G'],ll['G']));self.CL=np.vstack((ul['C'],ll['C']))
        self.Fin=block_diag(self.right_factors['upper'][1],self.right_factors['lower'][1])
        self.Fout=block_diag(self.right_factors['upper'][0],self.right_factors['lower'][0])
        self.Lin=block_diag(self.left_factors['upper'][1],self.left_factors['lower'][1])
        self.Lout=block_diag(self.left_factors['upper'][0],self.left_factors['lower'][0])
        difference=(self.left_zero-self.unshifted.conj().T).tocoo()
        axis=[]
        for name,c in (('outside',mu),('between',mu-1)):
            powers=[c-.5,c+.5] if c>0 else [-c+.5,-c-.5]
            axis.append(dict(segment=name,mu_minus_Aaxis=c,regular_powers=powers))
        self.assembly=dict(triangles=len(self.model.triangles),active_dof=len(self.active),
            opposite_adjoint_relative_error=float(np.max(abs(difference.data))/np.max(abs(self.unshifted.data))) if difference.nnz else 0.,
            minimum_cell_mass_eigenvalue=min(float(np.linalg.eigvalsh(c['mass'])[0]) for c in self.model.cells),
            axis_domain_equal=bool(np.array_equal(self.active,self.proxy.active)),
            axis_regular_powers=axis,axis_removed_dof=self.model.ndof-len(self.active),
            mouth_dimension_per_chirality=2*len(self.ns)+int(mu==.5),
            outer_dimension=len(self.config['outer_channels']),
            axis_scope='Both components vanish where |mu-Aaxis|>1/2. Higher regular powers are resolved only through mesh refinement.')
        self.probes=probe_geometry(.2*math.hypot(1.,rho_c),rho_c)

    def evaluate(self,kappa,save_matrices=False,lower_nmax=None):
        if not math.isfinite(float(kappa)) or kappa<=0:
            raise ValueError('Positive finite imaginary frequency required')
        kappa=float(kappa);started=time.perf_counter();cfg=self.config;N=cfg['nmax']
        Zraw,tail=tail_map(kappa,cfg['R'],cfg['L'],mu=cfg['mu'],nmax=cfg['outer_nmax'])
        ZR,FR,coordR=facet_map(self.outer,self.bases['outer'],Zraw)
        swapped=lambda points:self.bases['outer'](points)[::-1]
        ZL,FL,coordL=facet_map(self.left_outer,swapped,Zraw.conj().T)
        SR,resR=exterior_response(self.unshifted+kappa*self.mass,self.active,
            self.Gm,self.Cm,self.outer['G'],self.outer['C'],ZR)
        SL,resL=exterior_response(self.left_zero+kappa*self.mass,self.active,
            self.GL,self.CL,self.left_outer['G'],self.left_outer['C'],ZL)
        Rraw=np.linalg.solve(self.Fout,SR@self.Fin)
        Lraw=np.linalg.solve(self.Lout,SL@self.Lin)
        connected=block_diag(Rraw,Lraw)
        reference,capdiag=complete_reference(kappa,cfg['mu'],N,cfg['rho_c'])
        raw,pressure,_=raw_variation(kappa,N,cfg['rho_c'])
        shift=shift_derivatives(raw,kappa,N,self.probes)
        branches=[];matrices=[];coordinate_error=0.
        combined_in=block_diag(self.Fin,self.Lin);combined_out=block_diag(self.Fout,self.Lout)
        def outgoing_map(value):
            return np.linalg.solve(combined_in.T,(combined_out@value).T).T
        def incoming_map(value):
            return np.linalg.solve(combined_out.T,(combined_in@value).T).T
        for angle in (0.,math.pi):
            b=selected_map(raw,cfg['mu'],angle,N)
            dp=[selected_map(v,cfg['mu'],angle,N) for v in pressure]
            ds=[selected_map(v,cfg['mu'],angle,N) for v in shift]
            B=block_diag(b,b.conj().T)
            P=[block_diag(v,v.conj().T) for v in dp]
            Q=[block_diag(v,-v.conj().T) for v in ds]
            result=reference_relative(connected,reference,B,P,Q)
            if save_matrices:
                facet=reference_relative(block_diag(SR,SL),outgoing_map(reference),
                    incoming_map(B),[incoming_map(v) for v in P],[incoming_map(v) for v in Q])
                coordinate_error=max(coordinate_error,*[
                    float(np.max(abs(np.asarray(result[field+'_integrands'])-facet[field+'_integrands'])))
                    for field in FIELDS],abs(result['connected_logabs']-facet['connected_logabs']),
                    abs(result['reference_logabs']-facet['reference_logabs']))
                matrices.append(dict(B=_complex_json(B),pressure=_complex_json(P),shift=_complex_json(Q)))
            branches.append(dict(phase=angle,**result))
        output=dict(kappa=kappa,phases=branches,cap_diagnostics=capdiag,
            exterior_equation_relative_residuals=[resR,resL],tail_diagnostics=tail,
            outer_coordinate_diagnostics=[coordR,coordL],
            raw_facet_observable_error=coordinate_error if save_matrices else None,
            reference_phase_pressure_error=float(np.max(abs(np.array(branches[0]['reference_pressure_trace'])-branches[1]['reference_pressure_trace']))),
            reference_phase_shift_error=float(np.max(abs(np.array(branches[0]['reference_shift_trace'])-branches[1]['reference_shift_trace']))),
            seconds=time.perf_counter()-started)
        if save_matrices:
            output['matrices']=dict(S_connected=_complex_json(connected),S_reference=_complex_json(reference),
                phases=matrices,S_right_facet=_complex_json(SR),S_left_facet=_complex_json(SL),
                F_mouth_right_in=_complex_json(self.Fin),F_mouth_right_out=_complex_json(self.Fout),
                F_mouth_left_in=_complex_json(self.Lin),F_mouth_left_out=_complex_json(self.Lout),
                Z_raw_right=_complex_json(Zraw),Z_facet_right=_complex_json(ZR),Z_facet_left=_complex_json(ZL),
                F_outer_right=_complex_json(FR),F_outer_left=_complex_json(FL))
        if lower_nmax is not None and abs(cfg['m'])<=lower_nmax:
            low,info=reproject_environment(connected,cfg['mu'],N,lower_nmax,self.Fout,self.Lout)
            output['lower_cut_response']=source_response(low,kappa,cfg['mu'],lower_nmax,cfg['rho_c'],self.probes)
            output['lower_cut_response']['reprojection']=info
            output['lower_cut_response']['exterior_equation_relative_residuals']=[resR,resL]
        return output


def algebra_checks():
    """Count/selector checks plus independent free regular-axis spinors."""
    counts=[];selector_error=0.
    for N in (2,4,8):
        labels=[mouth_channels(m+.5,N) for m in range(-N,N+1)]
        counts.append(dict(N=N,paired_labels=sum(map(len,labels)),expected_paired=N*(N+2),
            single_weyl_mouth_dimension=sum(2*len(ns)+int(m==0) for m,ns in zip(range(-N,N+1),labels)),
            expected_dimension=2*N*(N+2)+1))
        full,_=reference_map(.4,N,.5)
        one=full[:2*N+1,:2*N+1]
        for m in range(-N,N+1):
            if m==0:continue
            rows,cols=selectors(m+.5,N)
            selected=one[np.ix_(cols,rows)]
            actual,_=complete_reference(.4,m+.5,N)
            selector_error=max(selector_error,float(np.max(abs(selected-actual[:len(rows),:len(rows)]))))
    E,kperp=1.3,.8;kz=math.sqrt(E*E-kperp*kperp);points=np.array([[.07,-.6],[.2,.2],[.9,.1]])
    rho,z=points.T;phase=np.exp(1j*kz*z);ratio=kperp/(E+kz);free=[]
    for m in range(-4,5):
        mu=m+.5
        u=phase[:,None]*np.column_stack((jv(m,kperp*rho),1j*ratio*jv(m+1,kperp*rho)))
        dr=kperp*phase[:,None]*np.column_stack((jvp(m,kperp*rho),1j*ratio*jvp(m+1,kperp*rho)))
        residual=(dr+u/(2*rho[:,None]))@SX.T+1j*kz*u@SZ.T+1j*mu/rho[:,None]*(u@SY.T)-1j*E*u
        free.append(dict(m=m,mu=mu,pde_error=float(np.max(abs(residual)))))
    rejected=[]
    for name,kwargs in [('non_half_mu',dict(mu=.3)),('empty_mouth',dict(mu=3.5,nmax=2)),
        ('empty_outer',dict(mu=4.5,nmax=4,outer_nmax=4)),('insufficient_tail_quadrature',dict(mu=.5,outer_nmax=64))]:
        try:configuration(**kwargs)
        except ValueError:rejected.append(name)
    return dict(counts=counts,reference_opposite_selector_error=selector_error,
        free_axis_spinors=free,invalid_configurations_rejected=rejected,
        free_formula='exp(i*kz*z)*(J_m(kperp*rho),i*kperp/(E+kz)*J_(m+1)(kperp*rho)), m=mu-1/2')


def legacy_comparison(row,saved):
    errors={}
    for name in ('S_connected','S_reference'):
        errors[name]=float(np.max(abs(_decode(row['matrices'][name])-_decode(saved['matrices'][name]))))
    for field in ('B','pressure','shift'):
        errors[field]=max(float(np.max(abs(_decode(a[field])-_decode(b[field]))))
            for a,b in zip(row['matrices']['phases'],saved['matrices']['phases']))
    errors['source_kernels']=max(float(np.max(abs(np.asarray(a[field+'_integrands'])-b[field+'_integrands'])))
        for a,b in zip(row['phases'],saved['phases']) for field in FIELDS)
    return errors


def run_pilot():
    started=time.perf_counter();inputs=_hashes();result=dict(stage='Selected-frequency pilot only; no new frequency integral',
        parameters=dict(q=1,b=.2,d=1.,R=2.,L=40.,rho_c=2.,phase_settings=[0.,math.pi],
            primary_h=.05,refined_h=.035,complete_outer_J=9),
        basis_checks=basis_checks(),algebra_checks=algebra_checks(),legacy_witnesses=[],
        complete_cuts=[],mesh_checks=[],input_sha256=inputs)
    old=json.loads(project_path('ads_reference_source_results.json').read_text())
    try:
        for mu in (-.5,.5,1.5):
            case=CompleteSector(mu,nmax=8,outer_nmax=8)
            row=case.evaluate(.4,save_matrices=True)
            group=next(g for g in old['mesh_refinement'] if g['config']['mu']==mu and g['config']['h']==.05)
            saved=next(r for r in group['samples'] if r['kappa']==.4)
            comparison=legacy_comparison(row,saved)
            result['legacy_witnesses'].append(dict(config=case.config,assembly=case.assembly,
                ports=case.port_checks,samples=[row],legacy_errors=comparison))
            print(json.dumps(dict(stage='legacy',mu=mu,max_error=max(comparison.values()),seconds=row['seconds'])),flush=True)
        for N,frequencies in ((2,(.01,.4,2.,12.)),(4,(2.,12.))):
            cut=dict(nmax=N,outer_nmax=9,h=.05,sectors=[],frequency_sums=[])
            result['complete_cuts'].append(cut)
            for m in range(-N,N+1):
                case=CompleteSector(m+.5,nmax=N,outer_nmax=9)
                samples=[case.evaluate(k,save_matrices=(k==2.)) for k in frequencies]
                cut['sectors'].append(dict(config=case.config,assembly=case.assembly,ports=case.port_checks,samples=samples))
                print(json.dumps(dict(stage='complete',N=N,m=m,seconds=sum(r['seconds'] for r in samples),
                    pi_pressure=[dict(kappa=r['kappa'],value=r['phases'][1]['pressure_integrands']) for r in samples])),flush=True)
            for i,kappa in enumerate(frequencies):
                phases={phase:{field:np.sum([s['samples'][i]['phases'][ip][field+'_integrands']
                    for s in cut['sectors']],axis=0).tolist() for field in FIELDS} for ip,phase in enumerate(PHASES)}
                cut['frequency_sums'].append(dict(kappa=kappa,phases=phases))
        for m in (-4,4):
            case=CompleteSector(m+.5,h=.035,nmax=4,outer_nmax=9)
            samples=[case.evaluate(k,save_matrices=(k==2.)) for k in (2.,12.)]
            base=next(s for s in result['complete_cuts'][1]['sectors'] if s['config']['m']==m)
            differences=[]
            for row,before in zip(samples,base['samples']):
                differences.append(dict(kappa=row['kappa'],phases={phase:{field:
                    (np.asarray(row['phases'][ip][field+'_integrands'])-before['phases'][ip][field+'_integrands']).tolist()
                    for field in FIELDS} for ip,phase in enumerate(PHASES)}))
            result['mesh_checks'].append(dict(config=case.config,assembly=case.assembly,ports=case.port_checks,
                samples=samples,refined_minus_primary=differences))
            print(json.dumps(dict(stage='mesh',m=m,seconds=sum(r['seconds'] for r in samples),
                pi_pressure_difference=[d['phases']['pi']['pressure'] for d in differences])),flush=True)
    except Exception as exc:
        result['execution_failure']=dict(type=type(exc).__name__,message=str(exc))
    groups=result['legacy_witnesses']+[s for c in result['complete_cuts'] for s in c['sectors']]+result['mesh_checks']
    samples=[row for g in groups for row in g['samples']]
    algebra=result['algebra_checks'];ports=[p for g in groups for p in g['ports']]
    def maximum(values):return max(values,default=0.)
    result['diagnostics']=dict(
        maximum_equation_residual=maximum(e for r in samples for e in r['exterior_equation_relative_residuals']),
        minimum_sewing_singular_value=min((v for r in samples for p in r['phases'] for v in p['minimum_singular_values']),default=0.),
        maximum_sewing_condition=maximum(v for r in samples for p in r['phases'] for v in p['condition_numbers']),
        maximum_port_gram_condition=maximum(v for p in ports for v in p['gram_condition_numbers']),
        maximum_resolvent_identity_error=maximum(p['resolvent_difference_identity_error'] for r in samples for p in r['phases']),
        maximum_coordinate_invariance_error=maximum(r['raw_facet_observable_error'] for r in samples if r['raw_facet_observable_error'] is not None),
        maximum_reference_phase_error=maximum(max(r['reference_phase_pressure_error'],r['reference_phase_shift_error']) for r in samples),
        maximum_reference_shift=maximum(abs(v) for r in samples for p in r['phases'] for v in p['reference_shift_trace']),
        maximum_tail_gram_error=maximum(r['tail_diagnostics']['angular_gram_error'] for r in samples),
        fresh_exterior_frequency_samples=len(samples),complete_cut_multiplicity='Each m and both Weyl operators actually solved; no degeneracy multiplier.')
    d=result['diagnostics']
    expected=3+5*4+9*2+2*2
    result['tests']=dict(
        execution_completed='execution_failure' not in result and len(samples)==expected,
        complete_basis=result['basis_checks']['passed'],
        complete_counts=all(c['paired_labels']==c['expected_paired'] and c['single_weyl_mouth_dimension']==c['expected_dimension'] for c in algebra['counts']),
        reference_selector_orientation=algebra['reference_opposite_selector_error']<1e-13,
        legacy_equal_cutoff_reproduction=len(result['legacy_witnesses'])==3 and maximum(v for g in result['legacy_witnesses'] for v in g['legacy_errors'].values())<2e-10,
        free_axis_equations=maximum(r['pde_error'] for r in algebra['free_axis_spinors'])<2e-13,
        input_validation=len(algebra['invalid_configurations_rejected'])==4,
        independent_opposite_assembly=maximum(g['assembly']['opposite_adjoint_relative_error'] for g in groups)<2e-12,
        positive_imaginary_mass=all(g['assembly']['minimum_cell_mass_eigenvalue']>0 and g['assembly']['axis_domain_equal'] for g in groups),
        physical_swapped_ports=maximum(max(p['basis_reproduction_error'],p['incoming_adjoint_error'],p['outgoing_adjoint_error'],p['factor_swap_error']) for p in ports)<2e-11,
        equation_residual=d['maximum_equation_residual']<2e-10,
        sampled_sewing_nonsingular=d['minimum_sewing_singular_value']>1e-8,
        resolvent_identity=d['maximum_resolvent_identity_error']<2e-10,
        coordinate_invariance=d['maximum_coordinate_invariance_error']<2e-10,
        reference_phase_and_flux=max(d['maximum_reference_phase_error'],d['maximum_reference_shift'])<2e-11,
        tail_basis=d['maximum_tail_gram_error']<2e-11,
        null_identities=maximum(float(np.max(abs(np.asarray(p['null_plus_integrands'])-np.asarray(p['pressure_integrands'])+2*np.asarray(p['shift_integrands'])))) for r in samples for p in r['phases'])<2e-13,
        frozen_inputs_preserved=_hashes()==inputs)
    result['passed']=all(result['tests'].values())
    result['probes']=probe_geometry(.2*math.sqrt(5),2.)
    result['observable_definitions']=dict(pressure='(1/pi) ReTr[(R_connected-R_AdS) dB_(rho+p)]',
        shift='(1/pi) ReTr[(R_connected-R_AdS) dB_shift]; equals minus the weighted lower-index energy flux',
        null_plus='pressure-2*shift',null_minus='pressure+2*shift',
        frequency_scope='Kernels evaluated at selected positive kappa; no integration or ultraviolet-tail bound in this artifact.')
    result['scope']=[
        'Reflected isolated-sphere spin/gauge sewing is a declared geometric branch. No physical finite collar or Einstein junction is solved.',
        'All m=-N,...,N are included at each complete mouth cutoff; higher n, outer J, frequency, spatial mesh and finite tail still require convergence.',
        'Outer J is independent of mouth N. The legacy witnesses deliberately use J=N=8; complete pilots use J=9.',
        'The two Weyl bulk operators and physical swapped Gram maps are assembled separately. No silent chirality or angular degeneracy multiplier is used.',
        'The AdS reference couples the two chiral zero modes only for m=0. Massive caps use their absolute n labels.',
        'Same-sign facet projection is a geometric approximation; port Gram conditioning and selected extreme-m mesh changes are reported.',
        'This is a finite-cut reference-relative source diagnostic, not a complete absolute stress, NEC result, material recipe or stable chronology-changing medium.']
    result['sources']=[dict(url='https://arxiv.org/pdf/hep-th/0206054',use='Section4.1 monopole Dirac spectrum; finite q=1 angular labels and degeneracies.'),
        dict(url='https://arxiv.org/pdf/2208.08252',use='AdS2 Dirac canonical operator and invariant normalizable reference domain, Eqs21--25 and SecIII.C--D.')]
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['seconds']=time.perf_counter()-started
    return result


def reuse_checks(pilot):
    """Saved, independently solved N2/N4 matrices test the exact reuse law."""
    rows=[]
    for low in pilot['complete_cuts'][0]['sectors']:
        m=low['config']['m'];high=next(s for s in pilot['complete_cuts'][1]['sectors'] if s['config']['m']==m)
        a=next(r for r in low['samples'] if r['kappa']==2.)
        b=next(r for r in high['samples'] if r['kappa']==2.)
        lv=a['matrices'];hv=b['matrices']
        value,info=reproject_environment(_decode(hv['S_connected']),m+.5,4,2,
            _decode(hv['F_mouth_right_out']),_decode(hv['F_mouth_left_out']))
        actual=_decode(lv['S_connected']);Wlow=block_diag(_decode(lv['F_mouth_right_out']),_decode(lv['F_mouth_left_out']))
        Whigh=block_diag(_decode(hv['F_mouth_right_out']),_decode(hv['F_mouth_left_out']))
        indices=info['outgoing_indices'];Whigh=Whigh.conj().T@Whigh
        gram_error=float(np.max(abs(Wlow.conj().T@Wlow-Whigh[np.ix_(indices,indices)])))
        kernel_error=0.
        for prior,matrices in zip(a['phases'],lv['phases']):
            test=reference_relative(value,_decode(lv['S_reference']),_decode(matrices['B']),
                _decode(matrices['pressure']),_decode(matrices['shift']))
            kernel_error=max(kernel_error,*[float(np.max(abs(np.asarray(prior[f+'_integrands'])-test[f+'_integrands']))) for f in FIELDS])
        rows.append(dict(m=m,kappa=2.,high_nmax=4,low_nmax=2,outer_nmax=9,h=.05,
            environment_error=float(np.max(abs(value-actual))),restricted_gram_error=gram_error,
            source_kernel_error=kernel_error,**info))
    return dict(rows=rows,passed=all(r['environment_error']<2e-12 and r['source_kernel_error']<2e-12 and r['restricted_gram_error']<2e-13 for r in rows),
        formula='S_low=W_low^-1 Pout W_high S_high Pin^T; W_low=Pout W_high Pout^T',
        assumptions='Identical m, mesh, bulk coefficients, outer cutoff J and tail; raw input lift is linear. Re-solve the smaller core/reference sewing. No reuse between different m or J.')


def integrate_cuts(case,order):
    """One exterior solve per node; both complete-cut source kernels share it."""
    stamp=time.perf_counter();nodes,weights=leggauss(order)
    cuts=(2,4) if abs(case.config['m'])<=2 else (4,)
    records={str(N):dict(config=dict(case.config,nmax=N,
        mouth_channels=mouth_channels(case.config['mu'],N).tolist(),
        reflection_cross_signs=reflection_parities(case.config['mu'],N).tolist()),
        order=order,edges=list(EDGES),panels=[]) for N in cuts}
    minima={str(N):np.full((2,2),np.inf) for N in cuts}
    conditions={str(N):np.zeros((2,2)) for N in cuts}
    errors={str(N):dict(equation=0.,identity=0.,reference_shift=0.,phase=0.) for N in cuts}
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        ks=(lo+hi)/2+(hi-lo)*nodes/2;ws=(hi-lo)*weights/2
        values={str(N):{p:{f:[] for f in FIELDS} for p in PHASES} for N in cuts}
        singular={str(N):[] for N in cuts};condition_samples={str(N):[] for N in cuts}
        for kappa in ks:
            high=case.evaluate(float(kappa),lower_nmax=2)
            for N in cuts:
                key=str(N);row=high if N==4 else high['lower_cut_response']
                errors[key]['equation']=max(errors[key]['equation'],*row['exterior_equation_relative_residuals'])
                sv=[p['minimum_singular_values'] for p in row['phases']]
                cc=[p['condition_numbers'] for p in row['phases']]
                minima[key]=np.minimum(minima[key],sv);conditions[key]=np.maximum(conditions[key],cc)
                singular[key].append(sv);condition_samples[key].append(cc)
                for ip,phase in enumerate(PHASES):
                    sample=row['phases'][ip]
                    for field in FIELDS:values[key][phase][field].append(sample[field+'_integrands'])
                    errors[key]['identity']=max(errors[key]['identity'],sample['resolvent_difference_identity_error'])
                    errors[key]['reference_shift']=max(errors[key]['reference_shift'],*abs(np.asarray(sample['reference_shift_trace'])))
                errors[key]['phase']=max(errors[key]['phase'],row['reference_phase_pressure_error'],row['reference_phase_shift_error'])
        for N in cuts:
            key=str(N);integrands=values[key]
            responses={p:{f:np.einsum('i,ij->j',ws,integrands[p][f]).tolist() for f in FIELDS} for p in PHASES}
            records[key]['panels'].append(dict(bounds=[lo,hi],kappa=ks.tolist(),weights=ws.tolist(),
                integrands=integrands,responses=responses,minimum_singular_values=singular[key],condition_numbers=condition_samples[key]))
        if hi in (1.,12.,96.):
            print(json.dumps(dict(stage='integration_panel',m=case.config['m'],h=case.config['h'],order=order,
                through_kappa=hi,seconds=time.perf_counter()-stamp)),flush=True)
    for N in cuts:
        key=str(N);r=records[key]
        r['totals']={p:{f:np.sum([panel['responses'][p][f] for panel in r['panels']],axis=0).tolist() for f in FIELDS} for p in PHASES}
        r['diagnostics']=dict(minimum_singular_values=minima[key].tolist(),maximum_condition_numbers=conditions[key].tolist(),
            maximum_equation_relative_residual=errors[key]['equation'],resolvent_identity_error=errors[key]['identity'],
            reference_shift_trace_error=errors[key]['reference_shift'],reference_phase_error=errors[key]['phase'])
        r['seconds_shared_exterior']=time.perf_counter()-stamp
        r['exterior_reuse']='Fresh N4 mouth-response columns' if N==4 else 'Exact weighted output-Gram reprojection of the same N4 exterior response; independent N2 core/reference sewing'
    return records


def outer_cutoff_checks(pilot):
    rows=[]
    for m in (-4,0,4):
        case=CompleteSector(m+.5,nmax=4,outer_nmax=11)
        old=next(s for s in pilot['complete_cuts'][1]['sectors'] if s['config']['m']==m)
        samples=[]
        for kappa in (2.,12.):
            row=case.evaluate(kappa);base=next(r for r in old['samples'] if r['kappa']==kappa)
            changes={phase:{field:(np.asarray(row['phases'][ip][field+'_integrands'])-base['phases'][ip][field+'_integrands']).tolist()
                for field in FIELDS} for ip,phase in enumerate(PHASES)}
            samples.append(dict(**row,J11_minus_J9=changes))
        rows.append(dict(config=case.config,assembly=case.assembly,ports=case.port_checks,samples=samples))
    return dict(rows=rows,scope='Selected-frequency independent outer-cutoff sensitivity, not an integrated or certified angular-tail error.')


def write_result(result):
    output=project_path('complete_angular_source_results.json')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(output.read_text(encoding='utf-8'))==result


def run_integration():
    """Retain the executed pilot, then complete the specifically authorized sum."""
    started=time.perf_counter();path=project_path('complete_angular_source_results.json')
    result=json.loads(path.read_text(encoding='utf-8')) if path.exists() else run_pilot()
    if 'pilot_provenance' not in result:
        text=Path(__file__).read_text(encoding='utf-8')
        result['pilot_provenance']=dict(execution_source_sha256=hashlib.sha256(text.encode('utf-8')).hexdigest(),
            execution_source_text=text,seconds=result['seconds'])
    provenance=result['pilot_provenance']
    if hashlib.sha256(provenance['execution_source_text'].encode('utf-8')).hexdigest()!=provenance['execution_source_sha256']:
        raise ValueError('Pilot source provenance does not match its retained text')
    if (result['input_sha256'].keys() != _hashes().keys() or not all(source_matches(name, digest) for name, digest in result['input_sha256'].items())):raise ValueError('A frozen pilot input has changed')
    pilot_tests=result.get('pilot_tests',result['tests'])
    if not all(pilot_tests.values()):raise ValueError('The retained pilot must pass before integration')
    result['pilot_tests']=pilot_tests;result['pilot_passed']=True
    result['pilot_diagnostics']=result.get('pilot_diagnostics',result['diagnostics'])
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    hashes={**_hashes(),'angular_reference_sum.py':hashlib.sha256(project_path('angular_reference_sum.py').read_bytes()).hexdigest()}
    result.update(passed=False,tests={},stage='Running complete finite-cut source integrals',
        source_sha256=source_hash,integration_input_sha256=hashes,integration_seconds=0.)
    result['cut_reuse_checks']=reuse_checks(result)
    if not result['cut_reuse_checks']['passed']:raise ValueError('Exact saved-matrix reuse witness failed')
    result['integration']=dict(parameters=dict(edges=list(EDGES),cutoffs=[2,4],outer_nmax=9,
        recipes=dict(coarse=dict(h=.05,order=8),fine=dict(h=.05,order=12),finer_mesh=dict(h=.035,order=8)),
        individual_observed_quadrature_target=2e-6,complete_N4_sum_target=1.8e-5,complete_N2_sum_target=1e-5,
        sign_is_not_an_acceptance_gate=True),cutoffs={'2':dict(sectors=[]),'4':dict(sectors=[])})
    write_result(result)
    try:
        result['outer_cutoff_checks']=outer_cutoff_checks(result)
        write_result(result)
        for m in range(-4,5):
            tick=time.perf_counter();case=CompleteSector(m+.5,nmax=4,outer_nmax=9)
            coarse=integrate_cuts(case,8);fine=integrate_cuts(case,12)
            small=CompleteSector(m+.5,h=.035,nmax=4,outer_nmax=9)
            finer=integrate_cuts(small,8)
            for key in coarse:
                record=dict(mu=m+.5,m=m,coarse=coarse[key],fine=fine[key],finer_mesh=finer[key])
                record['comparisons']=comparisons(record);record['controls']=record_controls(record)
                result['integration']['cutoffs'][key]['sectors'].append(record)
            result['stage']='Completed integration sector m='+str(m)
            result['integration_seconds']=time.perf_counter()-started;write_result(result)
            print(json.dumps(dict(stage='sector_integrated',m=m,seconds=time.perf_counter()-tick,
                N4_primary_pi=fine['4']['totals']['pi'],N4_finer_mesh_pi=finer['4']['totals']['pi'],
                controls=result['integration']['cutoffs']['4']['sectors'][-1]['controls'])),flush=True)
            del case,small
        for key,cut in result['integration']['cutoffs'].items():
            N=int(key);sectors=cut['sectors']
            if [s['m'] for s in sectors]!=list(range(-N,N+1)):raise ValueError('Incomplete angular multiplet')
            if any(s[r]['config']['outer_nmax']!=9 for s in sectors for r in RECIPES):raise ValueError('Mismatched outer J')
            sums={recipe:sum_records(sectors,recipe) for recipe in RECIPES}
            for record in sums.values():
                for name in ('m','mouth_channels','outer_channels','has_zero_mode'):
                    record['config'].pop(name,None)
                record['complete_mouth_nmax']=N;record['m_values']=list(range(-N,N+1))
            cut['partial_sums']=sums;cut['comparisons']=comparisons(sums)
            cut['sum_of_sector_quadrature_envelopes']={phase:{field:np.sum([
                r['comparisons']['quadrature_panel_envelope'][phase][field] for r in sectors],axis=0).tolist()
                for field in FIELDS} for phase in PHASES}
        cuts=result['integration']['cutoffs']
        result['integration']['N4_minus_N2']={recipe:{phase:{field:(
            np.asarray(cuts['4']['partial_sums'][recipe]['totals'][phase][field])-
            cuts['2']['partial_sums'][recipe]['totals'][phase][field]).tolist()
            for field in FIELDS} for phase in PHASES} for recipe in RECIPES}
        null_error=aggregate_error=0.
        for cut in cuts.values():
            for record in cut['partial_sums'].values():
                for phase in PHASES:
                    p=np.asarray(record['totals'][phase]['pressure']);s=np.asarray(record['totals'][phase]['shift'])
                    null_error=max(null_error,float(np.max(abs(np.asarray(record['totals'][phase]['null_plus'])-p+2*s))),
                        float(np.max(abs(np.asarray(record['totals'][phase]['null_minus'])-p-2*s))))
                    for field in FIELDS:
                        value=sum((np.einsum('i,ij->j',v['weights'],v['integrands'][phase][field]) for v in record['panels']),np.zeros(2))
                        aggregate_error=max(aggregate_error,float(np.max(abs(value-record['totals'][phase][field]))))
        def envelope(key):
            return max(v for phase in cuts[key]['sum_of_sector_quadrature_envelopes'].values() for field in phase.values() for v in field)
        tests=dict(pilot_passed=all(pilot_tests.values()),exact_cut_reuse=result['cut_reuse_checks']['passed'],
            all_complete_sector_controls=all(s['controls']['passed'] for cut in cuts.values() for s in cut['sectors']),
            N2_observed_quadrature=envelope('2')<5*2e-6,N4_observed_quadrature=envelope('4')<9*2e-6,
            summed_null_identities=null_error<1e-12,independent_integral_reaggregation=aggregate_error<1e-12,
            all_inputs_preserved=hashes=={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in hashes},
            execution_source_preserved=source_hash==hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        result.update(passed=all(tests.values()),tests={k:bool(v) for k,v in tests.items()},
            stage='Complete N2/N4 finite-cut integrals with actual paired Weyl responses',
            integration_seconds=time.perf_counter()-started)
        result['integration']['diagnostics']=dict(summed_null_identity_error=null_error,
            independent_reaggregation_error=aggregate_error,
            N2_maximum_summed_quadrature_envelope=envelope('2'),N4_maximum_summed_quadrature_envelope=envelope('4'),
            omitted_frequency_tail_bound=None,omitted_angular_tail_bound=None,continuum_error_bound=None,
            fresh_exterior_frequency_samples=9*8*(8+12+8)+6,
            N2_extra_exterior_solves=0)
        result['scope']=[v for v in result['scope'] if 'frequency integral' not in v and 'selected-frequency' not in v]
        result['scope']+=['The pilot used the exact retained source text in pilot_provenance. Integration and six J11 controls use the separately recorded final source SHA.',
            'The frequency integral covers [1e-6,96]. Observed quadrature, matched-order mesh, N4-minus-N2 and selected outer-J changes are separate diagnostics, not certified omitted-tail bounds.',
            'N2 environments are obtained by a saved-witness-validated weighted output-Gram reprojection of N4 at the same m, mesh and outer J. Their nonlinear throat/reference sewing is solved independently.',
            'Both phase settings compare connected and AdS-reference environments separately. No multiplication by a degeneracy or division by a volume supplies an absolute completed stress.',
            'The prescribed constant-radius core has zero geometric radial-null contraction. A nonzero completed source would require a changed self-consistent geometry.']
        result['observable_definitions']['frequency_scope']='Finite integrals over [1e-6,96], with every quadrature node/weight retained. No certified omitted tail.'
        write_result(result)
    except Exception as exc:
        result.update(passed=False,stage='Stopped with retained integration diagnostics',
            integration_failure=dict(type=type(exc).__name__,message=str(exc)),integration_seconds=time.perf_counter()-started)
        write_result(result)
        raise
    return result


if __name__=='__main__':
    data=run_integration()
    print(json.dumps(dict(passed=data['passed'],tests=data['tests'],integration_seconds=data['integration_seconds'],
        primary_N4=data['integration']['cutoffs'].get('4',{}).get('partial_sums',{}).get('fine',{}).get('totals'))),flush=True)
    if not data['passed']:raise SystemExit('A numerical gate failed; all calculated diagnostics were saved.')
