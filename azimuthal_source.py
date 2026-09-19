"""Explicit missing azimuthal sectors of the retained two-chirality source.

Only q=1 and mu=-1/2,+1/2,+3/2 are covered by the independently checked mouth
bases.  Both Weyl operators and their physical projected boundary maps are
assembled separately.  No angular degeneracy multiplier or absolute vacuum
source is inferred.  The synthetic geometry and its endpoint identification
remain specified inputs, not a solved Einstein--matter junction.

Run: .venv\\Scripts\\python.exe -B -W error azimuthal_source.py
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import block_diag as sparse_block_diag
from scipy.special import jv, jvp

from azimuthal_modes import mouth_characteristics, zero_mode, basis_checks, reflection_parities
from dirac_sewing_audit import (opposite_assembly, port_factors, convert_map,
    exterior_response, observables, integrate)
from exterior_boundary_map import ExteriorDG, ID, SX, SY, SZ, _complex_json
from exterior_tail_transport import facet_map
from imaginary_weyl_leads import tail_map
from spherical_weyl_leads import outer_characteristics
from throat_null_response import raw_variation, shift_derivatives, phase_map, probe_geometry
from two_mouth_mesh import mesh


def sew_matrix(value,mu,angle,nmax,reflected):
    rows=(np.arange(2*nmax+1) if mu==.5 else np.arange(1,2*nmax+1))
    cols=(np.arange(2*nmax+1) if mu==.5 else np.r_[np.arange(nmax),np.arange(nmax+1,2*nmax+1)])
    signs=reflection_parities(mu,nmax) if reflected else np.ones(nmax)
    result=phase_map(value,angle,nmax).copy();n=np.arange(1,nmax+1)
    result[n,nmax+n]*=signs;result[nmax+n,n-1]*=signs
    return result[np.ix_(rows,cols)]


def mouth_maps(kappa,nmax,probes,mu,rho_c=2.,reflected=False):
    """Remove an absent zero channel AFTER phases in the frozen mouth order.

    For mu!=1/2, incoming/outgoing order is [upper n1..N,lower n1..N].
    The mu=1/2 zero keeps its original index and sign.  Reflection parity
    changes both cross transmissions and their variations only.
    """
    raw,pressure,_=raw_variation(kappa,nmax,rho_c)
    shift=shift_derivatives(raw,kappa,nmax,probes)
    def trim(value,angle):
        return sew_matrix(value,mu,angle,nmax,reflected)
    Bs=[trim(raw,a) for a in (0.,math.pi)]
    Ps=[[trim(p,a) for p in pressure] for a in (0.,math.pi)]
    Qs=[[trim(p,a) for p in shift] for a in (0.,math.pi)]
    return Bs,Ps,Qs


class AzimuthalSector:
    """Minimal geometry adapter around the frozen DG and sewing utilities."""

    def __init__(self,mu,h=.08,R=2.,nmax=8,L=40.,rho_c=2.):
        if mu not in (-.5,.5,1.5):raise ValueError('Only the three checked q=1 azimuthal sectors are supported')
        if (not all(math.isfinite(float(v)) for v in (h,R,L,rho_c)) or h<=0
                or not L>R>.7 or rho_c<=0 or isinstance(nmax,bool)
                or int(nmax)!=nmax or nmax<2):
            raise ValueError('Require h/rho_c>0, L>R>.7, and integer nmax>=2')
        nmax=int(nmax)
        self.config=dict(mu=float(mu),h=float(h),R=float(R),nmax=nmax,L=float(L),rho_c=float(rho_c))
        self.config['reflection_cross_signs']=reflection_parities(mu,nmax).tolist()
        self.model=ExteriorDG(mesh(b=.2,d=1.,R=R,h=h),q=1,mu=mu,energy=0.)
        self.active=self.model.active;self.unshifted=self.model.matrix_full
        self.mass=sparse_block_diag([np.kron(c['mass'],ID) for c in self.model.cells],format='csc')
        self.left_zero,self.proxy=opposite_assembly(self.model)
        self.right_ports=[];self.left_ports=[];self.right_factors={};self.left_factors={}
        self.port_checks=[];self.bases={}
        for kind in ('upper','lower','outer'):
            def basis(points,kind=kind):
                if kind=='outer':
                    return outer_characteristics(mu,np.arctan2(points[:,0],points[:,1]),R,nmax)
                theta=np.arctan2(points[:,0],points[:,1]-(.5 if kind=='upper' else -.5))
                plus,minus=mouth_characteristics(mu,theta,.2,kind,nmax)
                zero=zero_mode(mu,theta,.2,kind)
                if (zero is not None)!=(mu==.5):raise AssertionError('Incorrect zero-channel count')
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
            self.port_checks.append(dict(kind=kind,incoming_count=right['G'].shape[1],outgoing_count=right['C'].shape[0],
                basis_reproduction_error=max(re,le),
                incoming_adjoint_error=float(np.max(abs(left['G']-right['C'].conj().T))),
                outgoing_adjoint_error=float(np.max(abs(left['C']-right['G'].conj().T))),
                factor_swap_error=max(float(np.max(abs(lf[j]-rf[1-j]))) for j in (0,1))))
        ur,lr,self.outer=self.right_ports;ul,ll,self.left_outer=self.left_ports
        self.Gm=np.column_stack((ur['G'],lr['G']));self.Cm=np.vstack((ur['C'],lr['C']))
        self.GL=np.column_stack((ul['G'],ll['G']));self.CL=np.vstack((ul['C'],ll['C']))
        self.Fin=block_diag(self.right_factors['upper'][1],self.right_factors['lower'][1])
        self.Fout=block_diag(self.right_factors['upper'][0],self.right_factors['lower'][0])
        self.Lin=block_diag(self.left_factors['upper'][1],self.left_factors['lower'][1])
        self.Lout=block_diag(self.left_factors['upper'][0],self.left_factors['lower'][0])
        diff=(self.left_zero-self.unshifted.conj().T).tocoo()
        minimum_mass=min(float(np.linalg.eigvalsh(c['mass'])[0]) for c in self.model.cells)
        axis=[]
        for label,c in (('outside',mu),('between',mu-1)):
            upper=(c-.5 if c>0 else -c+.5);lower=(c+.5 if c>0 else -c-.5)
            axis.append(dict(segment=label,c=c,regular_powers=[upper,lower]))
        self.assembly=dict(triangles=len(self.model.triangles),active_dof=len(self.active),
            opposite_adjoint_relative_error=float(np.max(abs(diff.data))/np.max(abs(self.unshifted.data))) if diff.nnz else 0.,
            minimum_cell_mass_eigenvalue=minimum_mass,axis_domain_equal=bool(np.array_equal(self.active,self.proxy.active)),
            axis_regular_powers=axis,axis_removed_dof=self.model.ndof-len(self.active),
            axis_scope='Both components vanish where |mu-Aaxis|>1/2; P1 resolves higher regular powers by mesh refinement, not exactly.')
        self.probes=probe_geometry(.2*math.hypot(1.,rho_c),rho_c);self.cache={}

    def evaluate(self,kappa,save_matrices=False,adjoint_control=True):
        kappa=float(kappa)
        if not math.isfinite(kappa) or kappa<=0:raise ValueError('Positive finite imaginary frequency required')
        key=(kappa,bool(save_matrices),bool(adjoint_control))
        if key in self.cache:return self.cache[key]
        stamp=time.perf_counter();N=self.config['nmax']
        rawZ,tail=tail_map(kappa,self.config['R'],self.config['L'],mu=self.config['mu'],nmax=N)
        ZR,FR,diagR=facet_map(self.outer,self.bases['outer'],rawZ)
        swapped=lambda points:self.bases['outer'](points)[::-1]
        ZL,FL,diagL=facet_map(self.left_outer,swapped,rawZ.conj().T)
        AR=self.unshifted+kappa*self.mass;AL=self.left_zero+kappa*self.mass
        SR,resR=exterior_response(AR,self.active,self.Gm,self.Cm,self.outer['G'],self.outer['C'],ZR)
        SL,resL=exterior_response(AL,self.active,self.GL,self.CL,self.left_outer['G'],self.left_outer['C'],ZL)
        if adjoint_control:
            SA,resA=exterior_response(AL,self.active,self.GL,self.CL,self.left_outer['G'],self.left_outer['C'],ZR.conj().T)
        else:SA,resA=None,None
        raw,dp,ds=mouth_maps(kappa,N,self.probes,self.config['mu'],self.config['rho_c'],True)
        BR=[convert_map(v,self.Fin,self.Fout) for v in raw]
        PR=[[convert_map(v,self.Fin,self.Fout) for v in p] for p in dp]
        QR=[[convert_map(v,self.Fin,self.Fout) for v in p] for p in ds]
        BL=[convert_map(v.conj().T,self.Lin,self.Lout) for v in raw]
        PL=[[convert_map(v.conj().T,self.Lin,self.Lout) for v in p] for p in dp]
        QL=[[-convert_map(v.conj().T,self.Lin,self.Lout) for v in p] for p in ds]
        right=observables(SR,BR,PR,QR);physical=observables(SL,BL,PL,QL)
        if adjoint_control:
            adjoint=observables(SA,[b.conj().T for b in BR],
                [[v.conj().T for v in p] for p in PR],[[-v.conj().T for v in p] for p in QR])
        else:adjoint=None
        fields=('energy_integrand','pressure_integrands','shift_integrands','null_plus_integrands','null_minus_integrands')
        paired={field:(np.asarray(right[field])+physical[field]).tolist() for field in fields}
        row=dict(kappa=kappa,right=right,left_physical=physical,left_adjoint_prescription=adjoint,
            actual_two_chirality_sum=paired,exterior_equation_relative_residuals=[resR,resL]+([resA] if adjoint_control else []),
            tail_channels=tail['channels'],tail_angular_gram_error=tail['angular_gram_error'],
            discrepancies=dict(physical_energy_minus_right=physical['energy_integrand']-right['energy_integrand'],
                physical_pressure_minus_right=(np.array(physical['pressure_integrands'])-right['pressure_integrands']).tolist(),
                physical_shift_plus_right=(np.array(physical['shift_integrands'])+right['shift_integrands']).tolist(),
                outer_gram_adjoint_defect=float(np.linalg.norm(ZL-ZR.conj().T)),
                mouth_gram_adjoint_defect=max(float(np.linalg.norm(a-b.conj().T)) for a,b in zip(BL,BR)),
                deliberate_adjoint_response_error=float(np.linalg.norm(SA-SR.conj().T)) if adjoint_control else None,
                deliberate_adjoint_energy_error=abs(adjoint['energy_integrand']-right['energy_integrand']) if adjoint_control else None,
                deliberate_adjoint_pressure_error=float(np.max(abs(np.array(adjoint['pressure_integrands'])-right['pressure_integrands']))) if adjoint_control else None,
                deliberate_adjoint_shift_error=float(np.max(abs(np.array(adjoint['shift_integrands'])+right['shift_integrands']))) if adjoint_control else None),
            seconds=time.perf_counter()-stamp)
        old_raw,old_dp,old_ds=mouth_maps(kappa,N,self.probes,self.config['mu'],self.config['rho_c'],False)
        oldR=observables(SR,[convert_map(v,self.Fin,self.Fout) for v in old_raw],
            [[convert_map(v,self.Fin,self.Fout) for v in p] for p in old_dp],
            [[convert_map(v,self.Fin,self.Fout) for v in p] for p in old_ds])
        oldL=observables(SL,[convert_map(v.conj().T,self.Lin,self.Lout) for v in old_raw],
            [[convert_map(v.conj().T,self.Lin,self.Lout) for v in p] for p in old_dp],
            [[-convert_map(v.conj().T,self.Lin,self.Lout) for v in p] for p in old_ds])
        old_pair={field:(np.asarray(oldR[field])+oldL[field]).tolist() for field in fields}
        row['inherited']=dict(right=oldR,left_physical=oldL,actual_two_chirality_sum=old_pair)
        if save_matrices:
            matrices=dict(S_right=SR,S_left_physical=SL,Z_raw_right=rawZ,Z_facet_right=ZR,Z_facet_left=ZL,
                F_outer_right_plus=FR[0],F_outer_right_minus=FR[1],F_outer_left_plus=FL[0],F_outer_left_minus=FL[1],
                F_mouth_right_in=self.Fin,F_mouth_right_out=self.Fout,F_mouth_left_in=self.Lin,F_mouth_left_out=self.Lout,
                B_raw_right=raw,pressure_raw_right=dp,shift_raw_right=ds,inherited_B_raw_right=old_raw,
                B_right=BR,B_left_physical=BL,pressure_right=PR,pressure_left_physical=PL,shift_right=QR,shift_left_physical=QL)
            if adjoint_control:matrices['S_left_adjoint_prescription']=SA
            row['matrices']={name:_complex_json(value) for name,value in matrices.items()}
        self.cache[key]=row;return row


def free_axis_checks():
    """Independent free Bessel spinors test the two new regular-axis domains."""
    energy=1.3;transverse=.8;kz=math.sqrt(energy**2-transverse**2);rows=[]
    for mu in (-.5,.5,1.5):
        m=int(mu-.5);ratio=transverse/(energy+kz)
        def exact(points):
            rho,z=points.T;phase=np.exp(1j*kz*z)
            return phase[:,None]*np.column_stack((jv(m,transverse*rho),1j*ratio*jv(m+1,transverse*rho)))
        pts=np.array([[.07,-.6],[.2,.2],[.9,.1]])
        rho=pts[:,0];u=exact(pts);phase=np.exp(1j*kz*pts[:,1])
        dr=transverse*phase[:,None]*np.column_stack((jvp(m,transverse*rho),1j*ratio*jvp(m+1,transverse*rho)))
        residual=(dr+u/(2*rho[:,None]))@SX.T+1j*kz*u@SZ.T+1j*mu/rho[:,None]*(u@SY.T)-1j*energy*u
        refinements=[]
        for h in (.16,.10):
            model=ExteriorDG(mesh(b=.2,d=1.,R=1.3,h=h),q=0,mu=mu,energy=energy)
            rhs,data,_=model.input(lambda kind,points,normal:exact(points))
            solution=model.solve(rhs);ledger=model.ledger(solution,data,rhs)
            refinements.append(dict(h=h,triangles=len(model.triangles),relative_error=model.l2_error(solution,exact),ledger=ledger))
        rows.append(dict(mu=mu,m=m,exact_pde_maximum_error=float(np.max(abs(residual))),refinements=refinements))
    return dict(formula='exp(i*kz*z)*(J_m(kperp*rho), i*kperp/(E+kz)*J_(m+1)(kperp*rho)), m=mu-1/2',
        energy=energy,k_perp=transverse,k_z=kz,rows=rows,
        scope='Free q=0 manufactured data on all boundaries validate the new axis sectors; they are not physical mouth injection or return probabilities.')


def pilot():
    start=time.perf_counter();controls=basis_checks();free=free_axis_checks();records=[]
    for mu in (-.5,.5,1.5):
        case=AzimuthalSector(mu)
        rows=[case.evaluate(k,save_matrices=True) for k in (.01,.4,2.)]
        records.append(dict(config=case.config,assembly=case.assembly,ports=case.port_checks,samples=rows))
        print(json.dumps(dict(mu=mu,assembly=case.assembly,samples=[dict(kappa=r['kappa'],paired=r['actual_two_chirality_sum'],
            inherited=r['inherited']['actual_two_chirality_sum'],minimum=min(r['left_physical']['minimum_singular_values']),seconds=r['seconds']) for r in rows])),flush=True)
    return dict(stage='Both specified endpoint sewings; selected-frequency pilot only.',basis_checks=controls,
        free_axis_checks=free,pilot=records,seconds=time.perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


class _BranchView:
    """Reuse the frozen quadrature on either result of the SAME Green solve."""
    def __init__(self,case,branch):
        self.case=case;self.branch=branch;self.config=dict(case.config,branch=branch)
    def evaluate(self,kappa,adjoint_control=False):
        row=self.case.evaluate(kappa,adjoint_control=adjoint_control)
        return row if self.branch=='reflection' else row['inherited']


def integrate_both(case,order,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
    return {branch:integrate(_BranchView(case,branch),order,edges) for branch in ('inherited','reflection')}


def replace_panel(fine,refined,index=2):
    panels=list(fine['panels']);panels[index]=refined['panels'][0]
    totals={name:{field:np.sum([p['responses'][name][field] for p in panels],axis=0).tolist()
                  for field in fine['totals'][name]} for name in fine['totals']}
    orders=[fine['order']]*len(panels);orders[index]=refined['order']
    return dict(fine,order=None,panel_orders=orders,panels=panels,totals=totals,
        scope='Explicit mixed16/24 quadrature; actual right and left physical responses, finite frequency window.')


def panel_envelope(before,after):
    fields=tuple(after['totals']['paired'])
    return {field:np.sum([abs(np.asarray(b['responses']['paired'][field])-a['responses']['paired'][field])
        for a,b in zip(before['panels'],after['panels'])],axis=0).tolist() for field in fields}


def final_envelope(coarse,fine,primary,index=2):
    fields=tuple(primary['totals']['paired'])
    changes=[]
    for i,(a,b,c) in enumerate(zip(coarse['panels'],fine['panels'],primary['panels'])):
        left,right=(b,c) if i==index else (a,b)
        changes.append({field:abs(np.asarray(right['responses']['paired'][field])-left['responses']['paired'][field]) for field in fields})
    return {field:np.sum([p[field] for p in changes],axis=0).tolist() for field in fields}


def compare_saved(current,saved):
    """Compare equal nodes/meshes, not an interpolation or an inferred factor."""
    for key in ('h','R','nmax','L','rho_c'):
        if current['config'][key]!=saved['config'][key]:raise AssertionError('Unmatched saved geometry: '+key)
    if current['edges']!=saved['edges']:raise AssertionError('Unmatched saved frequency panels')
    maximum=0.;node_error=0.
    for a,b in zip(current['panels'],saved['panels']):
        if len(a['kappa'])!=len(b['kappa']):raise AssertionError('Unmatched saved quadrature order')
        node_error=max(node_error,float(np.max(abs(np.array(a['kappa'])-b['kappa']))))
        for chirality in ('right','left_physical','paired'):
            for field in a['integrands'][chirality]:
                maximum=max(maximum,float(np.max(abs(np.array(a['integrands'][chirality][field])-b['integrands'][chirality][field]))))
    return dict(node_maximum_error=node_error,observable_maximum_error=maximum)


def _array(record):
    return np.asarray(record['real'])+1j*np.asarray(record['imag'])


def pressure_finite_difference(case,row):
    """Fresh finite-metric throat solves; exterior and endpoint bases are fixed."""
    N=case.config['nmax'];mu=case.config['mu'];kappa=row['kappa']
    SR=_array(row['matrices']['S_right']);SL=_array(row['matrices']['S_left_physical'])
    values={branch:[] for branch in ('inherited','reflection')}
    for step in (1e-3,5e-4):
        energies={branch:[] for branch in values}
        for eta in (-step,step):
            raw,_,_=raw_variation(kappa,N,case.config['rho_c'],eta,0,rtol=3e-13,atol=3e-15)
            for branch in values:
                logs=[]
                for phase in (0.,math.pi):
                    b=sew_matrix(raw,mu,phase,N,branch=='reflection')
                    br=convert_map(b,case.Fin,case.Fout);bl=convert_map(b.conj().T,case.Lin,case.Lout)
                    logs.append(np.linalg.slogdet(np.eye(len(br))-SR@br)[1]+np.linalg.slogdet(np.eye(len(bl))-SL@bl)[1])
                energies[branch].append(-(logs[1]-logs[0])/math.pi)
        for branch in values:values[branch].append(float((energies[branch][1]-energies[branch][0])/(2*step)))
    result={}
    for branch,estimates in values.items():
        source=row if branch=='reflection' else row['inherited']
        analytic=source['actual_two_chirality_sum']['pressure_integrands'][0]
        result[branch]=dict(analytic=analytic,steps=[1e-3,5e-4],centered_estimates=estimates,
            richardson_error=abs((4*estimates[1]-estimates[0])/3-analytic))
    return dict(mu=mu,kappa=kappa,probe='central',branches=result)


def sum_sectors(records,key):
    """Sum explicitly computed sectors at identical frequency recipe and mesh."""
    output={}
    for branch in ('inherited','reflection'):
        data=[r[key][branch] for r in records];first=data[0]
        for other in data[1:]:
            if other['edges']!=first['edges'] or other.get('panel_orders')!=first.get('panel_orders') or other['order']!=first['order']:
                raise AssertionError('Incompatible partial-sum quadratures')
            for name in ('h','R','nmax','L','rho_c'):
                if first['config'][name]!=other['config'][name]:raise AssertionError('Incompatible partial-sum geometry')
        config={name:first['config'][name] for name in ('h','R','nmax','L','rho_c','branch')}
        output[branch]=dict(config=config,sectors=[r['mu'] for r in records],order=first['order'],
            panel_orders=first.get('panel_orders'),edges=first['edges'],
            totals={chirality:{field:np.sum([r['totals'][chirality][field] for r in data],axis=0).tolist()
                    for field in first['totals'][chirality]} for chirality in first['totals']})
    return output


def checks():
    start=time.perf_counter();root=Path(__file__).resolve().parent
    inputs=('azimuthal_modes.py','dirac_sewing_audit.py','dirac_sewing_audit_results.json',
        'throat_null_response.py','throat_null_response_results.json','curved_throat_modes.py',
        'exterior_boundary_map.py','exterior_tail_transport.py','imaginary_weyl_leads.py',
        'magnetic_tail_map.py','spherical_weyl_leads.py','two_mouth_mesh.py')
    hashes={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs}
    saved=json.loads((root/'dirac_sewing_audit_results.json').read_text(encoding='utf-8'))
    basis=basis_checks();free=free_axis_checks();samples=[];integrations=[];fd=[];angular=[];legacy=[];all_rows=[]
    quadrature_target=2e-7;angular_target=2e-6
    for mu in (-.5,.5,1.5):
        small=AzimuthalSector(mu,h=.08)
        rows=[small.evaluate(k,save_matrices=True) for k in (.01,.4,2.)]
        samples.append(dict(config=small.config,assembly=small.assembly,ports=small.port_checks,samples=rows))
        all_rows.extend(small.cache.values());del small
        main=AzimuthalSector(mu,h=.05)
        rows=[main.evaluate(k,save_matrices=True) for k in (.01,.4,2.)]
        samples.append(dict(config=main.config,assembly=main.assembly,ports=main.port_checks,samples=rows))
        fd.append(pressure_finite_difference(main,rows[1]))
        stamp=time.perf_counter();coarse=integrate_both(main,8);fine=integrate_both(main,16)
        refined=integrate_both(main,24,(.02,.2))
        primary={branch:replace_panel(fine[branch],refined[branch]) for branch in fine}
        initial={branch:panel_envelope(coarse[branch],fine[branch]) for branch in fine}
        final={branch:final_envelope(coarse[branch],fine[branch],primary[branch]) for branch in fine}
        print(json.dumps(dict(mu=mu,h=.05,primary={b:primary[b]['totals']['paired'] for b in primary},
            initial_quadrature_changes=initial,final_quadrature_changes=final,seconds=time.perf_counter()-stamp)),flush=True)
        if mu==.5:
            for name,now,old in (('h05_order8',coarse['inherited'],saved['quadrature_refinement'][0]),
                ('h05_order16',fine['inherited'],saved['quadrature_refinement'][1]),
                ('h05_mixed_primary',primary['inherited'],saved['primary_quadrature'])):
                legacy.append(dict(record=name,**compare_saved(now,old)))
        all_rows.extend(main.cache.values())
        if mu!=.5:
            extra=AzimuthalSector(mu,h=.05,nmax=12);extra_row=extra.evaluate(.4,save_matrices=True)
            comparisons={}
            for branch in ('inherited','reflection'):
                one=rows[1] if branch=='reflection' else rows[1]['inherited']
                two=extra_row if branch=='reflection' else extra_row['inherited']
                comparisons[branch]={field:(np.asarray(two['actual_two_chirality_sum'][field])-one['actual_two_chirality_sum'][field]).tolist()
                    for field in one['actual_two_chirality_sum']}
            angular.append(dict(mu=mu,config=extra.config,sample=extra_row,change_N12_minus_N8=comparisons))
            all_rows.extend(extra.cache.values());del extra
        del main
        fine_mesh=AzimuthalSector(mu,h=.035)
        rows=[fine_mesh.evaluate(k,save_matrices=True) for k in (.01,.4,2.)]
        samples.append(dict(config=fine_mesh.config,assembly=fine_mesh.assembly,ports=fine_mesh.port_checks,samples=rows))
        stamp=time.perf_counter();fm=integrate_both(fine_mesh,8)
        if mu==.5:legacy.append(dict(record='h035_order8',**compare_saved(fm['inherited'],saved['integrated_mesh_refinement'][1])))
        print(json.dumps(dict(mu=mu,h=.035,order=8,totals={b:fm[b]['totals']['paired'] for b in fm},seconds=time.perf_counter()-stamp)),flush=True)
        integrations.append(dict(mu=mu,coarse=coarse,fine=fine,primary=primary,finer_mesh=fm,
            quadrature=dict(target=quadrature_target,initial_observed_panel_change_sum=initial,final_observed_panel_change_sum=final,
                initial_gate_passed=all(max(float(np.max(np.abs(v))) for v in values.values())<=quadrature_target for values in initial.values()),
                final_gate_passed=all(max(float(np.max(np.abs(v))) for v in values.values())<=quadrature_target for values in final.values()))))
        all_rows.extend(fine_mesh.cache.values());del fine_mesh
    sums={key:sum_sectors(integrations,key) for key in ('coarse','fine','primary','finer_mesh')}
    combined_changes={branch:{field:np.sum([r['quadrature']['final_observed_panel_change_sum'][branch][field] for r in integrations],axis=0).tolist()
        for field in integrations[0]['quadrature']['final_observed_panel_change_sum'][branch]} for branch in ('inherited','reflection')}
    root_targets={.08:(-.37788779014228824,-.1250938640349057),.05:(-.3803638300058992,-.12601622824121111),
                  .035:(-.38106699916232195,-.12627878704961049)}
    reflected_reference=[]
    for record in samples:
        if record['config']['mu']==.5:
            row=record['samples'][1]['actual_two_chirality_sum'];h=record['config']['h'];e,p=root_targets[h]
            reflected_reference.append(dict(h=h,energy_error=abs(row['energy_integrand']-e),central_rho_plus_pr_error=abs(row['pressure_integrands'][0]-p)))
    every_observable=[c for row in all_rows for branch in (row,row['inherited']) for c in (branch['right'],branch['left_physical'])]
    pilot_rows=[row for record in samples for row in record['samples']]
    angular_change=max(float(np.max(np.abs(value))) for record in angular for branch in record['change_N12_minus_N8'].values() for value in branch.values())
    tests=dict(monopole_basis_operator_and_spin_lift=basis['passed'],zero_modes_exactly_counted=basis['zero_mode_counts']=={'-0.5':0,'0.5':1,'1.5':0},
        free_Bessel_PDE=max(r['exact_pde_maximum_error'] for r in free['rows'])<1e-12,
        new_axis_mesh_refinement=all(r['refinements'][1]['relative_error']<.6*r['refinements'][0]['relative_error'] for r in free['rows']),
        free_current_identity=max(p['ledger']['identity_error'] for r in free['rows'] for p in r['refinements'])<1e-12,
        independent_opposite_assembly=max(c['assembly']['opposite_adjoint_relative_error'] for c in samples)<1e-12,
        same_axis_and_positive_imaginary_mass=all(c['assembly']['axis_domain_equal'] and c['assembly']['minimum_cell_mass_eigenvalue']>0 for c in samples),
        swapped_physical_ports=max(max(p[k] for k in ('basis_reproduction_error','incoming_adjoint_error','outgoing_adjoint_error','factor_swap_error')) for c in samples for p in c['ports'])<1e-12,
        exterior_equation_residual=max(max(r['exterior_equation_relative_residuals']) for r in all_rows)<1e-11,
        deliberately_adjoint_paired_control=max(r['discrepancies'][key] for r in pilot_rows for key in ('deliberate_adjoint_response_error','deliberate_adjoint_energy_error','deliberate_adjoint_pressure_error','deliberate_adjoint_shift_error'))<2e-9,
        real_sewing_determinants=max(abs(p['imag']) for r in every_observable for p in r['determinant_phases'])<1e-10,
        sampled_sewings_nonsingular=min(min(r['minimum_singular_values']) for r in every_observable)>1e-7,
        old_sector_inherited_nodes_and_integrals=max(max(r['node_maximum_error'],r['observable_maximum_error']) for r in legacy)<2e-11,
        independently_reconstructed_reflection_targets=max(max(r['energy_error'],r['central_rho_plus_pr_error']) for r in reflected_reference)<1e-12,
        finite_metric_throat_derivative=max(v['richardson_error'] for r in fd for v in r['branches'].values())<2e-8,
        observed_finite_window_quadrature=max(float(np.max(np.abs(v))) for branch in combined_changes.values() for v in branch.values())<=quadrature_target,
        selected_frequency_angular_change=angular_change<=angular_target,
        frozen_dependencies_preserved=hashes=={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs})
    tests={k:bool(v) for k,v in tests.items()}
    return dict(passed=all(tests.values()),tests=tests,stage='Explicit three-sector source contrast for two specified endpoint connections.',
        parameters=dict(q=1,sectors=[-.5,.5,1.5],b=.2,d=1.,R=2.,nmax=8,L=40.,rho_c=2.,units='hbar=c=d=1',phase_difference='pi minus0'),
        branches=dict(inherited='Original same-label modal sewing.',reflection='theta_lower=pi-theta_upper, same phi, spin lift -exp(i theta sigma_phi); paired cross parity (-1)^(n+|mu-.5|), zero unchanged.'),
        basis_checks=basis,free_axis_checks=free,mesh_refinement=samples,sector_integrations=integrations,partial_sums=sums,
        finite_metric_checks=fd,angular_refinement=angular,legacy_reproduction=legacy,independent_reflection_reference=reflected_reference,
        probes=probe_geometry(.2*math.sqrt(5),2.),
        diagnostics=dict(quadrature_target=quadrature_target,combined_observed_panel_change_sum=combined_changes,
            selected_angular_absolute_target=angular_target,maximum_selected_angular_change=angular_change,
            minimum_sampled_sewing_singular_value=min(min(r['minimum_singular_values']) for r in every_observable),
            maximum_exterior_residual=max(max(r['exterior_equation_relative_residuals']) for r in all_rows),
            frequency_geometry_evaluations=len(all_rows),omitted_frequency_tail_bound=None,full_angular_remainder_bound=None,continuum_error_bound=None),
        observable_definitions=dict(energy='-1/pi integral Delta log|det(I-S_ext B)| d kappa per complex Weyl sector, then sum actual right and left solves.',
            rho_plus_pr='Compact weighted energy-plus-radial-pressure PHASE DIFFERENCE, not radial pressure alone.',
            shift='Minus the weighted lower orthonormal radial energy-flux PHASE DIFFERENCE.',
            null_plus='rho_plus_pr - 2*shift',null_minus='rho_plus_pr + 2*shift'),
        scope=['Both chiralities and all three azimuthal sectors are explicitly solved; no degeneracy multiplier is used.',
            'Every sector retains coupled throat n1..8 and the existing exterior mixing. These three sectors complete only the n1 multiplet; higher multiplets remain incomplete.',
            'The reflected product-frame connection is a concrete specified isolated-sphere gluing, not a derivation of the physical finite magnetic/gravitational junction.',
            'The inherited and reflected branches are evaluated at identical local geometry and phase choices; they must not be mixed in a partial sum.',
            'A mere simultaneous raw basis rephasing of exterior S and throat B leaves the determinant invariant. The reflected branch changes B while holding S fixed and is a different specified connection.',
            'Physical raw spherical maps are converted by independent facet Gram factors for both chiralities. The deliberately adjoint-paired finite control is separate.',
            'Observed mesh, selected angular cutoff and frequency-panel changes are reported without a rigorous continuum, full-angular or missing-tail bound.',
            'The finite dissipative approximation does not establish a self-adjoint continuum determinant. The phase contrast is conditional on the common-reference symmetric-sea prescription.',
            'No absolute negative source, NEC violation, self-consistent geometry, stable medium, traversability or time travel is established.'],
        sources=['https://arxiv.org/pdf/hep-th/0206054','https://arxiv.org/pdf/1807.04726','https://arxiv.org/pdf/1606.07834'],
        input_source_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.perf_counter()-start)


if __name__=='__main__':
    result=checks();path=Path(__file__).with_name('azimuthal_source_results.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(path.read_text(encoding='utf-8'))==result
    print(json.dumps({k:result[k] for k in ('passed','tests','seconds','source_sha256')},indent=2),flush=True)
    if not result['passed']:raise SystemExit('An azimuthal source control failed; full diagnostics retained.')
