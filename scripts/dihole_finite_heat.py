"""Finite-box, finite-angular heat diagonal in the actual curved exterior.

Galerkin matrices use <H phi_i,H phi_j>, never (P H P)^2. Dirichlet sine
functions belong to the first-order graph domain after zero extension. This
is a pilot for the Friedrichs square on that box, not a continuum certificate
or a self-adjoint Dirac wall. Refinement changes are measured sensitivities.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh
from mpmath.ctx_mp import MPContext

from balanced_dihole import parameters, symbolic_fields
from dihole_null_scaling import physical_coefficients
from dihole_optical_heat import optical_fields
from null_saddle_leads import AU, AV, SY
from project_paths import project_path


def sine_table(x,count,radius):
    indices=np.arange(1,count+1)
    frequency=indices*np.pi/(2*radius)
    phase=(np.asarray(x)[:,None]+radius)*frequency
    return np.sin(phase)/np.sqrt(radius),np.cos(phase)*frequency/np.sqrt(radius)


def form_matrices(count,radius,order,coefficient):
    """M(delta)=M0+delta*M1+delta^2*M2 from an actual first-order action."""
    nodes,weights=np.polynomial.legendre.leggauss(order)
    nodes,weights=radius*nodes,radius*weights
    basis,derivative=sine_table(nodes,count,radius)
    values=np.einsum('pi,qj->pqij',basis,basis).reshape(order*order,count*count)
    du=np.einsum('pi,qj->pqij',derivative,basis).reshape(order*order,count*count)
    dv=np.einsum('pi,qj->pqij',basis,derivative).reshape(order*order,count*count)
    u,v=np.meshgrid(nodes,nodes,indexing='ij')
    w,cu,cv,mass,slope=coefficient(u.ravel(),v.ravel())
    def block(scalar,spin):
        return np.einsum('pk,ab->pakb',scalar,spin).reshape(2*order*order,2*count*count)
    K0=(-block(w[:,None]*du+cu[:,None]*values,AU.real)
        -block(w[:,None]*dv+cv[:,None]*values,AV.real)+block(mass[:,None]*values,SY.imag))
    K1=block(slope[:,None]*values,SY.imag)
    weight=np.outer(weights,weights).ravel()
    root=np.repeat(np.sqrt(weight),2)
    X,Y=root[:,None]*K0,root[:,None]*K1
    cross=X.T@Y
    matrices=(X.T@X,cross+cross.T,Y.T@Y)
    gram=values.T@(weight[:,None]*values)
    orthogonal_error=float(np.max(abs(gram-np.eye(count*count))))
    assert orthogonal_error<2e-12,orthogonal_error
    assert all(np.max(abs(M-M.T))<1e-10 for M in matrices)
    return matrices,dict(orthogonality_error=orthogonal_error,quadrature_order=order)


def point_evaluation(points,count,radius):
    result=[]
    for u,v in points:
        a=sine_table([u],count,radius)[0][0]
        b=sine_table([v],count,radius)[0][0]
        result.append(np.kron(np.outer(a,b).ravel()[None,:],np.eye(2)))
    return np.array(result)


def heat_sector(matrices,delta,evaluation,times):
    M=matrices[0]+delta*matrices[1]+delta*delta*matrices[2]
    cutoff=40/min(times)
    values,vectors=eigh(M,subset_by_value=(-np.inf,cutoff),driver='evr',check_finite=False)
    assert len(values)>0 and values[0]>-1e-8,(delta,values[:2])
    fields=evaluation@vectors
    amplitudes=np.sum(fields*fields,axis=1)
    heat=np.exp(-np.outer(times,values))@amplitudes.T
    # Completeness gives a conservative finite-matrix omitted-eigenvector bound.
    omitted=np.outer(np.exp(-np.asarray(times)*cutoff),np.sum(evaluation*evaluation,axis=(1,2)))
    return heat,omitted,float(values[0]),len(values)


def controls():
    n,A,order=6,3.,30
    w,slope,offset=.8,.3,.37
    constant=lambda u,v:(np.full_like(u,w),np.zeros_like(u),np.zeros_like(u),
                          np.full_like(u,-slope*offset),np.full_like(u,slope))
    matrices,quadrature=form_matrices(n,A,order,constant)
    points=[(0.,0.),(.4,-.2)];evaluation=point_evaluation(points,n,A)
    times=np.array([.5,1.,2.]); errors=[]
    kinetic=w*w*(np.pi/(2*A))**2*(np.arange(1,n+1)[:,None]**2+np.arange(1,n+1)[None,:]**2)
    for delta in (-3.5,.5,4.5):
        direct,omitted,_,_=heat_sector(matrices,delta,evaluation,times)
        M=matrices[0]+delta*matrices[1]+delta*delta*matrices[2]
        truth=np.repeat((kinetic+slope*slope*(delta-offset)**2).ravel(),2)
        assert np.max(abs(M-np.diag(truth)))<2e-12
        exact=np.empty_like(direct)
        for j,(u,v) in enumerate(points):
            bu=sine_table([u],n,A)[0][0]**2;bv=sine_table([v],n,A)[0][0]**2
            for i,tau in enumerate(times):
                exact[i,j]=2*np.sum(np.outer(bu,bv)*np.exp(-tau*(kinetic+slope*slope*(delta-offset)**2)))
        error=float(np.max(abs(direct-exact)))
        assert error<3e-13+float(np.max(omitted))
        errors.append(error)
    mp=MPContext();mp.dps=65
    z,a=mp.mpf('.7'),mp.mpf('.37')
    direct=mp.fsum(mp.exp(-z*(mp.mpf(j)+mp.mpf('.5')-a)**2) for j in range(-50,51))
    poisson=mp.sqrt(mp.pi/z)*(1+2*mp.fsum((-1)**k*mp.cos(2*mp.pi*k*a)*mp.exp(-mp.pi**2*k*k/z) for k in range(1,21)))
    assert abs(direct-poisson)<mp.mpf('1e-60')
    # Full angular/spin/coordinate normalization in a constant optical cylinder.
    alpha,beta,L,E=1.7,.9,.13,2.1
    speed,cslope=1/(alpha*L*E),1/(beta*E)
    tau=.8;s=tau/(E*E)
    plane=1/(4*np.pi*tau*speed*speed)
    angular_continuum=np.sqrt(np.pi/(tau*cslope*cslope))
    diagonal=4*plane*angular_continuum/(2*np.pi*L*L*alpha*alpha*beta)
    expected=4/(4*np.pi*s)**1.5
    assert abs(diagonal/expected-1)<1e-14
    return dict(passed=True,free_sine_matrix_max_error=2e-12,free_heat_errors=errors,
        noninteger_holonomy_poisson_error=mp.nstr(abs(direct-poisson),10),
        flat_volume_normalization_error=abs(diagonal/expected-1),quadrature=quadrature)


def physical_setup():
    saved=json.loads(project_path('dihole_probe_certificate_results.json').read_text())
    assert hashlib.sha256(project_path('dihole_probe_certificate_results.json').read_bytes()).hexdigest()=='d3c3e55a125c539baa8e9b5c41f97429e2759850a419c23280b70522e26f2843'
    constants=saved['target']['coordinate_constants']
    L,rs,vs=(float.fromhex(constants[key]) for key in ('L','rho_star','vstar'))
    E=vs/L;q=512000
    mp=MPContext();mp.dps=75;exact=optical_fields(mp)
    exact_rs=mp.mpf(rs.as_integer_ratio()[0])/rs.as_integer_ratio()[1]
    center_gauge=exact(exact_rs,mp.mpf(0))[2]
    mu0=math.floor(q*center_gauge)+.5
    p=parameters();function=symbolic_fields()[0]
    def coefficients(u,v):
        rho,z=rs+L*(u+v)/np.sqrt(2),L*(u-v)/np.sqrt(2)
        if np.any(rho<=0):raise ValueError('Box reaches the physical axis.')
        got=physical_coefficients(function,p,rho,z)
        w=got['v']/vs
        cu=L*(got['v_rho']+got['v_z'])/(2*np.sqrt(2)*vs)
        cv=L*(got['v_rho']-got['v_z'])/(2*np.sqrt(2)*vs)
        return w,cu,cv,got['b']*(mu0-q*got['a'])/E,got['b']/E
    points=[(0.,0.),(-1.,-1.),(1.,1.)]
    point_data=[];errors=[]
    for u,v in points+[(2.,-1.)]:
        rh=exact_rs+mp.mpf(L)*(mp.mpf(u)+v)/mp.sqrt(2);zz=mp.mpf(L)*(mp.mpf(u)-v)/mp.sqrt(2)
        aa,bb,gauge=exact(rh,zz)
        vfn=lambda x,y:mp.exp(-exact(x,y)[0])
        vr,vz=mp.diff(vfn,(rh,zz),(1,0)),mp.diff(vfn,(rh,zz),(0,1))
        truth=[mp.exp(-aa)/mp.mpf(vs),mp.mpf(L)*(vr+vz)/(2*mp.sqrt(2)*mp.mpf(vs)),
               mp.mpf(L)*(vr-vz)/(2*mp.sqrt(2)*mp.mpf(vs)),mp.exp(-bb)*(mp.mpf(mu0)-q*gauge)/(mp.mpf(vs)/mp.mpf(L)),
               mp.exp(-bb)/(mp.mpf(vs)/mp.mpf(L))]
        got=coefficients(np.array([u]),np.array([v]))
        errors.append(max(abs(float(raw[0])-float(value))/(1+abs(float(value))) for raw,value in zip(got,truth)))
        if (u,v) in points:
            point_data.append(dict(uv=[u,v],rho=float(rh),z=float(zz),optical_volume_density=float(mp.exp(2*aa+bb))))
    assert max(errors)<2e-9,errors
    return coefficients,dict(q=q,mu_center=mu0,center_gauge=float(center_gauge),length=L,energy=E,
        rho_center=rs,points=point_data,coefficient_control_max_normalized_error=max(errors),
        coefficients='Floating evaluation of exact-balanced formulas; checked against75-digit differentiation. No interval enclosure claimed.')


def calculate_case(count,radius,order,mode_cut,coefficient,target,times):
    started=time.perf_counter()
    matrices,quad=form_matrices(count,radius,order,coefficient)
    refined,_=form_matrices(count,radius,order+12,coefficient)
    changes=[float(np.linalg.norm(a-b,ord=2)) for a,b in zip(matrices,refined)]
    matrix_change=changes[0]+mode_cut*changes[1]+mode_cut**2*changes[2]
    assert matrix_change<1e-6,matrix_change
    points=[row['uv'] for row in target['points']]
    evaluation=point_evaluation(points,count,radius)
    heat=[];omitted=[];minimum=[];retained=[]
    for delta in range(-mode_cut,mode_cut+1):
        row,bound,lowest,number=heat_sector(matrices,delta,evaluation,times)
        heat.append(row);omitted.append(bound);minimum.append(lowest);retained.append(number)
        if (delta+mode_cut)%128==0:
            print(json.dumps(dict(stage='angular_sum',count=count,radius=radius,completed=delta+mode_cut+1,
                total=2*mode_cut+1,seconds=time.perf_counter()-started)),flush=True)
    heat=np.array(heat);omitted=np.array(omitted)
    E,L=target['energy'],target['length']
    normalization=np.array([2/(2*np.pi*L*L*p['optical_volume_density']) for p in target['points']])
    scale=(4*np.pi*np.asarray(times)/E**2)**1.5
    totals=[]
    for cutoff in sorted(set((mode_cut//3,2*mode_cut//3,mode_cut))):
        selected=abs(np.arange(-mode_cut,mode_cut+1))<=cutoff
        total=np.sum(heat[selected],axis=0)*normalization
        bound=np.sum(omitted[selected],axis=0)*normalization
        totals.append(dict(mode_cut=cutoff,angular_sectors=int(np.sum(selected)),heat_diagonal=total.tolist(),
            normalized_heat=(scale[:,None]*total).tolist(),
            finite_matrix_eigen_omission_bound=(scale[:,None]*bound).tolist()))
    assert min(minimum)>0,min(minimum)
    return dict(count=count,radius=radius,quadrature=quad,quadrature_refinement_order=order+12,
        component_matrix_changes=changes,all_retained_mu_matrix_change_estimate=matrix_change,
        mode_cut=mode_cut,real_matrix_dimension=2*count*count,
        minimum_retained_squared_eigenvalue=min(minimum),maximum_retained_eigenpairs=max(retained),
        totals=totals,sector_heat=heat.tolist(),seconds=time.perf_counter()-started,
        scope='Dirichlet/Friedrichs H^2 on this compact box. Quadrature changes and refinements are numerical sensitivities, not continuum or exterior error bounds.')


def run():
    started=time.perf_counter();checks=controls();coefficient,target=physical_setup()
    times=[.5,1.,2.]
    cases=[]
    for count,radius,order in ((10,4.,36),(14,4.,48),(18,4.,60),(22,5.,72)):
        cases.append(calculate_case(count,radius,order,384,coefficient,target,times))
    names=['balanced_dihole.py','dihole_null_scaling.py','dihole_return_geometry.py','dihole_optical_heat.py',
           'null_saddle_leads.py','dihole_probe_certificate_results.json','dihole_optical_heat_results.json','project_paths.py']
    result=dict(passed=True,controls=checks,target=target,scaled_proper_times=times,
        physical_proper_times=[tau/target['energy']**2 for tau in times],cases=cases,
        normalization='Fullrank4: Kopt=2/(2pi L^2 Jopt)*sum_mu tr2 Kq_mu, tau=s E^2. Jopt=alpha^2 beta.',
        form='M(delta)=<K0 phi,K0 phi>+delta*(<K0 phi,K1 phi>+transpose)+delta^2*<K1 phi,K1 phi>. Sines normalized analytically.',
        scope=['Finite-box and finite-basis heat diagonals, summed over all half-integer sectors in the displayed finite interval.',
            'Natural uncut exterior heat kernel and infinite angular sum are not certified by these refinement changes.',
            'Off-axis centered boxes cannot exhaust the whole exterior; boundary enlargement is only a sensitivity check.',
            'Finite-eigenvalue omission bound concerns the numerical matrix only; all other approximation errors remain separate.',
            'No proper-time energy integral, renormalized supporting stress, stable throat, material or backward passage.'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_finite_heat_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,seconds=result['runtime_seconds'],cases=[dict(count=row['count'],radius=row['radius'],totals=row['totals']) for row in cases]),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,5.2))
    tau=np.array(result['scaled_proper_times'])
    for row in result['cases']:
        values=np.array(row['totals'][-1]['normalized_heat'])
        left.plot(tau,values[:,0],'o-',label=f"N={row['count']}, A={row['radius']:g}")
    left.axhline(4,color='#888888',ls='--',label='Flat local leading term')
    left.set(xlabel='Scaled proper time',ylabel=r'$(4\pi s)^{3/2}\,K(s;x,x)$',title='Finite heat response at the magnetic null')
    left.legend(fontsize=9);left.grid(alpha=.15)
    latest=result['cases'][-1]
    for j,label in enumerate(('At null','c = -1','c = 1')):
        right.plot(tau,np.array(latest['totals'][-1]['normalized_heat'])[:,j],'o-',label=label)
    right.set(xlabel='Scaled proper time',ylabel='Normalized heat diagonal',title='Three locations in the same curved geometry')
    right.legend(fontsize=9);right.grid(alpha=.15)
    fig.suptitle('Finite angular sums toward the nonlocal quantum contribution',fontsize=14)
    fig.text(.5,.018,'Numerical finite-box pilot. Refinement changes are not error bounds; no renormalized stress or traversable medium is established.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.93))
    for extension in ('png','pdf'):
        fig.savefig(project_path('dihole_finite_heat_progress.'+extension),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
