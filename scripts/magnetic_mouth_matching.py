"""Finite source-free test-Maxwell matching for the prescribed sharp core.

Two flat exterior spheres, q=1, b=.2, d=1, are coupled by the static Maxwell
Dirichlet-to-Neumann map of the finite AdS2 x S2 core at rho_c=32. Decaying
zero-charge exterior multipoles repair the point-pair field's mismatch.
This does not solve a smooth collar, Einstein equations, quantum support or
stability. Finite multipole/ODE/quadrature calculations are not enclosures.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import eval_jacobi, eval_legendre
import sympy as sp


def angular(mu,nmax):
    """Y_l=(1-mu^2)P_l', derivative, P_l', and exact weighted norms."""
    mu=np.asarray(mu)
    degree=np.arange(1,nmax+1)
    pp=np.column_stack([(l+1)/2*eval_jacobi(l-1,1,1,mu) for l in degree])
    yy=(1-mu*mu)[:,None]*pp
    yp=np.column_stack([-l*(l+1)*eval_legendre(l,mu) for l in degree])
    gram=2*degree*(degree+1)/(2*degree+1)
    return yy,yp,pp,gram


def core_modes(nmax,b=.2,rho_c=32.,saved_bases=None):
    """Unit center even/odd ODEs; off-diagonal uses the exact unit Wronskian."""
    cut=math.atan(rho_c);ell=b*math.hypot(1.,rho_c)
    old={} if saved_bases is None else {row['degree']:row for row in saved_bases}
    x=sp.Symbol('x',real=True)
    rows=[]
    for degree in range(1,nmax+1):
        def rhs(at,state):
            e,ep,o,op=state
            potential=degree*(degree+1)/math.cos(at)**2
            return [ep,potential*e,op,potential*o]
        solution=solve_ivp(rhs,(0.,cut),(1.,0.,0.,1.),method='DOP853',
            rtol=2e-12,atol=2e-13,max_step=cut/128,dense_output=True)
        if not solution.success:raise RuntimeError(solution.message)
        e,ep,o,op=map(float,solution.y[:,-1])
        diagonal=(ep/e+op/o)/2
        off=-1/(2*e*o)
        row=dict(degree=degree,e_end=e,e_prime_end=ep,o_end=o,o_prime_end=op,
            physical_diagonal=diagonal/ell,physical_offdiagonal=off/ell,
            even_eigenvalue=(diagonal+off)/ell,odd_eigenvalue=(diagonal-off)/ell,
            lower_reflection_parity=(-1)**(degree+1),ode_evaluations=solution.nfev)
        if degree in old:
            even=sp.sympify(old[degree]['even'],locals={'x':x})
            odd=sp.sympify(old[degree]['odd'],locals={'x':x})
            exact=np.asarray(sp.lambdify(x,(even,sp.diff(even,x),odd,sp.diff(odd,x)),'numpy')(cut),float)
            row['saved_exact_basis_relative_error']=float(np.max(abs(solution.y[:,-1]/exact-1)))
            exact_diag=(exact[1]/exact[0]+exact[3]/exact[2])/2/ell
            exact_off=-1/(2*exact[0]*exact[2]*ell)
            row['saved_exact_dtn_absolute_error']=max(abs(exact_diag-diagonal/ell),abs(exact_off-off/ell))
        rows.append(row)
    return rows


def core_matrix(rows,nmax):
    result=np.zeros((2*nmax,2*nmax))
    for j,row in enumerate(rows[:nmax]):
        result[j,j]=result[nmax+j,nmax+j]=row['physical_diagonal']
        result[j,nmax+j]=result[nmax+j,j]=row['lower_reflection_parity']*row['physical_offdiagonal']
    return result


def incident_coefficients(nmax,q=1.,b=.2,d=1.):
    l=np.arange(1,nmax+1)
    upper=q/2*(-1.)**l*(b/d)**(l+1)/(l+1)
    lower=(-1.)**(l+1)*upper
    values=np.r_[upper,lower]
    return values,np.tile((l+1)/b,2)*values


def sphere_coordinates(mu,r,observer,source,d=1.):
    """Coordinates about source on the observer's varying-radius ray."""
    displacement=(source-observer)*d # indices upper=0, lower=1
    transverse=r*np.sqrt(np.maximum(0.,1-mu*mu))
    axial=displacement+r*mu
    radius=np.hypot(transverse,axial)
    cosine=axial/radius
    radial_derivative=(r+displacement*mu)/radius
    cosine_derivative=(mu-cosine*radial_derivative)/radius
    return radius,cosine,radial_derivative,cosine_derivative


def incident_field(mu,r,observer,q=1.,d=1.):
    _,cosine,_,derivative=sphere_coordinates(np.asarray(mu),r,observer,1-observer,d)
    sign=1 if observer==0 else -1
    return q/2*(sign*cosine-1),q/2*sign*derivative


def exterior_columns(mu,r,observer,nmax,b=.2,d=1.):
    """Value and outward-from-core radial derivative of every induced basis."""
    values=[];derivatives=[]
    degree=np.arange(1,nmax+1)
    for source in (0,1):
        radius,cosine,rp,cp=sphere_coordinates(np.asarray(mu),r,observer,source,d)
        yy,yp,_,_=angular(cosine,nmax)
        power=(b/radius[:,None])**degree
        values.append(power*yy)
        derivatives.append(power*(-degree*(rp/radius)[:,None]*yy+cp[:,None]*yp))
    return np.column_stack(values),np.column_stack(derivatives)


def projected_matrices(nmax,monitor,order,b=.2,d=1.,q=1.):
    mu,weights=np.polynomial.legendre.leggauss(order)
    _,_,pp,gram=angular(mu,monitor)
    projection=pp.T*weights[None,:]/gram[:,None]
    values=[];derivatives=[];incident=[];incident_derivatives=[]
    for observer in (0,1):
        vv,dd=exterior_columns(mu,b,observer,nmax,b,d)
        iv,idv=incident_field(mu,b,observer,q=q,d=d)
        values.append(projection@vv);derivatives.append(projection@dd)
        incident.append(projection@iv);incident_derivatives.append(projection@idv)
    return np.vstack(values),np.vstack(derivatives),np.concatenate(incident),np.concatenate(incident_derivatives)


def weighted_norm(coefficients,nmax):
    l=np.arange(1,nmax+1);gram=np.tile(2*l*(l+1)/(2*l+1),2)
    return float(np.sqrt(np.sum(gram*np.asarray(coefficients)**2)))


def solve_matching(nmax,order,core,b=.2,d=1.,q=1.):
    monitor=nmax+6
    evalue,eder,iprojected,dprojected=projected_matrices(nmax,monitor,order,b,d,q)
    inc,incder=incident_coefficients(monitor,q,b,d)
    k=core_matrix(core,monitor)
    low=np.r_[np.arange(nmax),monitor+np.arange(nmax)]
    high=np.setdiff1d(np.arange(2*monitor),low)
    operator=eder-k@evalue
    bare=incder-k@inc
    beta=np.linalg.solve(operator[low],-bare[low])
    boundary=inc+evalue@beta
    derivative=incder+eder@beta
    residual=derivative-k@boundary
    parity=(-1.)**(np.arange(1,nmax+1)+1)
    symmetry=float(np.max(abs(beta[nmax:]-parity*beta[:nmax])))
    self_value_error=max(float(np.max(abs(evalue[np.ix_(np.arange(nmax),np.arange(nmax))]-np.eye(nmax)))),
        float(np.max(abs(evalue[np.ix_(monitor+np.arange(nmax),nmax+np.arange(nmax))]-np.eye(nmax)))))
    self_derivative_error=max(float(np.max(abs(eder[np.ix_(np.arange(nmax),np.arange(nmax))]+np.diag(np.arange(1,nmax+1)/b)))),
        float(np.max(abs(eder[np.ix_(monitor+np.arange(nmax),nmax+np.arange(nmax))]+np.diag(np.arange(1,nmax+1)/b)))))
    modes=[]
    for j,row in enumerate(core[:monitor]):
        l=j+1;parity_l=(-1)**(l+1)
        left,right=parity_l*boundary[monitor+j],boundary[j]
        pair=np.array([left,right])
        dcore=np.array([[row['physical_diagonal'],row['physical_offdiagonal']],
                        [row['physical_offdiagonal'],row['physical_diagonal']]])
        quadratic=float(pair@dcore@pair)
        gram=2*l*(l+1)/(2*l+1)
        modes.append(dict(degree=l,core_left_value=float(left),core_right_value=float(right),
            center_value=float((left+right)/(2*row['e_end'])),
            center_chi_derivative=float((right-left)/(2*row['o_end'])),
            boundary_canonical_quadratic_form=quadratic,
            gauge_coupling_squared_core_energy=math.pi*gram*quadratic))
    full_norm=weighted_norm(residual,monitor)
    high_only=np.zeros_like(residual);high_only[high]=residual[high]
    return dict(induced_cutoff=nmax,monitor_degree=monitor,quadrature_order=order,
        beta_upper=beta[:nmax].tolist(),beta_lower=beta[nmax:].tolist(),
        incident_value_coefficients=inc.tolist(),incident_derivative_coefficients=incder.tolist(),
        repaired_value_coefficients=boundary.tolist(),repaired_derivative_coefficients=derivative.tolist(),
        bare_derivative_mismatch_coefficients=bare.tolist(),repaired_derivative_mismatch_coefficients=residual.tolist(),
        bare_mismatch_weighted_norm=weighted_norm(bare,monitor),
        retained_mismatch_maximum=float(np.max(abs(residual[low]))),
        higher_mismatch_maximum=float(np.max(abs(residual[high]))),
        higher_mismatch_weighted_norm=weighted_norm(high_only,monitor),
        monitored_total_mismatch_weighted_norm=full_norm,
        monitored_sheet_current_g2_L2=math.sqrt(2*math.pi)*full_norm,
        retained_matching_matrix_condition=float(np.linalg.cond(operator[low])),
        incident_value_projection_error=float(np.max(abs(iprojected-inc))),
        incident_derivative_projection_error=float(np.max(abs(dprojected-incder))),
        self_value_projection_error=self_value_error,self_derivative_projection_error=self_derivative_error,
        reflection_symmetry_error=symmetry,core_modes=modes,
        retained_core_energy_g2=sum(m['gauge_coupling_squared_core_energy'] for m in modes[:nmax]),
        monitored_core_energy_g2=sum(m['gauge_coupling_squared_core_energy'] for m in modes))


def local_controls(nmax=9,b=.2,d=1.):
    mu=np.array([-.73,.11,.68]);step=2e-6
    error=incident_error=charge_error=0.
    for observer in (0,1):
        analytic=exterior_columns(mu,b,observer,nmax,b,d)[1]
        values=[exterior_columns(mu,b+j*step,observer,nmax,b,d)[0] for j in (-2,-1,1,2)]
        finite=(values[0]-8*values[1]+8*values[2]-values[3])/(12*step)
        error=max(error,float(np.max(abs(finite-analytic))))
        exact_incident=incident_field(mu,b,observer,d=d)[1]
        vv=[incident_field(mu,b+j*step,observer,d=d)[0] for j in (-2,-1,1,2)]
        fd=(vv[0]-8*vv[1]+8*vv[2]-vv[3])/(12*step)
        incident_error=max(incident_error,float(np.max(abs(fd-exact_incident))))
        poles=exterior_columns(np.array([-1.,1.]),b,observer,nmax,b,d)[0]
        charge_error=max(charge_error,float(np.max(abs(2*math.pi*(poles[0]-poles[1])))))
    # Independent numerical angular identity in the source-free flat Maxwell ODE.
    theta=np.linspace(-.87,.91,43);h=1e-4
    yy,yp,_,_=angular(theta,nmax)
    forward=angular(theta+h,nmax)[1];backward=angular(theta-h,nmax)[1]
    second=(forward-backward)/(2*h)
    degree=np.arange(1,nmax+1)
    angular_residual=(1-theta*theta)[:,None]*second+degree*(degree+1)*yy
    reflected=angular(-theta,nmax)[0]
    reflection_error=float(np.max(abs(reflected-yy*((-1.)**(degree+1)))))
    cut=math.atan(32.);ell=b*math.hypot(1.,32.)
    proper_cut_factor=b/math.cos(cut)
    return dict(induced_normal_derivative_five_point_error=error,
        incident_normal_derivative_five_point_error=incident_error,
        extra_magnetic_flux_at_either_mouth=charge_error,
        angular_Maxwell_identity_finite_difference_relative_error=float(np.max(abs(angular_residual))/np.max(abs(degree*(degree+1)*yy))),
        lower_angular_reflection_error=reflection_error,
        core_cut_lapse=proper_cut_factor/ell,
        cut_proper_to_chi_derivative=1/proper_cut_factor,
        physical_dtn_clock_factor=1/ell,
        exact_flat_radial_Maxwell_identity='(r^(-l))_rr-l(l+1)r^(-l-2)=0; all included multipoles have l>=1.',
        charge_identity='Every induced Aphi vanishes on both axis endpoints; 2*pi*[Aphi(mu=-1)-Aphi(mu=+1)]=0.')


def run():
    inputs=('mmp_pulse.py','mmp_pulse_results.json','exterior_boundary_map.py')
    hashes={p:hashlib.sha256(project_path(p).read_bytes()).hexdigest() for p in inputs}
    saved=json.loads(project_path(inputs[1]).read_text(encoding='utf-8'))['magnetic_boundary_response']
    core=core_modes(15,saved_bases=saved['normalized_bases'])
    records=[solve_matching(n,order,core) for n in (3,6,9) for order in (32,64)]
    local=local_controls()
    quadrature=[];cutoff=[]
    for n in (3,6,9):
        a,b=[r for r in records if r['induced_cutoff']==n]
        quadrature.append(dict(induced_cutoff=n,
            maximum_beta_difference=float(np.max(abs(np.array(a['beta_upper'])-b['beta_upper']))),
            monitored_residual_norm_change=abs(a['monitored_total_mismatch_weighted_norm']-b['monitored_total_mismatch_weighted_norm']),
            retained_core_energy_change=abs(a['retained_core_energy_g2']-b['retained_core_energy_g2'])))
    fine=[r for r in records if r['quadrature_order']==64]
    for a,b in zip(fine[:-1],fine[1:]):
        n=a['induced_cutoff']
        cutoff.append(dict(lower_cutoff=n,upper_cutoff=b['induced_cutoff'],
            maximum_shared_beta_change=float(np.max(abs(np.array(a['beta_upper'])-b['beta_upper'][:n]))),
            leading_beta_change=b['beta_upper'][0]-a['beta_upper'][0],
            monitored_residual_reduction=a['monitored_total_mismatch_weighted_norm']/b['monitored_total_mismatch_weighted_norm']))
    tests=dict(
        saved_normalized_bases=max(r.get('saved_exact_basis_relative_error',0.) for r in core)<2e-10,
        positive_core_dtn=min(min(r['even_eigenvalue'],r['odd_eigenvalue']) for r in core)>0,
        incident_values_and_derivatives=max(max(r['incident_value_projection_error'],r['incident_derivative_projection_error']) for r in records)<2e-12,
        exact_self_projection=max(max(r['self_value_projection_error'],r['self_derivative_projection_error']) for r in records)<2e-11,
        retained_matching=max(r['retained_mismatch_maximum'] for r in records)<3e-14,
        reflected_two_mouth_symmetry=max(r['reflection_symmetry_error'] for r in records)<3e-14,
        physical_clock_normal_and_reflection=abs(local['core_cut_lapse']-1)<3e-14 and
            abs(local['cut_proper_to_chi_derivative']-local['physical_dtn_clock_factor'])<3e-14 and
            local['lower_angular_reflection_error']<3e-12,
        no_added_charge=local['extra_magnetic_flux_at_either_mouth']<2e-14,
        independent_field_derivatives=max(local['induced_normal_derivative_five_point_error'],local['incident_normal_derivative_five_point_error'])<3e-7,
        angular_Maxwell_control=local['angular_Maxwell_identity_finite_difference_relative_error']<3e-6,
        positive_core_energy=min(m['boundary_canonical_quadratic_form'] for r in records for m in r['core_modes'])>=0,
        quadrature_stability=max(r['maximum_beta_difference'] for r in quadrature)<2e-12,
        monitored_omitted_modes_decrease=all(r['monitored_residual_reduction']>100 for r in cutoff),
        inputs_preserved=all(hashlib.sha256(project_path(p).read_bytes()).hexdigest()==v for p,v in hashes.items()))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,
        parameters=dict(q=1,b=.2,d=1.,rho_c=32.,ell=.2*math.hypot(1.,32.),cut=math.atan(32.),
            cut_lapse=1.,induced_cutoffs=[3,6,9],monitor_extra_degrees=6,quadrature_orders=[32,64]),
        conventions=dict(geometry='Prescribed sharp AdS2 x S2 core with matched cut lapse one; flat test-field exterior outside two spheres.',
            field='Aphi=point_pair_Aphi+sum_(mouth,l) beta_mouth,l*(b/r_mouth)^l*(1-mu_mouth^2)*P_l_prime(mu_mouth).',
            ordering='Upper degrees1..M, then lower degrees1..M; every exterior radial derivative points away from its mouth center, from core into exterior.',
            core_orientation='chi=-xc is reflected lower mouth, chi=+xc is upper mouth. Lower core cosine=-global cosine; value and outward derivative get parity (-1)^(l+1).',
            physical_clock='Ncut=1 and dl/dchi=ell at either cut, hence physical core DtN=Dchi/ell; matching F_nphi then excludes an azimuthal sheet current in this prescribed metric.',
            projection='a_l=integral_-1^1 Aphi(mu)*P_l_prime(mu) dmu / G_l, G_l=2l(l+1)/(2l+1).',
            mismatch_norm='sqrt(sum_mouth,l G_l*(exterior_radial_derivative-core_DtN_value)_l^2). sqrt(2*pi) times this is g^2 times the surface-current L2 norm over both physical spheres, for monitored modes only.',
            energy='g^2 E_core,l=pi G_l * A_core^T (Dchi/ell) A_core. The background-monopole cross term integrates to zero. Saved sums are finite-mode values.',
            dimensions='d=1 length units; gauge coupling is left symbolic; no SI scale assigned.'),
        core_modes=core,local_controls=local,records=records,
        quadrature_sensitivity=quadrature,cutoff_sensitivity=cutoff,
        conclusions=dict(point_pair_alone_matches=False,finite_retained_magnetic_matching_solved=True,
            all_angular_matching_certified=False,smooth_Einstein_quantum_geometry_solved=False),
        scope=[
            'A bounded axisymmetric test-Maxwell magnetic matching calculation, not a smooth transition or a self-consistent supporting spacetime.',
            'The old point-pair field is source-free away from excluded charges but has a finite mouth derivative mismatch. Induced l>=1 multipoles preserve each magnetic flux.',
            'The core and flat exterior solve their source-free Maxwell equations in their respective pieces. Retained interface conditions are solved; residual through L+6 remains explicitly measured.',
            'The core means the full static Dirichlet extension of the actual exterior boundary trace, which contains arbitrarily high companion-induced degrees. Only its modes through L+6 are evaluated here; truncating the core field itself to L would also leave a value mismatch.',
            'Quadrature, angular-cutoff and ODE controls are floating-point consistency/sensitivity tests, not an infinite-mode or numerical error enclosure.',
            'The exact sphere charge is 2*pi*q at the upper mouth and -2*pi*q at the lower. No induced monopole, electric field or extra total charge is introduced.',
            'The sharp metric junction and its required stress remain. The corrected angular magnetic field changes the Dirac operator; previous quantum-source calculations are not automatically transferable.',
            'No material generation, topology creation, formation, stability or time-travel result follows. No Dirac, quantum or historical PDE calculation was rerun.'],
        input_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=saved['sources'])
    project_path('magnetic_mouth_matching_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary=dict(passed=result['passed'],tests=tests,local_controls=local,
        fine_records=[{k:r[k] for k in ('induced_cutoff','bare_mismatch_weighted_norm','retained_mismatch_maximum',
            'higher_mismatch_weighted_norm','monitored_total_mismatch_weighted_norm','beta_upper','retained_core_energy_g2')} for r in fine],
        quadrature_sensitivity=quadrature,cutoff_sensitivity=cutoff)
    print(json.dumps(summary,indent=2))
    if not result['passed']:raise AssertionError(tests)
    return result


if __name__=='__main__':
    run()
