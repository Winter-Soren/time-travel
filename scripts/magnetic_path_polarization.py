"""First local-log-action correction to fixed-flux, axisymmetric return paths.

The flat point-source exterior is a declared approximation. This is not a
finite-mouth matching solution or a complete nonlocal fermion determinant.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import solve_banded
from scipy.special import eval_legendre

from mmp_pulse import _magnetic_return_length


def flux_source(u,v):
    D=u*u-v*v; K=u*u*(1+3*v*v)+v**4-v*v
    return 12*u*(u*u-1)*(1-v*v)*(1+v*v)/(D*K)


def angular_basis(v,degrees):
    p=np.array([eval_legendre(int(l),v) for l in degrees])
    pm=np.array([eval_legendre(int(l-1),v) for l in degrees])
    Y=degrees[:,None]*(pm-v*p)
    Yv=-degrees[:,None]*(degrees[:,None]+1)*p
    return Y,Yv


def flux_solution(outer=20000.,radial=900,modes=24,angular=192,inner=120):
    """Separated angular flux harmonics; nonuniform radial differences in log u."""
    if not (outer>10 and radial>=100 and modes>=2 and angular>2*modes and inner>=30):
        raise ValueError('Insufficient outer radius or spectral/grid resolution')
    t=np.r_[0.,np.geomspace(1e-6,.1,inner),np.linspace(.1,math.log(outer),radial)[1:]]
    u=np.exp(t); degrees=np.arange(1,2*modes,2)
    nodes,weights=np.polynomial.legendre.leggauss(angular)
    Y,_=angular_basis(nodes,degrees)
    norm=2*degrees*(degrees+1)/(2*degrees+1)
    projections=flux_source(u[:,None],nodes)@(weights*Y/(1-nodes*nodes)).T/norm
    hm=t[1:-1]-t[:-2]; hp=t[2:]-t[1:-1]
    factor=-np.expm1(-2*t[1:-1])
    left=factor*(2/(hm*(hm+hp))+hp/(hm*(hm+hp)))
    diagonal=factor*(-2/(hm*hp)-(hp-hm)/(hm*hp))
    right=factor*(2/(hp*(hm+hp))-hm/(hp*(hm+hp)))
    solutions=[]; residuals=[]
    for j,l in enumerate(degrees):
        diag=diagonal-l*(l+1)
        bands=np.zeros((3,len(t)-2)); bands[1]=1.
        bands[0,1:]=right[:-1]/diag[:-1]
        bands[2,:-1]=left[1:]/diag[1:]
        interior=solve_banded((1,1),bands,projections[1:-1,j]/diag)
        row=np.r_[0.,interior,0.]
        residuals.append(float(np.max(abs(left*row[:-2]+diag*row[1:-1]+right*row[2:]-projections[1:-1,j]))))
        solutions.append(row)
    values=np.array(solutions).T
    spline=CubicSpline(t,values,axis=0)
    def evaluate(points,v):
        points,v=np.broadcast_arrays(points,v)
        if np.any(points<1) or np.any(points>outer): raise ValueError('Path leaves the solved prolate domain')
        basis,basis_v=angular_basis(v.ravel(),degrees)
        coefficients=spline(np.log(points.ravel())).T
        derivatives=spline(np.log(points.ravel()),1).T/points.ravel()
        return tuple(value.reshape(points.shape) for value in (
            np.sum(coefficients*basis,axis=0),np.sum(derivatives*basis,axis=0),
            np.sum(coefficients*basis_v,axis=0)))
    # Maximum principle and analytic barrier are checked away from focal corners.
    us=np.geomspace(1.0001,outer*.999,80); vs=np.linspace(-.995,.995,61)
    U,V=np.meshgrid(us,vs,indexing='ij'); psi,_,_=evaluate(U,V)
    psi0=2*U*(1-V*V)/(U*U-V*V); barrier=2*psi0*np.log(U)
    return evaluate,{'outer':outer,'radial_points':len(t),'modes':modes,'angular_order':angular,
        'inner_radial_points':inner,'psi1_at_u2_v0':float(evaluate(np.array([2.]),np.array([0.]))[0][0]),
        'max_discrete_residual':max(residuals),'maximum_positive_flux':float(np.max(psi)),
        'maximum_lower_barrier_violation':float(np.max(-barrier-psi)),
        'sample_u':us.tolist(),'sample_v':vs.tolist(),'sample_psi1':psi.tolist()}


def path_variation(nu,evaluate,order=192,outer=20000.):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    v=(nodes+1)/2; weights=weights/2
    h=1-v*v; root=np.hypot(h,nu*v); u=(h+root)/nu
    up=v*(nu-2*u)/root
    psi,psi_u,psi_v=evaluate(u,v)
    psi0_u=-2*h*(u*u+v*v)/(u*u-v*v)**2
    du=-psi/psi0_u
    log_derivative=-2*v/h+(2*u*up+2*v)/(u*u+v*v)-2*(2*u*up-2*v)/(u*u-v*v)
    dup=-(psi_u*up+psi_v)/psi0_u-du*log_derivative
    def arc(u,up):
        return .5*np.sqrt((u*u-v*v)*(up*up/(u*u-1)+1/h))
    base=arc(u,up); A=u*u-1; D=u*u-v*v; Q=up*up/A+1/h
    delta=base*(u/D-u*up*up/(A*A*Q))*du+base*up/(A*Q)*dup
    C=u*u*(1+3*v*v)-v*v+v**4
    kernel=3*u*(1+v*v)*(u*u-v*v)/(2*C**1.5)
    length=2*weights@base; direct=2*weights@delta
    variation=2*weights@(kernel*du)
    lower=-4*weights@(kernel*u*(u*u-v*v)/(u*u+v*v)*np.log(u))
    box_bound=4*outer*math.log(outer)/(outer*outer-1)**2*float(
        weights@(kernel*(u*u-1)*(u*u-v*v)**2/(u*u+v*v)))
    finite=[]
    for step in (.002,.001):
        finite.append(float(weights@(arc(u+step*du,up+step*dup)-arc(u-step*du,up-step*dup))/step))
    exact=_magnetic_return_length(float(nu))
    return {'nu':float(nu),'base_length':float(length),'delta_length_per_epsilon':float(variation),
        'relative_length_coefficient':float(variation/length),'base_length_error':abs(float(length)-exact),
        'direct_vs_positive_kernel':abs(float(direct-variation)),
        'global_lower_length_bound':float(lower),'outer_box_length_error_bound':box_bound,
        'finite_difference_derivatives':finite,
        'richardson_derivative_error':abs((4*finite[1]-finite[0])/3-float(variation)),
        'maximum_prolate_displacement_per_epsilon':float(np.max(abs(du))),
        'profile':{'v':v.tolist(),'u':u.tolist(),'delta_u':du.tolist()}}


def polarization_checks():
    prior=json.loads(project_path('dirac_return_closure_results.json').read_text())
    baseline=prior['records'][-1]['roots'][-1]
    eps=baseline['background']['parameters']['epsilon']; clock=baseline['background']['clock']
    distance=prior['parameters']['distance_over_re']; q=prior['parameters']['q']
    configs=[dict(outer=20000.,radial=900,modes=24,angular=192,inner=120),
             dict(outer=20000.,radial=1800,modes=36,angular=288,inner=240),
             dict(outer=40000.,radial=1927,modes=36,angular=288,inner=240)]
    records=[]; nus=np.array(baseline['nu']); weights=np.array(baseline['weights'])
    periods=np.array(baseline['periods']); frequency=(2*math.pi/periods)**2
    condition=prior['records'][-1]['conditioning_factor']
    for config in configs:
        solution,diag=flux_solution(**config)
        representative=[path_variation(nu,solution,256,config['outer']) for nu in (.02,.1,.5,1.,1.5,1.9)]
        paths=[path_variation(float(nu),solution,256,config['outer']) for nu in nus]
        df=np.array([p['delta_length_per_epsilon'] for p in paths])
        dperiod=eps*distance/clock*df
        contributions=-2*weights*frequency*dperiod/periods
        dF=float(np.sum(contributions)); dw=condition*dF
        box_dF=float(2*weights@(frequency/periods*eps*distance/clock*np.array([p['outer_box_length_error_bound'] for p in paths])))
        eta=6/np.sqrt(q*nus)*np.sqrt(1-nus*nus/4)
        records.append({'configuration':config,'field_diagnostics':diag,'representative_paths':representative,
            'state_variation_at_fixed_geometry':dF,'linear_closed_state_shift':dw,
            'formal_corrected_state':prior['corrected_root']+dw,
            'outer_box_state_shift_bound':condition*box_dF,
            'nu':nus.tolist(),'delta_length_per_epsilon':df.tolist(),
            'gradient_turning_screen':eta.tolist(),
            'formal_state_shift_by_channel':(condition*contributions).tolist(),
            'sampled_flux_failing_gradient_point1':float(weights[eta>.1].sum()),
            'sampled_formal_shift_from_failing_gradient_point1':float(condition*contributions[eta>.1].sum()),
            'sampled_baseline_state_from_failing_gradient_point1':float(weights[eta>.1]@frequency[eta>.1]),
            'max_base_length_relative_error':max(p['base_length_error']/p['base_length'] for p in paths),
            'max_kernel_disagreement_relative':max(p['direct_vs_positive_kernel']/(1+abs(p['delta_length_per_epsilon'])) for p in paths),
            'max_fd_disagreement_relative':max(p['richardson_derivative_error']/(1+abs(p['delta_length_per_epsilon'])) for p in paths)})
    coarse,fine,far=records
    representative_refinement=max(abs(a['relative_length_coefficient']-b['relative_length_coefficient']) for a,b in zip(coarse['representative_paths'],fine['representative_paths']))
    state_refinement=abs(coarse['linear_closed_state_shift']-fine['linear_closed_state_shift'])
    outer_change=abs(fine['linear_closed_state_shift']-far['linear_closed_state_shift'])
    tests={
        'flux_sign_and_barrier':all(r['field_diagnostics']['maximum_positive_flux']<=1e-9 and r['field_diagnostics']['maximum_lower_barrier_violation']<=1e-7 for r in records),
        'discrete_field_equations':max(r['field_diagnostics']['max_discrete_residual'] for r in records)<2e-8,
        'negative_length_variation':all(max(r['delta_length_per_epsilon'])<0 for r in records),
        'base_path_control':max(r['max_base_length_relative_error'] for r in records)<2e-7,
        'independent_length_forms':max(r['max_kernel_disagreement_relative'] for r in records)<2e-5,
        'length_finite_difference':max(r['max_fd_disagreement_relative'] for r in records)<2e-6,
        'grid_and_mode_refinement':state_refinement<2e-9 and representative_refinement<2e-3,
        'outer_domain_refinement':outer_change<2e-9,
    }
    return {'passed':all(tests.values()),'tests':tests,'epsilon':eps,'q':q,
        'baseline_state':prior['corrected_root'],'conditioning_factor':condition,
        'records':records,'representative_coefficient_refinement':representative_refinement,
        'state_shift_refinement':state_refinement,'outer_domain_state_shift_change':outer_change,
        'equations':{'field':'(u²−1)psi1_uu+(1−v²)psi1_vv=S; S=12u(u²−1)(1−v²)(1+v²)/[(u²−v²)(u²(1+3v²)+v4−v²)].',
            'boundary':'Zero flux correction on axes, segment and regular point foci; no added singular point multipoles. Decaying infinite solution approximated by zero outer boxes with a barrier bound.',
            'length':'delta f=2 integral K delta u dv, K>0; delta u=-psi1/psi0_u<=0. The fixed-point-mouth endpoint term vanishes.',
            'state':'delta w=[-2<Omega_nu² delta L_nu/L_nu>]/[1-Fprime]; delta L=epsilon*d*delta f/A. Only first-order path correction; no new full root solve.'},
        'scope':'Formal first-order local-log-action, flat fixed-flux point-source exterior. Full nonlocal current, curved metric, finite-mouth data and uncontrolled gradient tails are not solved.',
        'limitations':['The gradient screen is sampled at turning points and is not a uniform derivative-expansion error bound.',
            'Far-loop failure cannot be cured by discarding or renormalizing its flux; the full formal average and flagged contributions are both retained.',
            'The outer-box bound concerns the exact linear local PDE only; numerical quadrature/discretization and physical omissions are not included.',
            'Corrected Shapiro path integrals, global state preparation and finite-mouth flux multipoles remain omitted.',
            'This does not establish a stable physical medium, formation, new physics or travel to yesterday.'],
        'sources':['https://arxiv.org/html/1807.04726v3#S3.SS1','https://arxiv.org/pdf/hep-th/9804143','https://arxiv.org/pdf/hep-th/0406216'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result=polarization_checks()
    project_path('magnetic_path_polarization_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','tests','state_shift_refinement','outer_domain_state_shift_change')},indent=2))
    print(json.dumps({k:result['records'][-1][k] for k in ('state_variation_at_fixed_geometry','linear_closed_state_shift','formal_corrected_state','outer_box_state_shift_bound','sampled_flux_failing_gradient_point1','sampled_formal_shift_from_failing_gradient_point1','sampled_baseline_state_from_failing_gradient_point1','max_kernel_disagreement_relative','max_fd_disagreement_relative')},indent=2))
    if not result['passed']: raise SystemExit('A local exterior control failed; diagnostics retained.')
