"""Full-path derivative diagnostics and conditional distant-field sensitivity.

These are screens for a flat point-charge model, not physical error bounds
on its nonlocal fermion determinant or on finite wormhole entrances.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq, minimize_scalar

from magnetic_path_polarization import path_variation
from mmp_pulse import _magnetic_loop_source, _magnetic_return_path


def path_diagnostics(nu, v, q=100000):
    """Contract exact Cartesian monopole derivatives without far-field cancellation."""
    v = np.asarray(v, dtype=float)
    if not (0 < nu <= 2 and q > 0) or np.any((v < 0) | (v >= 1)):
        raise ValueError('Require 0<nu<=2, q>0 and 0<=v<1; point mouths use analytic limits')
    h = (1-v)*(1+v); root = np.hypot(h, nu*v); u = (h+root)/nu
    um = (2-nu)*h/(nu*(u+1)-2*h)
    near = um+1-v; far = u+v; D = near*far
    t = near/far; omt = 2*v/far
    omc = 2*h/D; opc = 2*um*(u+1)/D; c = 1-omc
    # t=r_near/r_far, c=cos(angle between radius vectors). Positive forms
    # avoid subtracting almost equal dipole fields at large distances.
    b2 = (omt*(1+t))**2+2*t*t*omc
    j2 = 6*((omt*(1+t+t*t))**2+3*t**3*omc*opc)
    k2 = 90*((omt*(1+t+t*t+t**3))**2+t**4*omc*(5*c*c+5*c+2))
    return np.stack((np.sqrt(j2/(q/2*b2**1.5)), np.sqrt(k2)/(q/2*b2),
                     (q/2)*np.sqrt(b2)/(.5*near)**2), axis=-1)


def turning_diagnostics(nu, q=100000):
    return np.array([6/math.sqrt(q*nu)*math.sqrt(1-nu*nu/4),
        3*math.sqrt(5*(48-30*nu*nu+5*nu**4))/(q*nu), q*nu**3/2])


def scan_path(nu, q=100000, resolution=256):
    # ponytail: numerical extrema, not an interval proof; analytic mouth limits
    # and endpoint-clustered refinement cover the singular coordinate corners.
    v = np.unique(np.r_[np.linspace(0, .999, resolution+1),
                        1-np.geomspace(1e-10, .1, resolution//2)])
    values = path_diagnostics(nu, v, q)
    mouth = [math.sqrt(12/q), 6*math.sqrt(10)/q]
    maxima = []; locations = []
    for column in (0, 1):
        i = int(np.argmax(values[:, column]))
        best, location = float(values[i, column]), float(v[i])
        if 0 < i < len(v)-1 and best > mouth[column]*(1+1e-9):
            result = minimize_scalar(lambda s: -path_diagnostics(nu, s, q)[column],
                bounds=(v[i-1], v[i+1]), method='bounded', options={'xatol': 1e-13})
            if -result.fun > best: best, location = float(-result.fun), float(result.x)
        if mouth[column] >= best*(1-2e-14): best, location = mouth[column], 1.
        maxima.append(best); locations.append(location)
    return {'nu': float(nu), 'max_eta1': maxima[0], 'max_eta2': maxima[1],
            'max_locations_v': locations, 'min_B_times_d2': float(values[:, 2].min()),
            'min_B_location_v': float(v[np.argmin(values[:, 2])])}


def eta1_endpoint_certificate():
    """Exact Descartes certificate: eta1 has no interior path maximum."""
    import sympy as s
    u,n,z = s.symbols('u n z', real=True)
    V = u*(2-n*u)/(2*u-n)
    K = u*u+(3*u*u-1)*V+V*V
    J = 3*u*u*(u*u-1)+(6*u**4+3)*V+(9*u*u-3)*V*V+V**3
    G = 2*K*s.diff(J,u)-3*J*s.diff(K,u)
    H = s.cancel(G*(2*u-n)**5/(12*u*(u*u-1)))
    T = s.Poly(s.cancel(H.subs(u,(n+2*z)/(n*(1+z)))*n**9*(1+z)**9),z)
    assert s.degree(H,u)==9 and T.degree()==9
    intervals = [(s.Rational(1),s.Rational(1))]+[(s.Rational(a,10000),s.Rational(b,10000))
        for a,b in ((11204,11205),(12218,12219),(13034,13035),(13671,13672),
                    (14155,14156),(14509,14510),(14751,14752),(14890,14891),(14932,14933))]
    rows = []
    for j in range(10):
        raw = T.nth(j); P = s.Poly(s.cancel(raw/(n-2)**3),n); npow = ppow = 0
        while P.eval(0)==0: P = P.exquo(s.Poly(n,n)); npow += 1
        while P.eval(-2)==0: P = P.exquo(s.Poly(n+2,n)); ppow += 1
        content,P = P.primitive()
        assert content>0 and P.LC()>0
        signs = [s.sign(c) for c in P.all_coeffs() if c]
        assert sum(a!=b for a,b in zip(signs,signs[1:]))==1
        lo,hi = intervals[j]
        assert (P.degree()==1 and P.eval(lo)==0) if lo==hi else P.eval(lo)<0<P.eval(hi)
        assert s.cancel(raw-content*n**npow*(n+2)**ppow*(n-2)**3*P.as_expr())==0
        if j: assert intervals[j-1][1]<lo
        rows.append({'power_z':j,'power_n':npow,'power_n_plus_2':ppow,'positive_factor':int(content),
                     'polynomial':str(P.as_expr()),'root_interval':[str(lo),str(hi)]})
    # Ordered positive roots make low-to-high coefficients negative then positive,
    # if they change sign at all. Since du/dz>0, a derivative zero is a minimum.
    return {'passed':True,'coefficient_certificate':rows,
        'identity':'q*eta1²=6J/K^(3/2); sign(d eta1²/du)=sign(T) for 1<u<2/nu.',
        'conclusion':'sup eta1=max(6/sqrt(qnu)*sqrt(1-nu²/4),sqrt(12/q)); no interior maximum.'}


def validity_checks():
    root = Path(__file__).resolve().parent
    prior = json.loads((project_path('dirac_return_closure_results.json')).read_text())
    source85 = project_path('magnetic_path_polarization_results.json')
    saved85 = json.loads(source85.read_text())
    state = prior['records'][-1]['roots'][-1]
    q = prior['parameters']['q']; distance = prior['parameters']['distance_over_re']
    clock = state['background']['clock']; eps = saved85['epsilon']
    nus = np.array(state['nu']); weights = np.array(state['weights']); periods = np.array(state['periods'])
    sample_nu = np.unique(np.r_[np.geomspace(1e-5, .1, 90), np.linspace(.1, 2., 130),
                                1.26, math.sqrt(8/5), nus])
    coarse = [scan_path(float(n), q, 128) for n in sample_nu]
    fine = [scan_path(float(n), q, 256) for n in sample_nu]
    refinement = max(abs(a[k]-b[k])/(1+b[k]) for a,b in zip(coarse, fine) for k in ('max_eta1','max_eta2'))
    endpoint_prediction_error = max(abs(r['max_eta1']-max(turning_diagnostics(r['nu'], q)[0],
        math.sqrt(12/q)))/(1+r['max_eta1']) for r in fine)
    counter_v = math.sqrt(28/833); counter = path_diagnostics(1.26,counter_v,q)[1]
    counter_exact = math.sqrt(945276640861045/2562028497522)/q
    certificate = eta1_endpoint_certificate()
    # Independent reference uses the existing Cartesian field/Jacobian helper;
    # its Jacobian's fourth-order finite difference checks the Hessian norm.
    references = []
    for nu,v in ((.02,0.),(.1,.4),(.5,.8),(1.,0.),(1.26,counter_v),(1.7,.6),(1.99,.2)):
        rho,z,_ = _magnetic_return_path(nu, v); xyz = np.array([rho,0.,z])
        row = _magnetic_loop_source(xyz, 0., circumference=1.)
        B = q/2*row['magnetic']; J = q/2*row['magnetic_derivative']
        analytic = path_diagnostics(nu, v, q)
        jerror = abs(np.linalg.norm(J)/np.linalg.norm(B)**1.5/analytic[0]-1)
        herrors = []
        radius = min(np.linalg.norm(xyz-[0,0,.5]), np.linalg.norm(xyz-[0,0,-.5]))
        for scale in (2e-4, 1e-4):
            step = scale*radius; derivatives = []
            for axis in range(3):
                offset = np.eye(3)[axis]*step
                samples = [_magnetic_loop_source(xyz+j*offset,0.,circumference=1.)['magnetic_derivative']
                           for j in (-2,-1,1,2)]
                derivatives.append(q/2*(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*step))
            eta2 = np.linalg.norm(derivatives)/np.linalg.norm(B)**2
            herrors.append(abs(eta2/analytic[1]-1))
        references.append({'nu':nu,'v':v,'jacobian_relative_error':float(jerror),
                           'hessian_relative_errors':herrors})
    turning_error = max(float(np.max(abs(path_diagnostics(float(n),0.,q)/turning_diagnostics(float(n),q)-1)))
                        for n in sample_nu if n < 2.)
    mouth_error = max(float(np.max(abs(path_diagnostics(float(n),1-1e-9,q)[:2]/
                    [math.sqrt(12/q),6*math.sqrt(10)/q]-1))) for n in (.001,.1,1.,1.9,2.))
    charge_error = float(np.max(abs(path_diagnostics(.1,.4,4*q)/path_diagnostics(.1,.4,q)-[.5,.25,4])))
    # A common nominal slow-variation screen uses eta1<h and eta2<h².
    # Choosing h=.1 supplies neither a determinant remainder nor a leakage bound.
    threshold = .1
    cut1 = brentq(lambda n: turning_diagnostics(n,q)[0]-threshold, .001, .2, xtol=1e-14)
    cut2 = brentq(lambda n: turning_diagnostics(n,q)[1]-threshold**2, .001, .2, xtol=1e-14)
    cut = max(cut1,cut2); U = 2/cut
    stored = [scan_path(float(n),q) for n in nus]
    flagged = np.array([r['max_eta1']>=threshold or r['max_eta2']>=threshold**2 for r in stored])
    Bmin = np.array([r['min_B_times_d2'] for r in stored])
    # Lowest positive NS frequency pi/L, not the state parameter 2pi/L.
    ns_gap = (math.pi*distance/(clock*periods))**2/(2*Bmin)
    # Exact homogeneous solution for arbitrary outer boundary perturbations.
    harmonic = lambda u,v: ((u*u-1)*(1-v*v), 2*u*(1-v*v), -2*v*(u*u-1))
    path = path_variation(1., harmonic, order=256, outer=U)
    J = path['delta_length_per_epsilon']; f = path['base_length']
    with mp.workdps(40):
        x = 1/mp.mpf(U)
        I1 = mp.atanh(x)-(mp.polylog(2,x)-mp.polylog(2,-x))/2
        tail_norm = 16*mp.pi/(3*mp.sqrt(3))
        tail_length_bound = float(tail_norm/2*I1*J)
    boundary_factor = 2*U/(U*U-1)**2
    base_shortening = -eps*next(p for p in saved85['records'][-1]['representative_paths'] if p['nu']==1.)['relative_length_coefficient']
    profiles = []
    for n in (.02,.05,1.,1.9):
        v = np.linspace(0,.999,250); d = path_diagnostics(n,v,q)
        profiles.append({'nu':n,'v':v.tolist(),'eta1':d[:,0].tolist(),'sqrt_eta2':np.sqrt(d[:,1]).tolist()})
    tests = {
        'cartesian_first_and_second_derivatives':max(r['jacobian_relative_error'] for r in references)<2e-11 and max(max(r['hessian_relative_errors']) for r in references)<2e-7,
        'turning_and_mouth_limits':turning_error<3e-11 and mouth_error<2e-7,
        'charge_scaling':charge_error<1e-13,
        'exact_eta1_endpoint_certificate':certificate['passed'],
        'full_path_extrema_refinement':refinement<2e-9 and endpoint_prediction_error<2e-8,
        'second_derivative_interior_counterexample':abs(counter/counter_exact-1)<2e-13 and counter>max(turning_diagnostics(1.26,q)[1],6*math.sqrt(10)/q),
        'field_minimum_at_turning':all(r['min_B_location_v']==0 for r in fine),
        'conditional_harmonic_length_response':path['direct_vs_positive_kernel']<2e-7 and path['richardson_derivative_error']<2e-7,
    }
    return {'passed':bool(all(tests.values())),'tests':{k:bool(v) for k,v in tests.items()},
        'q':q,'screen':{'eta1_limit':threshold,'eta2_limit':threshold**2,'nu_cut_eta1':cut1,'nu_cut_eta2':cut2,
            'joint_nu_cut':cut,'enclosing_prolate_u':U,'continuum_flux_below_cut':cut/2,
            'sampled_flagged_flux':float(weights[flagged].sum()),
            'sampled_flagged_baseline_state':float(weights[flagged]@(2*math.pi/periods[flagged])**2),
            'sampled_flagged_formal_shift85':float(np.array(saved85['records'][-1]['formal_state_shift_by_channel'])[flagged].sum())},
        'representative_paths':[scan_path(n,q) for n in (.02,.05,.1,.5,1.,1.5,1.9,2.)],
        'profiles':profiles,'extrema':fine,'references':references,'eta1_certificate':certificate,
        'refinement':refinement,'eta1_endpoint_max_prediction_error':endpoint_prediction_error,
        'eta2_endpoint_counterexample':{'nu':1.26,'u':14/9,'v':counter_v,'eta2':float(counter),
            'exact_q2_eta2_squared':'945276640861045/2562028497522',
            'eta2_turning':float(turning_diagnostics(1.26,q)[1]),'eta2_mouth':6*math.sqrt(10)/q},
        'turning_formula_error':turning_error,'mouth_limit_error':mouth_error,'charge_scaling_error':charge_error,
        'lowest_NS_gap':{'nu':nus.tolist(),'E2_over_2Bmin':ns_gap.tolist(),
            'maximum_on_screened_paths':float(ns_gap[~flagged].max()),
            'maximum_on_all_stored_paths':float(ns_gap.max()),
            'scope':'Lowest positive NS mode only, not the ultraviolet vacuum sum or an arbitrary probe.'},
        'conditional_tail':{'U':U,'nu1_harmonic_length_response':J,'nu1_base_length':f,
            'assumed_boundary_norm':'abs(delta_psi(U,v))<=beta*psi0(U,v); regular point-mouth data fixed and no extra interior source.',
            'relative_length_bound_per_beta':boundary_factor*J/f,
            'beta_bound_for_sign_in_linear_model':base_shortening/(boundary_factor*J/f),
            'illustrative_beta_point01_length_bound_ppm':.01*boundary_factor*J/f*1e6,
            'assumed_l1_source_norm':'abs(delta_s1(u))<=A1/u for u>=U; illustrative A1=16pi/(3sqrt3), not a physical bound.',
            'l1_source_norm':float(tail_norm),'l1_kernel_integral':float(I1),
            'l1_length_coefficient_bound':tail_length_bound,
            'l1_relative_length_bound_ppm_if_order_epsilon':eps*tail_length_bound/f*1e6},
        'scope':'Exact Cartesian derivative contractions and numerical extrema for flat opposite point magnetic charges. Thresholds are chosen screens, not physical error bounds. Tail bounds concern the linear flux equation and first length variation, conditional on stated norms.',
        'limitations':['No all-orders derivative expansion or massless nonlocal determinant is controlled by these diagnostics.',
            'The point-mouth limits are flat-space benchmarks; the actual curved finite-mouth collar needs matching.',
            'No bound on the unknown tail source/boundary amplitude has been derived.',
            'No stable medium, new law, or 13 Sep 2027 to 12 Sep 2027 passage is established.'],
        'sources':['https://arxiv.org/html/1807.04726v3#S3.SS1','https://arxiv.org/html/1807.04726v3','https://arxiv.org/abs/hep-th/9804143','https://arxiv.org/pdf/1503.06099'],
        'prior85_sha256':hashlib.sha256(source85.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result = validity_checks()
    project_path('magnetic_path_validity_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','tests','screen','refinement','conditional_tail')},indent=2))
    if not result['passed']: raise SystemExit('A validity diagnostic control failed; results retained.')
