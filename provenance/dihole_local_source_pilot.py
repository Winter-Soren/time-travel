"""Local action variations for the actual compact null-source probe.

These are coefficients of specified local action terms, not a calculation of
the finite renormalized Dirac stress. Couplings cannot be fitted to manufacture
support. Reuse the actual fields, profile, and independent curvature checker.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.special import expit

from balanced_dihole import symbolic_fields, curvature_check
from dihole_stretched_collar import exact_target_parameters, gauss, radius_profile, smooth_step, null_residual
from dihole_bridge_scattering import masked_geometric_target
from project_paths import project_path


def local_first(function, p, radius, theta):
    """Analytic local first derivatives; complex inputs retain complex-step jets."""
    radius, theta = np.broadcast_arrays(radius, theta)
    k = p['k']; distance = np.sqrt(radius**2+4*k*k+4*k*radius*np.cos(theta))
    dr = (radius+2*k*np.cos(theta))/distance
    dt = -2*k*radius*np.sin(theta)/distance
    raw = function(p['M']+(radius+distance)/2, (distance-radius)/(2*k), p['M'],p['a'],p['B'])
    raw = np.array([np.broadcast_to(x,radius.shape) for x in raw])
    first_r = raw[3:6]*(1+dr)/2+raw[6:9]*(dr-1)/(2*k)
    first_t = raw[3:6]*dt/2+raw[6:9]*dt/(2*k)
    first_r[:2] /= raw[:2]; first_t[:2] /= raw[:2]
    values = raw[:3].copy(); values[:2] = np.log(values[:2])
    return values, first_r, first_t


def mask_jets(u):
    t = u/.05
    h = -1/t+1/(1-t)
    hp = 1/t**2+1/(1-t)**2
    hpp = -2/t**3+2/(1-t)**3
    hppp = 6/t**4+6/(1-t)**4
    sigmoid = expit(h); product = expit(h)*expit(-h)
    f = smooth_step(t)[1]/2
    fs = product*((1-2*sigmoid)*hp**2+hpp)/(.1)
    fss = product*((1-6*sigmoid+6*sigmoid**2)*hp**3
                   +3*(1-2*sigmoid)*hp*hpp+hppp)/(.005)
    zero = np.zeros_like(f)
    return np.array([f,fs,zero,fss,zero,zero])


def prepare(function, p, order):
    nodes, weights = gauss(order)
    u, wu = .025*(nodes+1), .025*weights
    radius, chi, second = radius_profile(u,order=max(64,order))
    r, theta = np.broadcast_arrays(radius[:,None],np.arccos(nodes)[None,:])
    x, xr, xt = local_first(function,p,r,theta)
    tiny = 1e-25
    xrr = np.imag(local_first(function,p,r+1j*tiny,theta)[1])/tiny
    _, cross, xtt = local_first(function,p,r,theta+1j*tiny)
    cross, xtt = np.imag(cross)/tiny,np.imag(xtt)/tiny
    c, ss = chi[:,None],second[:,None]
    jets = np.array([x,xr*c,xt,xrr*c*c+xr*ss,cross*c,xtt])
    zero = np.zeros_like(r)
    log_r = np.array([np.log(r),np.broadcast_to(c/r,r.shape),zero,
                      np.broadcast_to(ss/r-c*c/r**2,r.shape),zero,zero])
    log_sin = np.array([np.log(np.sin(theta)),zero,1/np.tan(theta),zero,zero,-1/np.sin(theta)**2])
    log_metric = np.stack((jets[:,0],jets[:,1],jets[:,1]+log_r,
                           log_r+log_sin-jets[:,0]),axis=1)
    mask = np.broadcast_to(mask_jets(u)[:,None,:,None],(6,1,order,order)).copy()
    # Integrate both collars and phi; mu quadrature already supplies sin(theta).
    measure = 4*np.pi*wu[:,None]*weights[None,:]*np.exp(2*x[1])*r*r
    fields,_,_,extra = null_residual(function,p,r,c,ss,np.cos(theta))
    gphi = r*r*np.sin(theta)**2/fields['lapse']**2
    maxwell_null = 2*(xr[2]*c)**2/(fields['meridian_factor']**2*gphi)
    return dict(logs=log_metric.reshape(6,4,-1), mask=mask.reshape(6,1,-1),
                potential=jets[:,2].reshape(6,-1), measure=measure.ravel(),
                geometric_derivative=float(2*np.sum(measure*mask[0,0]*(extra+maxwell_null))),
                maxwell_derivative=float(8*np.sum(measure*mask[0,0]*(xr[2]*c)**2
                                                  /(fields['meridian_factor']**2*gphi))),
                u=u, nodes=nodes, order=order)


def metric_jets(prepared, epsilon=0.):
    logs = prepared['logs']+epsilon*prepared['mask']*np.array([1.,-1.,0.,0.])[None,:,None]
    g = np.exp(2*logs[0])*np.array([-1.,1.,1.,1.])[:,None]
    values = np.zeros((6,5,g.shape[-1]),dtype=g.dtype)
    values[0,:4] = g
    values[1:3,:4] = 2*g[None]*logs[1:3]
    for dest,a,b in ((3,1,1),(4,1,2),(5,2,2)):
        values[dest,:4] = 2*g*(logs[dest]+2*logs[a]*logs[b])
    values[:,4] = prepared['potential']
    return values


def invariants(values):
    """Batched version of balanced_dihole.curvature_check's diagonal contraction."""
    metric = values[0,:4]; inverse = 1/metric
    n = metric.shape[-1]
    dg = np.zeros((4,4,4,n),dtype=values.dtype)
    ddg = np.zeros((4,4,4,4,n),dtype=values.dtype)
    for i in range(4):
        dg[1,i,i],dg[2,i,i] = values[1,i],values[2,i]
        ddg[1,1,i,i],ddg[2,2,i,i] = values[3,i],values[5,i]
        ddg[1,2,i,i] = ddg[2,1,i,i] = values[4,i]
    connection = np.zeros((4,4,4,n),dtype=values.dtype)
    derivative = np.zeros((4,4,4,4,n),dtype=values.dtype)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                bracket = dg[b,a,c]+dg[c,a,b]-dg[a,b,c]
                connection[a,b,c] = inverse[a]*bracket/2
                for d in (1,2):
                    derivative[d,a,b,c] = (inverse[a]*(ddg[d,b,a,c]+ddg[d,c,a,b]-ddg[d,a,b,c])
                        -inverse[a]**2*dg[d,a,a]*bracket)/2
    ricci = np.zeros((4,4,n),dtype=values.dtype)
    for a in range(4):
        for b in range(4):
            ricci[a,b] = sum(derivative[c,c,a,b]-derivative[b,c,a,c]
                +sum(connection[c,c,d]*connection[d,a,b]-connection[c,b,d]*connection[d,a,c]
                     for d in range(4)) for c in range(4))
    scalar = sum(ricci[a,a]*inverse[a] for a in range(4))
    ricci_square = np.sum(ricci**2*inverse[:,None]*inverse[None,:],axis=(0,1))
    magnetic_square = 2*inverse[3]*(inverse[1]*values[1,4]**2+inverse[2]*values[2,4]**2)
    return np.array([scalar,scalar**2,ricci_square,magnetic_square]),ricci


def evaluate(prepared, epsilon=0.):
    values = metric_jets(prepared,epsilon)
    scalar,ricci = invariants(values)
    return scalar@prepared['measure'],values,ricci


def run():
    start = perf_counter()
    names=('balanced_dihole.py','dihole_stretched_collar.py','dihole_bridge_scattering.py','project_paths.py')
    hashes={str(project_path(n)):hashlib.sha256(project_path(n).read_bytes()).hexdigest() for n in names}
    function,_=symbolic_fields(); p=exact_target_parameters()
    records=[]; checks=[]
    for order in (32,64,96):
        prepared=prepare(function,p,order)
        baseline,values,ricci=evaluate(prepared)
        gradient=np.imag(evaluate(prepared,1e-25j)[0])/1e-25
        # A=-q/(2Q) times the classical Emparan potential; keep both conventions.
        charged_gradient=gradient[3]*p['physical_A_factor']**2
        records.append(dict(order=order,action_per_time=baseline.tolist(),
                            action_derivative=gradient.tolist(),energy_derivative=(-gradient).tolist(),
                            charge_included_F2_energy_derivative=float(-charged_gradient),
                            einstein_hilbert_identity_error=float(abs(gradient[0]-prepared['geometric_derivative'])),
                            maxwell_identity_error=float(abs(gradient[3]-prepared['maxwell_derivative']))))
        print(json.dumps(records[-1]),flush=True)
        if order==64:
            # Existing independently scalar-coded tensor routine at distinct points.
            for index in (17,1043,3009):
                one=values[:,:,index]
                checked=curvature_check(lambda *args:one.ravel(),p,0.,0.)
                checks.append(float(np.max(abs(np.asarray(checked['ricci_covariant'])-ricci[:,:,index]))))
            differences=[]
            for epsilon in (1e-4,5e-5):
                fd=(evaluate(prepared,epsilon)[0]-evaluate(prepared,-epsilon)[0])/(2*epsilon)
                differences.append(dict(epsilon=epsilon,relative_error=(abs(fd-gradient)/np.maximum(abs(gradient),1)).tolist()))
            fd_records=differences
    a,b=records[-2:]
    change=np.abs(np.array(a['action_derivative'])-b['action_derivative'])/np.maximum(abs(np.array(b['action_derivative'])),1)
    tests=dict(independent_ricci=max(checks)<1e-7,
               paired_EH_variation=max(r['einstein_hilbert_identity_error'] for r in records)<2e-7,
               paired_Maxwell_variation=max(r['maxwell_identity_error'] for r in records)<2e-9,
               quadrature_sensitivity=float(max(change))<1e-5,
               derivative_difference=max(fd_records[-1]['relative_error'])<2e-6,
               finite_coefficients=bool(np.all(np.isfinite(b['action_derivative']))),
               input_hashes=all(hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in hashes.items()))
    result=dict(passed=all(tests.values()),tests=tests,records=records,
                invariant_order=['R','R_squared','Ricci_ab_Ricci_up_ab','F_classical_squared'],
                convention='Lorentzian action S=sum c_i integral sqrt(-g) I_i; static energy E=-S/T; the paired volume measure is invariant.',
                units=['R action derivative: D','quadratic-curvature action derivatives: D^-1','classical F_squared derivative: D; charge-included F_squared derivative: D^-1'],
                independent_ricci_errors=checks,finite_difference_controls=fd_records,
                relative_quadrature_changes=change.tolist(),
                required_extra_source=masked_geometric_target(p,128)/(8*np.pi),
                scope=['Actual compact paired variation on both129 collars; unvaried bridge/exterior contribute zero to these local derivatives.',
                       'These are response coefficients of local action terms, not their measured couplings or the finite Dirac vacuum stress.',
                       'A common scheme shift changes quantum and gravitational/gauge bookkeeping together; it cannot create physical support.',
                       'Couplings and subtraction scale must be specified before an absolute source prediction.',
                       'Floating quadrature, complex-step and finite-difference checks are diagnostics, not interval enclosures.',
                       'No supporting current, stable material, formation mechanism or chronology change.'],
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashes,
                sources=['https://arxiv.org/abs/1210.4031','https://arxiv.org/abs/hep-th/0306138'],
                runtime_seconds=perf_counter()-start)
    project_path('dihole_local_source_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,runtime_seconds=result['runtime_seconds']),indent=2))
    if not result['passed']:raise RuntimeError('Local-source controls failed; result retained.')
    return result


if __name__=='__main__':
    run()
