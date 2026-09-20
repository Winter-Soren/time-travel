"""One actual conforming m=0 exterior trial and its differential residual.

The pilot is the right Weyl field, phase pi, upper compact pressure probe,
rho_c=32, x=0.44391506147705256. It recovers a continuous physical field
from one frozen-model DG solve, imposes an exact finite spherical trace,
and tapers to zero before the polygon outer boundary. Residual quadrature
is a numerical estimate, not a certified continuum error enclosure.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import splu

from angular_defect_bound import defects, norm
from complete_angular_source import CompleteSector, _decode, port_labels
from exterior_boundary_map import _complex_json, pair_potential
from spherical_weyl_leads import mouth_characteristics, zero_mode


def fade(radius, start, end):
    """C2 cutoff equal to one below start and zero above end."""
    s = np.clip((np.asarray(radius)-start)/(end-start),0.,1.)
    value = 1-10*s**3+15*s**4-6*s**5
    derivative = (-30*s**2+60*s**3-30*s**4)/(end-start)
    value = np.where(np.asarray(radius)>=end,0.,value)
    return value,derivative


def angular_basis(theta, port, b=.2, nmax=4):
    """Physical m0 (outgoing plus, incoming minus) and exact theta derivatives."""
    theta = np.asarray(theta)
    ep = np.stack((np.cos(theta/2),np.sin(theta/2)),axis=1)
    em = np.stack((-np.sin(theta/2),np.cos(theta/2)),axis=1)
    dep, dem = em/2,-ep/2
    cosine,sine = np.cos(theta),np.sin(theta)
    f,g,df,dg = [],[],[],[]
    for n in range(1,nmax+1):
        p = np.polynomial.legendre.Legendre.basis(n)
        first,second = p.deriv(),p.deriv(2)
        scale = math.sqrt((2*n+1)/2)/b
        lam = math.sqrt(n*(n+1))
        f.append(scale*p(cosine));df.append(-scale*sine*first(cosine))
        g.append(scale*sine*first(cosine)/lam)
        dg.append(scale*(cosine*first(cosine)-sine*sine*second(cosine))/lam)
    f,g,df,dg = (np.asarray(a).T for a in (f,g,df,dg))
    if port=='upper':
        plus = 1j*em[:,:,None]*g[:,None,:]
        dplus = 1j*(dem[:,:,None]*g[:,None,:]+em[:,:,None]*dg[:,None,:])
        minus = ep[:,:,None]*f[:,None,:]
        dminus = dep[:,:,None]*f[:,None,:]+ep[:,:,None]*df[:,None,:]
        minus = np.concatenate((ep[:,:,None]/(math.sqrt(2)*b),minus),axis=2)
        dminus = np.concatenate((dep[:,:,None]/(math.sqrt(2)*b),dminus),axis=2)
    elif port=='lower':
        plus = em[:,:,None]*f[:,None,:]
        dplus = dem[:,:,None]*f[:,None,:]+em[:,:,None]*df[:,None,:]
        minus = 1j*ep[:,:,None]*g[:,None,:]
        dminus = 1j*(dep[:,:,None]*g[:,None,:]+ep[:,:,None]*dg[:,None,:])
        plus = np.concatenate((em[:,:,None]/(math.sqrt(2)*b),plus),axis=2)
        dplus = np.concatenate((dem[:,:,None]/(math.sqrt(2)*b),dplus),axis=2)
    else:
        raise ValueError('Unknown mouth')
    return plus,minus,dplus,dminus


def stable_potential(points):
    """Return A and 1-A without subtracting their O(rho^2) axis corrections."""
    rho,z = np.asarray(points).T
    low,up = z+.5,z-.5
    rl,ru = np.hypot(rho,low),np.hypot(rho,up)
    sl,su = np.sign(low),np.sign(up)
    axis = (sl-su)/2
    correction = -.5*rho*rho*(sl/(rl*(rl+abs(low)))-su/(ru*(ru+abs(up))))
    return axis+correction,1-axis-correction


def differential(points,u,drho,dz,kappa):
    """i(H-i*kappa)u for q1,m0 in the actual pair potential; norm unchanged."""
    rho = points[:,0]
    if np.any(rho<=0):
        raise ValueError('Volume quadrature must avoid the measure-zero axis')
    a,one_minus_a = stable_potential(points)
    result = np.empty_like(u)
    result[:,0] = drho[:,1]+dz[:,0]+(one_minus_a/rho)[:,None]*u[:,1]+kappa*u[:,0]
    result[:,1] = drho[:,0]-dz[:,1]+(a/rho)[:,None]*u[:,0]+kappa*u[:,1]
    return result


class Trial:
    """Continuous P1 recovery, disjoint analytic mouth blends and outer taper."""

    def __init__(self,vertices,triangles,nodal,incoming,outgoing,b=.2,a=.30,
                 taper_start=1.5,taper_end=1.9):
        self.vertices,self.triangles = np.asarray(vertices),np.asarray(triangles)
        self.nodal = np.asarray(nodal)
        self.incoming,self.outgoing = np.asarray(incoming),np.asarray(outgoing)
        self.b,self.a = b,a
        self.taper_start,self.taper_end = taper_start,taper_end
        coordinates = self.vertices[self.triangles]
        inverse = np.linalg.inv(np.concatenate((np.ones((*coordinates.shape[:2],1)),coordinates),axis=2))
        self.affine = np.einsum('tij,tjsc->tisc',inverse,self.nodal[self.triangles])
        self.labels_in,self.labels_out = port_labels(.5,4)

    def mouth(self,points,port):
        local = np.array(points,copy=True)
        local[:,1] -= .5 if port=='upper' else -.5
        radius = np.linalg.norm(local,axis=1)
        theta = np.arctan2(local[:,0],local[:,1])
        plus,minus,dp,dm = angular_basis(theta,port,self.b)
        ip = [i for i,(p,_) in enumerate(self.labels_in) if p==port]
        op = [i for i,(p,_) in enumerate(self.labels_out) if p==port]
        boundary = plus@self.outgoing[op]+minus@self.incoming[ip]
        tangent = dp@self.outgoing[op]+dm@self.incoming[ip]
        value = (self.b/radius)[:,None,None]*boundary
        radial_derivative = -value/radius[:,None,None]
        theta_derivative = (self.b/radius)[:,None,None]*tangent
        sr,cr = local[:,0]/radius,local[:,1]/radius
        drho = sr[:,None,None]*radial_derivative+(cr/radius)[:,None,None]*theta_derivative
        dz = cr[:,None,None]*radial_derivative-(sr/radius)[:,None,None]*theta_derivative
        return value,drho,dz,radius,sr,cr

    def evaluate(self,points,triangle_ids):
        points = np.asarray(points)
        co = self.affine[np.asarray(triangle_ids)]
        value = co[:,0]+points[:,0,None,None]*co[:,1]+points[:,1,None,None]*co[:,2]
        drho,dz = co[:,1].copy(),co[:,2].copy()
        cg = value.copy()
        for port in ('upper','lower'):
            radius = np.hypot(points[:,0],points[:,1]-(.5 if port=='upper' else -.5))
            mask = radius<self.a
            if not mask.any():
                continue
            v,vr,vz,r,sr,cr = self.mouth(points[mask],port)
            chi,dchi = fade(r,self.b,self.a)
            difference = v-value[mask]
            drho[mask] += chi[:,None,None]*(vr-drho[mask])+(dchi*sr)[:,None,None]*difference
            dz[mask] += chi[:,None,None]*(vz-dz[mask])+(dchi*cr)[:,None,None]*difference
            value[mask] += chi[:,None,None]*difference
        before_taper = value.copy()
        radius = np.linalg.norm(points,axis=1)
        taper,dt = fade(radius,self.taper_start,self.taper_end)
        sr = np.divide(points[:,0],radius,out=np.zeros(len(radius)),where=radius>0)
        cr = np.divide(points[:,1],radius,out=np.zeros(len(radius)),where=radius>0)
        drho = taper[:,None,None]*drho+(dt*sr)[:,None,None]*value
        dz = taper[:,None,None]*dz+(dt*cr)[:,None,None]*value
        value = taper[:,None,None]*value
        return value,drho,dz,cg,before_taper


def recovered_vertices(model,solution):
    values = solution.reshape(len(model.triangles),3,2,-1)
    xy = model.vertices[model.triangles]
    first,second = xy[:,1]-xy[:,0],xy[:,2]-xy[:,0]
    area = abs(first[:,0]*second[:,1]-first[:,1]*second[:,0])/2
    nodal = np.zeros((len(model.vertices),2,solution.shape[1]),complex)
    weight = np.zeros(len(model.vertices))
    for vertex in range(3):
        np.add.at(nodal,model.triangles[:,vertex],area[:,None,None]*values[:,vertex])
        np.add.at(weight,model.triangles[:,vertex],area)
    used = weight>0
    nodal[used] /= weight[used,None,None]
    axis = np.flatnonzero(abs(model.vertices[:,0])<1e-13)
    wrong = []
    for i in axis:
        component = 0 if -.5<model.vertices[i,1]<.5 else 1
        wrong.append(float(np.max(abs(nodal[i,component]))))
        nodal[i,component] = 0.
    return nodal,dict(axis_vertices=len(axis),maximum_forbidden_axis_value_before_enforcing=max(wrong))


def controls(trial,model,kappa,source,bl,il,ol):
    theta = np.linspace(.08,math.pi-.08,51)
    basis_error = derivative_error = 0.
    for port in ('upper','lower'):
        p,m,dp,dm = angular_basis(theta,port)
        lp,lm = mouth_characteristics(theta,.2,port,4)
        zero = zero_mode(theta,.2,port)[:,:,None]
        if port=='upper':lm=np.concatenate((zero,lm),axis=2)
        else:lp=np.concatenate((zero,lp),axis=2)
        basis_error = max(basis_error,float(np.max(abs(p-np.moveaxis(lp,0,1)))),
                          float(np.max(abs(m-np.moveaxis(lm,0,1)))))
        step = 2e-6
        forward,backward = angular_basis(theta+step,port),angular_basis(theta-step,port)
        derivative_error = max(derivative_error,*[float(np.max(abs((forward[j]-backward[j])/(2*step)-der)))
                                                  for j,der in ((0,dp),(1,dm))])
    nodes,weights = np.polynomial.legendre.leggauss(96)
    th = np.arccos(nodes)
    recovered_in,recovered_out = np.zeros_like(trial.incoming),np.zeros_like(trial.outgoing)
    pole_error = high_projection = characteristic_error = 0.
    for port in ('upper','lower'):
        center = .5 if port=='upper' else -.5
        points = np.column_stack((.2*np.sin(th),center+.2*np.cos(th)))
        value = trial.mouth(points,port)[0]
        pp,mm = mouth_characteristics(th,.2,port,8)
        zero = zero_mode(th,.2,port)[:,:,None]
        if port=='upper':mm=np.concatenate((zero,mm),axis=2)
        else:pp=np.concatenate((zero,pp),axis=2)
        pp,mm = np.moveaxis(pp,0,1),np.moveaxis(mm,0,1)
        pin = .04*np.einsum('p,pai,pac->ic',weights,mm.conj(),value)
        pout = .04*np.einsum('p,pai,pac->ic',weights,pp.conj(),value)
        ip = [i for i,(p,_) in enumerate(trial.labels_in) if p==port]
        op = [i for i,(p,_) in enumerate(trial.labels_out) if p==port]
        recovered_in[ip],recovered_out[op] = pin[:len(ip)],pout[:len(op)]
        high_projection = max(high_projection,float(np.max(abs(pin[len(ip):]))),float(np.max(abs(pout[len(op):]))))
        gamma_value = -np.stack((np.cos(th)[:,None]*value[:,0]+np.sin(th)[:,None]*value[:,1],
                                np.sin(th)[:,None]*value[:,0]-np.cos(th)[:,None]*value[:,1]),axis=1)
        expected = pp[:,:,:len(op)]@trial.outgoing[op]-mm[:,:,:len(ip)]@trial.incoming[ip]
        characteristic_error = max(characteristic_error,float(np.max(abs(gamma_value-expected))))
    incoming_high = np.setdiff1d(np.arange(9),il)
    boundary_error = max(float(np.max(abs(recovered_in[il]-bl@recovered_out[ol]-source))),
                         float(np.max(abs(recovered_in[incoming_high]))))
    continuity_error = 0.
    for edge in model.interior:
        point = model.vertices[list(edge['vertices'])].mean(axis=0)[None,:]
        if point[0,0]<=0 or min(np.hypot(point[0,0],point[0,1]-.5),np.hypot(point[0,0],point[0,1]+.5))<.2:
            continue
        left = trial.evaluate(point,[edge['triangle']])[0]
        right = trial.evaluate(point,[edge['other']])[0]
        continuity_error = max(continuity_error,float(np.max(abs(left-right))))
    centroids = model.vertices[model.triangles].mean(axis=1)
    indices = np.flatnonzero((centroids[:,0]>.01)&(np.linalg.norm(centroids,axis=1)<1.88)&
        (np.hypot(centroids[:,0],centroids[:,1]-.5)>.205)&(np.hypot(centroids[:,0],centroids[:,1]+.5)>.205))
    indices = indices[np.linspace(0,len(indices)-1,min(160,len(indices)),dtype=int)]
    points = centroids[indices]
    value,drho,dz,*_ = trial.evaluate(points,indices)
    numerical = []
    step = 1e-6
    for axis in range(2):
        shift = np.eye(2)[axis]*step
        derivative = (-trial.evaluate(points+2*shift,indices)[0]+8*trial.evaluate(points+shift,indices)[0]
                      -8*trial.evaluate(points-shift,indices)[0]+trial.evaluate(points-2*shift,indices)[0])/(12*step)
        numerical.append(derivative)
    field_derivative_error = max(float(np.max(abs(a-b))) for a,b in zip((drho,dz),numerical))
    q = differential(points,value,drho,dz,kappa)
    qnumerical = differential(points,value,*numerical,kappa)
    potential_error = float(np.max(abs(stable_potential(points)[0]-pair_potential(points[:,0],points[:,1]))))
    # The minimum outer facet apothem is the distance to its supporting line.
    apothem = min(abs(float(np.dot(e['normal'],model.vertices[e['vertices'][0]])))
                  for e in model.boundary if e['kind']=='outer')
    axis_samples = []
    for edge in model.boundary:
        if edge['kind']!='axis':continue
        point = model.vertices[list(edge['vertices'])].mean(axis=0)
        if abs(point[1])>=trial.taper_end:continue
        ratios=[]
        forbidden = 0 if -.5<point[1]<.5 else 1
        for rho in (1e-6,1e-8):
            value = trial.evaluate(np.array([[rho,point[1]]]),[edge['triangle']])[0]
            ratios.append(float(np.linalg.norm(value[0,forbidden])/rho))
        axis_samples.append(dict(z=float(point[1]),wrong_component_over_rho=ratios))
    return dict(analytic_basis_error=basis_error,analytic_basis_theta_derivative_error=derivative_error,
        physical_round_boundary_closure_error=boundary_error,round_characteristic_error=characteristic_error,
        sampled_n5_to_n8_trace_projection=high_projection,shared_edge_continuity_error=continuity_error,
        field_cartesian_derivative_error=field_derivative_error,
        differential_residual_cartesian_check=float(np.max(abs(q-qnumerical))),
        stable_pair_potential_error=potential_error,minimum_outer_apothem=apothem,
        axis_samples=axis_samples)


def integrate_residual(trial,model,kappa,order):
    from physical_exterior_quadrature import rule
    quadrature = rule(model.vertices,model.triangles,order=order,b=.2,d=1.)
    if isinstance(quadrature,dict):
        points,weights,indices = (quadrature[k] for k in ('points','weights','triangle_indices'))
    else:
        points,weights,indices = quadrature
    points,weights,indices = np.asarray(points),np.asarray(weights),np.asarray(indices)
    positive = weights>0
    omitted_zero_weight_nodes = int(np.count_nonzero(~positive))
    points,weights,indices = points[positive],weights[positive],indices[positive]
    triangles = model.vertices[model.triangles[indices]]
    first,second,offset = triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0],points-triangles[:,0]
    determinant = first[:,0]*second[:,1]-first[:,1]*second[:,0]
    alpha = (offset[:,0]*second[:,1]-offset[:,1]*second[:,0])/determinant
    beta = (first[:,0]*offset[:,1]-first[:,1]*offset[:,0])/determinant
    minimum_barycentric = float(min(alpha.min(),beta.min(),(1-alpha-beta).min()))
    minimum_mouth_radius = float(np.minimum(np.hypot(points[:,0],points[:,1]-.5),np.hypot(points[:,0],points[:,1]+.5)).min())
    if minimum_barycentric < -1e-11 or minimum_mouth_radius < .2-1e-12:
        raise ArithmeticError('Physical quadrature ownership/domain failed')
    norm2 = residual2 = 0.
    regional = {k:0. for k in ('collars','interior','outer_taper')}
    columns = np.zeros(trial.incoming.shape[1])
    for start in range(0,len(points),12000):
        p,w,it = points[start:start+12000],weights[start:start+12000],indices[start:start+12000]
        u,ur,uz,*_ = trial.evaluate(p,it)
        residual = differential(p,u,ur,uz,kappa)
        density = np.sum(abs(residual)**2,axis=(1,2))
        norm2 += float(w@np.sum(abs(u)**2,axis=(1,2)))
        residual2 += float(w@density)
        columns += np.einsum('p,pac,pac->c',w,residual.conj(),residual).real
        in_collar = np.minimum(np.hypot(p[:,0],p[:,1]-.5),np.hypot(p[:,0],p[:,1]+.5))<trial.a
        in_taper = np.linalg.norm(p,axis=1)>trial.taper_start
        for name,mask in (('collars',in_collar),('outer_taper',in_taper),('interior',~(in_collar|in_taper))):
            regional[name] += float(w[mask]@density[mask])
    return dict(order=order,points=len(points),weighted_volume=float(weights.sum()),
        omitted_zero_weight_nodes=omitted_zero_weight_nodes,
        minimum_barycentric=minimum_barycentric,minimum_mouth_radius=minimum_mouth_radius,
        trial_l2=math.sqrt(norm2),residual_l2=math.sqrt(residual2),
        residual_columns=np.sqrt(columns).tolist(),
        regional_residual_l2={k:math.sqrt(v) for k,v in regional.items()},
        regional_sum_of_squares_error=abs(sum(regional.values())-residual2),
        interpretation='Numerical quadrature of the actual conforming trial differential residual; not an enclosure.')


def run():
    started = time.perf_counter()
    paths = ('long_throat_source_results.json','feedback_residual_bound_results.json',
        'complete_angular_source.py','complete_angular_modes.py','exterior_boundary_map.py',
        'dirac_sewing_audit.py','two_mouth_mesh.py','spherical_weyl_leads.py','physical_exterior_quadrature.py')
    hashes = {p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
    prior = json.loads(Path(paths[0]).read_text(encoding='utf-8'))
    group = next(g for g in prior['geometries'] if g['geometry']['rho_c']==32.)
    saved = next(s for s in group['sectors'] if s['m']==0)
    sample = saved['samples'][1]
    matrices = sample['matrices']
    kappa = sample['kappa']
    print(json.dumps(dict(stage='construct_same_finite_model',config=saved['config'])),flush=True)
    sector = CompleteSector(.5,h=.05,R=2.,nmax=4,outer_nmax=9,L=40.,rho_c=32.)
    model = sector.model
    active = sector.active
    matrix = (sector.unshifted+kappa*sector.mass)[active][:,active].tocsc()
    lu = splu(matrix)
    gm,go = sector.Gm[active]@sector.Fin,sector.outer['G'][active]
    cm,co = sector.Cm[:,active],sector.outer['C'][:,active]
    z = _decode(matrices['Z_facet_right'])
    lifts = lu.solve(np.column_stack((gm,go)))
    xm,xo = lifts[:,:9],lifts[:,9:]
    poisson = xm+xo@z@np.linalg.solve(np.eye(len(z))-co@xo@z,co@xm)
    s = np.linalg.solve(sector.Fout,cm@poisson)
    original = _decode(matrices['S_connected'])[:9,:9]
    b = _decode(matrices['phases'][1]['B'])[:9,:9]
    derivative = _decode(matrices['phases'][1]['pressure'])[1,:9,:9]
    labels_in,labels_out = port_labels(.5,4)
    il = np.array([i for i,(_,n) in enumerate(labels_in) if n<=2])
    ol = np.array([i for i,(_,n) in enumerate(labels_out) if n<=2])
    ih,oh = np.setdiff1d(np.arange(9),il),np.setdiff1d(np.arange(9),ol)
    bl = b[np.ix_(il,ol)]
    u,sv,vh = np.linalg.svd(derivative[np.ix_(il,ol)],full_matrices=False)
    source = u*np.sqrt(sv)
    right_data = np.linalg.solve(np.eye(len(il))-bl@s[np.ix_(ol,il)],source)
    coefficients = np.zeros((9,len(sv)),complex)
    coefficients[il] = right_data
    reduced_solution = poisson@coefficients
    full = np.zeros((model.ndof,len(sv)),complex)
    full[active] = reduced_solution
    outgoing = np.linalg.solve(sector.Fout,cm@reduced_solution)
    # Define the round trial trace algebraically from its outgoing coefficients.
    incoming = np.zeros_like(coefficients)
    incoming[il] = source+bl@outgoing[ol]
    discrete_residual = matrix@reduced_solution-gm@coefficients-go@z@(co@reduced_solution)
    nodal,recovery = recovered_vertices(model,full)
    trial = Trial(model.vertices,model.triangles,nodal,incoming,outgoing)
    print(json.dumps(dict(stage='conforming_trial_constructed',triangles=len(model.triangles),
        vertices=len(model.vertices),same_saved_response_error=norm(s-original),
        elapsed_seconds=time.perf_counter()-started)),flush=True)
    checks = controls(trial,model,kappa,source,bl,il,ol)
    quadratures = []
    for order in (6,10):
        value = integrate_residual(trial,model,kappa,order)
        quadratures.append(value)
        print(json.dumps(dict(stage='differential_residual',**value)),flush=True)
    rho = quadratures[-1]['residual_l2']
    budget_result = json.loads(Path(paths[1]).read_text(encoding='utf-8'))
    budget_row = next(r for r in budget_result['records'] if r['phase']=='pi' and r['probe']==1)
    target = next(r for r in budget_row['continuum_residual_budgets'] if r['signal_fraction']==.01)
    w,_ = defects(b[np.ix_(ih,oh)])
    high_trace = np.linalg.solve(w,outgoing[oh])
    tests = dict(
        same_saved_finite_response=norm(s-original)<3e-12,
        discrete_equation=float(np.linalg.norm(discrete_residual)/np.linalg.norm(gm@coefficients))<3e-12,
        discrete_low_closure=norm(incoming-coefficients)<3e-12,
        recovery_axis=recovery['maximum_forbidden_axis_value_before_enforcing']<1e-14,
        analytic_physical_basis=checks['analytic_basis_error']<2e-12 and checks['analytic_basis_theta_derivative_error']<2e-7,
        round_boundary_closure=checks['physical_round_boundary_closure_error']<3e-12,
        round_characteristics=checks['round_characteristic_error']<3e-12,
        absent_higher_trace=checks['sampled_n5_to_n8_trace_projection']<3e-12,
        global_continuity=checks['shared_edge_continuity_error']<3e-12,
        analytic_field_derivatives=checks['field_cartesian_derivative_error']<3e-7 and checks['differential_residual_cartesian_check']<3e-7,
        pair_potential=checks['stable_pair_potential_error']<2e-14,
        compact_support_inside_outer_polygon=trial.taper_end<checks['minimum_outer_apothem'],
        finite_differential_residual=math.isfinite(rho) and rho>0,
        inputs_preserved=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in hashes.items()))
    tests = {k:bool(v) for k,v in tests.items()}
    result = dict(passed=all(tests.values()),tests=tests,
        parameters=dict(m=0,mu=.5,q=1,weyl='right',phase='pi',probe='upper_pressure',
            rho_c=32.,x=sample['x'],kappa=kappa,retained_N=2,trace_N=4,outer_J=9,
            source_columns=len(sv),b=.2,collar_a=.30,taper_start=1.5,taper_end=1.9),
        recovery=recovery,controls=checks,quadratures=quadratures,
        numerical_residual_relative_quadrature_change=abs(quadratures[0]['residual_l2']/rho-1),
        finite_solver=dict(same_saved_response_error=norm(s-original),
            algebraic_equation_relative_residual=float(np.linalg.norm(discrete_residual)/np.linalg.norm(gm@coefficients)),
            canonical_low_boundary_residual=norm(incoming-coefficients),
            source_svd_factorization_error=norm(derivative[np.ix_(il,ol)]-source@(vh.conj().T*np.sqrt(sv)).conj().T)),
        comparison=dict(complete_trial_weighted_high_trace_frobenius=float(np.linalg.norm(high_trace)),
            previous_all_sector_two_weyl_residual_budget=target['equal_continuum_residual_norm_budget'],
            actual_single_piece_residual_over_that_budget=rho/target['equal_continuum_residual_norm_budget'],
            residual_to_weighted_leg_coefficient=budget_result['best_collar']['residual_to_weighted_leg'],
            estimated_residual_error_majorant=budget_result['best_collar']['residual_to_weighted_leg']*rho,
            necessary_single_piece_target_met=rho<=target['equal_continuum_residual_norm_budget']),
        trial_representation=dict(vertices=model.vertices.tolist(),triangles=model.triangles.tolist(),
            recovered_nodal_spinors=_complex_json(nodal),canonical_incoming=_complex_json(incoming),
            canonical_outgoing=_complex_json(outgoing),low_source_X=_complex_json(source),
            description='Trial(...) reconstructs this continuous compactly supported field without another exterior solve.'),
        conclusions=dict(actual_conforming_trial_constructed=True,
            continuum_residual_numerically_estimated=True,continuum_residual_certified=False,
            full_paired_source_certified=False),
        scope=[
            'This is one right-Weyl m0 source-weighted five-column pilot, not the two-leg all-sector certificate.',
            'Vertex averaging removes DG value jumps; piecewise gradient jumps require no distributional delta in a first-order operator.',
            'Round mouth traces are imposed exactly by the analytic full blend and contain only n<=4; the physical boundary relation uses the specified saved finite core B.',
            'Actual pair potential, collar blending derivatives and outer taper derivatives are included in the differential residual.',
            'The field has covariant H1 axis regularity, not a claim of smooth even radial Taylor expansions.',
            'Physical round holes are removed from the original polygon mesh quadrature. Quadrature refinement remains an error estimate, not a rigorous enclosure.',
            'The prior all-sector residual target is only a necessary budget for this individual piece; meeting it would not certify the missing dual, other Weyl or other sectors.',
            'The saved core/radial maps and arithmetic are not independently interval-certified. Numerical residual estimates must not be substituted for rigorous upper bounds.'],
        seconds=time.perf_counter()-started,input_sha256=hashes,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_name('conforming_exterior_trial_results.json').write_text(
        json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,comparison=result['comparison'],
        numerical_residual_relative_quadrature_change=result['numerical_residual_relative_quadrature_change'],
        seconds=result['seconds'])),flush=True)
    if not result['passed']:
        raise AssertionError(tests)
    return result


if __name__=='__main__':
    run()
