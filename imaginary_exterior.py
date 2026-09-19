"""Connected two-mouth supporting-field response at E=i*kappa, kappa>0.

The frozen real-energy DG assembly is used at E=0 and shifted by kappa
times its weighted mass.  The exact discrete identity includes bulk
attenuation 2*kappa*||u||_rho^2.  No imaginary-axis squared amplitude is
interpreted as a probability, escape fraction, or vacuum stress.
"""
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.sparse import block_diag
from scipy.sparse.linalg import splu

from exterior_boundary_map import ExteriorDG, ID, _complex_json
from exterior_tail_transport import facet_map
from spherical_weyl_leads import outer_characteristics, mouth_characteristics, zero_mode
from two_mouth_mesh import mesh


def shifted_model(kappa, h=.08, R=2., q=1):
    """Shift an instance, never the frozen class or its source file."""
    kappa = float(kappa)
    if not math.isfinite(kappa) or kappa <= 0:
        raise ValueError('kappa must be finite and positive')
    model = ExteriorDG(mesh(b=.2, d=1., R=R, h=h), q=q, energy=0.)
    mass = block_diag([np.kron(cell['mass'], ID) for cell in model.cells], format='csc')
    model.matrix_full = model.matrix_full+kappa*mass
    model.matrix = model.matrix_full[model.active][:, model.active].tocsc()
    model.lu = splu(model.matrix)
    return model, mass


def build_problem(kappa=.4, h=.08, R=2., nmax=8, L=40.):
    """Return the actual linear problem for independently declared sewings.

    Keys model, ports, feedback, source, receiver expose the existing solve:
    u=model.solve(source,feedback); T=receiver@u.  This API does not choose
    a throat propagation, spin phase, occupancy, or determinant prescription.
    """
    from imaginary_weyl_leads import mouth_decaying_ratios, tail_map
    started = time.perf_counter()
    model, mass = shifted_model(kappa, h, R)
    b, d = .2, 1.

    def mouth_basis(port, points):
        theta = np.arctan2(points[:, 0], points[:, 1]-(d/2 if port == 'upper' else -d/2))
        plus, minus = mouth_characteristics(theta, b, port, nmax)
        zero = zero_mode(theta, b, port)[:, :, None]
        if port == 'upper':
            minus = np.concatenate((zero, minus), axis=2)
        else:
            plus = np.concatenate((zero, plus), axis=2)
        return plus, minus

    upper = model.modal_port('upper', lambda points: mouth_basis('upper', points))
    lower = model.modal_port('lower', lambda points: mouth_basis('lower', points))
    outer_basis = lambda points: outer_characteristics(.5,
                            np.arctan2(points[:, 0], points[:, 1]), R, nmax)
    outer = model.modal_port('outer', outer_basis)
    spherical_Z, tail_diagnostics = tail_map(kappa, R, L, nmax=nmax)
    Z, factors, geometry = facet_map(outer, outer_basis, spherical_Z)
    left, singular, right = np.linalg.svd(Z)
    beta = mouth_decaying_ratios(kappa, b, nmax)
    G = np.column_stack((upper['G'][:, 1:], lower['G'], outer['G']@left))
    C = np.vstack((upper['C'], lower['C'][1:], right@outer['C']))
    ratios = np.r_[beta, beta, singular]
    return dict(model=model, mass=mass, ports=[upper, lower, outer],
                feedback=(G, ratios, C), source=upper['G'][:, 0], receiver=lower['C'][0],
                kappa=kappa, h=h, R=R, nmax=nmax, L=L, beta=beta,
                spherical_Z=spherical_Z, facet_Z=Z, gram_factors=factors,
                diagnostics=dict(tail=tail_diagnostics, facet_coordinates=geometry,
                     svd_reconstruction_error=float(np.linalg.norm((left*singular)@right-Z)),
                     build_seconds=time.perf_counter()-started))


def attenuation_ledger(model, mass, kappa, solution, incoming, rhs):
    """Replace the zero-energy balance by its exact dissipative identity."""
    base = model.ledger(solution, incoming, rhs)
    bulk = float(np.vdot(solution, mass@solution).real)
    quadrature_bulk = 0.
    for cell in model.cells:
        it = cell['triangle']
        values = cell['basis']@solution[6*it:6*it+6].reshape(3, 2)
        quadrature_bulk += float(np.einsum('p,ps,ps->', cell['weights'], values.conj(), values).real)
    attenuation = 2*kappa*bulk
    # The frozen routine's undamped residual is expected to be -attenuation.
    undamped = base.pop('balance_residual')
    base.pop('identity_error')
    base.update(bulk_weighted_norm=bulk, bulk_norm_quadrature_difference=abs(bulk-quadrature_bulk),
                bulk_attenuation=attenuation, undamped_balance_residual=undamped,
                balance_with_attenuation_residual=undamped+attenuation,
                attenuation_identity_error=abs(undamped+attenuation-base['algebraic_balance_term']))
    return base


def response_case(kappa=.4, h=.08, R=2., nmax=8, L=40.):
    started = time.perf_counter()
    problem = build_problem(kappa, h, R, nmax, L)
    model, ports = problem['model'], problem['ports']
    solution = model.solve(problem['source'], problem['feedback'])
    upper, lower, outer = ports
    beta, Z = problem['beta'], problem['facet_Z']
    outgoing = [port['C']@solution for port in ports]
    incoming = [np.r_[1., beta*outgoing[0]], beta*outgoing[1][1:], Z@outgoing[2]]
    data = model.modal_incoming(ports, incoming)
    rhs = sum((p['G']@v for p,v in zip(ports,incoming)), np.zeros(model.ndof, complex))
    ledger = attenuation_ledger(model, problem['mass'], kappa, solution, data, rhs)
    omitted = {p['kind']: float(ledger['outgoing'][p['kind']]-np.vdot(v,v).real)
               for p,v in zip(ports,outgoing)}
    transmission = problem['receiver']@solution
    outer_net = float(np.vdot(outgoing[2],outgoing[2]).real-np.vdot(incoming[2],incoming[2]).real)
    raw_out = np.linalg.solve(problem['gram_factors'][0], outgoing[2])
    raw_in = problem['spherical_Z']@raw_out
    spherical_net = float(np.vdot(raw_out,raw_out).real-np.vdot(raw_in,raw_in).real)
    gram_delta = [f.conj().T@f-np.eye(len(Z)) for f in problem['gram_factors']]
    geometry_net = float((np.vdot(raw_out,gram_delta[0]@raw_out)-np.vdot(raw_in,gram_delta[1]@raw_in)).real)
    mouth_net = float(np.sum((1-abs(beta)**2)*(abs(outgoing[0])**2+abs(outgoing[1][1:])**2)))
    balance = (abs(transmission)**2+outer_net+mouth_net+sum(omitted.values())
               +ledger['numerical_loss']+ledger['bulk_attenuation']-1)
    return dict(kappa=kappa, h=h, R=R, L=L, nmax=nmax, triangles=len(model.triangles),
                active_dof=len(model.active), transmission=_complex_json(transmission),
                response_magnitude_squared=float(abs(transmission)**2),
                log_response_magnitude=float(np.log(abs(transmission))),
                facet_outer_net=outer_net, spherical_outer_net=spherical_net,
                signed_outer_geometry_difference=geometry_net,
                outer_geometry_identity_error=abs(outer_net-spherical_net-geometry_net),
                nonzero_mouth_net=mouth_net, omitted_outgoing_norm=omitted, ledger=ledger,
                reduced_attenuation_balance_residual=float(balance),
                outgoing_coefficients={p['kind']:_complex_json(v) for p,v in zip(ports,outgoing)},
                incoming_coefficients={p['kind']:_complex_json(v) for p,v in zip(ports,incoming)},
                diagnostics=problem['diagnostics'], seconds=time.perf_counter()-started,
                interpretation='Imaginary-energy response and norm identity; magnitude squared is NOT a probability.')


def free_green_case(kappa=.4, h=.08, R=2., nmax=12):
    """Independent modified-Helmholtz Green spinor with q=0 everywhere."""
    from imaginary_weyl_leads import outer_decaying_ratios
    model, mass = shifted_model(kappa, h, R, q=0)

    def exact(points):
        rho, z = points.T
        z = z-.5
        radius = np.hypot(rho,z)
        green = np.exp(-kappa*radius)/(4*np.pi*radius)
        derivative = green*(-kappa-1/radius)
        return np.sqrt(2*np.pi)*np.column_stack((1j*kappa*green-1j*z/radius*derivative,
                                                 -1j*rho/radius*derivative))

    outer_basis = lambda points: outer_characteristics(.5,
                                np.arctan2(points[:,0],points[:,1]), R, nmax)
    outer = model.modal_port('outer', outer_basis)
    raw_Z = np.diag(outer_decaying_ratios(kappa,R,.5,nmax))
    Z, factors, geometry = facet_map(outer,outer_basis,raw_Z)
    left,singular,right = np.linalg.svd(Z)
    rhs,data,_ = model.input(lambda kind,points,normal: exact(points) if kind in ('upper','lower')
                              else np.zeros((len(points),2),complex))
    solution = model.solve(rhs,(outer['G']@left,singular,right@outer['C']))
    reflected = Z@(outer['C']@solution)
    extra = model.modal_incoming([outer],[reflected])
    data = [a+b for a,b in zip(data,extra)]
    rhs = rhs+outer['G']@reflected
    ledger = attenuation_ledger(model,mass,kappa,solution,data,rhs)
    # Compare the complex outgoing coefficient vector independently to the
    # exact spinor sampled on precisely the same polygon quadrature.
    sampled = np.concatenate([exact(e['points']) for e in outer['edges']])
    exact_coef = np.einsum('p,psn,ps->n',outer['weights'],outer['plus'].conj(),sampled)
    numerical_coef = outer['C']@solution
    return dict(kappa=kappa,h=h,R=R,nmax=nmax,triangles=len(model.triangles),
                relative_l2_error=model.l2_error(solution,exact),
                relative_complex_outer_coefficient_error=float(np.linalg.norm(numerical_coef-exact_coef)/np.linalg.norm(exact_coef)),
                ledger=ledger,facet_coordinates=geometry,
                exact_source='(i*kappa+H)[exp(-kappa*r)/(4*pi*r)] spin_up, point inside upper excision; exact data on BOTH cuts.')


def random_attenuation_identity():
    model,mass = shifted_model(.4,.16,1.3)
    rng = np.random.default_rng(9201)
    data = [np.zeros((len(e['weights']),2),complex) for e in model.boundary]
    errors=[]
    for _ in range(3):
        u=np.zeros(model.ndof,complex)
        u[model.active]=rng.normal(size=len(model.active))+1j*rng.normal(size=len(model.active))
        result=attenuation_ledger(model,mass,.4,u,data,np.zeros(model.ndof,complex))
        scale=sum(result['outgoing'].values())+result['numerical_loss']+result['bulk_attenuation']
        errors.append(result['attenuation_identity_error']/scale)
    return dict(maximum_relative_error=max(errors),samples=3,seed=9201)


def amplitude(row):
    return complex(row['transmission']['real'],row['transmission']['imag'])


def checks():
    started=time.perf_counter()
    cache={}

    def case(kappa=.4,h=.05,R=2.,nmax=12,L=40.):
        key=(kappa,h,R,nmax,L)
        if key not in cache:
            cache[key]=response_case(*key)
            row=cache[key]
            print(f"kappa={kappa:g}, h={h:g}, R={R:g}, N={nmax}, L={L:g}: T={amplitude(row)}",flush=True)
        return cache[key]

    free=[free_green_case(h=h) for h in (.12,.08,.05)]
    meshes=[case(h=h) for h in (.08,.05,.035)]
    radii=[case(R=R) for R in (2.,4.)]
    endpoints=[case(L=L) for L in (20.,40.,80.)]
    angular=[case(nmax=n) for n in (8,12)]
    samples=[case(kappa=k) for k in (.1,.4,1.)]
    representative=case(h=.035,R=4.)
    random=random_attenuation_identity()
    rows=list(cache.values())
    ledgers=[r['ledger'] for r in rows+free]
    tests=dict(
        sparse_linear_residual=max(r['linear_residual_relative'] for r in ledgers)<1e-10,
        exact_bulk_attenuation_identity=max(r['attenuation_identity_error'] for r in ledgers)<1e-10,
        reduced_connected_attenuation_balance=max(abs(r['reduced_attenuation_balance_residual']) for r in rows)<1e-10,
        bulk_mass_matches_quadrature=max(r['bulk_norm_quadrature_difference'] for r in ledgers)<1e-11,
        positive_bulk_attenuation=all(r['bulk_attenuation']>0 for r in ledgers),
        attenuation_is_not_omitted=all(abs(r['undamped_balance_residual'])>.1*r['bulk_attenuation'] for r in ledgers),
        free_green_l2_refines=all(a['relative_l2_error']>b['relative_l2_error'] for a,b in zip(free[:-1],free[1:])),
        free_green_complex_projection_refines=all(a['relative_complex_outer_coefficient_error']>b['relative_complex_outer_coefficient_error'] for a,b in zip(free[:-1],free[1:])),
        mesh_numerical_loss_refines=all(a['ledger']['numerical_loss']>b['ledger']['numerical_loss'] for a,b in zip(meshes[:-1],meshes[1:])),
        signed_geometry_identity=max(r['outer_geometry_identity_error'] for r in rows)<1e-12,
        random_field_attenuation_identity=random['maximum_relative_error']<1e-12,
        omitted_modal_norm_nonnegative=all(min(r['omitted_outgoing_norm'].values())>-1e-11 for r in rows),
        real_static_imaginary_response=max(abs(amplitude(r).imag) for r in rows)<1e-11,
    )
    assert all(tests.values()),tests
    root=Path(__file__).resolve().parent
    inputs=('exterior_boundary_map.py','exterior_tail_transport.py','imaginary_weyl_leads.py',
            'spherical_weyl_leads.py','two_mouth_mesh.py')
    return dict(passed=True,tests=tests,parameters=dict(q=1,mu=.5,b=.2,d=1.,energy='i*kappa'),
        representative=representative,free_green_controls=free,mesh_refinement=meshes,
        radius_refinement=radii,far_endpoint_refinement=endpoints,angular_refinement=angular,
        kappa_samples=samples,random_identity=random,
        sensitivities=dict(matching_radius_change=abs(amplitude(radii[1])-amplitude(radii[0])),
            endpoint_changes=[abs(amplitude(b)-amplitude(a)) for a,b in zip(endpoints[:-1],endpoints[1:])],
            angular_change=abs(amplitude(angular[1])-amplitude(angular[0])),
            fine_mesh_change=abs(amplitude(meshes[-1])-amplitude(meshes[-2])),
            infinity_or_continuum_error_certified=False),
        equations=dict(shift='L(i*kappa)=L(0)+kappa*M_rho',
            identity='incoming=outgoing+DG_jump+incoming_mismatch+2*kappa*bulk_weighted_norm',
            connected='1=abs(T)^2+outer_net+nonzero_mouth_net+omitted_modes+DG_loss+bulk_attenuation',
            interpretation='The identity is elliptic attenuation, not a real-energy probability ledger.'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs},
        seconds=time.perf_counter()-started,
        scope=['Actual connected finite-domain spinor equation at positive imaginary energy on the specified flat q=1 field.',
          'Semi-infinite product mouth leads and finite magnetic tail with free decaying endpoint remain stated model assumptions.',
          'The same-sign polygon Gram conversion is retained, including its signed current discrepancy; no singular values clipped.',
          'Input zero-mode and output zero-mode phases are fixed by the existing angular helper; raw complex responses are saved.',
          'Magnitude squared at imaginary energy is NOT a return, transmission, or escape probability.',
          'No spectral determinant, vacuum occupation, renormalized stress, curved-throat match, physical material or backward passage is established.'])


if __name__=='__main__':
    result=checks()
    target=Path(__file__).with_name('imaginary_exterior_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({key:result[key] for key in ('passed','tests','seconds','source_sha256')},indent=2))
