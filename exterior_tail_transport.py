"""Connect the frozen two-mouth DG benchmark to its magnetic radial tail.

This is supporting-field transport in a prescribed flat q=1 background, not
a material trajectory, self-consistent throat, or chronology-changing medium.
The polygonal same-sign characteristic projection remains a mesh approximation.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from exterior_boundary_map import ExteriorDG, _complex_json
from spherical_weyl_leads import (mouth_characteristics, mouth_outgoing_ratios,
                                  outer_characteristics, outgoing_ratios, zero_mode)
from two_mouth_mesh import mesh


def facet_map(port, basis_function, spherical_map):
    """Convert a spherical map to separately orthonormalized facet coordinates.

    If B_s is the same-sign projected raw basis and B_s^* W B_s=L_s L_s^*,
    U_s=B_s(L_s^*)^-1. Thus Z_facet=L_-^* Z (L_+^*)^-1.
    This does not include cross-sign facet projections or radial displacement
    of chord quadrature points; both geometry effects require h refinement.
    """
    factors = []
    for index, projector in ((0, 'pp'), (1, 'pm')):
        raw = np.concatenate([
            np.einsum('ab,pbc->pac', edge[projector],
                      np.moveaxis(basis_function(edge['points'])[index], 0, 1))
            for edge in port['edges']])
        gram = np.einsum('p,pai,paj->ij', port['weights'], raw.conj(), raw)
        factors.append(np.linalg.cholesky((gram+gram.conj().T)/2).conj().T)
    converted = np.linalg.solve(factors[0].T, (factors[1]@spherical_map).T).T
    # No singular-value clipping: geometry can introduce a passivity defect.
    return converted, factors, dict(
        coordinate_correction_norm=float(np.linalg.norm(converted-spherical_map, 2)),
        minimum_flux_eigenvalue=float(np.linalg.eigvalsh(
            np.eye(len(converted))-converted.conj().T@converted)[0]))


def transport_case(h=.05, R=2., nmax=12, energy=.4, L=160.,
                   free_outer=False, convert_coordinates=True):
    """One unit upper zero-mode input; save its full complex modal response."""
    from magnetic_tail_map import tail_map
    started = time.perf_counter()
    b, d = .2, 1.
    model = ExteriorDG(mesh(b=b, d=d, R=R, h=h), energy=energy)

    def mouth_basis(port, points):
        theta = np.arctan2(points[:, 0], points[:, 1]-(d/2 if port == 'upper' else -d/2))
        plus, minus = mouth_characteristics(theta, b, port, nmax)
        zero = zero_mode(theta, b, port)[:, :, None]
        if port == 'upper':
            minus = np.concatenate((zero, minus), axis=2)
        else:
            plus = np.concatenate((zero, plus), axis=2)
        return plus, minus

    outer_basis = lambda points: outer_characteristics(.5,
        np.arctan2(points[:, 0], points[:, 1]), R, nmax)
    upper = model.modal_port('upper', lambda points: mouth_basis('upper', points))
    lower = model.modal_port('lower', lambda points: mouth_basis('lower', points))
    outer = model.modal_port('outer', outer_basis)
    ports = [upper, lower, outer]
    beta = mouth_outgoing_ratios(energy, b, nmax)
    if free_outer:
        Z = np.diag(outgoing_ratios(energy, R, .5, nmax))
        tail_diagnostics = dict(omits_magnetic_field_beyond_R=True)
    else:
        Z, tail_diagnostics = tail_map(energy, R, L, nmax=nmax)
    spherical_Z = Z
    if convert_coordinates:
        Z, factors, coordinate_diagnostics = facet_map(outer, outer_basis, Z)
    else:
        factors = [np.eye(len(Z)), np.eye(len(Z))]
        coordinate_diagnostics = dict(coordinate_conversion=False)
    left, singular, right = np.linalg.svd(Z)
    svd_error = float(np.linalg.norm((left*singular)@right-Z))
    G = np.column_stack((upper['G'][:, 1:], lower['G'], outer['G']@left))
    C = np.vstack((upper['C'], lower['C'][1:], right@outer['C']))
    ratios = np.r_[beta, beta, singular]
    solution = model.solve(upper['G'][:, 0], (G, ratios, C))
    outgoing = [port['C']@solution for port in ports]
    incoming = [np.r_[1., beta*outgoing[0]], beta*outgoing[1][1:], Z@outgoing[2]]
    data = model.modal_incoming(ports, incoming)
    total_rhs = sum((port['G']@inc for port, inc in zip(ports, incoming)),
                    np.zeros(model.ndof, complex))
    ledger = model.ledger(solution, data, total_rhs)
    omitted = {port['kind']: float(ledger['outgoing'][port['kind']]-np.vdot(coef, coef).real)
               for port, coef in zip(ports, outgoing)}
    transmission = outgoing[1][0]
    outer_flux = float(np.vdot(outgoing[2], outgoing[2]).real-
                       np.vdot(incoming[2], incoming[2]).real)
    spherical_outgoing = np.linalg.solve(factors[0], outgoing[2])
    spherical_incoming = spherical_Z@spherical_outgoing
    spherical_flux = float(np.vdot(spherical_outgoing, spherical_outgoing).real-
                           np.vdot(spherical_incoming, spherical_incoming).real)
    gram_delta = [factor.conj().T@factor-np.eye(len(Z)) for factor in factors]
    geometry_flux = float((np.vdot(spherical_outgoing, gram_delta[0]@spherical_outgoing)-
                           np.vdot(spherical_incoming, gram_delta[1]@spherical_incoming)).real)
    mouth_flux = float(np.sum((1-abs(beta)**2)*(abs(outgoing[0])**2+abs(outgoing[1][1:])**2)))
    return dict(h=h, R=R, L=L, nmax=nmax, energy=energy,
        free_outer=free_outer, convert_coordinates=convert_coordinates,
        triangles=len(model.triangles), lower_zero_amplitude=_complex_json(transmission),
        lower_zero_flux=float(abs(transmission)**2), retained_tail_net_flux=outer_flux,
        spherical_tail_net_flux=spherical_flux, signed_tail_geometry_flux_difference=geometry_flux,
        tail_geometry_ledger_error=abs(outer_flux-spherical_flux-geometry_flux),
        mouth_nonzero_net_flux=mouth_flux, omitted_outgoing_norm=omitted, ledger=ledger,
        reduced_balance_residual=float(abs(transmission)**2+outer_flux+mouth_flux+
                                        sum(omitted.values())+ledger['numerical_loss']-1),
        outgoing_coefficients={p['kind']: _complex_json(v) for p, v in zip(ports, outgoing)},
        incoming_coefficients={p['kind']: _complex_json(v) for p, v in zip(ports, incoming)},
        tail=tail_diagnostics, facet_coordinates=coordinate_diagnostics,
        svd_reconstruction_error=svd_error, seconds=time.perf_counter()-started)


def amplitude(case):
    value = case['lower_zero_amplitude']
    return complex(value['real'], value['imag'])


def main():
    started = time.perf_counter()
    cache = {}

    def case(**kwargs):
        key = tuple(sorted(kwargs.items()))
        if key not in cache:
            cache[key] = transport_case(**kwargs)
            row = cache[key]
            print(f"R={row['R']:g} h={row['h']:g} N={row['nmax']} L={row['L']:g}: "
                  f"P={row['lower_zero_flux']:.10f}; loss={row['ledger']['numerical_loss']:.3g}",
                  flush=True)
        return cache[key]

    # Independent refinements, sharing identical cases rather than rerunning them.
    radius = [case(R=R, h=.05, nmax=12, L=160.) for R in (2., 4., 8.)]
    far_radius = [case(R=2., h=.05, nmax=12, L=L) for L in (40., 80., 160.)]
    angular = [case(R=2., h=.05, nmax=N, L=160.) for N in (8, 12)]
    meshes = {str(R): [case(R=R, h=h, nmax=12, L=160.) for h in (.08, .05, .035)]
              for R in (2., 4.)}
    reference = json.loads(Path('exterior_boundary_map_results.json').read_text(encoding='utf-8'))
    free_control = case(R=2., h=.05, nmax=12, L=160., free_outer=True,
                        convert_coordinates=False)
    old = reference['radius_refinement'][0]
    assert old['R'] == 2 and old['h'] == .05 and old['nmax'] == 12
    cases = list(cache.values())
    radius_change = [abs(amplitude(a)-amplitude(b)) for a, b in zip(radius, radius[1:])]
    old_radius = reference['radius_refinement']
    old_change = [abs(amplitude(a)-amplitude(b)) for a, b in zip(old_radius, old_radius[1:])]
    L_change = [abs(amplitude(a)-amplitude(b)) for a, b in zip(far_radius, far_radius[1:])]
    mesh_change = {R: abs(amplitude(rows[-1])-amplitude(rows[-2])) for R, rows in meshes.items()}
    tests = dict(
        frozen_free_map_reproduced=abs(amplitude(free_control)-amplitude(old)) < 2e-11,
        full_matrix_svd=max(row['svd_reconstruction_error'] for row in cases) < 1e-12,
        spherical_plus_geometry_current=max(row['tail_geometry_ledger_error'] for row in cases) < 1e-12,
        current_identity=max(row['ledger']['identity_error'] for row in cases) < 1e-10,
        reduced_flux_balance=max(abs(row['reduced_balance_residual']) for row in cases) < 1e-10,
        linear_residual=max(row['ledger']['linear_residual_relative'] for row in cases) < 1e-10,
        modal_remainders_nonnegative=min(min(row['omitted_outgoing_norm'].values()) for row in cases) > -1e-10,
        positive_excited_tail_flux=min(row['retained_tail_net_flux'] for row in cases) > 0,
        mesh_dissipation_decreases=all(all(a['ledger']['numerical_loss'] > b['ledger']['numerical_loss']
                                         for a, b in zip(rows, rows[1:])) for rows in meshes.values()))
    convergence = dict(complex_matching_radius_changes=radius_change,
        previous_free_matching_radius_changes=old_change, far_endpoint_changes=L_change,
        angular_change=abs(amplitude(angular[-1])-amplitude(angular[0])),
        mesh_changes=mesh_change, infinity_convergence_certified=False,
        note='Finite differences diagnose sensitivity; they are not rigorous truncation bounds.')
    root = Path(__file__).resolve().parent
    inputs = ('magnetic_tail_map.py', 'exterior_boundary_map.py', 'spherical_weyl_leads.py',
              'two_mouth_mesh.py', 'exterior_boundary_map_results.json')
    result = dict(passed=all(tests.values()), tests=tests, parameters=dict(q=1, mu=.5, b=.2, d=1., energy=.4),
        representative=meshes['4.0'][-1], radius_refinement=radius, far_endpoint_refinement=far_radius,
        angular_refinement=angular, mesh_refinement=meshes, frozen_free_control=free_control,
        convergence=convergence, seconds=time.perf_counter()-started,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs},
        scope='Prescribed flat q=1 field and product mouth leads. Magnetic tail kept only to finite L. '
              'Same-sign polygon projection with correct Gram coordinates; cross-sign projection and chord '
              'displacement are mesh approximations. Not the large-q target, quantum stress, stable medium, '
              'material traversal, new physical law, or backward time travel.')
    (root/'exterior_tail_transport_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, convergence=convergence,
                         seconds=result['seconds']), indent=2))
    assert result['passed'], tests


if __name__ == '__main__':
    main()
