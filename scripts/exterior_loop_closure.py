"""A declared zero-mode throat connection on the computed imaginary exterior.

The determinant is a finite DG boundary determinant ratio, not a renormalized
four-dimensional vacuum determinant. Throat length and endpoint phase are
inputs; curved nonzero throat modes, the other chirality and stress are absent.
"""
from project_paths import project_path
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from exterior_boundary_map import _complex_json


def closure_cases(kappa=.4, h=.08, R=2., nmax=8, L=40.):
    """Check direct boundary sewing against the scalar Schur complement."""
    from imaginary_exterior import build_problem
    problem = build_problem(kappa=kappa, h=h, R=R, nmax=nmax, L=L)
    model, source, receiver = (problem[key] for key in ('model', 'source', 'receiver'))
    G, ratios, C = problem['feedback']
    open_solution = model.solve(source, problem['feedback'])
    response = receiver@open_solution
    # Use the bare shifted LU to form the finite boundary determinant directly.
    bare_source = model.solve(source)
    lifted = model.solve(G*ratios[None, :])
    open_matrix = np.eye(len(ratios))-C@lifted
    sign_open, log_open = np.linalg.slogdet(open_matrix)
    if sign_open == 0:
        raise ArithmeticError('Singular reference boundary matrix')
    rows = []
    for length in (1., 4.):
        for phase in (0., np.pi):
            throat = np.exp(-kappa*length+1j*phase)
            determinant = 1-throat*response
            extended_G = np.column_stack((G, source))
            extended_C = np.vstack((C, receiver))
            extended_ratios = np.r_[ratios, throat]
            closed = model.solve(source, (extended_G, extended_ratios, extended_C))
            boundary_matrix = np.block([
                [open_matrix, -(C@bare_source)[:, None]*throat],
                [-(receiver@lifted)[None, :], np.array([[1-throat*(receiver@bare_source)]])]])
            sign_closed, log_closed = np.linalg.slogdet(boundary_matrix)
            direct_ratio = sign_closed/sign_open*np.exp(log_closed-log_open)
            residual = (model.matrix@closed[model.active]-source[model.active]
                        -(G@(ratios*(C@closed)))[model.active]
                        -throat*source[model.active]*(receiver@closed))
            rows.append(dict(optical_throat_length=length, endpoint_phase=phase,
                throat_factor=_complex_json(throat), determinant=_complex_json(determinant),
                direct_boundary_determinant_ratio=_complex_json(direct_ratio),
                determinant_identity_error=float(abs(direct_ratio-determinant)),
                closed_response=_complex_json(receiver@closed),
                scalar_response=_complex_json(response/determinant),
                response_identity_error=float(abs(receiver@closed-response/determinant)),
                full_field_identity_relative_error=float(np.linalg.norm(closed-open_solution/determinant)
                                                         /np.linalg.norm(closed)),
                sewn_equation_relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(source[model.active])),
                absolute_round_trip_factor=float(abs(throat*response)),
                log_absolute_determinant=float(np.log(abs(determinant)))))
    return dict(kappa=kappa, h=h, R=R, nmax=nmax, L=L,
                open_response=_complex_json(response), closures=rows)


def main():
    started = time.perf_counter()
    cases = [closure_cases(kappa=kappa) for kappa in (.1, .4, 1.)]
    rows = [row for case in cases for row in case['closures']]
    tests = dict(
        direct_boundary_determinant=max(r['determinant_identity_error'] for r in rows) < 1e-11,
        direct_sewn_response=max(r['response_identity_error'] for r in rows) < 1e-11,
        full_field_schur_complement=max(r['full_field_identity_relative_error'] for r in rows) < 1e-11,
        sewn_equation=max(r['sewn_equation_relative_residual'] for r in rows) < 1e-11,
        sampled_round_trip_below_one=max(r['absolute_round_trip_factor'] for r in rows) < 1,
        endpoint_phase_changes_determinant=all(abs(case['closures'][0]['log_absolute_determinant']-
                                                 case['closures'][1]['log_absolute_determinant']) > 1e-5
                                               for case in cases))
    root = Path(__file__).resolve().parent
    inputs = ('imaginary_exterior.py', 'imaginary_weyl_leads.py', 'exterior_tail_transport.py',
              'exterior_boundary_map.py', 'magnetic_tail_map.py', 'spherical_weyl_leads.py', 'two_mouth_mesh.py')
    result = dict(passed=all(tests.values()), tests=tests, cases=cases, seconds=time.perf_counter()-started,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs},
        equations=dict(response='T=c (A-G diag(beta) C)^(-1) g',
            throat='z=exp(i*phase-kappa*optical_length); upper_input=source+z*lower_output',
            determinant='det(M_sewn)/det(M_open)=1-z*T',
            driven_solution='u_sewn=u_open/(1-z*T)'),
        scope='Finite discretized boundary determinant ratio for one chosen zero-mode throat connection. '
              'Endpoint phases 0 and pi are explicit alternatives, not a derived physical spin/Wilson holonomy. '
              'The sampled round-trip check is not a global spectral stability proof. '
              'No four-dimensional determinant renormalization, second chirality, full curved throat, '
              'stress tensor, supporting material, or backward-time passage is established.',
        sources=['https://arxiv.org/pdf/1606.07834', 'https://arxiv.org/pdf/0711.2395',
                 'https://arxiv.org/html/1807.04726v3'])
    (project_path('exterior_loop_closure_results.json')).write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, seconds=result['seconds'],
        samples=[dict(kappa=c['kappa'], response=c['open_response'],
                      phase0_logdet=c['closures'][0]['log_absolute_determinant'],
                      phasepi_logdet=c['closures'][1]['log_absolute_determinant']) for c in cases]), indent=2))
    assert result['passed'], tests


if __name__ == '__main__':
    main()
