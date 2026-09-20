"""Curved finite-exterior return sewn to the actual bridge, at one frequency.

The remote boundary has zero incoming characteristics: it is not a transparent
Melvin infinity. Polygonal cuts and P1 DG errors are varied, not certified.
The log-determinant response below is not an absolute vacuum stress.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import block_diag

from balanced_dihole import symbolic_fields
from dihole_bridge_scattering import BridgeMatrices, complex_record, embed
from dihole_curved_exterior import curved_model, attenuation_ledger, graded_mesh
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path
from spherical_weyl_leads import mouth_characteristics, zero_mode


def cut_basis(kind, points, size):
    center = .5 if kind == 'upper' else -.5
    theta = np.arctan2(points[:, 0], points[:, 1]-center)
    plus, minus = mouth_characteristics(theta, .05, kind, size)
    zero = zero_mode(theta, .05, kind)[:, :, None]
    if kind == 'upper':
        minus = np.concatenate((zero, minus), axis=2)
    else:
        plus = np.concatenate((zero, plus), axis=2)
        # Explicit lower spin lift U_L=exp(+i theta_bridge sigma_y/2).
        plus *= (-1.)**(np.arange(size+1)+1)
        minus *= (-1.)**(np.arange(1, size+1)+1)
    return plus, minus


def cut_factors(port, size):
    factors = []
    for index, projector in ((0, 'pp'), (1, 'pm')):
        raw = np.concatenate([
            np.einsum('ab,pbc->pac', edge[projector],
                      np.moveaxis(cut_basis(port['kind'], edge['points'], size)[index], 0, 1))
            for edge in port['edges']])
        gram = np.einsum('p,pai,paj->ij', port['weights'], raw.conj(), raw)
        factors.append(np.linalg.cholesky((gram+gram.conj().T)/2).conj().T)
    return factors


def sew(scattering, phase, size, incoming_factor, outgoing_factor):
    """Bridge [aL,bR]->[bL,aR]; exterior [aR,bL]->[bR,aL]."""
    value = scattering.copy()
    p, m = size+1, size
    value[m:, :p] *= np.exp(1j*phase)
    value[:m, p:] *= np.exp(-1j*phase)
    rows = np.r_[m+np.arange(p), np.arange(m)]
    columns = np.r_[p+np.arange(m), np.arange(p)]
    raw = value[np.ix_(rows, columns)]
    return np.linalg.solve(outgoing_factor.T, (incoming_factor@raw).T).T


def run_case(bridge, varied, zero, h, outer_radius, size=4, kappa=.1, graded=False):
    started = perf_counter()
    domain = graded_mesh(Rtarget=outer_radius, h=h) if graded else None
    model, mass, geometry = curved_model(kappa=kappa, h=h, R=outer_radius,
                                        triangulation=domain)
    ports = [model.modal_port(kind, lambda points, kind=kind:cut_basis(kind, points, size))
             for kind in ('upper', 'lower')]
    lift = np.column_stack([p['G'] for p in ports])
    receiver = np.vstack([p['C'] for p in ports])
    solution = model.solve(lift)
    returned = receiver @ solution
    factors = [cut_factors(p, size) for p in ports]
    fi = block_diag(factors[0][1], factors[1][1])
    fo = block_diag(factors[0][0], factors[1][0])
    # Input is the true finite cut zero profile, with facet normalization explicit.
    raw_input = np.r_[zero, np.zeros(size)]
    projected_input = fi @ raw_input
    projected_input /= np.linalg.norm(projected_input)
    field = solution @ projected_input
    incoming = [projected_input[:size+1], projected_input[size+1:]]
    data = model.modal_incoming(ports, incoming)
    rhs = lift @ projected_input
    ledger = attenuation_ledger(model, mass, kappa, field, data, rhs)
    raw_return = np.linalg.solve(fo, returned @ (fi @ raw_input))
    zero_return = zero.conj() @ raw_return[size:]
    loops = []
    for phase in (0., np.pi):
        b = sew(bridge, phase, size, fi, fo)
        bp, bm = [sew(s, phase, size, fi, fo) for s in (varied[1], varied[0])]
        db = (bp-bm)/2e-4
        loop = b @ returned
        denominator = np.eye(len(loop))-loop
        logdet = float(np.linalg.slogdet(denominator)[1])
        derivative = -float(np.real(np.trace(np.linalg.solve(denominator, db @ returned))))
        finite_difference = float((np.linalg.slogdet(np.eye(len(loop))-bp@returned)[1]
                                  -np.linalg.slogdet(np.eye(len(loop))-bm@returned)[1])/2e-4)
        loops.append(dict(additional_sewing_phase=phase, real_logdet=logdet,
                          logdet_metric_derivative=derivative,
                          logdet_finite_difference=finite_difference,
                          derivative_identity_error=abs(derivative-finite_difference),
                          loop_spectral_radius=float(max(abs(np.linalg.eigvals(loop)))),
                          denominator_condition=float(np.linalg.cond(denominator)),
                          bridge_facet_max_singular=float(np.linalg.svd(b, compute_uv=False)[0]),
                          bridge_map=complex_record(b)))
    return dict(h=h, outer_radius=outer_radius, graded=graded, kappa=kappa, angular_pairs=size,
                geometry=geometry, active_dof=len(model.active),
                exterior_max_singular=float(np.linalg.svd(returned, compute_uv=False)[0]),
                cut_gram_error=max(p['gram_error'] for p in ports),
                cut_gram_condition=max(p['gram_condition'] for p in ports),
                zero_return=complex_record(zero_return),
                exterior_return=complex_record(returned), loops=loops,
                input_zero_ledger=ledger, runtime_seconds=perf_counter()-started)


def run():
    dependencies = ('dihole_curved_exterior.py', 'dihole_bridge_scattering.py',
                    'dihole_cylinder_transfer.py', 'dihole_stretched_collar.py',
                    'balanced_dihole.py', 'exterior_boundary_map.py', 'two_mouth_mesh.py',
                    'spherical_weyl_leads.py', 'project_paths.py')
    hashes = {str(project_path(n)):hashlib.sha256(project_path(n).read_bytes()).hexdigest()
              for n in dependencies}
    function, _ = symbolic_fields()
    matrices = BridgeMatrices(function, exact_target_parameters(), 4, 48)
    bridge, _ = embed(matrices.coefficients, 5, 4, .1j)
    varied = [embed(lambda s,e=e:matrices.coefficients(s,e),5,4,.1j)[0] for e in (-1e-4,1e-4)]
    zero = matrices.cut_zero()[0]
    records = []
    for h, radius, graded in ((.04, 1.1, False), (.028, 1.1, False),
                              (.028, 1.6, False), (.028, 8., True), (.028, 16., True)):
        row=run_case(bridge, varied, zero, h, radius, graded=graded)
        records.append(row)
        print(json.dumps(dict(h=h,R=radius,return00=row['zero_return'],
                              derivative=[x['logdet_metric_derivative'] for x in row['loops']],
                              seconds=row['runtime_seconds'])),flush=True)
    # Physical/convergence diagnostics remain values; they are not forced green.
    tests=dict(finite_outputs=all(np.isfinite(x['logdet_metric_derivative']) for r in records for x in r['loops']),
               passive_finite_exterior=max(r['exterior_max_singular'] for r in records)<1+1e-8,
               passive_facet_bridge=max(x['bridge_facet_max_singular'] for r in records for x in r['loops'])<1+1e-8,
               facet_orthonormality=max(r['cut_gram_error'] for r in records)<2e-11,
               curved_flux_identity=max(r['input_zero_ledger']['attenuation_identity_error'] for r in records)<2e-10,
               physical_norm_quadrature=max(r['input_zero_ledger']['bulk_norm_quadrature_difference'] for r in records)<2e-10,
               actual_linear_residual=max(r['input_zero_ledger']['linear_residual_relative'] for r in records)<2e-10,
               determinant_derivative_identity=max(x['derivative_identity_error'] for r in records for x in r['loops'])<2e-6,
               frozen_inputs=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()))
    result=dict(passed=all(tests.values()),tests=tests,records=records,
                scope=['Actual curved coefficients, one q1 mu.5 Weyl block, one imaginary frequency.',
                       'Finite outer zero-incoming characteristic boundary is not transparent Melvin infinity.',
                       'Radial mesh grading extends the far domain while keeping the mouth grid fixed; distant cells are coarser and require their own sensitivity check.',
                       'Polygonal cuts use separate projected Gram factors; cross-sign facet and displacement errors remain.',
                       'Finite logdet derivative is an interaction diagnostic, not an absolute vacuum-energy or null-stress result.',
                       'Neither a quantum-support sign nor a continuum/geometric convergence claim follows from passing algebraic controls.',
                       'Additional phase0/pi is kept separate from the explicitly fixed lower spin-frame parity.'],
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashes)
    render(result)
    project_path('dihole_return_pilot_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests),indent=2))
    if not result['passed']: raise RuntimeError('Curved return pilot algebraic controls failed.')
    return result


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.6,4.4),layout='constrained')
    records=result['records']; x=np.arange(len(records))
    for j,label in enumerate(('additional phase 0','additional phase pi')):
        ax.plot(x,[r['loops'][j]['logdet_metric_derivative'] for r in records],'o-',label=label)
    ax.axhline(0,color='gray',lw=.7)
    ax.set_xticks(x,[f"h={r['h']}\nR={r['outer_radius']}"+ ('\ngraded' if r['graded'] else '') for r in records])
    ax.set_ylabel('Metric derivative of Re log det(I - bridge x return)')
    ax.set_title('Connected curved coefficients: finite-boundary sensitivity\nq=1, kappa D=0.1; no physical stress conclusion',fontsize=11)
    ax.legend();ax.grid(alpha=.2)
    for extension in ('png','pdf'):fig.savefig(project_path('dihole_return_pilot_progress.'+extension),dpi=180)
    plt.close(fig)


if __name__=='__main__':
    run()
