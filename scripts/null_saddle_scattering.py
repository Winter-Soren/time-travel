"""Two-dimensional local Dirac saddle with frozen Hermite lead matching.

P1 discontinuous Galerkin uses exact straight faces and reports numerical
jump loss. Finite-box/frozen leads are not the actual dihole exterior.
"""
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from exterior_boundary_map import ExteriorDG, SX, SY, SZ, ID, _quadrature
from project_paths import project_path

AU, AV = (SX+SZ)/math.sqrt(2), (SX-SZ)/math.sqrt(2)
FACES = ('u+', 'u-', 'v+', 'v-')


def square_mesh(radius, spacing):
    if not np.isfinite([radius, spacing]).all() or not 0 < spacing <= radius:
        raise ValueError('Require finite radius>=spacing>0.')
    count = 2*math.ceil(radius/spacing)
    axis = np.linspace(-radius, radius, count+1)
    uu, vv = np.meshgrid(axis, axis)
    vertices = np.column_stack([uu.ravel(), vv.ravel()])
    triangles = []
    for j in range(count):
        for i in range(count):
            bl = j*(count+1)+i
            br, tl, tr = bl+1, bl+count+1, bl+count+2
            if (i+j) % 2:
                triangles.extend([(bl, br, tl), (br, tr, tl)])
            else:
                triangles.extend([(bl, br, tr), (bl, tr, tl)])
    return vertices, np.asarray(triangles), 2*radius/count


class SaddleDG(ExteriorDG):
    """Reuse established lifts/solves/traces; replace cylindrical assembly."""
    def __init__(self, radius=2.5, spacing=.25, eta=0., energy=0., constant_mass=None):
        if not np.isfinite([eta, energy]).all() or (constant_mass is not None and not math.isfinite(constant_mass)):
            raise ValueError('Finite real coefficients required.')
        self.vertices, self.triangles, self.spacing = square_mesh(radius, spacing)
        self.radius, self.eta, self.energy = radius, eta, energy
        self.ndof = 6*len(self.triangles)
        self.boundary, self.interior, self.cells = [], [], []
        bary, qw, ex, ew = _quadrature(5)
        rows, cols, values, owners = [], [], [], {}

        def block(left, right, matrix):
            ii, jj = 6*left+np.arange(6), 6*right+np.arange(6)
            rows.extend(np.repeat(ii, 6)); cols.extend(np.tile(jj, 6))
            values.extend(matrix.ravel())

        for it, tri in enumerate(self.triangles):
            xy = self.vertices[tri]
            jac = np.linalg.det(np.stack([xy[1]-xy[0], xy[2]-xy[0]]))
            if jac <= 0:
                raise ValueError('Mesh orientation failed.')
            grad = np.linalg.inv(np.column_stack([np.ones(3), xy]))[1:].T
            points, weights = bary@xy, qw*jac
            mass = np.einsum('p,pa,pb->ab', weights, bary, bary)
            moment = np.einsum('p,pa->a', weights, bary)
            potential = eta-2*points[:, 0]*points[:, 1] if constant_mass is None else np.full(len(points), constant_mass)
            angular = np.einsum('p,pa,pb->ab', weights*potential, bary, bary)
            local = np.empty((3, 2, 3, 2), complex)
            for a in range(3):
                for b in range(3):
                    local[a, :, b, :] = moment[a]*(grad[b, 0]*AU+grad[b, 1]*AV)+1j*angular[a, b]*SY-1j*energy*mass[a, b]*ID
            block(it, it, local.reshape(6, 6))
            self.cells.append(dict(triangle=it, points=points, weights=weights, basis=bary, mass=mass))
            for a, b in ((0, 1), (1, 2), (2, 0)):
                owners.setdefault(tuple(sorted((int(tri[a]), int(tri[b])))), []).append((it, a, b))

        for key, cells in owners.items():
            if len(cells) not in (1, 2):
                raise ValueError('Nonmanifold mesh.')
            it, a, b = cells[0]
            start, end = self.vertices[self.triangles[it, [a, b]]]
            tangent = end-start
            length = np.linalg.norm(tangent)
            normal = np.array([tangent[1], -tangent[0]])/length
            points = (1-ex[:, None])*start+ex[:, None]*end
            weights = ew*length
            basis = np.zeros((len(ex), 3)); basis[:, a], basis[:, b] = 1-ex, ex
            an = normal[0]*AU+normal[1]*AV
            pm, pp = (ID-an)/2, (ID+an)/2
            face_mass = np.einsum('p,pa,pb->ab', weights, basis, basis)
            block(it, it, np.kron(face_mass, pm))
            edge = dict(triangle=it, points=points, weights=weights, basis=basis,
                        normal=normal, pm=pm, pp=pp, vertices=key)
            if len(cells) == 2:
                jt = cells[1][0]
                other = np.zeros_like(basis)
                for local_index, vertex in enumerate(self.triangles[it]):
                    found = np.where(self.triangles[jt] == vertex)[0]
                    if len(found):
                        other[:, found[0]] = basis[:, local_index]
                mixed = np.einsum('p,pa,pb->ab', weights, basis, other)
                other_mass = np.einsum('p,pa,pb->ab', weights, other, other)
                block(it, jt, -np.kron(mixed, pm))
                block(jt, jt, np.kron(other_mass, pp))
                block(jt, it, -np.kron(mixed.T, pp))
                edge.update(other=jt, other_basis=other)
                self.interior.append(edge)
            else:
                axis = int(np.argmax(abs(normal)))
                edge['kind'] = ('u' if axis == 0 else 'v')+('+' if normal[axis] > 0 else '-')
                assert np.max(abs(points[:, axis]-math.copysign(radius, normal[axis]))) < 1e-12
                self.boundary.append(edge)
        self.matrix_full = coo_matrix((values, (rows, cols)), shape=(self.ndof, self.ndof)).tocsc()
        self.active = np.arange(self.ndof)
        self.matrix = self.matrix_full
        self.lu = splu(self.matrix)

    def ledger(self, solution, incoming_data, rhs):
        outgoing = dict.fromkeys(FACES, 0.)
        total_in, mismatch, jump = 0., 0., 0.
        for edge, incoming in zip(self.boundary, incoming_data):
            trace = self.trace(solution, edge)
            norm = lambda v: float(np.einsum('p,ps,ps->', edge['weights'], v.conj(), v).real)
            outgoing[edge['kind']] += norm(trace@edge['pp'].T)
            total_in += norm(incoming)
            mismatch += norm(trace@edge['pm'].T-incoming)
        for edge in self.interior:
            j = edge['other']
            right = edge['other_basis']@solution[6*j:6*j+6].reshape(3, 2)
            difference = self.trace(solution, edge)-right
            jump += float(np.einsum('p,ps,ps->', edge['weights'], difference.conj(), difference).real)
        residual = self.matrix@solution-rhs
        balance = sum(outgoing.values())+jump+mismatch-total_in
        algebraic = float(2*np.vdot(solution, residual).real)
        return dict(outgoing=outgoing, full_incoming_norm=total_in, jump_loss=jump,
                    incoming_mismatch=mismatch, numerical_loss=jump+mismatch,
                    balance_residual=balance, identity_error=abs(balance-algebraic),
                    linear_residual_relative=float(np.linalg.norm(residual)/np.linalg.norm(rhs)))


def modal_port(model, lead):
    """Distinct measured Gram transforms; analytic lead map stays unchanged."""
    face = lead.metadata()['face']
    edges = [edge for edge in model.boundary if edge['kind'] == face]
    points = np.concatenate([edge['points'] for edge in edges])
    weights = np.concatenate([edge['weights'] for edge in edges])
    raw = lead.basis_function(points)
    bases, factors, grams = [], [], []
    for functions in raw:
        basis = np.moveaxis(functions, 0, 1)
        gram = np.einsum('p,pai,paj->ij', weights, basis.conj(), basis)
        gram = (gram+gram.conj().T)/2
        factor = np.linalg.cholesky(gram).conj().T
        basis = basis@np.linalg.inv(factor)
        bases.append(basis); factors.append(factor); grams.append(gram)
    plus, minus = bases
    rp, rm = factors
    c = np.zeros((plus.shape[2], model.ndof), complex)
    g = np.zeros((model.ndof, minus.shape[2]), complex)
    index = 0
    for edge in edges:
        stop = index+len(edge['weights']); dof = 6*edge['triangle']+np.arange(6)
        c[:, dof] += np.einsum('p,pa,psn->nas', edge['weights'], edge['basis'], plus[index:stop].conj()).reshape(len(c), 6)
        g[dof] += np.einsum('p,pa,psn->asn', edge['weights'], edge['basis'], minus[index:stop]).reshape(6, -1)
        index = stop
    return dict(kind=face, edges=edges, plus=plus, minus=minus, weights=weights, G=g, C=c,
                Z=rm@lead.Z@np.linalg.inv(rp), Rplus=rp, Rminus=rm,
                Gramplus=grams[0], Gramminus=grams[1],
                gram_deviation=max(float(np.linalg.norm(x-np.eye(len(x)), 2)) for x in grams),
                gram_condition=max(float(np.linalg.cond(x)) for x in grams))


def complex_json(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def solve_case(eta=0., radius=2.5, spacing=.25, pairs=3, save_field=False):
    from null_saddle_leads import frozen_lead
    started = time.perf_counter()
    model = SaddleDG(radius, spacing, eta)
    leads = [frozen_lead(face, radius, eta, pairs=pairs, energy=0.) for face in FACES]
    ports = [modal_port(model, lead) for lead in leads]
    g, c, z = np.column_stack([p['G'] for p in ports]), np.vstack([p['C'] for p in ports]), block_diag(*[p['Z'] for p in ports])
    incoming0 = np.zeros((g.shape[1], 2), complex)
    offset = 0
    for port in ports:
        if port['kind'] in ('v+', 'v-'):
            column = ('v+', 'v-').index(port['kind'])
            incoming0[offset:offset+len(port['Rminus']), column] = port['Rminus'][:, 0]
        offset += port['G'].shape[1]
    rhs = g@incoming0
    direct, lifted = model.solve(rhs), model.solve(g@z)
    feedback = np.eye(len(c))-c@lifted
    solution = direct+lifted@np.linalg.solve(feedback, c@direct)
    outgoing, incoming = c@solution, incoming0+z@(c@solution)
    total_rhs = g@incoming
    outputs, ledgers, omitted = [], [], []
    gram_correction = np.zeros(2)
    plus_offset = minus_offset = 0
    incoming_by_port = []
    for port in ports:
        count_plus, count_minus = port['C'].shape[0], port['G'].shape[1]
        analytic = np.linalg.solve(port['Rplus'], outgoing[plus_offset:plus_offset+count_plus])
        analytic_in = np.linalg.solve(port['Rminus'], incoming[minus_offset:minus_offset+count_minus])
        for sign, coef, gram in ((1, analytic, port['Gramplus']), (-1, analytic_in, port['Gramminus'])):
            gram_correction += sign*np.einsum('ij,ik,kj->j', coef.conj(), gram-np.eye(len(gram)), coef).real
        if port['kind'] in ('u+', 'u-'):
            outputs.append(analytic[0])
        incoming_by_port.append(incoming[minus_offset:minus_offset+count_minus])
        plus_offset += count_plus; minus_offset += count_minus
    scattering = np.asarray(outputs)
    for column in range(2):
        data = model.modal_incoming(ports, [x[:, column] for x in incoming_by_port])
        ledger = model.ledger(solution[:, column], data, total_rhs[:, column])
        unresolved = sum(ledger['outgoing'].values())-float(np.vdot(outgoing[:, column], outgoing[:, column]).real)
        ledgers.append(ledger); omitted.append(unresolved)
    result = dict(eta=eta, radius=radius, requested_spacing=spacing, actual_spacing=model.spacing,
                  pairs=pairs, energy=0., elements=len(model.triangles), dofs=model.ndof,
                  scattering=complex_json(scattering), outgoing_zero_powers=(abs(scattering)**2).tolist(),
                  zero_power_sum=np.sum(abs(scattering)**2, axis=0).tolist(),
                  numerical_unitarity_defect=float(np.linalg.norm(scattering.conj().T@scattering-np.eye(2), 2)),
                  ledger=ledgers, omitted_outgoing_norm=omitted,
                  signed_finite_face_Gram_term=gram_correction.tolist(),
                  zero_power_ledger_error=[float(abs(1-np.sum(abs(scattering[:, j])**2)-ledgers[j]['numerical_loss']-omitted[j]-gram_correction[j])) for j in range(2)],
                  maximum_face_Gram_deviation=max(p['gram_deviation'] for p in ports),
                  feedback_condition=float(np.linalg.cond(feedback)),
                  seconds=time.perf_counter()-started,
                  scope='Finite box and frozen linear-mass Hermite leads; raw amplitudes retain DG and omitted-mode losses. Not the full dihole return matrix.')
    if save_field:
        np.savez_compressed(project_path('null_saddle_fields.npz'), vertices=model.vertices,
                            triangles=model.triangles, solution=solution, eta=eta, radius=radius, spacing=model.spacing)
    return result


def plane_wave_control(spacing):
    mass, ku, kv = .7, .8, .3
    energy = math.sqrt(mass*mass+ku*ku+kv*kv)
    values, vectors = np.linalg.eigh(ku*AU+kv*AV+mass*SY)
    spinor = vectors[:, -1]
    assert abs(values[-1]-energy) < 1e-14
    exact = lambda points: np.exp(1j*(ku*points[:, 0]+kv*points[:, 1]))[:, None]*spinor
    model = SaddleDG(1., spacing, energy=energy, constant_mass=mass)
    rhs, data, _ = model.input(lambda kind, points, normal: exact(points))
    solution = model.solve(rhs)
    return dict(spacing=model.spacing, l2_error=model.l2_error(solution, exact),
                ledger=model.ledger(solution, data, rhs))


def checks():
    from null_saddle_leads import self_check
    from scipy import __version__ as scipy_version
    lead_checks = self_check()
    plane = [plane_wave_control(h) for h in (.25, .125)]
    cases = {}
    # ponytail: bounded mesh/face/pair comparisons, not a continuum certificate.
    settings = [('coarse', 1., 2.5, .25, 3), ('base', 1., 2.5, .125, 3),
                ('fine', 1., 2.5, .0625, 3), ('reflected', -1., 2.5, .125, 3),
                ('symmetric', 0., 2.5, .125, 3), ('larger_box', 1., 3., .0625, 3),
                ('extra_pair', 1., 2.5, .0625, 4)]
    for name, eta, radius, spacing, pairs in settings:
        cases[name] = solve_case(eta, radius, spacing, pairs, save_field=name == 'fine')
        print(json.dumps(dict(case=name, powers=cases[name]['outgoing_zero_powers'],
                              seconds=cases[name]['seconds'])), flush=True)
    matrix = lambda name: np.asarray(cases[name]['scattering']['real'])+1j*np.asarray(cases[name]['scattering']['imag'])
    power = lambda name: np.asarray(cases[name]['outgoing_zero_powers'])
    symmetry = float(np.max(abs(matrix('reflected')-SX@matrix('base')@np.diag([-1., 1.]))))
    equal_split = float(np.max(abs(power('symmetric')[0]-power('symmetric')[1])))
    rows = list(cases.values())
    ledger_error = max(x['identity_error'] for row in rows for x in row['ledger'])
    decomposition = max(max(row['zero_power_ledger_error']) for row in rows)
    linear = max(x['linear_residual_relative'] for row in rows for x in row['ledger'])
    losses = [cases[name]['ledger'][0]['numerical_loss'] for name in ('coarse', 'base', 'fine')]
    changes = {name: float(np.max(abs(power(name)-power('fine')))) for name in ('base', 'larger_box', 'extra_pair')}
    tests = dict(lead_operator=lead_checks['passed'], manufactured_refinement=plane[1]['l2_error'] < plane[0]['l2_error']/3,
                 flux_identity=ledger_error < 2e-11, raw_zero_power_ledger=decomposition < 2e-11,
                 solve_residual=linear < 2e-11, exact_reflection=symmetry < 2e-11,
                 equal_split=equal_split < 2e-11, decreasing_DG_loss=losses[2] < losses[1] < losses[0],
                 outgoing_projection_positive=all(min(row['omitted_outgoing_norm']) > -2e-11 for row in rows))
    tests = {name: bool(passed) for name, passed in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests, lead_checks=lead_checks, plane_wave=plane,
                  cases=cases, comparisons=dict(reflection_amplitude_error=symmetry, equal_power_error=equal_split,
                  maximum_flux_identity_error=ledger_error, maximum_zero_power_ledger_error=decomposition,
                  maximum_linear_residual=linear, mesh_DG_losses=losses, maximum_power_changes_from_fine=changes),
                  convention=dict(equation='H=-i AU d_u-i AV d_v+SY(eta-2uv), E=0',
                                  incoming=['v+', 'v-'], outgoing=['u+', 'u-'], full_line_incoming_amplitude=1),
                  limits=['Finite square with frozen linear-mass leads, not the bilinear infinite-plane scattering matrix.',
                          'Mesh, radius and pair changes are measured sensitivities, not certified error bars.',
                          'No extrapolation to q=1, global curved-mouth return, vacuum stress, backreaction or time travel.'],
                  versions=dict(numpy=np.__version__, scipy=scipy_version),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in
                                ('null_saddle_leads.py', 'exterior_boundary_map.py', 'two_mouth_mesh.py', 'project_paths.py',
                                 'dihole_null_scaling.py', 'dihole_null_scaling_results.json')},
                  field_sha256=hashlib.sha256(project_path('null_saddle_fields.npz').read_bytes()).hexdigest())
    if not result['passed']:
        raise AssertionError(json.dumps(dict(tests=tests, comparisons=result['comparisons']), indent=2))
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from matplotlib.colors import PowerNorm
    with np.load(project_path('null_saddle_fields.npz')) as field:
        xy, triangles = field['vertices'], field['triangles']
        psi = field['solution'][:, 0].reshape(-1, 3, 2).mean(axis=1)
    centers = xy[triangles].mean(axis=1)
    density = np.sum(abs(psi)**2, axis=1)
    current = np.column_stack([np.einsum('pa,ab,pb->p', psi.conj(), a, psi).real for a in (AU, AV)])
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.8), layout='constrained')
    ax = axes[0]
    im = ax.tripcolor(Triangulation(xy[:, 0], xy[:, 1], triangles), facecolors=density,
                      cmap='magma', norm=PowerNorm(.45), rasterized=True)
    take = np.arange(0, len(psi), 36)
    take = take[density[take] > .035*np.max(density)]
    direction = current[take]/np.maximum(np.linalg.norm(current[take], axis=1)[:, None], 1e-30)
    ax.quiver(centers[take, 0], centers[take, 1], direction[:, 0], direction[:, 1], color='white',
              scale=23, width=.005, alpha=.65)
    ax.set(xlabel='u (local units)', ylabel='v (local units)', aspect='equal', title=r'Computed field: $\eta=1$, input from top')
    fig.colorbar(im, ax=ax, shrink=.72, label=r'$|\psi|^2$ (unit input flux)')
    ax = axes[1]
    eta = [-1., 0., 1.]
    for index, label, color in ((0, 'Right output', '#257e80'), (1, 'Left output', '#b37a25')):
        powers = [result['cases'][name]['outgoing_zero_powers'][index][0] for name in ('reflected', 'symmetric', 'base')]
        ax.plot(eta, powers, 'o--', color=color, label=label)
    ax.set(xlabel=r'Local channel offset $\eta$', ylabel='Raw outgoing zero-channel power', ylim=(0, 1.05),
           xticks=eta, title='Quantum routing in the local model')
    ax.legend(fontsize=9)
    ax.text(.01, -.25, 'Three calculated points; dashed lines guide the eye.\nFinite mesh and frozen outer leads.', transform=ax.transAxes, fontsize=9)
    ax = axes[2]
    rows = [result['cases'][name] for name in ('coarse', 'base', 'fine')]
    h = [row['actual_spacing'] for row in rows]
    ax.loglog(h, [row['ledger'][0]['numerical_loss'] for row in rows], 'o-', label='DG numerical loss', color='#536ab3')
    ax.loglog(h, [row['omitted_outgoing_norm'][0] for row in rows], 's-', label='Omitted outgoing modes', color='#b37a25')
    ax.set(xlabel='Mesh spacing h (smaller is finer)', ylabel='Flux deficit contribution', title='Errors measured, not normalized away')
    ax.invert_xaxis(); ax.legend(fontsize=9)
    ax.text(.01, -.25, 'Separate box-size and lead-pair checks saved.\nNo full-geometry or time-travel claim.', transform=ax.transAxes, fontsize=9)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Local quantum mixing near a magnetic null | Two-dimensional Dirac calculation', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('null_saddle_scattering_progress.'+ext), dpi=175)
    plt.close(fig)


if __name__ == '__main__':
    result = checks()
    project_path('null_saddle_scattering_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result)
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'comparisons')}, indent=2))
