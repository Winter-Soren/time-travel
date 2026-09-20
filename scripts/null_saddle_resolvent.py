"""Whole-plane imaginary-frequency response of the local Dirac saddle.

One fixed Gaussian source, global Schwartz trial fields, and the entire
continuum residual. This is not a channel transmission experiment or a
quantum stress calculation. Width changes change the physical source.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.special import eval_hermite, roots_hermite

from null_saddle_leads import AU, AV, SY, _hermite_functions
from project_paths import project_path


def real_generator(count, eta=1., width=1.):
    """H=iK; order (Hermite u index, Hermite v index, global spin).

    The returned finite matrix is P_N K P_N. For a full continuum action
    on N-by-N coefficients, pad into N+1 per coordinate before applying K.
    """
    if (isinstance(count, bool) or int(count) != count or count < 2
            or not np.isfinite([eta, width]).all() or width <= 0):
        raise ValueError('Require integer count>=2, finite eta and width>0.')
    count = int(count)
    off = np.sqrt(np.arange(1, count)/2)
    position = sparse.diags([width*off, width*off], [-1, 1])
    derivative = sparse.diags([-off/width, off/width], [-1, 1])
    eye = sparse.eye(count)
    return (-sparse.kron(sparse.kron(derivative, eye), AU.real)
            -sparse.kron(sparse.kron(eye, derivative), AV.real)
            +sparse.kron(eta*sparse.eye(count*count)
                         -2*sparse.kron(position, position), SY.imag)).tocsc()


def full_residual(amplitude, kappa, eta=1., width=1., adjoint=False):
    """Real rho=f-(H-i*kappa)(i*a), or its adjoint counterpart.

    No artificial spatial boundary: all omitted Hermite coefficients,
    including the raised corner, are present in the returned array.
    """
    count = amplitude.shape[0]
    padded = np.zeros((count+1, count+1, 2))
    padded[:count, :count] = amplitude
    source = np.zeros(padded.size)
    source[0] = 1.
    action = real_generator(count+1, eta, width)@padded.ravel()
    sign = 1. if adjoint else -1.
    return (source+action+sign*kappa*padded.ravel()).reshape(padded.shape)


def solve_probe(count, kappa, eta=1., width=1.):
    """Solve a fixed normalized ground-Gaussian spin-up source and its dual."""
    if not math.isfinite(kappa) or kappa <= 0:
        raise ValueError('Require finite kappa>0; no real-axis limit is claimed.')
    started = time.perf_counter()
    generator = real_generator(count, eta, width)
    defect = generator+generator.T
    symmetry_error = float(np.max(abs(defect.data), initial=0.))
    source = np.zeros(2*count*count)
    source[0] = 1.
    factor = splu(generator-kappa*sparse.eye(len(source)))
    # psi=i*a and chi=i*b solve the original and adjoint equations.
    a = factor.solve(-source).reshape(count, count, 2)
    b = factor.solve(source, trans='T').reshape(count, count, 2)
    rho = full_residual(a, kappa, eta, width)
    sigma = full_residual(b, kappa, eta, width, adjoint=True)
    retained = rho[:count, :count]
    shell2 = np.sum(rho[count, :, :]**2)+np.sum(rho[:count, count, :]**2)
    rho_norm, sigma_norm = float(np.linalg.norm(rho)), float(np.linalg.norm(sigma))
    norm = float(np.linalg.norm(a))
    correction = -float(np.sum(b*retained))
    row = dict(count=count, kappa=kappa, eta=eta, width=width,
        real_degrees_of_freedom=len(source), seconds=time.perf_counter()-started,
        source='f=phi_0(u; width)*phi_0(v; width)*(1,0), normalized on R^2',
        solution='psi=i*a; adjoint solution chi=i*b; fixed global spin coordinates',
        generator_skew_error=symmetry_error, solution_norm=norm,
        imaginary_G=float(a[0, 0, 0]), imaginary_G_dual_correction=correction,
        corrected_imaginary_G=float(a[0, 0, 0])+correction,
        retained_residual_norm=float(np.linalg.norm(retained)),
        omitted_shell_residual_norm=float(math.sqrt(shell2)),
        full_residual_norm=rho_norm, adjoint_residual_norm=sigma_norm,
        numerical_solution_error_bound=rho_norm/kappa,
        numerical_G_disk_radius=rho_norm*sigma_norm/kappa,
        absorption_identity_error=abs(float(a[0, 0, 0])-kappa*norm*norm),
        residual_split_error=abs(rho_norm**2-float(np.linalg.norm(retained))**2-shell2))
    print(json.dumps({key: row[key] for key in ('count', 'kappa', 'imaginary_G',
        'full_residual_norm', 'numerical_G_disk_radius', 'seconds')}), flush=True)
    return row, a, b


def independent_control():
    """Direct polynomial differentiation and exact Gaussian quadrature control."""
    count, width, eta, kappa = 5, .8, -.7, .6
    rng = np.random.default_rng(421)
    amplitude = rng.normal(size=(count, count, 2))
    # Physicists' polynomial evaluation and derivative, independently of the
    # coefficient ladder and the existing normalized-function recurrence.
    nodes, weights = roots_hermite(count+3)
    basis, derivative = [], []
    for n in range(count+1):
        poly = eval_hermite(n, nodes)
        dp = 2*n*eval_hermite(n-1, nodes) if n else np.zeros_like(nodes)
        scale = math.sqrt(width*math.sqrt(math.pi)*2**n*math.factorial(n))
        basis.append(np.exp(-nodes*nodes/2)*poly/scale)
        derivative.append(np.exp(-nodes*nodes/2)*(dp-nodes*poly)/(width*scale))
    basis, derivative = np.array(basis).T, np.array(derivative).T
    field = np.einsum('pi,qj,ijs->pqs', basis[:, :count], basis[:, :count], amplitude)
    du = np.einsum('pi,qj,ijs->pqs', derivative[:, :count], basis[:, :count], amplitude)
    dv = np.einsum('pi,qj,ijs->pqs', basis[:, :count], derivative[:, :count], amplitude)
    mass = eta-2*width**2*np.outer(nodes, nodes)
    source = np.zeros_like(field)
    source[:, :, 0] = np.outer(basis[:, 0], basis[:, 0])
    direct = (source-du@AU.real.T-dv@AV.real.T
              +mass[:, :, None]*(field@SY.imag.T)-kappa*field)
    coefficients = full_residual(amplitude, kappa, eta, width)
    reconstructed = np.einsum('pi,qj,ijs->pqs', basis, basis, coefficients)
    values_error = float(np.max(abs(direct-reconstructed)))
    quadrature = width**2*np.outer(weights*np.exp(nodes**2), weights*np.exp(nodes**2))
    direct_norm2 = float(np.sum(quadrature*np.sum(direct**2, axis=2)))
    coefficient_norm2 = float(np.sum(coefficients**2))
    norm_error = abs(direct_norm2-coefficient_norm2)/coefficient_norm2
    assert values_error < 2e-13 and norm_error < 2e-13
    # A manufactured finite Schwartz field and its exact polynomial source.
    n = 8
    exact = np.zeros((n, n, 2)); exact[:3, :3] = rng.normal(size=(3, 3, 2))
    operator = real_generator(n, eta, width)-kappa*sparse.eye(2*n*n)
    recovered = splu(operator).solve(operator@exact.ravel())
    recovery_error = float(np.linalg.norm(recovered-exact.ravel())/np.linalg.norm(exact))
    assert recovery_error < 2e-13
    return dict(coordinate_operator_maximum_error=values_error,
                Gaussian_quadrature_relative_error=norm_error,
                manufactured_relative_error=recovery_error,
                raised_corner_norm=float(np.linalg.norm(coefficients[-1, -1])))


def draw(result, fields):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.1), layout='constrained')
    row = result['cases']['k05_n256']
    amplitude = fields['k05_n256_a']
    points = np.linspace(-6, 6, 241)
    basis = _hermite_functions(points, amplitude.shape[0]-1)
    values = np.einsum('pi,ijs,qj->pqs', basis, amplitude, basis, optimize=True)
    density = np.sum(values**2, axis=2)
    view = axes[0].imshow(density.T, extent=(-6, 6, -6, 6), origin='lower', cmap='magma')
    fig.colorbar(view, ax=axes[0], label=r'$|\psi(u,v)|^2$')
    axes[0].set(xlabel='u', ylabel='v', title='One smooth field, including the corners')
    for kappa, color in ((1., '#347f8f'), (.5, '#9b5bb5')):
        rows = [x for x in result['cases'].values() if x['kappa'] == kappa]
        axes[1].semilogy([x['count'] for x in rows], [x['numerical_solution_error_bound'] for x in rows],
                         'o-', color=color, label=fr'$\kappa={kappa:g}$')
    axes[1].set(xlabel='Hermite functions per coordinate', ylabel=r'Full residual bound $\|r\|/\kappa$',
                title='Continuum error, including omitted modes')
    axes[1].legend()
    axes[2].axis('off')
    axes[2].text(0, .96, 'What this establishes', weight='bold', fontsize=13, va='top')
    axes[2].text(0, .85, 'A whole-plane response to one Gaussian probe.\nNo overlapping exterior pieces or remote caps.\nThe full equation residual controls the error.', fontsize=11, va='top', linespacing=1.6)
    radius = row['interval_certificate']['G_disk_radius_upper']
    axes[2].text(0, .59, fr'$\kappa=0.5$: Im G = {row["imaginary_G"]:.9f}'+'\n'
                 +fr'Certified scalar error radius: {radius:.2e}'+'\nFor this fixed probe and ideal saddle model.', fontsize=11, va='top', linespacing=1.6)
    axes[2].text(0, .33, 'Still missing', weight='bold', fontsize=13, va='top')
    axes[2].text(0, .23, 'Curved-geometry matching and local-model error.\nThe full quantum stress and gravitational support.\nA stable throat or backward passage.', fontsize=11, va='top', linespacing=1.6)
    fig.suptitle('Whole-plane local Dirac model | A controlled response, not a time-travel medium', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('null_saddle_resolvent_progress.'+ext), dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    from null_saddle_residual_interval import enclose_residual, self_check
    control = independent_control()
    interval_control = self_check()
    cases, fields = {}, {}
    for name, kappa, count in [('k1_n32', 1., 32), ('k1_n64', 1., 64), ('k1_n128', 1., 128),
                              ('k05_n64', .5, 64), ('k05_n128', .5, 128), ('k05_n256', .5, 256),
                              ('k025_n128', .25, 128)]:
        cases[name], fields[name+'_a'], fields[name+'_b'] = solve_probe(count, kappa)
    for name in ('k1_n128', 'k05_n256'):
        started = time.perf_counter()
        row = cases[name]
        row['interval_certificate'] = enclose_residual(fields[name+'_a'], fields[name+'_b'], row['kappa'])
        print(json.dumps(dict(certified=name, seconds=time.perf_counter()-started,
            G_interval=row['interval_certificate']['imaginary_G_interval'])), flush=True)
    comparisons = []
    for first, second in (('k1_n32', 'k1_n64'), ('k1_n64', 'k1_n128'),
                          ('k05_n64', 'k05_n128'), ('k05_n128', 'k05_n256')):
        a, difference = fields[first+'_a'], fields[second+'_a'].copy()
        difference[:len(a), :len(a)] -= a
        comparisons.append(dict(first=first, second=second,
            field_L2_difference=float(np.linalg.norm(difference)),
            imaginary_G_difference=abs(cases[first]['imaginary_G']-cases[second]['imaginary_G']),
            summed_numerical_solution_error_bounds=sum(cases[key]['numerical_solution_error_bound'] for key in (first, second))))
    tests = dict(independent_coordinate_operator=control['coordinate_operator_maximum_error'] < 2e-13,
        independent_Gaussian_quadrature=control['Gaussian_quadrature_relative_error'] < 2e-13,
        manufactured_solution=control['manufactured_relative_error'] < 2e-13,
        real_generator_skew=all(row['generator_skew_error'] < 2e-13 for row in cases.values()),
        absorption_sign=all(row['absorption_identity_error'] < 2e-12 for row in cases.values()),
        full_residual_split=all(row['residual_split_error'] < 2e-12 for row in cases.values()),
        field_refinement_inside_residual_bounds=all(row['field_L2_difference'] <= row['summed_numerical_solution_error_bounds'] for row in comparisons),
        interval_helper=interval_control['passed'],
        kappa_one_certified_G_radius=cases['k1_n128']['interval_certificate']['G_disk_radius_upper'] < 4.4e-8,
        kappa_half_certified_G_radius=cases['k05_n256']['interval_certificate']['G_disk_radius_upper'] < 4.1e-6)
    assert all(tests.values()), tests
    result = dict(passed=True, tests=tests, control=control, interval_control=interval_control,
        cases=cases, comparisons=comparisons,
        model='Whole-plane H=-i AU d_u-i AV d_v+SY*(eta-2uv), eta=1, one fixed width-1 Gaussian spin-up source.',
        proof=dict(domain='Unique self-adjoint closure on R^2; smooth Hermitian mass, complete constant propagation metric.',
                   reference='https://arxiv.org/abs/1803.08106, Theorem 2.1(i)',
                   full_residual='rho=f-(H-i*kappa)psi_N, sigma=f-(H+i*kappa)chi_N; both evaluated on the entire plane.',
                   solution_bound='||psi-psi_N||<=||rho||/kappa.',
                   dual_bound='|G-G_N-<chi_N,rho>|<=||rho||||sigma||/kappa; no assumed exact Galerkin solve.',
                   real_symmetry='Entrywise complex conjugation sends H to -H; for real f uniqueness makes psi and G purely imaginary. H is self-adjoint.'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                      for name in ('null_saddle_leads.py', 'null_saddle_residual_interval.py', 'project_paths.py')},
        limits=['Local saddle extended to its entire plane, not the exact curved dihole exterior.',
                'Fixed Gaussian probe; no channel probability, pointwise Green kernel or vacuum stress.',
                'Kappa is strictly positive. This supplies no uniform real-axis or low-frequency limit.',
                'Kappa=.25, N128 has a very large residual bound; the low-frequency calculation is unresolved.'])
    np.savez_compressed(project_path('null_saddle_resolvent_fields.npz'), **fields)
    result['field_sha256'] = hashlib.sha256(project_path('null_saddle_resolvent_fields.npz').read_bytes()).hexdigest()
    project_path('null_saddle_resolvent_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result, fields)
    print(json.dumps(dict(passed=result['passed'], tests=tests, case_count=len(cases)), indent=2))
