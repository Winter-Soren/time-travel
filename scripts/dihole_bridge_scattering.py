"""Finite-mode scattering through the actual prescribed two-collar bridge.

One q=1, mu=1/2 Weyl sector is retained.  The canonical flux variable obeys
[-i J d_s+T(s)]w=E W(s)w in a fixed regular angular basis.  Both collars use
the actual dihole fields; no exterior return or vacuum stress is supplied.
The four-block invariant embedding never forms a growing transfer matrix.
"""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import solve_ivp

from balanced_dihole import symbolic_fields
from curved_throat_modes import _constant_exact
from dihole_cylinder_transfer import angular_matrices
from dihole_stretched_collar import (exact_target_parameters, gauss,
                                     log_speed_derivatives, radius_profile,
                                     smooth_step)
from project_paths import project_path


def complex_record(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def collar_mask(s):
    """Compact paired lapse/radial variation, identical on the two collars."""
    u = abs(float(s))-.05
    if u <= 0 or u >= .05:
        return 0.
    return float(smooth_step(u/.05)[1]/2)


class BridgeMatrices:
    def __init__(self, function, parameters, size, order):
        self.function, self.parameters = function, parameters
        self.size, self.order = size, order
        self.plus, self.minus = size+1, size
        self.gram_error = 0.

    @lru_cache(maxsize=4096)
    def at_radius(self, radius):
        transverse, clock, chirality, gram = angular_matrices(
            self.function, self.parameters, self.size, self.order, radius)
        self.gram_error = max(self.gram_error, gram)
        p = self.plus
        return clock[:p, :p], clock[p:, p:], transverse[:p, p:]

    def coefficients(self, s, epsilon=0.):
        u = min(.05, max(0., abs(float(s))-.05))
        radius = float(radius_profile(u)[0])
        wp, wm, cross = self.at_radius(radius)
        mask = collar_mask(s)
        # Exact fixed-angular-metric variation: N_e=N exp(e f),
        # sqrt(h_ss,e)=sqrt(h_ss) exp(-e f), A fixed.
        return (np.exp(-2*epsilon*mask)*wp,
                np.exp(-2*epsilon*mask)*wm,
                np.exp(-epsilon*mask)*cross)

    def cut_zero(self):
        cross = self.at_radius(.05)[2]
        vectors, singular, _ = np.linalg.svd(cross, full_matrices=True)
        zero = vectors[:, -1]
        if zero[0] < 0:
            zero = -zero
        return zero, float(np.linalg.norm(cross.T @ zero)), float(singular[-1])


def embed(coefficients, plus, minus, energy, interval=(-.1, .1),
          rtol=2e-11, atol=2e-13):
    """Inputs [a_left,b_right], outputs [b_left,a_right]; unequal port sizes."""
    energy = complex(energy)
    if not np.isfinite(energy) or energy.imag < 0:
        raise ValueError('Energy must be finite in the closed upper half-plane.')
    sizes = (plus*minus, plus*plus, minus*minus, minus*plus)
    ends = np.cumsum((0,)+sizes)
    shapes = ((plus, minus), (plus, plus), (minus, minus), (minus, plus))

    def unpack(state):
        return tuple(state[ends[j]:ends[j+1]].reshape(shapes[j]) for j in range(4))

    initial = np.concatenate((np.zeros(sizes[0]), np.eye(plus).ravel(),
                              np.eye(minus).ravel(), np.zeros(sizes[3]))).astype(complex)

    def rhs(s, state):
        right, forward, reverse, left = unpack(state)
        wp, wm, cross = coefficients(s)
        a, b = 1j*energy*wp, -1j*cross
        c, d = 1j*cross.conj().T, -1j*energy*wm
        ac = a-right @ c
        cd = c @ right+d
        derivatives = (ac @ right+b-right @ d,
                       ac @ forward, -reverse @ cd,
                       -reverse @ c @ forward)
        return np.concatenate([value.ravel() for value in derivatives])

    start = perf_counter()
    # Split at the two flat collar endpoints; this is not a geometry approximation.
    boundaries = [interval[0]]+[s for s in (-.05, .05) if interval[0] < s < interval[1]]+[interval[1]]
    state, evaluations, steps = initial, 0, 0
    for left_endpoint, right_endpoint in zip(boundaries[:-1], boundaries[1:]):
        solution = solve_ivp(rhs, (left_endpoint, right_endpoint), state,
                             method='DOP853', rtol=rtol, atol=atol)
        if not solution.success:
            raise RuntimeError('Stable bridge embedding failed: '+solution.message)
        state = solution.y[:, -1]
        evaluations += solution.nfev
        steps += len(solution.t)-1
    right, forward, reverse, left = unpack(state)
    scattering = np.block([[left, reverse], [forward, right]])
    if not np.all(np.isfinite(scattering)):
        raise ArithmeticError('Nonfinite bridge scattering matrix.')
    defect = np.eye(plus+minus)-scattering.conj().T @ scattering
    defect = (defect+defect.conj().T)/2
    diagnostics = dict(energy=complex_record(energy), rtol=rtol, atol=atol,
                       nfev=int(evaluations), accepted_steps=int(steps),
                       runtime_seconds=perf_counter()-start,
                       unitarity_error=float(np.max(abs(defect))),
                       minimum_passivity_defect=float(np.linalg.eigvalsh(defect)[0]),
                       maximum_singular_value=float(np.linalg.svd(scattering, compute_uv=False)[0]))
    return scattering, diagnostics


def zero_response(scattering, zero, minus):
    plus = len(zero)
    reflected = scattering[:minus, :plus] @ zero
    transmitted = scattering[minus:, :plus] @ zero
    amplitude = zero.conj() @ transmitted
    nonzero_transmitted = transmitted-zero*amplitude
    return dict(zero_transmission=complex_record(amplitude),
                transmission_magnitude=float(abs(amplitude)),
                transmission_phase_radians=float(np.angle(amplitude)),
                reflected_norm=float(np.linalg.norm(reflected)),
                transmitted_nonzero_norm=float(np.linalg.norm(nonzero_transmitted)),
                full_output_norm_squared=float(np.vdot(reflected, reflected).real
                                              +np.vdot(transmitted, transmitted).real))


def masked_geometric_target(parameters, order):
    nodes, weights = gauss(order)
    u, wu = .025*(nodes+1), .025*weights
    radius, _, second = radius_profile(u)
    mu, wm = gauss(order)
    lr = log_speed_derivatives(radius[:, None], mu[None, :], parameters)['radial'] @ wm/2
    source = -4*math.pi*second*(2*radius-radius*radius*lr)
    mask = smooth_step(u/.05)[1]/2
    return float(2*wu @ (mask*source))


def run():
    started = perf_counter()
    input_names = ('balanced_dihole.py', 'dihole_stretched_collar.py',
                   'dihole_cylinder_transfer.py', 'spherical_weyl_leads.py',
                   'curved_throat_modes.py', 'project_paths.py')
    hashes = {str(project_path(name)):hashlib.sha256(project_path(name).read_bytes()).hexdigest()
              for name in input_names}
    function, _ = symbolic_fields()
    parameters = exact_target_parameters()
    records, matrices, models = [], {}, {}
    for size, order in ((4, 48), (8, 64)):
        model = BridgeMatrices(function, parameters, size, order)
        models[size] = model
        zero, residual, gap = model.cut_zero()
        for kappa in (0., .1, 1.):
            scattering, diagnostics = embed(model.coefficients, model.plus, model.minus, 1j*kappa)
            matrices[size, kappa] = scattering
            row = dict(size=size, angular_order=order, kappa=kappa,
                       cut_zero_coefficients=zero.tolist(), cut_zero_residual=residual,
                       finite_cut_transverse_gap=gap, diagnostics=diagnostics,
                       response=zero_response(scattering, zero, size),
                       scattering=complex_record(scattering))
            records.append(row)
            print(f'nmax={size} kappa={kappa:g}: |t00|={row["response"]["transmission_magnitude"]:.12g}', flush=True)

    model = models[8]
    zero = model.cut_zero()[0]
    real_scattering, real_diagnostics = embed(model.coefficients, model.plus, model.minus, .1)
    tighter, tighter_diagnostics = embed(model.coefficients, model.plus, model.minus, .1j,
                                         rtol=4e-12, atol=4e-14)
    ode_change = float(np.max(abs(tighter-matrices[8, .1])))

    # Independent existing scalar benchmark, with the same common-clock scaling.
    scalar_errors = []
    for energy in (0., .2, .3j):
        scalar, _ = embed(lambda s: (np.array([[2.3]]), np.array([[2.3]]), np.array([[1.7]])),
                          1, 1, energy, (0., .4))
        scalar_errors.append(float(np.max(abs(scalar-_constant_exact(2.3*energy, 1.7, .4)))))

    # One finite-difference response to the actual compact null-source direction.
    epsilon = 1e-4
    varied = []
    for sign in (-1, 1):
        scattering, diagnostics = embed(lambda s: model.coefficients(s, sign*epsilon),
                                         model.plus, model.minus, .1j)
        varied.append(scattering)
    derivative = (varied[1]-varied[0])/(2*epsilon)
    derivative_t00 = zero.conj() @ derivative[model.minus:, :model.plus] @ zero
    target_coarse = masked_geometric_target(parameters, 64)
    target_fine = masked_geometric_target(parameters, 128)
    variation_signal = float(np.max(abs(varied[1]-varied[0])))

    refinement = []
    for kappa in (0., .1, 1.):
        low, high = [r for r in records if r['kappa'] == kappa]
        # Shared characteristic ordering, including the one plus zero coordinate.
        indices_in = np.r_[np.arange(5), 9+np.arange(4)]
        indices_out = np.r_[np.arange(4), 8+np.arange(5)]
        common = matrices[8, kappa][np.ix_(indices_out, indices_in)]
        refinement.append(dict(kappa=kappa,
                               common_scattering_max_change=float(np.max(abs(common-matrices[4, kappa]))),
                               zero_transmission_magnitude_change=abs(high['response']['transmission_magnitude']
                                                                     -low['response']['transmission_magnitude']),
                               reflected_norm_change=abs(high['response']['reflected_norm']-low['response']['reflected_norm'])))

    generator_error = 0.
    for s in (-.08, 0., .08):
        wp, wm, cross = model.coefficients(s)
        j = np.diag(np.r_[np.ones(model.plus), -np.ones(model.minus)])
        w = np.block([[wp, np.zeros((model.plus, model.minus))],
                      [np.zeros((model.minus, model.plus)), wm]])
        for energy in (.1, .1j):
            k = np.block([[1j*energy*wp, -1j*cross],
                          [1j*cross.T, -1j*energy*wm]])
            generator_error = max(generator_error, float(np.max(abs(k.conj().T@j+j@k+2*np.imag(energy)*w))))
    tests = dict(
        fixed_regular_basis=max(m.gram_error for m in models.values()) < 2e-12,
        actual_cut_zero=max(r['cut_zero_residual'] for r in records) < 2e-10,
        generator_current_identity=generator_error < 2e-11,
        constant_scalar_scattering=max(scalar_errors) < 2e-10,
        real_energy_unitarity=real_diagnostics['unitarity_error'] < 2e-9,
        imaginary_energy_passivity=min(r['diagnostics']['minimum_passivity_defect'] for r in records) > -2e-9,
        ode_refinement=ode_change < 2e-9,
        finite_zero_response_refinement=max(r['zero_transmission_magnitude_change'] for r in refinement) < 2e-7
            and max(r['reflected_norm_change'] for r in refinement) < 2e-7,
        variation_signal_resolved=variation_signal > 1000*max(ode_change, 1e-15),
        matched_mask_target_refinement=abs(target_coarse-target_fine) < 2e-10,
        inputs_preserved=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()))
    tests = {key:bool(value) for key, value in tests.items()}
    result = dict(
        passed=all(tests.values()), tests=tests, parameters=parameters,
        sector=dict(q=1, mu=.5, Weyl_blocks=1, plus_count='nmax+1', minus_count='nmax'),
        geometry=dict(cylinder_radius=.025, cut_radius=.05, half_collar_width=.05,
                      central_coordinate_length=.1, interval=[-.1, .1]),
        records=records, refinement=refinement,
        controls=dict(real_energy=real_diagnostics,
                      real_zero_response=zero_response(real_scattering, zero, model.minus),
                      real_scattering=complex_record(real_scattering),
                      constant_scalar_errors=scalar_errors,
                      current_identity_max_error=generator_error,
                      tighter_ode=tighter_diagnostics, ode_max_change=ode_change,
                      basis_gram_error={str(n):m.gram_error for n,m in models.items()}),
        paired_variation=dict(epsilon=epsilon, size=8, kappa=.1,
                      mask='f(s)=S_prime((abs(s)-.05)/.05)/2 on collars; zero central/outside.',
                      metric='N_e=N exp(e f), h_ss,e=h_ss exp(-2e f); angular metric and A fixed.',
                      coefficient_identity='W_e=exp(-2e f) W, T_e=exp(-e f) T.',
                      minus_scattering=complex_record(varied[0]), plus_scattering=complex_record(varied[1]),
                      scattering_derivative=complex_record(derivative),
                      zero_transmission_derivative=complex_record(derivative_t00),
                      reflected_vector_derivative_norm=float(np.linalg.norm(derivative[:model.minus, :model.plus]@zero)),
                      max_matrix_derivative=float(np.max(abs(derivative))),
                      finite_difference_signal=variation_signal,
                      geometric_weighted_null_integral=target_fine,
                      geometric_target_coarse=target_coarse,
                      geometric_target_units='D; divide by 8pi G for the required compact energy variation.',
                      limitation='One epsilon finite difference; no rigorous ODE, angular or finite-difference error enclosure.'),
        conventions=dict(
            input='[a_plus(left), b_minus(right)]', output='[b_minus(left), a_plus(right)]',
            normalization='Canonical current w_dagger J w; raw fixed orthonormal angular coordinates.',
            zero_phase='Real null vector of X_cut^T, unit norm with first component positive, same fixed core angular frame at both cuts.',
            embedding='Rr_prime=A Rr+B-Rr C Rr-Rr D; Tlr_prime=(A-Rr C)Tlr; Trl_prime=-Trl(C Rr+D); Rl_prime=-Trl C Tlr.',
            initial='Rr=Rl=0; Tlr=I_plus; Trl=I_minus.'),
        scope=[
            'Actual prescribed nonspherical bridge coefficients and both fixed collars, one q=1 mu=1/2 Weyl sector only.',
            'Finite angular characteristic scattering; truncation/refinement and ODE checks are numerical diagnostics, not continuum certificates.',
            'A single physical chiral zero input enters from the left in this Weyl convention; no extra right incoming zero is manufactured.',
            'Real-energy unitarity tests finite-cut flux conservation, not asymptotic scattering probabilities in the exterior. The opposite Weyl block is separate.',
            'At positive imaginary frequency the displayed amplitudes and reflected norms are not probabilities.',
            'The paired variation probes this same bridge and fixed cuts; its matrix derivative is not a vacuum-energy derivative or renormalized null stress.',
            'No exterior return, global spin holonomy, selected sewn quantum state, physical source/current, stability or backward-time path is supplied.',
            'The angular basis is fixed along s. Half-density and flux normalization remove longitudinal connection derivatives; no moving-basis connection is omitted.',
            'No singular-value clipping, forced unitarization, transfer-matrix overflow cancellation or interpolation of coefficients is used.'
        ],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=hashes, runtime_seconds=perf_counter()-started)
    result['figures'] = render(result)
    output = project_path('dihole_bridge_scattering_results.json')
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, refinement=refinement,
                         variation_derivative=complex_record(derivative_t00),
                         required_masked_target=target_fine, runtime_seconds=result['runtime_seconds']), indent=2))
    if not result['passed']:
        raise RuntimeError('A bounded bridge scattering control failed; inspect the saved evidence.')
    return result


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), layout='constrained')
    for size, marker in ((4, 'o'), (8, 'x')):
        rows = [r for r in result['records'] if r['size'] == size]
        x = [r['kappa'] for r in rows]
        axes[0].semilogy(x, [r['response']['transmission_magnitude'] for r in rows],
                         marker+'-', label=f'nmax {size}')
        axes[1].plot(x, [r['response']['reflected_norm'] for r in rows], marker+'-', label=f'nmax {size}')
    axes[0].set_ylabel('Cut zero-to-zero transmission magnitude')
    axes[1].set_ylabel('Reflected characteristic norm')
    for axis in axes:
        axis.set_xlabel('Imaginary frequency kappa D')
        axis.grid(alpha=.25)
        axis.legend()
    fig.suptitle('Scattering through both actual collars and the central bridge\n'
                 'q=1 finite angular pilot; imaginary-frequency amplitudes are not probabilities', fontsize=10)
    files = [project_path('dihole_bridge_scattering_progress.'+ext) for ext in ('png', 'pdf')]
    for path in files:
        fig.savefig(path, dpi=180)
    plt.close(fig)
    return [str(path) for path in files]


if __name__ == '__main__':
    run()
