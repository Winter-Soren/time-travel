"""All-radius central-slab coercivity and an intrinsic quantum transmission bound.

The intrinsic map prescribes incoming data on the two slab planes. Exterior
reflections and throat feedback are not included in this map.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm

from project_paths import project_path


def certify_slab():
    """Interval boxes plus analytic inner/outer estimates; no sampled gap claim."""
    iv, old = mp.iv, mp.iv.dps
    iv.dps = 40
    try:
        charge, k = iv.mpf(1)/5, iv.mpf(1)/2
        mass = (iv.sqrt(k*k+4*charge*k)-k)/2
        a = iv.sqrt(k*k-mass*mass)
        rp = mass+k
        field = mass/((k+a)*rp)
        target, height = iv.mpf('.18'), iv.mpf('.25')

        def box(lo, hi, ylo, yhi):
            rho = iv.mpf([str(lo), str(hi)])
            y = iv.mpf([str(ylo), str(yhi)])
            delta = rho*rho/y
            r = mass+iv.sqrt(k*k+delta)
            sigma = r*r-a*a*(1-y)
            w = (r*r-a*a)**2+delta*a*a*y
            j = delta+a*a*y+2*field*mass*r*a*y+field*field*w*y/4
            lapse = j/sigma
            potential = (2*mass*r*a+field*w/2)*y/(2*charge*j)
            # U*rho is nonsingular at rho=0. z^2=(k^2+rho^2/y)(1-y).
            return lapse*lapse*(potential-iv.mpf('.5')), iv.sqrt((k*k+delta)*(1-y))

        stack = [(.05, 8., .75, 1., 0)]
        accepted, discarded, max_depth = [], 0, 0
        worst = iv.mpf(100)
        while stack:
            lo, hi, yl, yh, depth = stack.pop()
            coefficient, z = box(lo, hi, yl, yh)
            if z.a > height:
                discarded += 1
                continue
            lower = coefficient/iv.mpf(str(hi))
            if lower.a > target:
                accepted.append([lo, hi, yl, yh])
                if lower.a < worst.a:
                    worst = lower.a
                max_depth = max(max_depth, depth)
                continue
            if depth >= 28:
                raise AssertionError('Interval subdivision did not certify the slab.')
            if (hi-lo)/max(lo, .1) > 4*(yh-yl):
                middle = (lo+hi)/2
                stack.extend([(lo, middle, yl, yh, depth+1), (middle, hi, yl, yh, depth+1)])
            else:
                middle = (yl+yh)/2
                stack.extend([(lo, hi, yl, middle, depth+1), (lo, hi, middle, yh, depth+1)])

        inner = [box(0, .05, .75+j/128, .75+(j+1)/128)[0]/iv.mpf('.05') for j in range(32)]
        assert all(bound.a > target for bound in inner)
        # For rho>=R, r>=R, y>=3/4, Sigma<=r^2 and rho<=r.
        # J>=j0*r^4 and K=J(A-1/2)>=k0*r^4, so U>=j0*k0*R^3.
        y, radius = iv.mpf(3)/4, iv.mpf(8)
        w0 = (1-a*a/(radius*radius))**2
        j0 = field*field*y*w0/4
        k0 = y*field*(2-field*charge)*w0/(8*charge)-1/(2*radius*radius)-a*a/(2*radius**4)
        assert 1/charge-field > 0 and k0.a > 0
        outer = j0*k0*radius**3
        assert outer.a > target and worst.a > target
        return dict(passed=True, interval_digits=40, slab_half_width='.25 D',
                    certified_U_lower='.18/D', compact_rho=['.05', '8'],
                    compact_boxes=len(accepted), outside_slab_boxes=discarded,
                    maximum_subdivision_depth=max_depth, compact_minimum_lower_interval=str(worst),
                    inner_minimum_lower=min(float(x.a) for x in inner), outer_lower_interval=str(outer),
                    accepted_partition_sha256=hashlib.sha256(json.dumps(accepted).encode()).hexdigest(),
                    coordinate_coverage='rho>0, |z|<=.25; y=sin²theta>=3/4; exact Delta=rho²/y',
                    speed_bound='v=h²/Sigma²<1 since Sigma-h=2Mr-M²(1-cos²theta)>0')
    finally:
        iv.dps = old


def algebra_checks():
    """Noncommuting finite-channel controls of the exact Stokes identities."""
    sx = np.array([[0., 1.], [1., 0.]], complex)
    sy = np.array([[0., -1j], [1j, 0.]])
    sz = np.diag([1., -1.])
    rng = np.random.default_rng(112)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    momentum = (raw+raw.conj().T)/6
    velocity = np.array([.45, .6, .8, 1.])
    u = np.array([.19, .31, .55, .8])
    sigma_x, current = np.kron(np.eye(4), sx), np.kron(np.eye(4), sz)
    plus, minus = np.arange(0, 8, 2), np.arange(1, 8, 2)
    rows = []
    for real, imaginary in ((0., 0.), (0., .2), (.12, 0.), (.12, .2), (0., 2.)):
        generator = np.kron(momentum, sy)+np.kron(np.diag(u/velocity), sx)
        generator += np.kron(np.diag((1j*real-imaginary)/velocity), sz)
        target_s = 2*np.kron(np.diag(u/velocity), np.eye(2))+2*real*np.kron(np.diag(1/velocity), sy)
        target_j = -2*imaginary*np.kron(np.diag(1/velocity), np.eye(2))
        identity_error = max(np.linalg.norm(generator.conj().T@sigma_x+sigma_x@generator-target_s),
                             np.linalg.norm(generator.conj().T@current+current@generator-target_j))
        transfer = expm(.5*generator)
        aa, bb = transfer[np.ix_(plus, plus)], transfer[np.ix_(plus, minus)]
        cc, dd = transfer[np.ix_(minus, plus)], transfer[np.ix_(minus, minus)]
        reflection = -np.linalg.solve(dd, cc)
        transmission = aa+bb@reflection
        norm = float(np.linalg.norm(transmission, 2))
        bound = 1/math.cosh(.5*(.18-abs(real)))
        assert identity_error < 2e-14 and norm <= bound+2e-14
        rows.append(dict(real_energy=real, imaginary_energy=imaginary,
                         transmission_norm=norm, rigorous_intrinsic_bound=bound,
                         identity_error=float(identity_error)))
    # A constant scalar barrier at real E=0 attains the sech bound exactly.
    uniform = expm(.5*.18*sx)
    amplitude = float((1/uniform[1, 1]).real)
    assert abs(amplitude-1/math.cosh(.09)) < 3e-16
    return rows


def run():
    certificate, controls = certify_slab(), algebra_checks()
    amplitude = 1/math.cosh(.09)
    # The cutoff estimate concerns a regular raw exterior solution with passive
    # infinity. These are floating evaluations of the analytic inequalities.
    gap_low = math.hypot(.18, .2)
    every_cutoff_lower = math.exp(-gap_low*.25)/(.2*.25**2*gap_low)
    high = lambda m: math.exp(-48*(.25-.025-.025))/(4*64*.025**2*(math.hypot(m, 64)-48))
    assert every_cutoff_lower > 277 and high(.18) < 2.646e-5
    result = dict(passed=True, parameters=dict(Q_over_D='1/5', D=1, q=1, mu='1/2', slab_half_width_over_D='1/4'),
                  certificate=certificate, finite_channel_controls=controls,
                  intrinsic_imaginary_frequency_amplitude_bound=amplitude,
                  intrinsic_imaginary_frequency_squared_norm_bound=amplitude**2,
                  raw_exterior_cutoff_diagnostics=dict(
                      low_kappa=.2, lower_bound_on_every_cutoff_candidate=every_cutoff_lower,
                      resulting_low_kappa_norm_bound=1., high_kappa=64.,
                      high_cutoffs=dict(d=.025, h=.025, alpha=48.),
                      high_kappa_norm_bound=high(.18), zero_magnetic_gap_comparator=high(0.),
                      interpretation='Low-frequency bound is vacuous; high-frequency suppression is almost entirely imaginary-frequency absorption.'),
                  theorem='In the q=1, mu=1/2 sector, for E=xi+i*kappa, kappa>=0 and |xi|<.18/D, ||T_slab||<=sech[.5*(.18-|xi|D)].',
                  boundary_condition='Only one slab face has prescribed incoming data; the other face has zero incoming characteristic data. Canonical current trace w=sqrt(v)u.',
                  scope=['The all-radius slab gap and intrinsic norm estimate use the full two-dimensional meridional operator at fixed q=1, mu=1/2; no radial-mode decoupling is assumed. Other azimuthal sectors are not covered.',
                    'Transmission estimates are a priori for regular solutions with the stated boundary data; this calculation does not construct or prove existence of the full infinite-channel scattering map.',
                    'Finite-channel matrices check algebra; they are not the dihole scattering simulation.',
                    'Exterior reflections can feed both slab faces. This is not an upper-mouth to lower-mouth return bound after feedback.',
                    'The bound is weak: squared norm below about .992 does not show an opaque barrier or determine the actual return amplitude.',
                    'No vacuum stress, supported throat, stability or chronology-changing medium has been obtained.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={'project_paths.py':hashlib.sha256(project_path('project_paths.py').read_bytes()).hexdigest()},
                  sources=['https://arxiv.org/abs/hep-th/9906160'])
    project_path('dihole_slab_bound_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=True, compact_boxes=certificate['compact_boxes'],
                         slab_gap=certificate['certified_U_lower'], amplitude_bound=amplitude,
                         squared_norm_bound=amplitude**2, largest_identity_error=max(x['identity_error'] for x in controls)), indent=2))
    return result


if __name__ == '__main__':
    run()
