"""Actual q=1 cylinder propagation; no exterior return or vacuum stress yet.

The physical observable here is the imaginary-frequency attenuation exponent
over the displayed central length.  The spherical basis is only a complete
coordinate basis; N, P and A are the actual matched-cylinder fields from129.
Finite truncation comparisons are numerical diagnostics, not enclosures.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from balanced_dihole import symbolic_fields
from dihole_stretched_collar import exact_target_parameters, field_values, gauss
from project_paths import project_path
from spherical_weyl_leads import _mouth_arguments, mouth_characteristics, zero_mode


def cylinder_fields(function, parameters, theta, radius=.025):
    data = field_values(function, parameters, radius, np.cos(theta))
    n, p, a = (data[key] for key in ('lapse', 'meridian_factor', 'unit_A_phi'))
    # The difference is regular at the poles, unlike either summand alone.
    delta = (n*p*(.5-a)-.5*np.cos(theta))/np.sin(theta)
    return p/n, delta


def zero_moments(function, parameters, order=96, inner_order=48):
    """Integrate the exact first-order zero profile after factoring sqrt(sin)."""
    cosine, weights = gauss(order)
    theta = np.arccos(cosine)
    nodes, iw = gauss(inner_order)
    half = (theta-np.pi/2)/2
    sample = np.pi/2+half[:, None]*(nodes+1)
    _, correction = cylinder_fields(function, parameters, sample)
    log_profile = half*(correction @ iw)
    density = np.exp(2*(log_profile-np.max(log_profile)))
    density /= weights @ density
    optical, _ = cylinder_fields(function, parameters, theta)
    mean = float(weights @ (density*optical))
    variance = float(weights @ (density*(optical-mean)**2))
    return dict(mean_W=mean, variance_W=variance, rms_leakage_per_kappa=float(np.sqrt(variance)),
                cosine=cosine, weights=weights, density=density, W=optical)


def angular_matrices(function, parameters, size, order, radius=.025):
    """Regular monopole basis, q=1, mu=1/2; lower physical profile is i*g."""
    cosine, weights = gauss(order)
    theta = np.arccos(cosine)
    _, _, er_plus, er_minus = _mouth_arguments(theta, 1., 'upper')
    plus, minus = mouth_characteristics(theta, 1., 'upper', size)
    upper = np.sum(er_plus[:, :, None]*minus, axis=0).real
    lower = np.sum(er_minus[:, :, None]*plus, axis=0).imag
    zeroth = np.sum(er_plus*zero_mode(theta, 1., 'upper'), axis=0)
    upper = np.column_stack((zeroth, upper))
    optical, delta = cylinder_fields(function, parameters, theta, radius)
    # Both basis profiles omit sqrt(sin theta); dcos(theta) absorbs it.
    cross = upper.T @ ((weights*delta/radius)[:, None]*lower)
    cross[1:, :] += np.diag(np.sqrt(np.arange(1, size+1)*np.arange(2, size+2))/radius)
    dimension = 2*size+1
    transverse, clock = np.zeros((dimension, dimension)), np.zeros((dimension, dimension))
    transverse[:size+1, size+1:] = cross
    transverse[size+1:, :size+1] = cross.T
    clock[:size+1, :size+1] = upper.T @ ((weights*optical)[:, None]*upper)
    clock[size+1:, size+1:] = lower.T @ ((weights*optical)[:, None]*lower)
    chirality = np.r_[np.ones(size+1), -np.ones(size)]
    gram_error = max(np.max(abs(upper.T @ (weights[:, None]*upper)-np.eye(size+1))),
                     np.max(abs(lower.T @ (weights[:, None]*lower)-np.eye(size))))
    return transverse, clock, chirality, float(gram_error)


def attenuation(transverse, clock, chirality, kappa, reference):
    generator = -kappa*chirality[:, None]*clock-1j*chirality[:, None]*transverse
    values, vectors = eigh(generator)
    # At kappa=0 this is the unique zero. Follow by overlap, not eigenvalue order.
    overlaps = abs(reference.conj() @ vectors)**2
    index = int(np.argmax(overlaps))
    vector, value = vectors[:, index], float(values[index])
    return dict(kappa=kappa, exponent_per_coordinate_length=-value,
                overlap_with_zero=float(overlaps[index]),
                projected_residual=float(np.linalg.norm(generator@vector-value*vector)),
                hermiticity_error=float(np.max(abs(generator-generator.conj().T)))), vector


def run():
    inputs = [project_path(name) for name in ('balanced_dihole.py', 'dihole_stretched_collar.py',
              'spherical_weyl_leads.py', 'project_paths.py')]
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    function, _ = symbolic_fields()
    parameters = exact_target_parameters()
    moments = zero_moments(function, parameters)
    finer = zero_moments(function, parameters, 160, 80)
    frequencies = [0., .01, .03, .1, .3, 1., 3.]
    runs = []
    for size, order in ((8, 64), (12, 96), (20, 160)):
        transverse, clock, chirality, gram_error = angular_matrices(function, parameters, size, order)
        tvalues, tvectors = eigh(transverse)
        zero = tvectors[:, np.argmin(abs(tvalues))].astype(complex)
        rows = []
        for kappa in frequencies:
            row, _ = attenuation(transverse, clock, chirality, kappa, zero)
            scalar = kappa*finer['mean_W']
            row.update(scalar_exponent=scalar, exponent_correction=row['exponent_per_coordinate_length']-scalar,
                       log_amplitude_correction_length_point1=-.1*(row['exponent_per_coordinate_length']-scalar))
            rows.append(row)
        zero_clock = float(np.real(zero.conj() @ clock @ zero))
        zero_var = float(np.linalg.norm(clock@zero-zero_clock*zero)**2)
        runs.append(dict(size=size, quadrature_order=order, rows=rows,
                         gram_error=gram_error, projected_zero_W=zero_clock,
                         projected_zero_variance=zero_var,
                         projected_first_nonzero_transverse_abs=float(np.sort(abs(tvalues))[1])))
    mean_change = abs(moments['mean_W']-finer['mean_W'])
    variance_change = abs(moments['variance_W']-finer['variance_W'])
    spectral_change = max(abs(a['exponent_per_coordinate_length']-b['exponent_per_coordinate_length'])
                          for a,b in zip(runs[-2]['rows'], runs[-1]['rows']))
    tests = dict(normalized_zero=abs(finer['weights']@finer['density']-1)<2e-14,
                 zero_moment_refinement=mean_change<1e-9 and variance_change<1e-11,
                 regular_basis_gram=max(r['gram_error'] for r in runs)<2e-12,
                 hermitian_generator=max(x['hermiticity_error'] for r in runs for x in r['rows'])<1e-11,
                 zero_frequency_mode=abs(runs[-1]['rows'][0]['exponent_per_coordinate_length'])<1e-10,
                 exact_profile_vs_projection=abs(runs[-1]['projected_zero_W']-finer['mean_W'])<1e-8
                   and abs(runs[-1]['projected_zero_variance']-finer['variance_W'])<1e-9,
                 attenuation_refinement=spectral_change<2e-9,
                 inputs_preserved=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()))
    result = dict(passed=all(tests.values()), tests={k:bool(v) for k,v in tests.items()},
                  q=1, mu=.5, cylinder_radius=.025, displayed_coordinate_length=.1,
                  moments={k:v for k,v in finer.items() if np.isscalar(v)}, runs=runs,
                  refinement=dict(mean_W=mean_change, variance_W=variance_change, exponent=spectral_change),
                  scope=['Actual central-cylinder coefficients; q=1,mu=1/2 only.',
                         'Zero moments integrate the exact first-order profile; errors are numerical sensitivity.',
                         'Nonzero-frequency eigenvalues use finite regular angular projections, not continuum enclosures.',
                         'One attenuation branch is not the finite-segment scattering matrix: endpoint mode matching is missing.',
                         'No collars, physical exterior return, renormalized stress, current, stability or chronology solved.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), input_sha256=hashes)
    render(result, finer)
    project_path('dihole_cylinder_transfer_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('passed','tests','moments','refinement')},indent=2))
    if not result['passed']:
        raise RuntimeError('Cylinder propagation checks failed.')
    return result


def render(result, moments):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout='constrained')
    axes[0].plot(moments['cosine'], moments['W']-moments['mean_W'])
    axes[0].set(xlabel='cos(theta)', ylabel='Optical factor W minus zero-mode mean',
                title='Actual variation across the bridge')
    rows=result['runs'][-1]['rows'][1:]
    axes[1].semilogx([r['kappa'] for r in rows], [r['log_amplitude_correction_length_point1'] for r in rows], 'o-')
    axes[1].set(xlabel='Imaginary frequency kappa [1/D]', ylabel='Change in log attenuation amplitude',
                title='Coupled branch versus scalar approximation')
    for ax in axes: ax.grid(alpha=.25)
    fig.suptitle('Actual q=1 central cylinder; finite angular calculation, no supporting source', fontsize=10)
    for extension in ('png','pdf'):
        fig.savefig(project_path('dihole_cylinder_transfer_progress.'+extension), dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    run()
