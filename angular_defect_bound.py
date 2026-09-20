"""Check a defect-weighted angular remainder inequality on frozen matrices.

This is our operator estimate, not a certified infinite-mode error bar.
The Julia-unitary feedback identity removes an unweighted high-block
inverse from the estimate. Weighted exterior mixing remains unbounded.
Run from this directory; no exterior solve or old artifact is rerun.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from complete_angular_source import _decode, port_labels
from curved_throat_modes import throat_scattering


def norm(a):
    return float(np.linalg.norm(a, 2)) if a.size else 0.


def trace_norm(a):
    return float(np.linalg.svd(a, compute_uv=False).sum()) if a.size else 0.


def defects(b):
    """Strict defects; never clip a nonpositive eigenvalue into acceptance."""
    answer = []
    for gram in (b.conj().T @ b, b @ b.conj().T):
        value, vectors = np.linalg.eigh(np.eye(len(b)) - gram)
        if value.min() <= 0:
            raise ValueError('The finite core must be strictly contractive')
        answer.append((vectors * np.sqrt(value)) @ vectors.conj().T)
    return answer


def response(s, b):
    return np.linalg.solve(np.eye(len(b)) - s @ b, s)


def estimate(s, reference, b, derivatives, incoming_low, outgoing_low):
    """Input/output spaces are distinct; their retained selectors can differ."""
    if max(norm(s), norm(reference)) > 1 + 2e-12:
        raise ValueError('This inequality requires passive finite maps')
    di, do = defects(b)  # Domain of B (exterior output), codomain of B.
    f, fr = response(s, b), response(reference, b)
    k = f - fr
    il, ol = np.asarray(incoming_low, int), np.asarray(outgoing_low, int)
    ih = np.setdiff1d(np.arange(len(b)), il)
    oh = np.setdiff1d(np.arange(len(b)), ol)
    bh, sh = b[np.ix_(ih, oh)], s[np.ix_(oh, ih)]
    dhi, dho = defects(bh)
    inverse_dhi, inverse_dho = np.linalg.inv(dhi), np.linalg.inv(dho)
    high_inverse = np.linalg.solve(np.eye(len(bh)) - sh @ bh, np.eye(len(bh)))
    weighted_inverse = dhi @ high_inverse @ dhi
    low_feedback_bound = 0.
    compressed = np.empty((0, 0), complex)
    effective_passivity = 0.
    if len(il):
        bl, sl, rl = b[np.ix_(il, ol)], s[np.ix_(ol, il)], reference[np.ix_(ol, il)]
        slh, shl = s[np.ix_(ol, ih)], s[np.ix_(oh, il)]
        uw, vw = norm(slh @ bh @ inverse_dhi), norm(inverse_dhi @ shl)
        gl = 1 / (1 - norm(bl))
        low_feedback_bound = gl * gl * uw * vw
        compressed = response(sl, bl) - response(rl, bl)
        effective = sl + slh @ bh @ high_inverse @ shl
        effective_passivity = norm(effective)
    else:
        uw = vw = gl = 0.
    rows = []
    for derivative in derivatives:
        dl, dh = derivative[np.ix_(il, ol)], derivative[np.ix_(ih, oh)]
        cross = max(norm(derivative[np.ix_(il, oh)]), norm(derivative[np.ix_(ih, ol)]))
        if cross > 2e-12:
            raise ValueError('The metric derivative must preserve the angular split')
        full = float(np.trace(k @ derivative).real / np.pi)
        retained = float(np.trace(k[np.ix_(ol, il)] @ dl).real / np.pi)
        cut = float(np.trace(compressed @ dl).real / np.pi)
        direct_bound = 2 * trace_norm(inverse_dho @ dh @ inverse_dhi) / np.pi
        feedback_bound = 2 * low_feedback_bound * trace_norm(dl) / np.pi
        rows.append(dict(full=full, raw_compressed=cut, direct_high=full-retained,
            low_feedback=retained-cut, error=full-cut,
            direct_bound=direct_bound, feedback_bound=feedback_bound,
            total_bound=direct_bound+feedback_bound))
    return dict(connected_norm=norm(s), reference_norm=norm(reference), core_norm=norm(b),
        connected_julia_norm=norm(di @ f @ do-b.conj().T),
        reference_julia_norm=norm(di @ fr @ do-b.conj().T),
        weighted_relative_norm=norm(di @ k @ do),
        high_weighted_inverse_norm=norm(weighted_inverse),
        high_unweighted_inverse_norm=norm(high_inverse),
        effective_return_norm=effective_passivity, u_weighted=uw, v_weighted=vw, g_low=gl,
        reference_cross_norm=max(norm(reference[np.ix_(ol, ih)]), norm(reference[np.ix_(oh, il)])),
        core_cross_norm=max(norm(b[np.ix_(il, oh)]), norm(b[np.ix_(ih, ol)])),
        probes=rows)


def self_check():
    """Independent nonnormal complex contractions exercise the Julia identity."""
    rng = np.random.default_rng(102)
    worst_julia = worst_relative = worst_inverse = 0.
    for size in (2, 5, 9):
        for scale in (.2, .8, .999):
            matrices = []
            for _ in range(3):
                a = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
                matrices.append(a * scale / norm(a))
            b, s, r = matrices
            di, do = defects(b)
            fs, fr = response(s, b), response(r, b)
            c = di @ fs @ do - b.conj().T
            inverse = np.linalg.solve(np.eye(size)-s@b, np.eye(size))
            assert norm(di @ inverse @ di-(np.eye(size)+c@b)) < 2e-13
            worst_julia = max(worst_julia, norm(c))
            worst_relative = max(worst_relative, norm(di@(fs-fr)@do))
            worst_inverse = max(worst_inverse, norm(di@inverse@di))
    assert worst_julia <= 1+2e-12 and worst_relative <= 2+2e-12 and worst_inverse <= 2+2e-12
    return dict(max_julia_norm=worst_julia, max_relative_norm=worst_relative,
                max_weighted_inverse_norm=worst_inverse, nonnormal_cases=9)


def absorption_checks():
    """Positive defect ODE versus the independent complex scattering solver."""
    cut = np.arctan(2.)
    length, ell = 2*cut, .2*np.sqrt(5.)
    rows = []
    for n in (3, 8, 16):
        lam, mstar = np.sqrt(n*(n+1)), np.sqrt(n*(n+1))/np.cos(cut)
        for kappa in (1e-4, .4, 2., 12., 96.):
            k = kappa*ell
            def rhs(sigma, state):
                y, logt, delta = state
                mass = lam/np.cos(sigma)
                return [mass*(1-y*y)-2*k*y, -k-mass*y,
                        -2*mass*y*delta+2*k*(2*y*y+np.exp(2*logt))]
            sol = solve_ivp(rhs, (-cut,cut), [0.,0.,0.], method='DOP853',
                            rtol=5e-12, atol=5e-14, max_step=length/24)
            if not sol.success:
                raise RuntimeError(sol.message)
            y, logt, delta = sol.y[:,-1]
            omega = np.hypot(lam,k)
            th = np.tanh(omega*length/2)
            y0 = lam*th/(omega+k*th)
            lower = 2*k/mstar*y0*y0*(-np.expm1(-mstar*length))
            scatter, _ = throat_scattering(1j*kappa,n,ell,2.,rtol=5e-12,atol=5e-14)
            defect = np.eye(2)-scatter.conj().T@scatter
            rows.append(dict(n=n,kappa=kappa,defect_positive_ode=float(delta),
                explicit_lower_bound=float(lower),
                scattering_defect_error=float(np.max(abs(defect-delta*np.eye(2)))),
                subtraction_error=float(abs(delta-(1-y*y-np.exp(2*logt))))))
    assert all(0 < row['explicit_lower_bound'] <= row['defect_positive_ode']+2e-12 for row in rows)
    assert max(row['scattering_defect_error'] for row in rows) < 2e-10
    assert max(row['subtraction_error'] for row in rows) < 2e-10
    return rows


def run():
    inputs = ('complete_angular_source_results.json', 'complete_angular_source.py',
              'complete_angular_modes.py', 'azimuthal_modes.py', 'curved_throat_modes.py')
    hashes = {name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in inputs}
    saved = json.loads(Path(inputs[0]).read_text(encoding='utf-8'))
    records = []
    for sector in saved['complete_cuts'][1]['sectors']:
        mu, m = sector['config']['mu'], sector['config']['m']
        sample = next(row for row in sector['samples'] if row['kappa'] == 2.)
        matrices = sample['matrices']
        s, reference = (_decode(matrices[key]) for key in ('S_connected', 'S_reference'))
        labels_in, labels_out = port_labels(mu, 4)
        count = len(labels_in)
        pin = [i for i, (_, n) in enumerate(labels_in) if n <= 2]
        pout = [i for i, (_, n) in enumerate(labels_out) if n <= 2]
        il, ol = np.r_[pin, np.asarray(pout)+count], np.r_[pout, np.asarray(pin)+count]
        for phase, matrix, saved_phase in zip(('0', 'pi'), matrices['phases'], sample['phases']):
            record = estimate(s, reference, _decode(matrix['B']), _decode(matrix['pressure']), il, ol)
            error = float(np.max(abs(np.array([p['full'] for p in record['probes']])-saved_phase['pressure_integrands'])))
            records.append(dict(m=m, mu=mu, phase=phase, saved_kernel_error=error, **record))
    controls = self_check()
    absorption = absorption_checks()
    tests = dict(
        complete_sectors=len(records)==18 and sorted({r['m'] for r in records})==list(range(-4,5)),
        saved_kernel_reconstruction=max(r['saved_kernel_error'] for r in records)<2e-12,
        julia_contractions=max(max(r['connected_julia_norm'],r['reference_julia_norm']) for r in records)<=1+2e-12,
        weighted_relative_bound=max(r['weighted_relative_norm'] for r in records)<=2+2e-12,
        weighted_inverse_bound=max(r['high_weighted_inverse_norm'] for r in records)<=2+2e-12,
        effective_passivity=max(r['effective_return_norm'] for r in records)<=1+2e-12,
        block_structure=max(max(r['reference_cross_norm'],r['core_cross_norm']) for r in records)<2e-12,
        source_bounds=all(abs(p[key])<=p[bound]+2e-12 for r in records for p in r['probes']
            for key,bound in (('direct_high','direct_bound'),('low_feedback','feedback_bound'),('error','total_bound'))),
        positive_absorption_and_lower_bound=all(0<r['explicit_lower_bound']<=r['defect_positive_ode']+2e-12 for r in absorption),
        independent_scattering_absorption=max(r['scattering_defect_error'] for r in absorption)<2e-10)
    sums = {}
    for phase in ('0','pi'):
        subset = [r for r in records if r['phase']==phase]
        sums[phase] = {key:[sum(abs(r['probes'][j][key]) if key=='error' else r['probes'][j][key]
            for r in subset) for j in range(2)] for key in ('error','total_bound','direct_bound','feedback_bound')}
    result = dict(passed=all(tests.values()), tests=tests, self_check=controls, absorption_checks=absorption,
        parameters=dict(kappa=2., high_nmax=4, low_nmax=2, outer_nmax=9, h=.05, rho_c=2.),
        input_sha256=hashes, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        records=records, sums_of_sector_absolute_errors_and_bounds=sums,
        scope='Finite saved-matrix checks against exact raw compression; actual N2 Gram correction is separate.',
        unresolved=['Infinite defect-weighted exterior mixing','Low-frequency integrability',
                    'Continuum and exterior-tail numerical errors','Physical geometry and supporting source'])
    Path(__file__).with_name('angular_defect_bound_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,sums=sums),indent=2))
    assert result['passed']
    return result


if __name__ == '__main__':
    run()
