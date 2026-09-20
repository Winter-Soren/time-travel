"""Source-weighted two-leg feedback checks and a continuum residual budget.

Only frozen rho=32 matrices are used; no exterior equation is solved here.
The finite N2-to-N4 witness is distinct from an infinite-angular enclosure.
The latter needs conforming low-core-closed trial fields and their actual
continuum L2 residuals. Those fields/residuals are NOT present in this file.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from angular_defect_bound import defects, norm, response
from buffered_angular_bound import graph_constant, partner_bounds
from complete_angular_source import _decode, port_labels
from radial_angular_envelope import coefficients


def two_leg(s, b, derivative, incoming_low, outgoing_low, masses_in=None,
            masses_out=None):
    """Distinct input/output spaces; B: outgoing -> incoming, S the reverse.

    The derivative argument is its retained low block. D=X Y^*, where X
    lies in incoming-low space and Y in outgoing-low space. No singular
    values are discarded. W_H is the DOMAIN defect of B_H.
    """
    il, ol = np.asarray(incoming_low, int), np.asarray(outgoing_low, int)
    ih = np.setdiff1d(np.arange(len(b)), il)
    oh = np.setdiff1d(np.arange(len(b)), ol)
    if not len(il) or not len(ih) or len(il) != len(ol):
        raise ValueError('Nonempty paired low and high port spaces required')
    if norm(s) > 1+2e-12:
        raise ValueError('The canonical exterior return must be contractive')
    cross = max(norm(b[np.ix_(il, oh)]), norm(b[np.ix_(ih, ol)]))
    if cross > 2e-12:
        raise ValueError('Core must preserve the chosen angular split')
    bl, bh = b[np.ix_(il, ol)], b[np.ix_(ih, oh)]
    sll, slh = s[np.ix_(ol, il)], s[np.ix_(ol, ih)]
    shl, shh = s[np.ix_(oh, il)], s[np.ix_(oh, ih)]
    gl = np.linalg.solve(np.eye(len(ol))-sll@bl, np.eye(len(ol)))
    gr = np.linalg.solve(np.eye(len(il))-bl@sll, np.eye(len(il)))
    th = shh+shl@bl@gl@slh
    w, v = defects(bh)
    iw = np.linalg.inv(w)
    inverse = np.linalg.solve(np.eye(len(oh))-th@bh, np.eye(len(oh)))
    middle = w@inverse@w
    u, singular, vh = np.linalg.svd(derivative, full_matrices=False)
    x, y = u*np.sqrt(singular), vh.conj().T*np.sqrt(singular)
    left_input, right_input = gl.conj().T@y, gr@x
    left = left_input.conj().T@slh@bh@iw
    right = iw@shl@right_input
    b0 = np.zeros_like(b)
    b0[np.ix_(il, ol)] = bl
    exact = (response(s, b)-response(s, b0))[np.ix_(ol, il)]
    schur = gl@slh@bh@inverse@shl@gr
    exact_trace = np.trace(exact@derivative)
    factored_trace = np.trace(left@middle@right)
    ln, rn = float(np.linalg.norm(left)), float(np.linalg.norm(right))
    phi = w@response(th, bh)@v-bh.conj().T
    result = dict(
        incoming_low=len(il), outgoing_low=len(ol), high_ports=len(ih),
        source_columns=len(singular), connected_norm=norm(s), core_cross_norm=cross,
        closed_high_return_norm=norm(th), weighted_middle_norm=norm(middle),
        high_julia_norm=norm(phi), left_frobenius=ln, right_frobenius=rn,
        signed_finite_feedback=float(exact_trace.real/math.pi),
        finite_absolute_complex_trace=float(abs(exact_trace)/math.pi),
        universal_two_leg_bound=2*ln*rn/math.pi,
        finite_middle_bound=norm(middle)*ln*rn/math.pi,
        schur_identity_error=norm(exact-schur),
        trace_identity_error=float(abs(exact_trace-factored_trace)),
        factorization_error=norm(derivative-x@y.conj().T),
        julia_inverse_identity_error=norm(middle-np.eye(len(oh))-phi@bh),
        nonnormal_high_core_commutator=norm(bh.conj().T@bh-bh@bh.conj().T))
    if masses_in is not None and masses_out is not None:
        result.update(left_input_U=float(np.linalg.norm(left_input)),
            left_input_V=float(np.linalg.norm(np.asarray(masses_out)[:, None]*left_input)),
            right_input_U=float(np.linalg.norm(right_input)),
            right_input_V=float(np.linalg.norm(np.asarray(masses_in)[:, None]*right_input)))
    return result


def residual_coefficient(kappa, c_delta, a, b=.2):
    """Weighted high trace error <= coefficient * continuum volume residual.

    Both the exact and trial fields satisfy f_L-B_L g_L=X, f_H=0 (dual:
    B_L^* and swapped characteristics). The error's passive low boundary
    closure gives ||e|| <= rho/kappa. For u=chi P_H e, incoming trace is
    zero and H0 commutes P_H. With h=1/(a-b),

      ||H0 u|| <= (kappa+v+h)*rho/kappa + ||chi P_H residual||;
      || |D|^(1/2) e_H(b)||^2 <= a/sqrt(c0) * ||H0 u||^2.

    Here ||W_H^-1||^2 <= c_delta*|D_H|, c0=1-1/sqrt(2).
    This uses rho_local <= rho; no residual derivative is needed.
    """
    if not (math.isfinite(kappa) and kappa>0 and math.isfinite(c_delta)
            and c_delta>0 and 0<b<a<.5):
        raise ValueError('Positive frequency, positive defect bound and disjoint collars required')
    c0 = 1-1/math.sqrt(2)
    v = partner_bounds(a)['potential']
    h = 1/(a-b)
    prefactor = math.sqrt(c_delta*a/math.sqrt(c0))
    value = prefactor*((kappa+v+h)/kappa+1)
    return dict(a=a, b=b, kappa=kappa, c_delta=c_delta, c0=c0,
        potential_bound=v, cutoff_slope=h, trace_prefactor=prefactor,
        residual_to_weighted_leg=value)


def symmetric_budget(left, right, target):
    """Largest equal nonnegative leg error with 2/pi*(L+e)*(R+e)<=target."""
    spare = math.pi*target/2-left*right
    if spare <= 0:
        return None
    # Rationalization avoids subtraction near the already achieved bound.
    return 2*spare/(math.sqrt((left+right)**2+4*spare)+left+right)


def identity_controls():
    """Independent nonnormal B and distinct input/output permutations."""
    rng = np.random.default_rng(105)
    rows = []
    for size, low in ((5,2), (9,5), (15,9)):
        z = rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
        s = .96*z/norm(z)
        b = np.zeros_like(s)
        for cut, scale in ((slice(0,low),.83), (slice(low,size),.97)):
            n = len(b[cut])
            z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            b[cut,cut] = scale*z/norm(z)
        d = rng.normal(size=(low,low))+1j*rng.normal(size=(low,low))
        pin, pout = rng.permutation(size), rng.permutation(size)
        il, ol = np.flatnonzero(pin<low), np.flatnonzero(pout<low)
        new_d = d[np.ix_(pin[il],pout[ol])]
        rows.append(two_leg(s[np.ix_(pout,pin)], b[np.ix_(pin,pout)], new_d, il, ol))
    return rows


def run():
    paths = ('long_throat_source_results.json', 'angular_defect_bound.py',
             'buffered_angular_bound.py', 'radial_angular_envelope.py',
             'complete_angular_source.py', 'complete_angular_modes.py', 'azimuthal_modes.py')
    hashes = {p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
    saved = json.loads(Path(paths[0]).read_text(encoding='utf-8'))
    group = next(g for g in saved['geometries'] if g['geometry']['rho_c']==32.)
    x, ell = group['sums']['4'][1]['x'], group['geometry']['ell']
    kappa = x/ell
    cd = coefficients(x,32.,0,cutoff=2)['delta_inverse_coefficient']
    collars = [residual_coefficient(kappa,cd,a) for a in (.3,.4,.49)]
    best = min(collars,key=lambda r:r['residual_to_weighted_leg'])
    controls, records = identity_controls(), []
    for ip, phase in enumerate(('0','pi')):
        for probe in range(2):
            sectors = []
            full_source = 0.
            for sector in group['sectors']:
                mat = sector['samples'][1]['matrices']
                s, ref = (_decode(mat[k]) for k in ('S_connected','S_reference'))
                b = _decode(mat['phases'][ip]['B'])
                derivative = _decode(mat['phases'][ip]['pressure'])[probe]
                full_source += float(np.trace((response(s,b)-response(ref,b))@derivative).real/math.pi)
                if abs(sector['m'])>2:
                    continue
                li, lo = port_labels(sector['mu'],4)
                ni, no = np.array([n for _,n in li+lo]), np.array([n for _,n in lo+li])
                il, ol = np.flatnonzero(ni<=2), np.flatnonzero(no<=2)
                ih, oh = np.flatnonzero(ni>2), np.flatnonzero(no>2)
                row = two_leg(s,b,derivative[np.ix_(il,ol)],il,ol,
                    np.sqrt(ni[il]*(ni[il]+1)),np.sqrt(no[ol]*(no[ol]+1)))
                row.update(m=sector['m'],reference_cross_norm=max(
                    norm(ref[np.ix_(ol,ih)]),norm(ref[np.ix_(oh,il)])))
                sectors.append(row)
            left = math.sqrt(sum(r['left_frobenius']**2 for r in sectors))
            right = math.sqrt(sum(r['right_frobenius']**2 for r in sectors))
            feedback = sum(r['signed_finite_feedback'] for r in sectors)
            bound = 2*left*right/math.pi
            expected = group['sums']['4'][1]['kernels'][phase]['pressure'][probe]
            graph = []
            for collar in collars:
                inputs = {}
                for side in ('left','right'):
                    U = math.sqrt(sum(r[side+'_input_U']**2 for r in sectors))
                    V = math.sqrt(sum(r[side+'_input_V']**2 for r in sectors))
                    inputs[side] = graph_constant(kappa,U,V,collar['a'])
                graph.append(dict(a=collar['a'],inputs=inputs,
                    conditional_graph_two_leg_bound=2/math.pi*cd/math.sqrt(12)*
                    inputs['left']['C_X']*inputs['right']['C_X']))
            budgets = []
            for fraction in (.01,.1,1.):
                target = fraction*abs(expected)
                epsilon = symmetric_budget(left,right,target)
                rho = None if epsilon is None else epsilon/best['residual_to_weighted_leg']
                achieved = None if epsilon is None else 2/math.pi*(left+epsilon)*(right+epsilon)
                budgets.append(dict(signal_fraction=fraction,target_feedback_budget=target,
                    equal_weighted_leg_error_budget=epsilon,
                    equal_continuum_residual_norm_budget=rho,
                    reconstructed_feedback_budget=achieved))
            records.append(dict(phase=phase,probe=probe,sectors=sectors,
                full_saved_N4_source=expected,reconstructed_source_error=abs(full_source-expected),
                signed_finite_N2_to_N4_low_feedback=feedback,
                left_frobenius=left,right_frobenius=right,
                universal_two_leg_bound=bound,
                bound_fraction_of_full_saved_source=bound/abs(expected),
                sum_sector_two_leg_bounds=sum(r['universal_two_leg_bound'] for r in sectors),
                finite_middle_bound=sum(r['finite_middle_bound'] for r in sectors),
                maximum_finite_middle_norm=max(r['weighted_middle_norm'] for r in sectors),
                continuum_residual_budgets=budgets,graph_only_comparison=graph))
    all_rows = controls+[s for r in records for s in r['sectors']]
    tests = dict(
        complete_retained_multiplets=all(len(r['sectors'])==5 and
            sum(s['source_columns'] for s in r['sectors'])==34 for r in records),
        nonnormal_controls=min(r['nonnormal_high_core_commutator'] for r in controls)>1e-3,
        identities=max(max(r[k] for k in ('schur_identity_error','trace_identity_error',
            'factorization_error','julia_inverse_identity_error')) for r in all_rows)<3e-12,
        passive_high_closures=max(r['closed_high_return_norm'] for r in all_rows)<=1+2e-12,
        julia_bounds=max(r['high_julia_norm'] for r in all_rows)<=1+2e-12 and
            max(r['weighted_middle_norm'] for r in all_rows)<=2+2e-12,
        finite_source_bounds=all(r['finite_absolute_complex_trace']<=r['finite_middle_bound']+2e-12
            and r['finite_middle_bound']<=r['universal_two_leg_bound']+2e-12 for r in all_rows),
        reference_low_feedback_cancels=max(s['reference_cross_norm'] for r in records for s in r['sectors'])<2e-12,
        saved_source_reconstruction=max(r['reconstructed_source_error'] for r in records)<2e-12,
        global_direct_sum_bound=all(abs(r['signed_finite_N2_to_N4_low_feedback'])<=r['sum_sector_two_leg_bounds']+2e-12
            and r['sum_sector_two_leg_bounds']<=r['universal_two_leg_bound']+2e-12 for r in records),
        residual_budget_equation=all(b['reconstructed_feedback_budget'] is None or
            abs(b['reconstructed_feedback_budget']/b['target_feedback_budget']-1)<3e-15
            for r in records for b in r['continuum_residual_budgets']),
        inputs_preserved=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in hashes.items()))
    tests = {k:bool(v) for k,v in tests.items()}
    result = dict(passed=all(tests.values()),tests=tests,
        parameters=dict(q=1,b=.2,d=1.,rho_c=32.,x=x,kappa=kappa,ell=ell,
            retained_N=2,saved_high_N=4,source='pressure',normalization='unscaled d-kappa kernel'),
        collars=collars,best_collar=best,identity_controls=controls,records=records,
        conclusions=dict(continuum_source_sign_certified=False,continuum_residuals_available=False,
            saved_finite_feedback_bound_below_one_percent_signal=all(
                r['bound_fraction_of_full_saved_source']<.01 for r in records)),
        residual_trial_requirements=[
            'Conforming physical spinor in the actual exterior operator domain, with the monopole gauge/spin transitions and physical round-mouth flux measure.',
            'Right trial satisfies incoming_L-B_L*outgoing_L=X and incoming_H=0 exactly; dual trial uses B_L-adjoint, Y and swapped characteristics at -i*kappa.',
            'Continuum L2 residuals are for these source-weighted columns, aggregated by Frobenius norm over all retained m, both mouths and both Weyls once.',
            'The trial belongs to the L2 exterior domain at infinity, or obeys an exactly passive outer continuation. A conforming compactly supported global trial is admissible.',
            'The complete high angular trace must be computed or enclosed, including modes above saved N4; a truncated matrix trace alone is not the trial trace.',
            'If a boundary relation, gauge/interface matching, outer condition or core map is approximate, its mismatch must be lifted or separately bounded; it is not part of an algebraic LU residual.'],
        scope=[
            'Exact operator identities and analytic inequalities are evaluated on frozen finite matrices in floating point, not outward-rounded enclosures.',
            'The universal factor 2 applies to every passive high closure; the much smaller finite middle norm is recorded only as a finite witness.',
            'Residual budgets assume conforming dressed trials whose full weighted high trace norms are no larger than the recorded illustrative finite leg norms, or those norms plus separately enclosed corrections.',
            'The low-core-closed residual argument avoids a separate true-S_LL enclosure, only if the trial low feedback boundary relation is exact.',
            'Budgets allocate 1%, 10% or 100% of a saved finite kernel to this feedback error alone. No direct angular, frequency, continuum discretization or physical-support error is included.',
            'The graph-only two-leg comparison is also conditional on true continuum data and is not substituted for measured continuum residuals.',
            'All source and feedback values are pointwise positive-frequency d-kappa kernels; no frequency integral or new exterior solve is performed.'],
        input_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=['https://arxiv.org/pdf/1803.09329','https://arxiv.org/pdf/hep-th/0206054'])
    if not result['passed']:
        raise AssertionError(tests)
    Path(__file__).with_name('feedback_residual_bound_results.json').write_text(
        json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,best_collar=best,
        phase_pi=[{k:v for k,v in r.items() if k not in ('sectors','graph_only_comparison')}
                  for r in records if r['phase']=='pi']),indent=2))
    return result


if __name__=='__main__':
    run()
