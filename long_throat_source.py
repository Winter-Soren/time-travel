"""Selected-frequency source sensitivity to a longer prescribed product throat.

Only rho_c changes: b=.2, d=1, q=1 and the exterior mesh/cutoffs are fixed.
ell=b*sqrt(1+rho_c**2) keeps the endpoint lapse one.  This changes the
redshift and physical core length, not just coordinates.  The same fractional
sigma probes change their proper support; their full geometry is recorded.

The rho_c=2 kernels are reused exactly from frozen101 fine quadrature nodes.
Fresh exterior responses at matched x=kappa*ell are computed for rho_c=8,32.
These are three frequency samples, NOT an integral or a support/NEC verdict.
The sharp throat/exterior join still requires an unresolved material source.

Run: .venv\\Scripts\\python.exe -B -W error long_throat_source.py
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

from ads_reference_source import cap_reflection, reference_relative
from complete_angular_source import (CompleteSector, FIELDS, PHASES, _decode,
    complete_reference, port_labels, reproject_environment, selected_map,
    source_response)
from exterior_boundary_map import _complex_json
from throat_null_response import (mode_variation, probe_geometry, raw_variation,
    shift_derivatives)


BASELINE = 'complete_angular_source_results.json'
OUTPUT = Path(__file__).with_name('long_throat_source_results.json')
B = .2
N = 4
J = 9
RHO_VALUES = (2., 8., 32.)


def geometry(rho_c):
    ell = B * math.hypot(1., rho_c)
    cut = math.atan(rho_c)
    probes = probe_geometry(ell, rho_c)
    for probe in probes:
        probe['proper_center'] = B*math.asinh(math.tan(probe['center']))
        probe['proper_support'] = [B*math.asinh(math.tan(s)) for s in probe['support']]
        probe['proper_width'] = probe['proper_support'][1]-probe['proper_support'][0]
        probe['proper_coordinate'] = 'l=b*asinh(tan(sigma)), with l=0 at the throat center'
    return dict(rho_c=rho_c, b=B, ell=ell, sigma_cut=cut,
        proper_length=2*B*math.asinh(rho_c), center_lapse=1/math.hypot(1., rho_c),
        endpoint_lapse=B/(ell*math.cos(cut)), optical_length=2*ell*cut,
        reference_global_optical_length=math.pi*ell,
        metric='r=b; N=b/(ell*cos(sigma)); A=b/cos(sigma); ell=b*sqrt(1+rho_c^2)',
        probes=probes)


def kernels(phases):
    return {phase: {field: row[field+'_integrands'] for field in FIELDS}
        for phase, row in zip(PHASES, phases)}


def annotate(row, x, ell):
    row['x'] = x
    row['kernels'] = kernels(row['phases'])
    row['scaled_kernels'] = {p: {f: (np.asarray(v)/ell).tolist()
        for f, v in fields.items()} for p, fields in row['kernels'].items()}
    row['scaled_frequency_identity_error'] = abs(row['kappa']*ell-x)
    if 'lower_cut_response' in row:
        annotate(row['lower_cut_response'], x, ell)
    return row


def baseline_nodes(prior):
    """Select exact saved nodes, rather than interpolate an old kernel."""
    panels = prior['integration']['cutoffs']['4']['partial_sums']['fine']['panels']
    nodes = [(float(k), ip, ik) for ip, panel in enumerate(panels)
        for ik, k in enumerate(panel['kappa'])]
    ell = geometry(2.)['ell']
    chosen = [min(nodes, key=lambda row: abs(row[0]-target)) for target in (.1, 1., 10.)]
    return [dict(kappa_rho2=k, x=k*ell, panel_index=ip, node_index=ik,
        source_recipe='integration.cutoffs[cutoff].sectors[m].fine.panels')
        for k, ip, ik in chosen]


def baseline_group(prior, nodes):
    geom = geometry(2.)
    group = dict(geometry=geom, h=.05, outer_nmax=J, nmax=N,
        origin='Exact frozen101 kernels; no exterior or radial solve here', sectors=[])
    for high in prior['integration']['cutoffs']['4']['sectors']:
        m = high['m']
        low = (next(s for s in prior['integration']['cutoffs']['2']['sectors'] if s['m']==m)
            if abs(m)<=2 else None)
        sector = dict(m=m, mu=high['mu'], config=high['fine']['config'], samples=[])
        for node in nodes:
            def extract(record):
                ip, ik = node['panel_index'], node['node_index']
                panel = record['fine']['panels'][ip]
                if panel['kappa'][ik] != node['kappa_rho2']:
                    raise ValueError('Baseline nodes are not exactly aligned')
                phases = [dict(phase=angle,
                    **{f+'_integrands':panel['integrands'][p][f][ik] for f in FIELDS},
                    minimum_singular_values=panel['minimum_singular_values'][ik][jp],
                    condition_numbers=panel['condition_numbers'][ik][jp])
                    for jp, (p, angle) in enumerate(zip(PHASES, (0., math.pi)))]
                return dict(kappa=node['kappa_rho2'], phases=phases,
                    reused_from=dict(file=BASELINE, cutoff=record['fine']['config']['nmax'],
                        m=m, recipe='fine', panel_index=ip, node_index=ik),
                    matrices_available=False)
            row = extract(high)
            if low is not None:
                row['lower_cut_response'] = extract(low)
            sector['samples'].append(annotate(row, node['x'], geom['ell']))
        group['sectors'].append(sector)
    group['matrix_scope'] = ('The chosen saved quadrature nodes contain scalar kernels, not matrices. '
        'Existing same-mesh rho2,kappa2 matrices remain in '+BASELINE+
        ' at complete_cuts[1].sectors; they are a different frequency and are not substituted here.')
    group['sums'] = sum_group(group)
    return group


def sum_group(group):
    sums = {}
    for nmax in (2, 4):
        sectors = [s for s in group['sectors'] if abs(s['m'])<=nmax]
        if [s['m'] for s in sectors] != list(range(-nmax, nmax+1)):
            raise ValueError('Missing member of a complete finite multiplet')
        rows = []
        for index in range(len(sectors[0]['samples'])):
            selected = [s['samples'][index] if nmax==4 else s['samples'][index]['lower_cut_response']
                for s in sectors]
            rows.append(dict(x=selected[0]['x'], kappa=selected[0]['kappa'],
                m_values=list(range(-nmax, nmax+1)), nmax=nmax,
                **{kind: {p: {f: np.sum([r[kind][p][f] for r in selected], axis=0).tolist()
                    for f in FIELDS} for p in PHASES} for kind in ('kernels', 'scaled_kernels')}))
        sums[str(nmax)] = rows
    return sums


def weighted_reprojection_check(row, mu):
    """Check the Gram-weighted restriction using independent label matrices."""
    values = row['matrices']
    high = _decode(values['S_connected'])
    fout = block_diag(_decode(values['F_mouth_right_out']), _decode(values['F_mouth_left_out']))
    w = fout.conj().T @ fout
    hi, ho = port_labels(mu, 4)
    li, lo = port_labels(mu, 2)
    pi = [hi.index(v) for v in li]
    po = [ho.index(v) for v in lo]
    count = len(hi)
    ins = np.r_[pi, count+np.asarray(po)]
    outs = np.r_[po, count+np.asarray(pi)]
    P = np.eye(len(high))[outs]
    Q = np.eye(len(high))[ins]
    low, info = reproject_environment(high, mu, 4, 2,
        _decode(values['F_mouth_right_out']), _decode(values['F_mouth_left_out']))
    residual = P@w@P.T@low-P@w@high@Q.T
    recomputed = source_response(low, row['kappa'], mu, 2, row['rho_c'],
        geometry(row['rho_c'])['probes'])
    kernel_error = max(float(np.max(abs(np.asarray(a[f+'_integrands'])-b[f+'_integrands'])))
        for a,b in zip(recomputed['phases'], row['lower_cut_response']['phases']) for f in FIELDS)
    return dict(weighted_equation_error=float(np.max(abs(residual))), source_kernel_error=kernel_error,
        naive_submatrix_difference=info['naive_submatrix_difference'],
        scope='No second exterior solve; independent low-cut witness was already established in frozen101.')


def tight_radial_checks(row):
    """Refine only local radial/cap ODEs while holding the saved exterior fixed."""
    rho = 32.
    geom = geometry(rho)
    ell, kappa = geom['ell'], row['kappa']
    modes, caps = [], []
    for n in (0, 1, 4):
        base, db, diag = mode_variation(1j*kappa, n, ell, rho)
        fine, df, fd = mode_variation(1j*kappa, n, ell, rho, rtol=3e-13, atol=3e-15)
        entry = dict(n=n, scattering_change=float(np.max(abs(base-fine))),
            tangent_change=float(np.max(abs(db-df))),
            log_transmission_change=float(abs(_decode(diag['log_transmission'])-_decode(fd['log_transmission']))),
            baseline=diag, tightened=fd)
        if n==0:
            t = math.exp(-kappa*geom['optical_length'])
            exact = np.array([[0.,t],[t,0.]], complex)
            exact_d = np.array([2*kappa*p['I_f']*exact for p in geom['probes']])
            entry.update(exact_scattering_error=float(np.max(abs(fine-exact))),
                exact_tangent_error=float(np.max(abs(df-exact_d))))
        modes.append(entry)
        lam = math.sqrt(n*(n+1))
        value, cd = cap_reflection(1j*kappa, lam, ell, rho)
        better, cf = cap_reflection(1j*kappa, lam, ell, rho,
            start=5e-6, rtol=3e-13, atol=3e-15)
        cap = dict(n=n, cutoff_tolerance_change=float(abs(value-better)),
            baseline=cd, tightened=cf, reflection=_complex_json(value))
        if n:
            ground, _ = cap_reflection((lam+.5)/ell, lam, ell, rho,
                start=5e-6, rtol=3e-13, atol=3e-15)
            cap['ground_mode_error'] = float(abs(ground-np.exp(-1j*geom['sigma_cut'])))
            # Independent normalizable linear column, removing delta**lambda.
            start, end = 5e-6, cf['delta_end']
            scale = kappa*ell
            den = 1+2*lam
            y = 1-2*scale/den*start+2*scale**2/den**2*start**2
            y += (2*lam*scale/(3*den)-4*scale**3/den**3)/(3+2*lam)*start**3
            s3 = np.diag([1.,-1.])
            s2 = np.array([[0.,-1j],[1j,0.]])
            def rhs(delta, column):
                return (scale*s3-lam/math.sin(delta)*s2-lam/delta*np.eye(2))@column
            solved = solve_ivp(rhs, (start,end), np.array([1.,-1j*y]),
                method='DOP853',rtol=3e-13,atol=3e-15)
            if not solved.success:
                raise RuntimeError('Independent cap-column control failed: '+solved.message)
            direct = solved.y[1,-1]/solved.y[0,-1]
            cap['linear_column_error'] = float(abs(better-direct))
        else:
            cap['exact_massless_error'] = float(abs(better+1j*math.exp(-2*kappa*ell*cf['delta_end'])))
        caps.append(cap)

    connected = _decode(row['matrices']['S_connected'])
    raw, pressure, _ = raw_variation(kappa, N, rho, rtol=3e-13, atol=3e-15)
    shift = shift_derivatives(raw, kappa, N, geom['probes'])
    reference, _ = complete_reference(kappa, .5, N, rho)
    reference = reference.copy()
    count = 2*N+1
    for n in range(1,N+1):
        value, _ = cap_reflection(1j*kappa, math.sqrt(n*(n+1)), ell, rho,
            start=5e-6,rtol=3e-13,atol=3e-15)
        reference[n-1,n] = value
        reference[N+n,N+n] = value
        reference[count+n,count+n-1] = value.conjugate()
        reference[count+N+n,count+N+n] = value.conjugate()
    fine_phases = []
    for angle in (0.,math.pi):
        b = selected_map(raw,.5,angle,N)
        dp = [selected_map(v,.5,angle,N) for v in pressure]
        ds = [selected_map(v,.5,angle,N) for v in shift]
        fine_phases.append(reference_relative(connected,reference,block_diag(b,b.conj().T),
            [block_diag(v,v.conj().T) for v in dp],
            [block_diag(v,-v.conj().T) for v in ds]))
    changes = {p:{f:(np.asarray(fine_phases[j][f+'_integrands'])-
        row['phases'][j][f+'_integrands']).tolist() for f in FIELDS} for j,p in enumerate(PHASES)}
    return dict(rho_c=rho,kappa=kappa,x=row['x'],mu=.5,h=.05,nmax=N,
        modes=modes,caps=caps,tightened_source_phases=fine_phases,
        tightened_minus_baseline=changes,
        maximum_source_change=max(abs(v) for fields in changes.values() for values in fields.values() for v in values),
        exterior_rebuilt=False,scope='Observed local ODE/Frobenius refinement and independent analytic/linear controls; not interval bounds.')


def differences(left, right):
    return {key:[dict(x=a['x'],**{kind:{p:{f:(np.asarray(a[kind][p][f])-b[kind][p][f]).tolist()
        for f in FIELDS} for p in PHASES} for kind in ('kernels','scaled_kernels')})
        for a,b in zip(left['sums'][key],right['sums'][key])] for key in ('2','4')}


def write_result(result):
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if json.loads(OUTPUT.read_text(encoding='utf-8'))!=result:
        raise ValueError('Saved result did not round-trip exactly')


def validate(result):
    """Numerical consistency gates; no source sign is an acceptance test."""
    all_groups = result['geometries']+[result['mesh_control']]
    equation = coordinate = null_error = resolvent = reference_flow = 0.
    min_singular = float('inf')
    max_condition = 0.
    port_error = assembly_error = scaled_identity = 0.
    min_mass = float('inf')
    gram_condition = 0.
    for group in all_groups:
        ell = group['geometry']['ell']
        for sector in group['sectors']:
            if 'assembly' in sector:
                assembly_error=max(assembly_error,sector['assembly']['opposite_adjoint_relative_error'])
                min_mass=min(min_mass,sector['assembly']['minimum_cell_mass_eigenvalue'])
                for p in sector['ports']:
                    port_error=max(port_error,*[p[k] for k in ('basis_reproduction_error',
                        'incoming_adjoint_error','outgoing_adjoint_error','factor_swap_error')])
                    gram_condition=max(gram_condition,*p['gram_condition_numbers'])
            for original in sector['samples']:
                for row in [original]+([original['lower_cut_response']] if 'lower_cut_response' in original else []):
                    equation=max(equation,*row.get('exterior_equation_relative_residuals',[0.]))
                    coordinate=max(coordinate,row.get('raw_facet_observable_error') or 0.)
                    scaled_identity=max(scaled_identity,row['scaled_frequency_identity_error'])
                    for phase in row['phases']:
                        p=np.asarray(phase['pressure_integrands']);s=np.asarray(phase['shift_integrands'])
                        null_error=max(null_error,float(np.max(abs(np.asarray(phase['null_plus_integrands'])-p+2*s))),
                            float(np.max(abs(np.asarray(phase['null_minus_integrands'])-p-2*s))))
                        resolvent=max(resolvent,phase.get('resolvent_difference_identity_error',0.))
                        reference_flow=max(reference_flow,max(abs(v) for v in phase.get('reference_shift_trace',[0.])))
                        min_singular=min(min_singular,*phase['minimum_singular_values'])
                        max_condition=max(max_condition,*phase['condition_numbers'])
        for rows in group['sums'].values():
            for row in rows:
                for p in PHASES:
                    for f in FIELDS:
                        scaled_identity=max(scaled_identity,float(np.max(abs(
                            np.asarray(row['scaled_kernels'][p][f])-np.asarray(row['kernels'][p][f])/ell))))
    checks=result['weighted_reprojection_checks']
    weight_error=max(r['weighted_equation_error'] for r in checks)
    weight_kernel=max(r['source_kernel_error'] for r in checks)
    radial=result['radial_checks']
    mode_change=max(max(r['scattering_change'],r['tangent_change'],r['log_transmission_change']) for r in radial['modes'])
    cap_change=max(max(r.get(k,0.) for k in ('cutoff_tolerance_change','ground_mode_error','linear_column_error','exact_massless_error'))
        for r in radial['caps'])
    zero=radial['modes'][0]
    tests=dict(baseline101_passed=result['baseline_controls']['passed'],
        baseline_exact_reuse=result['baseline_controls']['exact_reuse_error']==0.,
        complete_multiplets=all(len(g['sectors'])==9 for g in all_groups),
        exterior_equations=equation<1e-8,opposite_weyl_assembly=assembly_error<1e-11,
        positive_imaginary_mass=min_mass>0.,ports_and_gram=port_error<1e-10 and gram_condition<1e8,
        raw_facet_coordinates=coordinate<1e-8,sewing_nonsingular=min_singular>1e-10,
        resolvent_identity=resolvent<1e-8,reference_flow_zero=reference_flow<1e-8,
        null_identities=null_error<1e-11,scaled_frequency_and_jacobian=scaled_identity<1e-11,
        exact_weighted_reprojection=weight_error<1e-10 and weight_kernel<1e-10,
        radial_ode_refinement=mode_change<1e-8,normalizable_cap_controls=cap_change<1e-9,
        exact_zero_mode=max(zero['exact_scattering_error'],zero['exact_tangent_error'])<1e-9,
        radial_source_refinement=radial['maximum_source_change']<1e-8,
        inputs_preserved=result['input_sha256']==hashes(result['input_sha256']),
        execution_source_preserved=result['source_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result['tests']={k:bool(v) for k,v in tests.items()}
    result['passed']=all(tests.values())
    result['diagnostics']=dict(maximum_exterior_equation_residual=equation,maximum_coordinate_error=coordinate,
        maximum_null_identity_error=null_error,maximum_resolvent_identity_error=resolvent,
        minimum_sewing_singular_value=min_singular,maximum_sewing_condition=max_condition,
        maximum_port_error=port_error,maximum_gram_condition=gram_condition,
        maximum_assembly_adjoint_error=assembly_error,minimum_cell_mass_eigenvalue=min_mass,
        maximum_scaled_frequency_or_jacobian_error=scaled_identity,
        maximum_weighted_reprojection_error=weight_error,maximum_reprojected_kernel_error=weight_kernel,
        maximum_radial_mode_change=mode_change,maximum_cap_control_error=cap_change,
        maximum_radial_source_change=radial['maximum_source_change'],
        fresh_exterior_frequency_samples=63,N2_extra_exterior_solves=0,
        frequency_integral_computed=False,omitted_frequency_tail_bound=None,
        omitted_angular_tail_bound=None,continuum_error_bound=None)


def hashes(names):
    return {name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in names}


def run_checks():
    started=time.perf_counter()
    prior=json.loads(Path(BASELINE).read_text(encoding='utf-8'))
    names=sorted(set(prior['integration_input_sha256'])|{'complete_angular_source.py',BASELINE})
    result=dict(passed=False,tests={},stage='Preparing exact reused baseline',seconds=0.,
        input_sha256=hashes(names),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parameters=dict(q=1,b=B,d=1.,R=2.,L=40.,nmax=4,lower_nmax=2,outer_nmax=J,
            primary_h=.05,refined_h=.035,rho_values=list(RHO_VALUES),primary_phase=math.pi,
            phase_settings=[0.,math.pi],units='hbar=c=d=1',sign_is_not_an_acceptance_gate=True),
        geometries=[],weighted_reprojection_checks=[],
        scope=[
            'A source-only sensitivity test on prescribed sharp-join geometries; no solved collar, shell material, formation, stability, or physical support is inferred.',
            'q=1 is not the large-q MMP regime. Changing rho changes physical redshift, proper/optical length and the AdS reference with the endpoint lapse held one.',
            'Both actual Weyl responses and all m=-N..N for N2/N4 are retained; these finite cutoffs do not complete the angular sum.',
            'The reflected spin/gauge endpoint identification is a declared matching branch, not a derived physical junction.',
            'Probe profiles are fixed fractions of sigma_cut: their proper support and geometric weights change with rho. All such data are retained.',
            'Each phase compares the connected environment with its same-local-geometry invariant AdS reference. This is not a phase-pi-minus-zero result.',
            'Three matched x=kappa*ell samples do not determine a frequency integral or its sign. Scaled kernels include the eventual d-kappa=dx/ell Jacobian; no quadrature weights are used.',
            'The unchanged constant-radius core has G_kk=0. A nonzero completed quantum source would require geometry response, not support the unchanged cylinder.',
            'New numerical sources are still finite dissipative DG approximations. Observed mesh, angular-cut and ODE changes are not certified continuum or tail bounds.'],
        observable_definitions=dict(pressure='Smeared connected-minus-AdS-reference rho+p_r kernel per d-kappa',
            shift='Negative smeared lower orthonormal radial flux kernel per d-kappa',
            null_plus='pressure-2*shift',null_minus='pressure+2*shift',
            scaled_kernels='kernels/ell, for a future integral in x=kappa*ell; no such integral is computed here'),
        sources=prior.get('sources',[]))
    try:
        if not prior['passed'] or not all(prior['tests'].values()):
            raise ValueError('The frozen101 baseline is not passing')
        if prior['source_sha256']!=result['input_sha256']['complete_angular_source.py']:
            raise ValueError('Frozen101 source provenance changed')
        nodes=baseline_nodes(prior)
        result['parameters']['common_scaled_nodes']=nodes
        base=baseline_group(prior,nodes)
        result['geometries'].append(base)
        reuse_error=0.
        for key in ('2','4'):
            panels=prior['integration']['cutoffs'][key]['partial_sums']['fine']['panels']
            for row,node in zip(base['sums'][key],nodes):
                panel=panels[node['panel_index']]
                for p in PHASES:
                    for f in FIELDS:
                        reuse_error=max(reuse_error,float(np.max(abs(np.asarray(row['kernels'][p][f])-
                            panel['integrands'][p][f][node['node_index']]))))
        result['baseline_controls']=dict(passed=prior['passed'],tests=prior['tests'],
            exact_reuse_error=reuse_error,exact_weighted_reuse_witness=prior['cut_reuse_checks'],
            current_baseline_source_hash=prior['source_sha256'],old_matrix_witness_path=
                BASELINE+'::complete_cuts[1].sectors[*].samples[kappa=2].matrices',
            absent_middle_node_matrix_scope='Frozen fine-node kernels contain no matrices; no extra rho2 exterior solve was authorized or performed.')
        write_result(result)
        for rho in (8.,32.):
            geom=geometry(rho)
            group=dict(geometry=geom,h=.05,nmax=N,outer_nmax=J,sectors=[],origin='Fresh paired-Weyl exterior responses at matched x')
            result['geometries'].append(group)
            for m in range(-4,5):
                tick=time.perf_counter()
                case=CompleteSector(m+.5,nmax=N,outer_nmax=J,rho_c=rho)
                sector=dict(m=m,mu=m+.5,config=case.config,assembly=case.assembly,ports=case.port_checks,samples=[])
                group['sectors'].append(sector)
                for index,node in enumerate(nodes):
                    row=case.evaluate(node['x']/geom['ell'],save_matrices=index==1,lower_nmax=2)
                    row['rho_c']=rho
                    sector['samples'].append(annotate(row,node['x'],geom['ell']))
                    if index==1 and abs(m)<=2:
                        result['weighted_reprojection_checks'].append(dict(rho_c=rho,h=.05,m=m,
                            **weighted_reprojection_check(row,m+.5)))
                result.update(stage=f'Completed rho={rho:g}, m={m}',seconds=time.perf_counter()-started)
                write_result(result)
                print(json.dumps(dict(stage='sector_complete',rho_c=rho,m=m,seconds=time.perf_counter()-tick,
                    middle_pi=sector['samples'][1]['kernels']['pi']['pressure'])),flush=True)
                del case
            group['sums']=sum_group(group)
            print(json.dumps(dict(stage='rho_complete',rho_c=rho,source_pi=[r['kernels']['pi']['pressure'] for r in group['sums']['4']],
                scaled_pi=[r['scaled_kernels']['pi']['pressure'] for r in group['sums']['4']])),flush=True)
        middle=nodes[1]
        geom=geometry(32.)
        small=dict(geometry=geom,h=.035,nmax=N,outer_nmax=J,sectors=[],origin='Matched finer mesh at the middle scaled frequency')
        result['mesh_control']=small
        for m in range(-4,5):
            case=CompleteSector(m+.5,h=.035,nmax=N,outer_nmax=J,rho_c=32.)
            row=case.evaluate(middle['x']/geom['ell'],save_matrices=True,lower_nmax=2)
            row['rho_c']=32.
            annotate(row,middle['x'],geom['ell'])
            small['sectors'].append(dict(m=m,mu=m+.5,config=case.config,assembly=case.assembly,ports=case.port_checks,samples=[row]))
            if abs(m)<=2:
                result['weighted_reprojection_checks'].append(dict(rho_c=32.,h=.035,m=m,
                    **weighted_reprojection_check(row,m+.5)))
            result.update(stage=f'Completed finer mesh m={m}',seconds=time.perf_counter()-started)
            write_result(result)
            print(json.dumps(dict(stage='finer_mesh_complete',m=m,seconds=row['seconds'])),flush=True)
            del case
        small['sums']=sum_group(small)
        primary=result['geometries'][2]
        midpoint=dict(sums={key:[rows[1]] for key,rows in primary['sums'].items()})
        result['mesh_control']['finer_minus_primary']=differences(small,midpoint)
        mzero=next(s for s in primary['sectors'] if s['m']==0)['samples'][1]
        result['radial_checks']=tight_radial_checks(mzero)
        result['rho_changes_from_2']={str(int(g['geometry']['rho_c'])):differences(g,base) for g in result['geometries'][1:]}
        for group in result['geometries']+[small]:
            group['N4_minus_N2']=[dict(x=a['x'],**{kind:{p:{f:(np.asarray(a[kind][p][f])-b[kind][p][f]).tolist()
                for f in FIELDS} for p in PHASES} for kind in ('kernels','scaled_kernels')})
                for a,b in zip(group['sums']['4'],group['sums']['2'])]
        validate(result)
        result.update(stage='Completed selected-frequency long-throat source pilot',seconds=time.perf_counter()-started)
        write_result(result)
    except Exception as exc:
        result.update(passed=False,stage='Stopped with saved diagnostics',failure=dict(type=type(exc).__name__,message=str(exc)),
            seconds=time.perf_counter()-started)
        write_result(result)
        raise
    return result


if __name__=='__main__':
    data=run_checks()
    print(json.dumps(dict(passed=data['passed'],tests=data['tests'],seconds=data['seconds'],diagnostics=data.get('diagnostics'))),flush=True)
    if not data['passed']:
        raise SystemExit('A numerical gate failed; all calculated diagnostics were saved.')
