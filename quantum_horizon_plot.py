"""Render the horizon geometry, observer data, and scalar source test."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parent
    record = json.loads((root / 'coupled_ori_results.json').read_text())['quantum_horizon_geometry']
    assert record['observer']['unit_norm_and_four_geodesic_equations']
    a, w = .02, 0.
    mu = -1-2*a*a/9
    B = np.diag([-.25-a/6, -.25+a/6])
    q = 1/np.sqrt(2)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12, 6.3))
    left = fig.add_subplot(121, projection='3d')
    right = fig.add_subplot(122)
    fig.subplots_adjust(left=.05, right=.97, top=.76, bottom=.24, wspace=.15)
    fig.suptitle('A quantum obstacle at the first time loops', fontsize=21, fontweight='bold', y=.96)
    fig.text(.5, .87, 'Complete ideal geometry • selected negative-source branch • untwisted example', ha='center', color='#475569')
    az = np.linspace(0, 2*np.pi, 160)
    left.plot(3*np.cos(az), 3*np.sin(az), az*0, color='#a96b00', lw=2.4, label='Central null circle')
    angle, section = np.meshgrid(np.linspace(0, 2*np.pi, 40), np.linspace(0, 2*np.pi, 12))
    left.plot_wireframe((3+.9*np.cos(section))*np.cos(angle),
                        (3+.9*np.cos(section))*np.sin(angle), .9*np.sin(section),
                        color='#718096', alpha=.12, lw=.45)
    past = np.linspace(0, 4*np.pi, 400)
    for z0, color in zip((0., 2.1, 4.2), ('#007e87', '#326bb3', '#7951b0')):
        r = np.exp(2*np.diag(B)[:, None]*past)*np.array([[.7], [.5]])
        zz = z0-past
        curve = np.array([(3+r[0])*np.cos(zz), (3+r[0])*np.sin(zz), r[1]])
        left.plot(*curve, color=color, lw=1.7)
        i = 20; delta = curve[:, i+4]-curve[:, i]
        left.quiver(*curve[:, i], *delta, color=color, arrow_length_ratio=.6, linewidth=1.8)
    left.set_box_aspect((1, 1, .42)); left.view_init(elev=30, azim=-62); left.set_axis_off()
    left.set_title('Trace boundary light paths into the past:\nthey approach the same closed circle', pad=8)
    left.legend(loc='lower center', frameon=False, fontsize=9)
    lam = np.linspace(-2, 0, 250)
    density = -4*mu/(4*q-lam)**2
    tidal = np.sqrt(8*(mu*mu+a*a))/(4*q-lam)**2
    assert np.all(np.isfinite(density)) and np.all(np.isfinite(tidal))
    right.plot(lam, density, color='#7951b0', lw=2, label=r'$|\rho|\,\kappa$')
    right.plot(lam, tidal, color='#007e87', lw=2, label='Tidal curvature norm')
    right.axvline(0, color='#9f3345', ls='--', lw=1.2, label='First null circle')
    right.set(xlim=(-2, .06), ylim=(0, .58), xlabel='Observer proper time λ (geometry units)',
              ylabel='Normalized classical quantities', title='Classical measurements stay finite')
    right.spines[['top', 'right']].set_visible(False); right.grid(alpha=.18)
    right.legend(loc='upper left', frameon=False, fontsize=9)
    fig.text(.07, .15, 'Left: mathematical embedding of a bounded horizon patch; arrows point into the past. No physical tunnel is shown.', fontsize=9, color='#475569')
    fig.text(.07, .105, 'KRW theorem: initially regular free-scalar quantum correlations cannot remain Hadamard at the central circle.', fontsize=10, color='#9f3345')
    fig.text(.07, .065, 'Finite classical curves do not establish a regular quantum crossing. The theorem does not give an injury prediction.', fontsize=9, color='#475569')
    fig.text(.07, .025, '13 September 2027 → 12 September 2027: no demonstrated route.', fontsize=10, fontweight='bold', color='#9f3345')
    for suffix in ('png', 'pdf'):
        fig.savefig(root / f'quantum_horizon_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root / 'quantum_horizon_progress.png')


def source_progress():
    root = Path(__file__).resolve().parent
    record = json.loads((root/'ori_quantum_source_results.json').read_text())
    a = float(record['TT_coefficient'])
    conformal = float(record['nonminimal_coupling']['conformal_TT_coefficient'])
    assert a > 0 and conformal > a
    fig = plt.figure(figsize=(12, 6.3))
    left = fig.add_subplot(121, projection='3d')
    right = fig.add_subplot(122)
    fig.subplots_adjust(left=.04, right=.97, top=.77, bottom=.25, wspace=.2)
    fig.suptitle('A field state is possible here. Its stress does not fit.',
                 fontsize=20, fontweight='bold', y=.96)
    fig.text(.5, .87, 'Fixed ideal geometry • one periodic massless scalar • before the first time loops',
             ha='center', color='#475569')
    az, cross = np.meshgrid(np.linspace(0, 2*np.pi, 65), np.linspace(0, 2*np.pi, 25))
    left.plot_wireframe((3+.7*np.cos(cross))*np.cos(az),
                        (3+.7*np.cos(cross))*np.sin(az), .7*np.sin(cross),
                        color='#007e87', alpha=.3, lw=.65)
    circle = np.linspace(0, 2*np.pi, 160)
    left.plot(3*np.cos(circle), 3*np.sin(circle), np.zeros_like(circle),
              color='#9f3345', lw=2.5, label='Central circle at the horizon')
    left.set_box_aspect((1, 1, .4)); left.view_init(elev=29, azim=-61)
    left.set_axis_off()
    left.set_title('Visual map of the assumed geometry\nThe field has not been shown to support it', pad=8)
    left.legend(loc='lower center', frameon=False, fontsize=9)
    T = np.linspace(-1, -.1, 250)
    right.axhline(0, color='#007e87', lw=2, label='Required by this geometry: zero')
    right.plot(T, 1/T**3, color='#9f3345', lw=2, label='Calculated scalar: ξ = 0')
    right.plot(T, conformal/a/T**3, color='#7951b0', lw=2,
               label='Zero-trace choice: ξ = 1/6')
    right.axvline(0, color='#475569', ls='--', lw=1)
    right.set_yscale('symlog', linthresh=1)
    right.set(xlim=(-1.02, .045), ylim=(-1e4, 2), xlabel='Axis coordinate T (model units)',
              ylabel=r'Null stress component $T_{TT}/A_0$',
              title='One necessary support condition fails')
    right.spines[['top', 'right']].set_visible(False); right.grid(alpha=.18)
    right.legend(loc='lower left', frameon=False, fontsize=8.7)
    right.text(.98, .98, 'Horizon: T = 0', transform=right.transAxes,
               ha='right', va='top', fontsize=9, color='#475569')
    fig.text(.06, .16, 'Left: a bounded patch drawn with the periodic coordinate wrapped into a circle; not a material tube or an entrance.',
             fontsize=9, color='#475569')
    fig.text(.06, .113, 'Result: no curvature-coupling value makes this particular state supply the full required stress tensor.',
             fontsize=10, color='#9f3345')
    fig.text(.06, .07, 'The graph is a fixed-background calculation, not a collapse simulation or a prediction of what a person would feel.',
             fontsize=9, color='#475569')
    fig.text(.06, .025, '13 September 2027 → 12 September 2027: still no demonstrated route.',
             fontsize=10, fontweight='bold', color='#9f3345')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'ori_quantum_source_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root/'ori_quantum_source_progress.png')


def response_progress():
    root = Path(__file__).resolve().parent
    data = json.loads((root/'ori_quantum_source_results.json').read_text())
    response = data['first_order_response']
    eta, q, mu = 1e-12, 1/np.sqrt(2), -11251/11250
    reference = float(response['observer_cases']['minimal']['events'][2]['delta'])
    scaled = np.logspace(2,-2,350)
    delta = reference*scaled
    s = q*delta+delta**2/4
    fig, axes = plt.subplots(1,2,figsize=(12,6.5))
    fig.subplots_adjust(left=.08,right=.97,top=.75,bottom=.25,wspace=.3)
    fig.suptitle('A tiny light-bundle change does not establish stability',
                 fontsize=19,fontweight='bold',y=.96)
    fig.text(.5,.875,'Calculated quantum energy grows faster than this first-order area response',
             ha='center',color='#475569',fontsize=11)
    fig.text(.5,.825,'Model: η = 10⁻¹², observer q = 1/√2; initial beam area = 1 and expansion = 0 at T = −1',
             ha='center',color='#475569',fontsize=9)
    for name,color,label,coefficients in (
        ('minimal','#9f3345','Minimal coupling',data['full_axis_coefficients']),
        ('conformal','#7951b0','Conformal coupling',data['NEC_compensation']['conformal_axis_coefficients'])):
        a,b,d = [float(coefficients[k]) for k in ('TT','Tz','zz')]
        p = a*(2*q+delta)**2-2*b*delta*(2*q+delta)+d*delta**2
        ratio = 4*eta*p/(abs(mu)*delta**3*(4*q+delta))
        area_change = eta*a*(1-s)**2/(2*s)
        marker = response['observer_cases'][name]['events'][2]
        assert float(marker['eta'])==eta and np.all((s>0)&(s<=1))
        xmarker = float(marker['delta'])/reference
        assert np.all(area_change>0) and np.all(ratio>0)
        axes[0].loglog(scaled,ratio,color=color,lw=2,label=label)
        axes[1].loglog(scaled,area_change,color=color,lw=2,label=label)
        axes[0].scatter([xmarker],[1],color=color,zorder=5,s=30)
        axes[1].scatter([xmarker],[float(marker['null_area_first_order_change'])],
                        color=color,zorder=5,s=30)
    axes[0].axhline(1,color='#007e87',ls='--',lw=1.5,label='Equal source magnitudes')
    axes[0].set(title='Freely falling observer: energy comparison',
                ylabel='|Quantum energy| / |prescribed classical energy|')
    axes[1].set(title='Null bundle: necessary optical response',
                ylabel='Fractional area change (first order)')
    for ax in axes:
        ax.set(xlim=(100,.01),xlabel='Remaining proper time / minimal-coupling equality time')
        ax.spines[['top','right']].set_visible(False)
        ax.grid(alpha=.18,which='both'); ax.legend(frameon=False,fontsize=9)
    fig.text(.06,.16,'Dots mark equal quantum and classical energy magnitudes. The first-order area change there is only about 10⁻¹².',
             fontsize=9,color='#475569')
    fig.text(.06,.112,'Both calculations start from the prescribed state and geometry. This is not a simulation of a stable medium.',
             fontsize=10,color='#9f3345')
    fig.text(.06,.067,'Right: the z-invariant null congruence is sampled along the observer’s path; the observer does not travel with one beam.',
             fontsize=9,color='#475569')
    fig.text(.06,.025,'13 September 2027 → 12 September 2027: no demonstrated crossing or calendar mapping.',
             fontsize=10,fontweight='bold',color='#9f3345')
    for suffix in ('png','pdf'):
        fig.savefig(root/f'ori_quantum_response_progress.{suffix}',dpi=170,facecolor='white')
    plt.close(fig)
    print(root/'ori_quantum_response_progress.png')


def probe_future_progress():
    root = Path(__file__).resolve().parent
    probe = json.loads((root/'ori_quantum_source_results.json').read_text())['abiotic_test_probe']
    trip = json.loads((root/'future_time_results.json').read_text())['cases']['tenfold_one_lab_day']
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.8))
    fig.subplots_adjust(left=.075, right=.97, top=.76, bottom=.29, wspace=.31)
    fig.suptitle('Tiny payload first: what the calculations establish',
                 fontsize=19, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Past: a test on an assumed geometry     •     Future: a route using ordinary spacetime',
             ha='center', fontsize=11, color='#475569')
    lam, q = np.linspace(-1, 0, 300), 1/np.sqrt(2)
    for name, color in (('x', '#007e87'), ('y', '#7951b0')):
        record = probe['directions'][name]
        nu = float(record['nu']); r = (4*q-lam)/(4*q+1)
        f = ((1+nu)*r**((1-nu)/2)+(nu-1)*r**((1+nu)/2))/(2*nu)
        assert abs(f[-1]-float(record['classical_boundary_separation_factor']))<1e-12
        axes[0].plot(lam, 100*(f-1), lw=2, color=color,
                     label=f'{name} direction: {100*(f[-1]-1):.2f}%')
    axes[0].axvline(0, color='#9f3345', ls='--', lw=1.5)
    axes[0].set(title='Unbound tracer pair in the past-travel model',
                xlabel='Proper-time parameter λ (classical boundary at 0)',
                ylabel='Change in tracer separation (%)', xlim=(-1.02, .035))
    axes[0].legend(frameon=False, loc='upper left', fontsize=9)
    c, au = 299792458., 149597870700.
    h, a = float(trip['peak_rapidity']), float(trip['proper_acceleration_m_per_s2'])
    t0, x0 = 0., 0.
    for theta0, sign in ((0, 1), (h, -1), (0, -1), (-h, 1)):
        theta = theta0+sign*h*np.linspace(0, 1, 100)
        t = t0+c/(sign*a)*(np.sinh(theta)-np.sinh(theta0))
        x = x0+c*c/(sign*a)*(np.cosh(theta)-np.cosh(theta0))
        axes[1].plot(x/au, t/3600, color='#326bb3', lw=2.5)
        t0, x0 = t[-1], x[-1]
    assert abs(t0-86400)<1e-7 and abs(x0/au)<1e-10
    axes[1].plot([0, 0], [0, 24], ls='--', color='#475569', lw=1.5)
    axes[1].scatter([0, 0], [0, 24], color='#007e87', s=35, zorder=4)
    axes[1].text(4, -.65, '13 Sep: depart', fontsize=9, color='#475569')
    axes[1].text(4, 24.25, '14 Sep: reunite at rest', fontsize=9, color='#475569')
    axes[1].text(28, 12, 'Probe clock\n2.4 hours', ha='center', va='center',
                 fontsize=12, color='#326bb3', fontweight='bold')
    axes[1].set(title='Forward trip: station clock advances 24 hours',
                xlabel='Distance from station (AU; 1 AU ≈ Earth–Sun distance)',
                ylabel='Time at station (hours)', ylim=(-1, 25), xlim=(-4, 100))
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.17)
    fig.text(.06, .21, 'Left: formal 1 µm rescaling gives +86 nm / +83 nm. Valid as a small-separation approximation only; no bonded-grain model.',
             fontsize=9, color='#475569')
    fig.text(.06, .158, 'Right: prescribed acceleration ≈ 63,700g; peak speed ≈ 0.999753c. Forces, propulsion and payload survival remain unmodelled.',
             fontsize=9, color='#475569')
    fig.text(.06, .104, 'Past route: no validated crossing. Forward route: exact ideal motion using experimentally supported time dilation.',
             fontsize=10, color='#475569')
    fig.text(.06, .041, '13 → 12 September 2027: unproved.     13 → 14 September: less probe time is physically allowed; this flight is not demonstrated.',
             fontsize=10, color='#007e87', fontweight='bold')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'probe_future_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root/'probe_future_progress.png')


def proton_ring_progress():
    root = Path(__file__).resolve().parent
    data = json.loads((root/'future_time_results.json').read_text())['proton_ring']
    nominal, radiation = data['nominal'], data['radiation']
    radius, gamma = nominal['bending_radius_metres'], nominal['gamma']
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.7))
    fig.subplots_adjust(left=.075, right=.97, top=.76, bottom=.32, wspace=.32)
    fig.suptitle('A proton can take the forward route in a ring',
                 fontsize=20, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Calculated ideal orbit using magnetic force • not a complete accelerator design',
             ha='center', color='#475569', fontsize=11)
    angle = np.linspace(0, 2*np.pi, 401)
    axes[0].plot(radius*np.cos(angle), -radius*np.sin(angle), color='#007e87', lw=3)
    axes[0].scatter([radius], [0], color='#326bb3', s=65, zorder=5)
    for phi in (np.pi/3, np.pi, 5*np.pi/3):
        start = np.array([radius*np.cos(phi), -radius*np.sin(phi)])
        end = np.array([radius*np.cos(phi+.15), -radius*np.sin(phi+.15)])
        axes[0].annotate('', xy=end, xytext=start,
                         arrowprops=dict(arrowstyle='->', color='#007e87', lw=2))
    axes[0].annotate('', xy=(radius, 0), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='<->', color='#475569'))
    axes[0].text(radius/2, 3, f'{radius:.2f} m', ha='center', fontsize=11)
    axes[0].text(-6, -12, 'B = 1 T\nout of page ⊙', ha='center', color='#475569', fontsize=10)
    axes[0].set(title='Uniform-field orbit, viewed from above', xlabel='x (metres)',
                ylabel='y (metres)', xlim=(-40, 40), ylim=(-40, 40), aspect='equal')
    hours = np.linspace(0, 24, 100)
    axes[1].plot(hours, hours, color='#475569', lw=2, label='Lab clock')
    axes[1].plot(hours, hours/gamma, color='#326bb3', lw=2.5, label='Proton proper time')
    axes[1].scatter([24, 24], [24, 24/gamma], color=['#475569', '#326bb3'], s=40, zorder=4)
    axes[1].annotate('24 hours', xy=(24, 24), xytext=(16, 24), color='#475569')
    axes[1].annotate('2.4 hours', xy=(24, 2.4), xytext=(17, 4.5), color='#326bb3', fontweight='bold')
    axes[1].set(title='During constant-energy circulation (γ = 10)',
                xlabel='Time elapsed in lab (hours)', ylabel='Accumulated time (hours)',
                xlim=(0, 25), ylim=(0, 26))
    axes[1].legend(frameon=False, loc='upper left', fontsize=10)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.06, .20, f'Speed: {100*nominal["speed_over_c"]:.4f}% of light speed. Proton kinetic energy: {nominal["payload_kinetic_energy_GeV"]:.2f} GeV; this excludes apparatus energy.',
             fontsize=10, color='#475569')
    fig.text(.06, .146, f'Free-space radiation estimate: {radiation["energy_per_day_eV_if_gamma_held"]/1e6:.3f} MeV per proton per day. Ideal energy replacement holds the time ratio fixed.',
             fontsize=9, color='#475569')
    fig.text(.06, .092, 'Uniform B bends the orbit but does not focus vertical errors. No one-day confinement or complete injection/extraction path is demonstrated.',
             fontsize=9, color='#9f3345')
    fig.text(.06, .037, '13 → 14 September 2027: forward time dilation.     13 → 12 September: still no demonstrated route.',
             fontsize=11, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'proton_ring_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root/'proton_ring_progress.png')


def focusing_feedback_progress():
    root = Path(__file__).resolve().parent
    focus = json.loads((root/'future_time_results.json').read_text())['proton_focusing']
    focused, uniform = focus['cases']['0.4'], focus['cases']['0.0']
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2))
    fig.subplots_adjust(left=.085, right=.97, top=.75, bottom=.33, wspace=.33)
    fig.suptitle('Future route: the focusing field contains a disturbed proton',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Calculated ideal field: 1 T at a 31.14 m orbit radius • focusing index n = 0.4',
             ha='center', fontsize=10, color='#475569')
    for case, label, color in ((uniform, 'Uniform field: vertical drift', '#9f3345'),
                               (focused, 'Focusing field: bounded motion', '#007e87')):
        axes[0].plot(case['plot']['reference_turn'], case['plot']['vertical_mm'],
                     label=label, color=color, lw=2)
    axes[0].axhline(10, ls='--', color='#64748b', lw=1.2, label='Chosen tube limit: +10 mm')
    axes[0].set(title='Same starting offset and small momentum kicks',
                xlabel='Reference revolutions', ylabel='Vertical displacement (mm)', xlim=(0, 10))
    axes[0].legend(frameon=False, fontsize=8, loc='upper left')
    axes[1].plot(focused['plot']['radial_mm'], focused['plot']['vertical_mm'], color='#007e87', lw=1.2)
    axes[1].scatter([1], [.5], color='#326bb3', s=35, zorder=3, label='Start: +1 mm / +0.5 mm')
    axes[1].set(title='Magnified cross-section of the focused orbit',
                xlabel='Radial displacement (mm)', ylabel='Vertical displacement (mm)',
                xlim=(-1.3, 1.3), ylim=(-1.05, 1.05))
    axes[1].legend(frameon=False, fontsize=8, loc='lower left')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.06, .19, '100-turn check: focused motion stays within 1.08 mm radially and 0.71 mm vertically; uniform-field vertical drift reaches 196 mm.',
             fontsize=9, color='#475569')
    fig.text(.06, .13, 'An energy bound also prevents escape in the ideal static model. Radiation, RF, field errors and injection are outside this bound.',
             fontsize=9, color='#475569')
    fig.text(.06, .065, '100 turns ≈ 65.6 microseconds. This is progress toward confinement, not a demonstrated one-day trip.',
             fontsize=10, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'proton_focusing_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)

    data = json.loads((root/'mmp_pulse_results.json').read_text())['boundary_feedback']
    response = data['cases']['0.1']['plot']
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2))
    fig.subplots_adjust(left=.065, right=.97, top=.74, bottom=.34, wspace=.25)
    fig.suptitle('Wormhole model: the supporting quantum stress responds',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'First-interval calculation with moving mouth clocks and a finite exterior delay',
             ha='center', fontsize=11, color='#475569')
    xx = np.linspace(-1, 1, 200); yy = .25+.7*xx**2
    axes[0].plot(xx, yy, color='#326bb3', lw=2)
    axes[0].plot(xx, -yy, color='#326bb3', lw=2)
    axes[0].fill_between(xx, -yy, yy, color='#326bb3', alpha=.08)
    for x, label in ((-1, 'Mouth L'), (1, 'Mouth R')):
        axes[0].plot([x, x], [-.95, .95], color='#007e87', lw=4)
        axes[0].text(x, -1.13, label, ha='center', fontsize=10, color='#007e87')
    axes[0].annotate('', xy=(-.8, 0), xytext=(.8, 0),
                     arrowprops=dict(arrowstyle='->', color='#b7791f', lw=2.5))
    axes[0].text(0, .11, 'Injected signal', ha='center', fontsize=10, color='#b7791f')
    axes[0].annotate('', xy=(1, 1), xytext=(-1, 1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-.5',
                                     ls='--', color='#7951b0', lw=1.5))
    axes[0].text(0, 1.65, 'Transparent exterior return: D = 0.1', ha='center', fontsize=9, color='#7951b0',
                 bbox=dict(facecolor='white', edgecolor='none', pad=2))
    axes[0].text(0, -.7, 'Assumed throat', ha='center', fontsize=10, color='#326bb3')
    axes[0].text(0, -1.55, 'Schematic cross-section; not a construction plan', ha='center', fontsize=9, color='#475569')
    axes[0].set(xlim=(-1.4, 1.4), ylim=(-1.8, 1.9)); axes[0].axis('off')
    xx = np.asarray(response['x']); after = xx >= .4
    axes[1].plot(xx[after], np.asarray(response['j_plus_over_alpha'])[after],
                 color='#9f3345', lw=2, label='Forward-moving stress change')
    axes[1].plot(xx[after], np.asarray(response['j_minus_over_alpha'])[after],
                 color='#007e87', lw=2, label='Opposite-moving stress change')
    axes[1].axhline(0, color='#64748b', lw=1, ls='--', label='Frozen-source prediction after pulse')
    axes[1].set(title='Response continues after injection ends', xlabel='Dimensionless throat null time x',
                ylabel='Null stress perturbation j / α', xlim=(.4, np.pi), ylim=(-.009, .013))
    axes[1].spines[['top', 'right']].set_visible(False); axes[1].grid(alpha=.16)
    axes[1].legend(frameon=False, fontsize=8, loc='upper right')
    fig.text(.06, .20, 'The field response, mouth clocks and throat-area change now satisfy the same reduced equations. Independent calculations agree.',
             fontsize=9, color='#475569')
    fig.text(.06, .14, 'Assumptions: fixed exterior geometry, one common delay and linear quantum mean-field response. No full-return stability test yet.',
             fontsize=9, color='#475569')
    fig.text(.06, .067, '13 September 2027 → 12 September 2027: still unproved. This improves the model, not evidence of an operating time machine.',
             fontsize=10, fontweight='bold', color='#9f3345')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'mmp_feedback_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root/'proton_focusing_progress.png')
    print(root/'mmp_feedback_progress.png')


def return_loss_progress():
    root = Path(__file__).resolve().parent
    loops = json.loads((root/'mmp_pulse_results.json').read_text())['loop_return']['cases']
    loss = json.loads((root/'future_time_results.json').read_text())['proton_slow_loss']
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.3))
    fig.subplots_adjust(left=.085, right=.97, top=.74, bottom=.33, wspace=.31)
    fig.suptitle('Next checks: complete returns and radiation losses', fontsize=19, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Our reduced calculations • finite intervals and stated approximations',
             ha='center', fontsize=11, color='#475569')
    for mode, label, color in (('one_pass', 'One-pass carrier', '#326bb3'),
                               ('circulating', 'Once-injected circulating carrier', '#7951b0')):
        data = loops[f'0.1_{mode}']['plot']
        axes[0].plot(data['loop_number'], 100*(np.asarray(data['right_clock_rate'])-1),
                     color=color, lw=1.5, label=label)
    axes[0].set(title='Wormhole model: finite exterior delay D = 0.1',
                xlabel='Elapsed time / background loop time', ylabel='Right mouth-clock rate change (%)',
                xlim=(0, 8), ylim=(-1.85, 2.5))
    axes[0].legend(frameon=False, fontsize=8, loc='upper left')
    data = loss['plot']; hours = np.asarray(data['lab_hours'])
    center = np.asarray(data['radius_shift_mm'])
    amplitude = np.asarray(data['radial_envelope_mm'])-np.abs(center)
    axes[1].fill_between(hours, center-amplitude, center+amplitude,
                         color='#007e87', alpha=.2, label='Estimated radial oscillation envelope')
    axes[1].plot(hours, center, color='#007e87', lw=2, label='Slowly drifting circular orbit')
    axes[1].axhline(0, color='#64748b', lw=1, ls='--')
    axes[1].set(title='Proton: one-day adiabatic radiation estimate',
                xlabel='Lab time (hours)', ylabel='Radial offset from original orbit (mm)',
                xlim=(0, 24), ylim=(-3, 2))
    axes[1].legend(frameon=False, fontsize=8, loc='lower left')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.06, .21, 'Left: all sampled clock/ray rates and return-time gaps remain positive. This is not a proof of nonlinear stability or past travel.',
             fontsize=9, color='#475569')
    fig.text(.06, .151, 'Right: orbit centre shifts inward 1.40 mm; estimated maximum radial excursion is 2.48 mm. Chosen tube limits are ±10 mm.',
             fontsize=9, color='#475569')
    fig.text(.06, .091, 'Radiation changes the conserved quantities used in the earlier confinement proof. These envelopes require adiabatic and linear assumptions.',
             fontsize=9, color='#475569')
    fig.text(.06, .035, '13 → 12 Sep 2027: unproved.     Forward example: 24 lab hours ≈ 2 h 24 min 0.116 s of proton proper time in this loss model.',
             fontsize=10, color='#007e87', fontweight='bold')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'return_loss_progress.{suffix}', dpi=170, facecolor='white')
    plt.close(fig)
    print(root/'return_loss_progress.png')


def offset_rf_progress():
    root = Path(__file__).resolve().parent
    offset = json.loads((root/'mmp_pulse_results.json').read_text())['time_offset']
    rf = json.loads((root/'future_time_results.json').read_text())['proton_rf']
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 6.5))
    fig.subplots_adjust(left=.045, right=.98, top=.74, bottom=.32, wspace=.42)
    fig.suptitle('A time offset changes the support; RF restores proton energy',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Our mathematical checks • assumed wormhole geometry and an ideal accelerator model',
             ha='center', fontsize=11, color='#475569')
    x = np.linspace(-1, 1, 150); width = .2+.45*x*x
    axes[0].fill_between(x, -width, width, color='#326bb3', alpha=.10)
    axes[0].plot(x, width, x, -width, color='#326bb3', lw=2)
    for mouth, label in ((-1, 'Mouth A'), (1, 'Mouth B')):
        axes[0].plot([mouth, mouth], [-.65, .65], color='#007e87', lw=4)
        axes[0].text(mouth, -.92, label, ha='center', fontsize=10)
    for start, end, height, line_width, color in ((-.8, .8, .09, 1.5, '#7951b0'),
                                                (.8, -.8, -.1, 4, '#9f3345')):
        axes[0].annotate('', xy=(end, height), xytext=(start, height),
                         arrowprops=dict(arrowstyle='->', lw=line_width, color=color))
    axes[0].annotate('', xy=(1, .7), xytext=(-1, .7),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-.45',
                                     color='#64748b', ls='--'))
    axes[0].text(0, 1.5, 'Exterior return', ha='center', fontsize=10, color='#64748b')
    axes[0].text(0, -1.35, 'Unequal quantum stresses\nforce the mouth data to change',
                 ha='center', fontsize=10, color='#9f3345')
    axes[0].set(title='Assumed passage with a physical offset', xlim=(-1.5, 1.5), ylim=(-1.8, 1.9))
    axes[0].axis('off')
    a = np.linspace(0, .99, 300)
    axes[1].semilogy(a, 1/(1-a)**2, color='#9f3345', lw=2, label='Short-direction stress magnitude')
    axes[1].semilogy(a, 1/(1-a*a), color='#7951b0', lw=2, label='Rest-frame energy magnitude')
    for row in offset['samples']:
        aa = row['offset_over_loop_length']
        assert np.isclose(row['short_chiral_vacuum_magnitude_ratio'], 1/(1-aa)**2)
    axes[1].set(title='Prescribed vacuum at fixed loop length', xlabel='Time offset / loop length',
                ylabel='Magnitude / value at zero offset', xlim=(0, 1))
    axes[1].legend(frameon=False, fontsize=8, loc='upper left')
    data = rf['plot']
    axes[2].plot(data['phase_error_degrees'], data['energy_deviation_eV'], color='#007e87', lw=2)
    axes[2].scatter(data['phase_error_degrees'][0], data['energy_deviation_eV'][0],
                    color='#326bb3', s=35, label='Starting disturbance', zorder=3)
    axes[2].axhline(0, color='#94a3b8', lw=.7); axes[2].axvline(0, color='#94a3b8', lw=.7)
    axes[2].set(title='Proton: restoring phase–energy motion', xlabel='Arrival-phase error (degrees)',
                ylabel='Energy deviation (eV)')
    axes[2].legend(frameon=False, fontsize=8, loc='upper right')
    for ax in axes[1:]:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.045, .215, 'Past: this shifted vacuum cannot support the unchanged static throat. The curves are a family of states, not an evolution or a material recipe.',
             fontsize=9, color='#475569')
    fig.text(.045, .15, 'Future: a 1 mV effective RF gap balances mean loss and gives a restoring cycle of about 3.9 seconds. Full cavity and transverse motion remain open.',
             fontsize=9, color='#475569')
    fig.text(.045, .085, 'Next: couple physical mouth driving to its quantum source; connect proton RF to injection and extraction. Schematic is not to scale.',
             fontsize=9, color='#475569')
    fig.text(.045, .032, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is supported; the complete proposed apparatus is not demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'offset_rf_progress.{suffix}', dpi=180, facecolor='white')
    plt.close(fig)
    print(root/'offset_rf_progress.png')


def balance_kicker_progress():
    root = Path(__file__).resolve().parent
    balance = json.loads((root/'mmp_pulse_results.json').read_text())['stationary_balance']
    kicker = json.loads((root/'future_time_results.json').read_text())['proton_kicker']
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 6.5))
    fig.subplots_adjust(left=.06, right=.975, top=.74, bottom=.32, wspace=.39)
    fig.suptitle('Balanced throat stresses still give forward returns; a finite kick starts an exit',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Checked mathematical models • stationary quantum source and a short proton trajectory',
             ha='center', fontsize=11, color='#475569')
    length = balance['fixed_length_parameters']['L']
    critical = balance['fixed_length_parameters']['support_threshold_delta']
    delta = np.linspace(0, length, 401)
    opening = ((2*np.pi/(length+delta))**2-1)/3
    axes[0].plot(delta, opening, lw=2.5, color='#007e87')
    axes[0].fill_between(delta, 0, np.maximum(opening, 0), color='#007e87', alpha=.1)
    axes[0].axhline(0, color='#64748b', lw=1)
    axes[0].axvline(critical, color='#9f3345', ls='--', lw=1.5)
    axes[0].text(.86, .9, 'Fixed throat clock scale\nL = π + 0.1', transform=axes[0].transAxes,
                 ha='right', va='top', fontsize=9, color='#475569')
    axes[0].set(title='Best balanced stationary support', xlabel='Dimensionless mouth time offset |Δ|',
                ylabel='Maximum opening coefficient β / α', xlim=(0, length), ylim=(-.075, 1.02))
    rows = balance['matched_common_delay_samples'][:5]
    axes[1].plot([r['physical_offset_over_ell_v'] for r in rows],
                 [r['shorter_return_over_ell_v'] for r in rows], 'o-', color='#7951b0', lw=2,
                 label='Matched common-delay samples')
    axes[1].axhline(0, color='#9f3345', lw=1.5, ls='--', label='Zero time: closed-loop threshold')
    axes[1].set_xscale('symlog', linthresh=.1)
    axes[1].set(title='Allowing the physical throat to lengthen', xlabel='Physical offset / reference length ℓᵥ',
                ylabel='Complete return time / ℓᵥ  (c = 1)', xlim=(-.01, 150), ylim=(-.5, 11.5))
    axes[1].legend(frameon=False, fontsize=8, loc='upper left')
    data = kicker['plot']; t = 1000*np.asarray(data['lab_microseconds'])
    radial = np.asarray(data['radial_mm'])
    crossing = kicker['cases']['nominal_strength']['departure']['lab_seconds']*1e9
    before, after = t<=crossing, t>crossing
    axes[2].plot(t[before], radial[before], color='#326bb3', lw=2.5, label='Finite-field proton trajectory')
    axes[2].plot(t[after], radial[after], color='#326bb3', lw=1.5, ls='--')
    axes[2].plot(t, data['no_kick_radial_mm'], color='#64748b', lw=1.2, label='No-kick control')
    axes[2].axhline(10, color='#b7791f', ls='--', lw=1.5, label='Chosen radial boundary: 10 mm')
    axes[2].scatter([crossing], [10], color='#b7791f', s=40, zorder=4)
    axes[2].annotate(f'{crossing:.1f} ns', xy=(crossing, 10), xytext=(45, 11.7),
                     arrowprops=dict(arrowstyle='->', color='#b7791f'), fontsize=9, color='#b7791f')
    axes[2].set(title='Finite magnetic kick: outward departure', xlabel='Lab time from entrance (nanoseconds)',
                ylabel='Radial displacement (mm)', xlim=(0, t[-1]), ylim=(-.6, 13.5))
    axes[2].legend(frameon=False, fontsize=7.5, loc='lower right', bbox_to_anchor=(1., .08))
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.055, .215, 'Past: neutral fermion excitations can balance the mean stress, but weaken its support. Changing the matched throat length keeps the return positive.',
             fontsize=9, color='#475569')
    fig.text(.055, .15, 'Future: the finite fringe field sends a nominal proton across 1 cm. Dashed continuation still uses the ideal ring field; a real outgoing region is missing.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Next: earliest-arrival and total-stress checks for a driven throat; a septum and matched outgoing field for the proton. No construction is demonstrated.',
             fontsize=9, color='#475569')
    fig.text(.055, .032, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: supported time dilation; complete entry, storage and exit remain to be joined.',
             fontsize=10, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'balance_kicker_progress.{suffix}', dpi=180, facecolor='white')
    plt.close(fig)
    print(root/'balance_kicker_progress.png')


def driven_outgoing_progress():
    from matplotlib.patches import Circle, FancyArrowPatch
    root = Path(__file__).resolve().parent
    drive = json.loads((root/'mmp_pulse_results.json').read_text())['driven_clock']
    outgoing = json.loads((root/'future_time_results.json').read_text())['proton_outgoing_candidate']
    fig, axes = plt.subplots(1, 3, figsize=(14, 7))
    fig.subplots_adjust(left=.065, right=.975, top=.74, bottom=.32, wspace=.36)
    fig.suptitle('A driven throat needs extra forcing; the proton has a candidate outgoing field',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Verified equations and ideal trajectories — neither drawing is a construction blueprint',
             ha='center', fontsize=11, color='#475569')
    ax = axes[0]
    ax.set(xlim=(0, 10), ylim=(0, 8), title='The proposed passage being tested')
    ax.axis('off')
    for x, label in ((1.5, 'Left mouth'), (8.5, 'Right mouth')):
        ax.add_patch(Circle((x, 4), .62, edgecolor='#007e87', facecolor='#e0f2f1', lw=2))
        ax.text(x, 2.8, label, ha='center', fontsize=10)
    ax.add_patch(FancyArrowPatch((2.12, 4.2), (7.88, 4.2), connectionstyle='arc3,rad=-.6',
                                arrowstyle='->', mutation_scale=14, lw=1.8, color='#64748b'))
    ax.text(5, 6.4, 'Exterior route: D = 0.1', ha='center', fontsize=9, color='#475569')
    ax.add_patch(FancyArrowPatch((2.12, 3.95), (7.88, 3.95), connectionstyle='arc3,rad=.28',
                                arrowstyle='->', mutation_scale=16, lw=9, alpha=.2, color='#7951b0'))
    ax.add_patch(FancyArrowPatch((2.12, 3.95), (7.88, 3.95), connectionstyle='arc3,rad=.28',
                                arrowstyle='->', mutation_scale=13, lw=1.8, color='#7951b0'))
    ax.text(5, 3.7, 'Assumed throat', ha='center', fontsize=10, color='#7951b0')
    ax.text(5, 1.75, 'Right clock ramps to 1.05 × left rate\nQuantum stress follows repeated returns',
            ha='center', fontsize=9, linespacing=1.6, color='#475569')
    ax.text(5, .45, 'Creation and physical driving remain unsolved',
            ha='center', fontsize=9, fontweight='bold', color='#9f3345')

    ax = axes[1]; data = drive['plot']; thresholds = drive['thresholds']
    t = np.asarray(data['entry_times'])
    ax.plot(t, data['one_way_times'], color='#7951b0', lw=2, label='One throat passage')
    ax.plot(t, data['complete_return_times'], color='#007e87', lw=2, label='Throat + outside return')
    ax.axhline(.1, color='#64748b', lw=1, ls=':', label='Outside one-way time')
    ax.axhline(0, color='#9f3345', lw=1)
    for key, color in (('left_entry_time_first_shortcut', '#b7791f'),
                       ('left_entry_time_first_closed_null_loop', '#9f3345')):
        ax.axvline(thresholds[key], color=color, lw=1, ls='--')
    ax.text(thresholds['left_entry_time_first_shortcut']-.15, .45, 'Shortcut',
            rotation=90, va='top', ha='right', color='#b7791f', fontsize=8)
    ax.text(thresholds['left_entry_time_first_closed_null_loop']-.15, .45, 'Closed-loop threshold',
            rotation=90, va='top', ha='right', color='#9f3345', fontsize=8)
    ax.set(title='Two different timing thresholds', xlabel='Left entry time (dimensionless)',
           ylabel='Elapsed exterior time (dimensionless)', xlim=(56, 66), ylim=(-.15, .5))
    ax.legend(frameon=False, fontsize=8, loc='lower left')

    ax = axes[2]; data = outgoing['plot']
    ax.plot(data['s_metres'], data['candidate_radius_above_nominal_mm'], color='#326bb3',
            lw=2.5, label='Candidate outgoing field')
    ax.plot(data['s_metres'], data['control_radius_above_nominal_mm'], color='#64748b',
            lw=1.5, ls='--', label='Original-ring control')
    ax.set(title='Proton in a smooth outgoing field', xlabel='Distance along local tangent (m)',
           ylabel='Distance outside nominal ring radius (mm)', xlim=(0, 5), ylim=(0, 315))
    ax.text(.04, .74, 'Same initial proton state\nDifferent field models', transform=ax.transAxes,
            fontsize=9, color='#475569', linespacing=1.4)
    ax.text(.04, .50, '5 m checkpoint: 16.77 ns lab\n1.677 ns proton clock', transform=ax.transAxes,
            fontsize=9, color='#326bb3', linespacing=1.4)
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    for ax in axes[1:]:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.055, .215, 'Past: the imposed drive fails its gravity balance. At the assumed shortcut, one chiral stress is about 200 × its initial magnitude.',
             fontsize=9, color='#475569')
    fig.text(.055, .15, 'Future: Maxwell equations and the local outgoing trajectory pass. Matching one entry point does not join this field to the earlier storage ring.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Next: specify a physical boundary source and solve its feedback; design one common storage-to-exit magnetic field. Diagram is schematic.',
             fontsize=9, color='#475569')
    fig.text(.055, .032, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is supported; the complete proposed apparatus remains unbuilt.',
             fontsize=10, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'driven_outgoing_progress.{suffix}', dpi=180, facecolor='white')
    plt.close(fig)
    print(root/'driven_outgoing_progress.png')


def coupled_coil_progress():
    from matplotlib.patches import Ellipse
    root = Path(__file__).resolve().parent
    pulse = json.loads((root/'mmp_pulse_results.json').read_text())['nonlinear_pulse']
    coil = json.loads((root/'future_time_results.json').read_text())['proton_finite_coil']
    fig, axes = plt.subplots(1, 3, figsize=(14, 7))
    fig.subplots_adjust(left=.065, right=.975, top=.74, bottom=.32, wspace=.34)
    fig.suptitle('The throat response now balances; one coil field confines the model proton',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Two verified model advances — an assumed throat with dynamical clocks, and ideal filament-coil storage',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]; data = pulse['cases']['0.1']['plot']
    loops = np.asarray(data['exterior_time'])/(np.pi+pulse['parameters']['exterior_delay'])
    ax.plot(loops, data['right_clock_rate'], color='#7951b0', lw=2, label='Right mouth')
    ax.plot(loops, data['left_clock_rate'], color='#007e87', lw=2, label='Left mouth')
    ax.axhline(1, color='#64748b', lw=1, ls=':')
    ax.axvspan(0, pulse['parameters']['pulse_width']/(np.pi+.1), color='#b7791f', alpha=.2)
    ax.set(title='Clocks respond to the specified pulse', xlabel='Elapsed background loop times',
           ylabel='Throat clock rate / exterior clock rate', xlim=(0, 4), ylim=(.78, 1.2))
    ax.legend(frameon=False, loc='upper left', fontsize=9)
    ax.text(.97, .92, 'Stronger input: ε = 0.1', transform=ax.transAxes,
            ha='right', fontsize=9, color='#475569')
    ax.text(.04, .06, 'Shortest sampled return: 3.064 > 0\nDimensionless time; no past return',
            transform=ax.transAxes, fontsize=9, color='#475569', linespacing=1.4)

    ax = axes[1]
    ax.set(xlim=(-65, 65), ylim=(-58, 58), title='Two source coils (schematic)')
    ax.axis('off')
    for z, label in ((32, 'z = +40 m'), (-32, 'z = −40 m')):
        ax.add_patch(Ellipse((0, z), 100, 14, edgecolor='#b7791f', facecolor='#fef3c7', lw=2))
        ax.text(0, z, label, ha='center', va='center', fontsize=9, color='#92400e')
    ax.add_patch(Ellipse((0, 0), 62.28, 10, edgecolor='#326bb3', facecolor='none', lw=2.3))
    ax.scatter([31.14], [0], s=30, color='#326bb3', zorder=3)
    ax.text(0, 13, 'Proton orbit: R = 31.14 m', ha='center', fontsize=9, color='#326bb3')
    ax.text(0, -15, '1 T at the orbit; γ = 10', ha='center', fontsize=9, color='#326bb3')
    ax.text(0, 49, 'Each coil radius: 50 m', ha='center', fontsize=9, color='#475569')
    ax.text(0, -50, 'About 1.06 × 10⁸ ampere-turns per coil\nIdeal thin currents; hardware is unresolved',
            ha='center', fontsize=8.5, color='#92400e', linespacing=1.4)

    ax = axes[2]; data = coil['plot']
    ax.plot(data['reference_turn'], data['radial_mm'], color='#326bb3', lw=1.4, label='Radial displacement')
    ax.plot(data['reference_turn'], data['vertical_mm'], color='#007e87', lw=1.3, label='Vertical displacement')
    ax.set(title='A displaced proton stays near its orbit', xlabel='Reference revolutions',
           ylabel='Displacement from nominal orbit (mm)', xlim=(0, 20), ylim=(-1.6, 2.2))
    ax.legend(frameon=False, loc='upper right', fontsize=8)
    ax.text(.02, .04, 'Conserved-invariant bound: inside ±10 mm\nApplies to this ideal conservative state',
            transform=ax.transAxes, fontsize=8.5, color='#475569', linespacing=1.4)
    for ax in (axes[0], axes[2]):
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.055, .215, 'Past: gravity matching, returned quantum stress and the energy ledger agree. All tested pulses still give forward returns.',
             fontsize=9, color='#475569')
    fig.text(.055, .15, 'Future: one field from explicit coil locations supports the orbit and its confinement bound. Leaving it requires a separately modeled control.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Next: stress-test the coupled throat within its validity limits; add extraction control to this same coil field. No physical throat has been created.',
             fontsize=9, color='#475569')
    fig.text(.055, .032, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: supported time dilation; one-day storage and retrieval are not demonstrated here.',
             fontsize=10, fontweight='bold', color='#007e87')
    for suffix in ('png', 'pdf'):
        fig.savefig(root/f'coupled_coil_progress.{suffix}', dpi=180, facecolor='white')
    plt.close(fig)
    print(root/'coupled_coil_progress.png')


def interior_escape_progress():
    """Plot only the final interior and same-coil escape records."""
    root = Path(__file__).resolve().parent
    mmp = json.loads((root/'mmp_pulse_results.json').read_text())
    cases = mmp['nonlinear_interior']['cases']
    required = {'0.0_4', '0.1_16', '0.3_4', '0.5_4'}
    if not required.issubset(cases):
        raise ValueError('Final nonlinear_interior cases are not all available; no figure written')
    escape = json.loads((root/'future_time_results.json').read_text())['proton_finite_coil_escape']
    loop_time = np.pi+mmp['nonlinear_pulse']['parameters']['exterior_delay']
    minimum_return = min(cases[key]['minimum_sampled_complete_return_time'] for key in required)
    assert minimum_return > 0
    for key in required:
        case = cases[key]; interior = case['interior']
        values = np.asarray([r['minimum_phi_over_alpha'] for r in interior['slice_summary']])
        assert np.all(np.isfinite(values)) and interior['alpha'] > 0
        assert np.isclose(1+interior['alpha']*values.min(),
                          interior['minimum_sampled_full_area_over_extremal_area'],
                          rtol=0, atol=2e-15)
        assert interior['minimum_sampled_full_area_over_extremal_area'] > 0
        assert case['minimum_sampled_clock_rate'] > 0
    gamma = escape['source_parameters']['reference_gamma']
    for case in escape['cases'].values():
        assert np.isclose(case['elapsed_lab_seconds']/case['elapsed_proton_proper_seconds'],
                          gamma, rtol=1e-10)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.08, right=.975, top=.835, bottom=.20, hspace=.60, wspace=.27)
    fig.suptitle('Throat interior tested; a proton leaves its storage orbit',
                 fontsize=19, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Two ideal model calculations • finite throat samples and conditional proton initial data',
             ha='center', fontsize=10.5, color='#475569')
    colors = {'0.1_16': '#007e87', '0.5_4': '#7951b0'}
    labels = {'0.1_16': 'ε = 0.1; 16 loop times', '0.5_4': 'ε = 0.5; 4 loop times'}

    ax = axes[0, 0]
    x = np.linspace(-1, 1, 180); radius = .26+.42*x*x
    ax.fill_between(x, -radius, radius, color='#326bb3', alpha=.10)
    ax.plot(x, radius, x, -radius, color='#326bb3', lw=2)
    for mouth, label in ((-1, 'Left mouth'), (1, 'Right mouth')):
        ax.plot([mouth, mouth], [-.68, .68], color='#007e87', lw=3)
        ax.text(mouth, -.90, label, ha='center', fontsize=10)
    ax.annotate('', xy=(.97, .76), xytext=(-.97, .76),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-.20',
                                color='#64748b', ls='--', lw=1.3))
    ax.text(0, 1.36, 'Outside return included', ha='center', color='#475569', fontsize=9)
    ax.text(0, 0, 'Assumed throat', ha='center', va='center', color='#326bb3', fontsize=11)
    ax.text(0, -1.22, f'Smallest sampled complete return: {minimum_return:.3f} > 0',
            ha='center', fontsize=9.5, color='#9f3345')
    ax.text(0, -1.51, 'Model time units • no past return in these samples',
            ha='center', fontsize=9, color='#475569')
    ax.set(title='A. Assumed geometry — schematic, not a blueprint',
           xlim=(-1.65, 1.65), ylim=(-1.67, 1.55))
    ax.axis('off')

    ax = axes[0, 1]; maximum_excess = 0.
    for key in ('0.0_4', '0.1_16', '0.5_4'):
        interior = cases[key]['interior']; rows = interior['slice_summary']
        times = np.asarray([r['exterior_time'] for r in rows])/loop_time
        excess = np.asarray([r['minimum_phi_over_alpha'] for r in rows])
        maximum_excess = max(maximum_excess, float(excess.max()))
        assert np.all(np.diff(times) > 0)
        ax.plot(times, excess, lw=1.8, color=colors.get(key, '#64748b'),
                ls='--' if key == '0.0_4' else '-',
                label=labels.get(key, 'Unperturbed; 4 loop times'))
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set(title='B. Minimum area excess in the sampled interior',
           xlabel='Elapsed background loop times', ylabel='Minimum area excess / α',
           xlim=(0, 16), ylim=(-.035, 1.20*maximum_excess))
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    alpha = cases['0.1_16']['interior']['alpha']
    strip = cases['0.1_16']['interior']['chi_interval']
    assert all(cases[key]['interior']['alpha'] == alpha and
               cases[key]['interior']['chi_interval'] == strip for key in required)
    ax.text(.98, .97, 'Area / reference area = 1 + α × plotted value\n'
            f'α = {alpha:.0e}; sampled χ ∈ [{strip[0]:g}, {strip[1]:g}]',
            transform=ax.transAxes, ha='right', va='top', fontsize=8.4, color='#475569')

    ax = axes[1, 0]; maximum_rate = 0.
    for key in colors:
        data = cases[key]['plot']
        times = np.asarray(data['exterior_time'])/loop_time
        rate = np.minimum(data['right_clock_rate'], data['left_clock_rate'])
        maximum_rate = max(maximum_rate, float(rate.max()))
        assert len(times) == len(rate) and np.all(np.isfinite(rate)) and np.all(rate > 0)
        ax.plot(times, rate, color=colors[key], lw=1.8, label=labels[key])
    ax.axhline(1, color='#64748b', ls=':', lw=1)
    ax.set(title='C. Both mouth clocks keep advancing',
           xlabel='Elapsed background loop times', ylabel='Lower mouth rate / exterior rate',
           xlim=(0, 16), ylim=(0, 1.12*maximum_rate))
    ax.legend(frameon=False, fontsize=8, loc='lower right')

    ax = axes[1, 1]
    for key, color, label in (('below_threshold_5_degrees', '#7951b0', '5°: outward, then returns'),
                               ('above_threshold_10_degrees', '#326bb3', '10°: reaches 1 km')):
        data = escape['plot'][key]
        t, radius = np.asarray(data['lab_microseconds']), np.asarray(data['radius_metres'])
        assert len(t) == len(radius) and np.all(np.isfinite(radius)) and np.all(radius > 0)
        assert np.isclose(radius[-1], escape['cases'][key]['endpoint_radius_metres'])
        ax.plot(t, radius, color=color, lw=2.2, label=label)
        ax.scatter([t[-1]], [radius[-1]], color=color, s=23, zorder=3)
    ax.set_yscale('log')
    ax.set_yticks([30, 100, 300, 1000], labels=['30', '100', '300', '1000'])
    ax.minorticks_off()
    ax.set(title='D. Same coils; two prescribed initial directions',
           xlabel='Elapsed lab time (microseconds)', ylabel='Distance from coil axis (m; log scale)',
           xlim=(0, 4.05), ylim=(29, 1350))
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    outgoing = escape['cases']['above_threshold_10_degrees']
    ax.text(.96, .13, f'At 1 km: lab {outgoing["elapsed_lab_seconds"]*1e6:.3f} μs\n'
            f'Proton proper time {outgoing["elapsed_proton_proper_seconds"]*1e6:.3f} μs (γ = {gamma:g})',
            transform=ax.transAxes, ha='right', fontsize=8.8, color='#326bb3', linespacing=1.5)
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16)
    fig.text(.06, .13, 'Throat: leading reconstruction on a finite strip and time interval. Positive sampled area does not prove stability or safe passage.',
             fontsize=9, color='#475569')
    fig.text(.06, .085, 'Proton: the changed direction is prescribed. A physical kick, capture, readout and one-day storage remain unmodeled.',
             fontsize=9, color='#475569')
    fig.text(.06, .033, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is supported; this complete journey is not demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'interior_escape_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def angular_control_progress():
    """Compare one static angular correction with the final retarded-control run."""
    from matplotlib.patches import Ellipse
    root = Path(__file__).resolve().parent
    angular = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['angular_metric']
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'proton_retarded_control' not in future:
        raise ValueError('Final retarded-control results are unavailable; no figure written')
    control = future['proton_retarded_control']
    parameters = control['parameters']
    pulse, route = control['plot']['pulse'], control['plot']['route']
    clear_ns = control['pulse_clear']['lab_seconds']*1e9
    width_ns = parameters['source_pulse_duration_seconds']*1e9
    ns, current = np.asarray(pulse['lab_nanoseconds']), np.asarray(pulse['source_current_MA'])
    gamma = np.asarray(pulse['gamma'])
    us, radius = np.asarray(route['lab_microseconds']), np.asarray(route['radius_metres'])
    proper = np.asarray(route['proton_proper_microseconds'])
    assert len(ns) == len(current) == len(gamma) and np.all(np.diff(ns) > 0)
    assert len(us) == len(radius) == len(proper) and np.all(np.diff(us) > 0)
    assert np.all(np.isfinite(current)) and np.all(np.isfinite(gamma)) and np.all(gamma >= 1)
    assert np.all(np.isfinite(radius)) and np.all(radius > 0) and us[-1] > proper[-1] > 0
    assert 0 < width_ns < clear_ns < ns[-1]
    assert np.isclose(radius[-1], control['endpoint']['radius_metres'])
    assert np.isclose(us[-1]*1e-6, control['endpoint']['lab_seconds'], rtol=1e-10, atol=1e-15)
    rows = angular['compact_points']
    before = np.abs([row['uncorrected_angular_residual'] for row in rows])
    after = np.abs([row['corrected_angular_residual'] for row in rows])
    assert angular['symbolic_verification']['zero_and_first_order_angular_residuals_vanish']
    assert np.all(np.isfinite(before)) and np.all(np.isfinite(after))
    assert np.all((after > 0) & (after < before))
    source = parameters['storage_source']
    orbit = source['reference_orbit_radius_metres']

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.08, right=.92, top=.835, bottom=.20, hspace=.61, wspace=.38)
    fig.suptitle('Gravity correction checked; a finite pulse steers the proton',
                 fontsize=19, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Ideal model calculations • continuous field-driven motion • no constructed apparatus',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    for z in source['coil_z_metres']:
        ax.add_patch(Ellipse((0, z), 2*source['coil_radius_metres'], 10,
                            edgecolor='#b7791f', facecolor='#fef3c7', lw=1.8))
    ax.add_patch(Ellipse((0, 0), 2*orbit, 5, edgecolor='#326bb3', facecolor='none', lw=1.8))
    for x, _, z in parameters['control_loop_centers_metres']:
        ax.add_patch(Ellipse((x, z), 2*parameters['control_loop_radius_metres'], 3,
                            edgecolor='#7951b0', facecolor='#ede9fe', lw=2))
    phase = parameters['starting_azimuth_radians']
    ax.scatter([orbit*np.cos(phase)], [2.5*np.sin(phase)], s=23, color='#326bb3', zorder=4)
    ax.text(0, 52, 'Original storage coils stay on', ha='center', fontsize=9.5, color='#92400e')
    ax.text(-13, -13, 'Initial proton orbit', ha='center', fontsize=9, color='#326bb3')
    ax.annotate('Two localized\npulsed loops', xy=(orbit+3, 6), xytext=(44, 20),
                ha='center', fontsize=9, color='#7951b0',
                arrowprops=dict(arrowstyle='->', color='#7951b0', lw=1))
    ax.text(0, -53, 'Field changes propagate at c\nPrescribed currents; no hardware model',
            ha='center', va='top', fontsize=9, color='#475569', linespacing=1.5)
    ax.set(title='A. Coil locations — schematic, not to scale', xlim=(-63, 66), ylim=(-78, 62))
    ax.axis('off')

    ax = axes[0, 1]; energy_axis = ax.twinx()
    current_line, = ax.plot(ns, current, color='#7951b0', lw=2, label='Current in each control loop')
    energy_line, = energy_axis.plot(ns, gamma, color='#326bb3', lw=2, label='Proton energy / rest energy')
    off_line = ax.axvline(width_ns, color='#64748b', ls=':', lw=1.2,
                         label=f'Source current ends: {width_ns:g} ns')
    clear_line = ax.axvline(clear_ns, color='#b7791f', ls='--', lw=1.2,
                           label=f'Control field clears proton: {clear_ns:.1f} ns')
    current_span = max(float(np.ptp(current)), 1.)
    energy_pad = .12*max(float(np.ptp(gamma)), .001)
    ax.set(title='B. The retarded field changes proton energy',
           xlabel='Lab time (nanoseconds)', ylabel='Impressed current per loop (MA)',
           xlim=(ns[0], ns[-1]), ylim=(current.min()-.10*current_span, current.max()+.15*current_span))
    energy_axis.set(ylabel='Proton energy / rest energy, γ',
                    ylim=(gamma.min()-energy_pad, gamma.max()+energy_pad))
    ax.yaxis.label.set_color('#7951b0'); energy_axis.yaxis.label.set_color('#326bb3')
    ax.legend(handles=[current_line, energy_line, off_line, clear_line],
              framealpha=.95, edgecolor='none', fontsize=7.4, loc='lower right')
    energy_axis.spines[['top', 'left']].set_visible(False)

    ax = axes[1, 0]
    ax.plot(us, radius, color='#326bb3', lw=2.3, label='Finite current pulse')
    ax.plot(us, route['zero_current_radius_metres'], color='#64748b', ls=':', lw=1.4,
            label='Zero-current control')
    ax.scatter([us[-1]], [radius[-1]], color='#326bb3', s=27, zorder=4)
    ax.axvline(clear_ns/1000, color='#b7791f', ls='--', lw=1.1)
    ax.set_yscale('log')
    ax.set_yticks([30, 100, 1000, 5000], labels=['30', '100', '1000', '5000'])
    ax.minorticks_off()
    ax.set(title='C. One continuous Lorentz-force trajectory',
           xlabel='Elapsed lab time (microseconds)', ylabel='Radius from storage axis (m; log scale)',
           xlim=(0, us[-1]*1.05), ylim=(radius.min()*.9, radius.max()*1.35))
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax.text(.96, .12, f'At {radius[-1]/1000:g} km: lab {us[-1]:.3f} μs\n'
            f'Proton proper time {proper[-1]:.3f} μs\nNo prescribed momentum turn',
            transform=ax.transAxes, ha='right', fontsize=8.8, color='#326bb3', linespacing=1.5)

    ax = axes[1, 1]; x = np.arange(len(rows))
    ax.bar(x-.17, before, width=.32, color='#b7791f', label='Before metric correction')
    ax.bar(x+.17, after, width=.32, color='#007e87', label='After first-order correction')
    ax.set_yscale('log')
    ax.set_xticks(x, labels=[f'ρ = {row["rho"]:g}' for row in rows])
    ax.set(title='D. One static angular gravity equation',
           xlabel='Fixed compact sample positions', ylabel='Absolute residual (dimensionless)',
           ylim=(after.min()/5, before.max()*100))
    ax.legend(framealpha=.95, edgecolor='none', fontsize=8, loc='center left')
    ax.text(.98, .96, 'First-order benchmark only\nNot the complete gravity equations',
            transform=ax.transAxes, ha='right', va='top', fontsize=8.2, color='#475569')
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16, axis='y')
        ax.set_axisbelow(True)
    fig.text(.06, .13, 'Gravity: this static correction removes one leading residual. The full quantum source, dynamics and exterior matching remain open.',
             fontsize=9, color='#475569')
    fig.text(.06, .085, 'Proton: retarded fields come from ideal prescribed currents. No finite-conductor circuit, capture, reunion or one-day storage is modeled.',
             fontsize=9, color='#475569')
    fig.text(.06, .033, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is supported; this complete journey is not demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'angular_control_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def polyakov_source_progress():
    """Render the final static Polyakov solution and coherent source radiation."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['polyakov_static']
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'source_feasibility' not in future:
        raise ValueError('Final source-feasibility results are unavailable; no figure written')
    source = future['source_feasibility']; radiation = source['radiation']
    case = gravity['cases'][0]; alpha = case['alpha']
    solution, leading = case['plot'], gravity['leading_comparison']
    x = np.asarray(solution['x'])
    y = np.asarray(solution['area_excess_over_alpha'])
    first = np.asarray(leading['area_excess_over_alpha'])
    assert alpha == gravity['alpha'] and np.array_equal(x, leading['x'])
    assert len(x) == len(y) == len(first) and np.all(np.diff(x) > 0)
    assert np.all(np.isfinite(y)) and np.all(1+alpha*y > 0)
    assert np.isclose(first[0], gravity['Omega_squared']/3, rtol=1e-13)
    area_difference = alpha*(y-first)
    assert np.isclose(np.max(np.abs(area_difference)),
                      alpha*case['maximum_scaled_area_difference_from_first_order'], rtol=1e-7, atol=0)
    exact = radiation['energy_joules']; dipole = radiation['dipole_energy_joules']
    assert 0 < exact < dipole and np.isclose(exact/dipole, radiation['exact_to_dipole_ratio'])
    scan = source['pulse_duration_scan']
    duration = np.asarray([row['duration_seconds'] for row in scan])*1e9
    ratio = np.asarray([row['exact_to_dipole_ratio'] for row in scan])
    assert np.all(np.diff(duration) > 0) and np.all(np.isfinite(ratio)) and np.all(ratio > 0)
    for row in scan:
        assert np.isclose(row['energy_joules']/row['dipole_energy_joules'], row['exact_to_dipole_ratio'])
    control = future['proton_retarded_control']
    width_ns = control['parameters']['source_pulse_duration_seconds']*1e9
    peak_MA = control['parameters']['peak_impressed_current_per_loop_amperes']/1e6
    assert np.isclose(duration[0], width_ns) and np.isclose(scan[0]['energy_joules'], exact)
    assert source['parameters']['peak_current_per_loop_amperes'] == peak_MA*1e6
    assert source['parameters']['source_pulse_duration_seconds'] == width_ns*1e-9
    assert source['parameters']['loop_radius_metres'] == control['parameters']['control_loop_radius_metres']
    assert source['parameters']['loop_centers_metres'] == control['parameters']['control_loop_centers_metres']

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.64, wspace=.30)
    fig.suptitle('A specified quantum source supports a static model; its apparatus is unresolved',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Finite throat equations solved together • electromagnetic source cost checked separately',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    position = np.linspace(-1, 1, 150); width = .30+.10*position**2
    ax.fill_between(position, -width, width, color='#007e87', alpha=.13)
    ax.plot(position, width, position, -width, color='#007e87', lw=2)
    for sign in (-1, 1):
        ax.plot([sign, sign], [-.54, .54], color='#326bb3', lw=2)
        ax.plot([sign, 1.38*sign], [0, 0], color='#94a3b8', ls='--', lw=1.2)
        ax.text(sign, -.81, f'x = {sign*x[-1]:.3f}', ha='center', fontsize=9, color='#475569')
        ax.text(1.50*sign, .24, 'Exterior\nunsolved', ha='center', fontsize=8.5, color='#64748b')
    ax.text(0, 0, 'Static throat patch', ha='center', va='center', color='#007e87', fontsize=11)
    ax.text(0, 1.13, 'Ideal 2D Polyakov mean source\nPrescribed loop state; radius-independent action',
            ha='center', fontsize=9, color='#475569', linespacing=1.5)
    ax.text(0, -1.15, 'All static equations checked in this source model\nGeometry sketch, not a physical material or an entrance',
            ha='center', va='top', fontsize=9, color='#475569', linespacing=1.5)
    ax.set(title='A. Assumed finite throat region — schematic',
           xlim=(-1.88, 1.88), ylim=(-1.85, 1.75))
    ax.axis('off')

    ax = axes[0, 1]
    ax.plot(x, area_difference, color='#007e87', lw=2.3)
    ax.axhline(0, color='#94a3b8', ls=':', lw=1)
    span = max(float(np.ptp(area_difference)), np.finfo(float).tiny)
    ax.set(title='B. Full solution minus the first-order area',
           xlabel='Position x from center toward one cutoff',
           ylabel='Difference in area / (4π rₑ²)', xlim=(0, x[-1]),
           ylim=(area_difference.min()-.12*span, area_difference.max()+.30*span))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    ax.text(.96, .95, f'α = {alpha:.0e}; reflection symmetry\nDifference only; total area stays positive',
            transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color='#475569')

    ax = axes[1, 0]
    energies = np.array([exact, dipole])/1e9
    ax.bar([0, 1], energies, width=.52, color=['#326bb3', '#b7791f'])
    for i, value in enumerate(energies):
        ax.text(i, value+.025*energies.max(), f'{value:.3f} GJ', ha='center',
                fontsize=10, fontweight='bold', color='#326bb3' if i == 0 else '#92400e')
    ax.set_xticks([0, 1], labels=['Full retardation\nand interference', 'Point-dipole\napproximation'])
    ax.set(title=f'C. Radiation from the same {width_ns:g} ns pulse',
           ylabel='Energy radiated by both control loops (GJ)', ylim=(0, energies.max()*1.22))
    ax.text(.03, .95, f'{peak_MA:g} MA peak impressed current per loop',
            transform=ax.transAxes, va='top', fontsize=8.5, color='#475569')

    ax = axes[1, 1]
    ax.semilogx(duration, ratio, 'o-', color='#7951b0', lw=2, ms=5)
    ax.scatter([duration[0]], [ratio[0]], s=55, color='#326bb3', zorder=4)
    ax.axhline(1, color='#64748b', ls=':', lw=1.2, label='Dipole agreement')
    ax.set_xticks(duration, labels=[f'{value:g}' for value in duration])
    ax.minorticks_off()
    ax.set(title='D. Longer pulses approach the dipole limit',
           xlabel='Pulse duration (ns; logarithmic scale)',
           ylabel='Radiated energy: full / dipole',
           xlim=(duration[0]/1.17, duration[-1]*1.18), ylim=(0, max(1., ratio.max())*1.15))
    ax.text(.97, .07, 'Geometry and peak current held fixed\n'
            'Radiation scan only; longer-pulse trajectories untested\n'
            'Current test: full / dipole = '+f'{exact/dipole:.3f}',
            transform=ax.transAxes, ha='right', fontsize=8.6, color='#475569', linespacing=1.5)
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    gain = control['electric_work']['gain_MeV']
    fig.text(.055, .13, 'Gravity: specified reduced quantum action; full 4D material, evolving stability, preparation and exterior matching remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, f'{width_ns:g} ns pulse: {exact/1e9:.3f} GJ radiated; one proton gained {gain:.2f} MeV. Hardware and complete circuit energy accounting remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is supported; the proposed apparatus remains unbuilt.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'polyakov_source_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def junction_pulse_progress():
    """Compare the restricted shell obstruction and a verified lower-radiation pulse."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'polyakov_junction' not in gravity or 'lower_radiation_control' not in future:
        raise ValueError('Final junction and selected-pulse results are required; no figure written')
    junction = gravity['polyakov_junction']
    selected = future['lower_radiation_control']
    control = selected['selected_control']
    radiation = selected['selected_source_feasibility']
    baseline_control = future['proton_retarded_control']
    baseline_radiation = future['source_feasibility']
    cases = junction['cases']
    assert len(cases) == 3 and np.isclose(cases[-1]['x'], np.arctan(4), atol=1e-14)
    for case in cases:
        pressure = case['all_nonnegative_energy_shells_pressure_lower_bound']
        energy = case['all_nonnegative_energy_shells_energy_upper_bound']
        assert np.isfinite(pressure) and np.isfinite(energy) and pressure > energy > 0
        assert np.isclose(pressure, case['pressure_from_normal_force_identity'], rtol=1e-12)
        assert case['derivative_matched_surface_energy_times_8piG_re'] == 0
    boundary = cases[-1]
    p_min = boundary['all_nonnegative_energy_shells_pressure_lower_bound']
    sigma_max = boundary['all_nonnegative_energy_shells_energy_upper_bound']
    route = control['plot']['route']
    us = np.asarray(route['lab_microseconds'], dtype=float)
    radius = np.asarray(route['radius_metres'], dtype=float)
    proper = np.asarray(route['proton_proper_microseconds'], dtype=float)
    no_pulse = np.asarray(route['zero_current_radius_metres'], dtype=float)
    assert len(us) == len(radius) == len(proper) == len(no_pulse)
    assert all(np.all(np.isfinite(a)) for a in (us, radius, proper, no_pulse))
    assert np.all(np.diff(us) > 0) and np.all(np.diff(proper) > 0)
    assert np.all(radius > 0) and np.all(no_pulse > 0) and us[-1] > proper[-1] > 0
    assert np.isclose(us[-1]*1e-6, control['endpoint']['lab_seconds'], rtol=1e-10, atol=1e-15)
    assert np.isclose(proper[-1]*1e-6, control['endpoint']['proton_proper_seconds'], rtol=1e-10, atol=1e-15)
    assert np.isclose(radius[-1], control['endpoint']['radius_metres'], rtol=1e-10)
    energy_joules = np.asarray([baseline_radiation['radiation']['energy_joules'],
                               radiation['radiation']['energy_joules']])
    assert np.all(np.isfinite(energy_joules)) and energy_joules[0] > energy_joules[1] > 0
    widths = []
    currents = []
    for trajectory, audit in ((baseline_control, baseline_radiation), (control, radiation)):
        parameters = trajectory['parameters']; source_parameters = audit['parameters']
        width = parameters['source_pulse_duration_seconds']
        current = parameters['peak_impressed_current_per_loop_amperes']
        assert width == source_parameters['source_pulse_duration_seconds']
        assert current == source_parameters['peak_current_per_loop_amperes']
        assert parameters['control_loop_radius_metres'] == source_parameters['loop_radius_metres']
        assert parameters['control_loop_centers_metres'] == source_parameters['loop_centers_metres']
        widths.append(width); currents.append(current/1e6)
    assert control['parameters']['storage_source'] == baseline_control['parameters']['storage_source']
    clear_us = control['pulse_clear']['lab_seconds']*1e6
    assert 0 < widths[1]*1e6 < clear_us < us[-1]
    reduction = energy_joules[0]/energy_joules[1]

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.66, wspace=.30)
    fig.suptitle('The exterior join remains unresolved; a revised pulse radiates less',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate gravity and proton-source calculations • ideal models • no complete apparatus',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    inside = np.linspace(-1.6, 0, 120)
    throat_width = .24+.10*((inside+1.6)/1.6)**2
    ax.fill_between(inside, -throat_width, throat_width, color='#007e87', alpha=.15)
    ax.plot(inside, throat_width, inside, -throat_width, color='#007e87', lw=2)
    exterior = np.linspace(0, 1.6, 100)
    exterior_width = throat_width[-1]+.18*exterior
    ax.plot(exterior, exterior_width, exterior, -exterior_width,
            color='#326bb3', lw=1.7, ls='--')
    ax.axvspan(-.045, .045, ymin=.22, ymax=.69, color='#b7791f', alpha=.30)
    ax.plot([0, 0], [-.65, .65], color='#b7791f', lw=2)
    ax.text(-.85, 0, 'Static throat\nwith quantum source', ha='center', va='center',
            fontsize=9.5, color='#007e87', linespacing=1.4)
    ax.text(.92, 0, 'Proposed same-charge\nRN exterior', ha='center', va='center',
            fontsize=9.5, color='#326bb3', linespacing=1.4)
    ax.annotate('Required boundary stress', xy=(0, .56), xytext=(0, 1.03),
                ha='center', fontsize=9.5, color='#92400e',
                arrowprops=dict(arrowstyle='->', color='#b7791f', lw=1.2))
    ax.text(0, -.93, 'Abrupt source removal is the tested assumption\n'
            'Exterior geometry and quantum state remain unresolved',
            ha='center', va='top', fontsize=8.8, color='#475569', linespacing=1.5)
    ax.set(title='A. Proposed junction — geometry schematic', xlim=(-1.83, 1.83), ylim=(-1.55, 1.55))
    ax.axis('off')

    ax = axes[0, 1]
    bound_values = np.asarray([p_min, sigma_max])
    ax.bar([0, 1], bound_values, width=.47, color=['#b7791f', '#007e87'])
    ax.set_yscale('log')
    ax.set_xticks([0, 1], labels=['Pressure lower bound', 'Energy upper bound'])
    ax.set_yticks([1e-9, 1e-6, 1e-3, 1], labels=['10⁻⁹', '10⁻⁶', '10⁻³', '1'])
    ax.minorticks_off()
    ax.set(title='B. Shell energy cannot cover the pressure',
           ylabel='Surface stress (shared normalization)',
           ylim=(sigma_max/8, p_min*90))
    ax.text(0, p_min*1.8, f'≈ {p_min:.5f}', ha='center', fontsize=10, color='#92400e', fontweight='bold')
    ax.text(1, sigma_max*1.8, f'≈ {sigma_max/1e-9:.3f} × 10⁻⁹', ha='center', fontsize=10,
            color='#007e87', fontweight='bold')
    ax.text(.5, .96, 'Outermost tested cutoff: x = arctan(4)',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569')
    ax.text(.5, .49, f'Bound ratio ≈ {p_min/sigma_max/1e6:.1f} million\n'
            'Fails the condition Σ ≥ |p| (surface DEC)', transform=ax.transAxes,
            ha='center', va='center', fontsize=8.7, color='#475569', linespacing=1.5,
            bbox=dict(facecolor='white', alpha=.94, edgecolor='none', pad=3))

    ax = axes[1, 0]
    ax.plot(us, radius, color='#326bb3', lw=2.2, label='Selected finite current pulse')
    ax.plot(us, no_pulse, color='#64748b', ls=':', lw=1.4, label='Zero-current control')
    ax.axvline(clear_us, color='#b7791f', ls='--', lw=1.1,
               label=f'Control field clears proton: {clear_us:.2f} μs')
    ax.scatter([us[-1]], [radius[-1]], s=25, color='#326bb3', zorder=4)
    ax.set_yscale('log')
    ticks = np.asarray([10, 30, 100, 1000, 5000], dtype=float)
    ticks = ticks[(ticks >= radius.min()*.85) & (ticks <= radius.max()*1.3)]
    ax.set_yticks(ticks, labels=[f'{value:g}' for value in ticks]); ax.minorticks_off()
    ax.set(title='C. The revised pulse has a continuous route',
           xlabel='Elapsed lab time (microseconds)', ylabel='Radius from storage axis (m; log scale)',
           xlim=(0, us[-1]*1.05), ylim=(min(radius.min(), no_pulse.min())*.85, radius.max()*1.4))
    ax.legend(framealpha=.95, edgecolor='none', fontsize=7.8, loc='center right')
    ax.text(.97, .12, f'At {radius[-1]/1000:g} km: lab {us[-1]:.3f} μs\n'
            f'Proton proper time {proper[-1]:.3f} μs\nSource ends at {widths[1]*1e6:g} μs',
            transform=ax.transAxes, ha='right', fontsize=8.8, color='#326bb3', linespacing=1.5)

    ax = axes[1, 1]
    energy_MJ = energy_joules/1e6
    ax.bar([0, 1], energy_MJ, width=.5, color=['#94a3b8', '#7951b0'])
    for i, value in enumerate(energy_MJ):
        ax.text(i, value*1.25, f'{value:,.3f} MJ', ha='center', fontsize=10,
                color='#475569' if i == 0 else '#7951b0', fontweight='bold')
    width_labels = [f'{width*1e9:g} ns' if width < 1e-6 else f'{width*1e6:g} μs' for width in widths]
    ax.set_xticks([0, 1], labels=[f'Baseline\n{width_labels[0]}, {currents[0]:g} MA',
                                 f'Selected test\n{width_labels[1]}, {currents[1]:g} MA'])
    ax.set_yscale('log'); ax.minorticks_off()
    ax.set(title='D. Less radiation for the rechecked pulse',
           ylabel='Radiated energy (MJ; log scale)',
           ylim=(energy_MJ.min()/5, energy_MJ.max()*12))
    ax.text(.97, .95, f'{reduction:,.0f}× lower radiated energy\n'
            'Both loops; full retardation and interference\n'
            'Width, current and initial phase changed', transform=ax.transAxes,
            ha='right', va='top', fontsize=8.7, color='#475569', linespacing=1.5)
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    fig.text(.055, .13, 'Junction: restricted static, outward RN family; all three tested cutoffs fail this condition. Other continuations remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Pulse: ideal imposed currents; radiation is only part of the source budget. Finite conductors, capture and one-day storage are unmodeled.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a working device remains unbuilt.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'junction_pulse_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def finite_paths_conductors_progress():
    """Plot final finite return-path and isolated control-conductor diagnostics."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'magnetic_return_paths' not in gravity or 'finite_conductor' not in future:
        raise ValueError('Final return-path and conductor results are required; no figure written')
    paths = gravity['magnetic_return_paths']
    conductor = future['finite_conductor']
    selected = future['lower_radiation_control']
    selected_control = selected['selected_control']
    radiated_MJ = selected['selected_source_feasibility']['radiation']['energy_joules']/1e6
    flux_labels = (.5, 1., 1.5, 1.9)
    samples = []
    for label in flux_labels:
        matches = [row for row in paths['field_line_samples'] if np.isclose(row['nu'], label)]
        assert len(matches) == 1
        row = matches[0]
        rho = np.asarray(row['rho_over_d'], dtype=float)
        z = np.asarray(row['z_over_d'], dtype=float)
        assert len(rho) == len(z) and len(rho) > 2
        assert np.all(np.isfinite(rho)) and np.all(np.isfinite(z)) and np.all(rho >= 0)
        assert np.allclose(rho[[0, -1]], 0, atol=1e-14)
        assert np.allclose(z[[0, -1]], [-.5, .5], atol=1e-14)
        assert row['length_over_separation'] >= 2/label
        samples.append((label, rho, z))
    support_cases = [row for row in paths['fixed_scale_cases']
                     if np.isclose(row['separation_over_ell'], .1)]
    assert len(support_cases) == 1
    support = support_cases[0]
    old_beta = support['common_shortest_path_beta_over_alpha']
    new_beta = support['mean_beta_over_alpha']
    assert old_beta > new_beta > 0
    assert np.isclose(new_beta, (support['mean_omega_squared']-1)/3, rtol=1e-13)
    assert np.isclose(old_beta, (support['common_shortest_path_omega_squared']-1)/3, rtol=1e-13)
    reduction_percent = 100*(1-new_beta/old_beta)
    rows = conductor['rows']
    wire_mm = np.asarray([row['wire_radius_metres'] for row in rows], dtype=float)*1000
    stored_MJ = np.asarray([row['control_pair_peak_quasistatic_magnetic_energy_joules']
                            for row in rows], dtype=float)/1e6
    emf_GV = np.asarray([row['quasistatic_peak_inductive_emf_per_loop_volts']
                         for row in rows], dtype=float)/1e9
    assert len(wire_mm) == len(stored_MJ) == len(emf_GV) and len(wire_mm) > 1
    assert all(np.all(np.isfinite(a)) and np.all(a > 0) for a in (wire_mm, stored_MJ, emf_GV))
    assert np.all(np.diff(wire_mm) > 0) and np.all(np.diff(stored_MJ) < 0)
    assert np.all(np.diff(emf_GV) < 0) and 0 < radiated_MJ < stored_MJ.min()
    width_us = selected_control['parameters']['source_pulse_duration_seconds']*1e6
    peak_MA = selected_control['parameters']['peak_impressed_current_per_loop_amperes']/1e6
    source_parameters = conductor['parameters']
    assert source_parameters['source_pulse_duration_seconds'] == selected_control['parameters']['source_pulse_duration_seconds']
    assert source_parameters['signed_peak_current_per_loop_amperes'] == [peak_MA*1e6, peak_MA*1e6]
    assert source_parameters['loop_major_radius_metres'] == selected_control['parameters']['control_loop_radius_metres']
    assert source_parameters['loop_centers_metres'] == selected_control['parameters']['control_loop_centers_metres']
    assert conductor['verification']['positive_inductance_matrix_all_rows']
    for row in rows:
        inductance = row['equal_current_quasistatic_inductance_per_loop_henries']
        assert np.isclose(row['control_pair_peak_quasistatic_magnetic_energy_joules'],
                          inductance*(peak_MA*1e6)**2, rtol=1e-12)
        assert np.isclose(row['quasistatic_peak_inductive_emf_per_loop_volts'],
                          inductance*source_parameters['maximum_abs_dI_dt_amperes_per_second'], rtol=1e-12)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Finite paths lower average support; conductor demands remain large',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Leading exterior-field model and control-coil estimates • separate calculations • no built medium',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    colors = ['#326bb3', '#007e87', '#7951b0', '#b7791f']
    for (label, rho, z), color in zip(samples, colors):
        ax.plot(rho, z, color=color, lw=1.8, label=f'ν = {label:g}')
    ax.scatter([0, 0], [-.5, .5], s=28, color='#1e293b', zorder=5, label='Ideal mouths')
    rho_max = max(float(rho.max()) for _, rho, _ in samples)
    z_max = max(float(np.abs(z).max()) for _, _, z in samples)
    ax.set(title='A. Finite exterior magnetic field paths',
           xlabel='Radial position ρ / d', ylabel='Axial position z / d',
           xlim=(-.05*rho_max, rho_max*1.05), ylim=(-z_max*1.12, z_max*1.12))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([0, 1, 2]); ax.set_yticks([-.5, 0, .5])
    ax.legend(frameon=False, fontsize=8, loc='center left', bbox_to_anchor=(1.01, .5))
    ax.text(.5, -.30, 'd: mouth separation • ν: channel label', transform=ax.transAxes,
            ha='center', fontsize=8.2, color='#475569')

    ax = axes[0, 1]
    support_values = np.asarray([old_beta, new_beta])
    ax.bar([0, 1], support_values, width=.5, color=['#94a3b8', '#007e87'])
    for i, value in enumerate(support_values):
        ax.text(i, value+.026, f'{value:.4f}', ha='center', fontweight='bold',
                fontsize=10, color='#475569' if i == 0 else '#007e87')
    ax.set_xticks([0, 1], labels=['Shared shortest\nexterior path', 'Full distribution\nof exterior paths'])
    ax.set(title='B. Average support at the same throat scale',
           ylabel='Leading support coefficient β / α', ylim=(0, old_beta*1.30))
    ax.text(.98, .95, f'{reduction_percent:.1f}% lower\nd / ℓ = {support["separation_over_ell"]:g}',
            transform=ax.transAxes, ha='right', va='top', fontsize=9.5, color='#007e87', linespacing=1.4)
    ax.text(.5, .09, 'Angular variation still requires a geometry solution',
            transform=ax.transAxes, ha='center', fontsize=8.5, color='#475569',
            bbox=dict(facecolor='white', alpha=.94, edgecolor='none', pad=3))

    ax = axes[1, 0]
    ax.plot(wire_mm, stored_MJ, 'o-', color='#7951b0', lw=2, ms=5,
            label='Peak quasistatic stored energy')
    ax.axhline(radiated_MJ, color='#326bb3', ls='--', lw=1.5,
               label=f'Selected pulse radiation: {radiated_MJ:.3f} MJ')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xticks(wire_mm, labels=[f'{value:g}' for value in wire_mm])
    ax.minorticks_off()
    ax.set(title='C. Energy scales for the control pair only',
           xlabel='Wire radius (mm; log scale)', ylabel='Energy (MJ; log scale)',
           xlim=(wire_mm[0]/1.28, wire_mm[-1]*1.28),
           ylim=(radiated_MJ/3, stored_MJ.max()*3.3))
    ax.legend(framealpha=.95, edgecolor='none', fontsize=8, loc='center right')
    ax.text(.5, .96, 'Stored energy is not dissipation or a full apparatus budget',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569')

    ax = axes[1, 1]
    ax.semilogx(wire_mm, emf_GV, 'o-', color='#b7791f', lw=2, ms=5)
    for x, value in zip(wire_mm, emf_GV):
        ax.text(x, value+.035*emf_GV.max(), f'{value:.2f}', ha='center',
                fontsize=9, color='#92400e')
    ax.set_xticks(wire_mm, labels=[f'{value:g}' for value in wire_mm]); ax.minorticks_off()
    ax.set(title='D. Quasistatic voltage scale — estimate only',
           xlabel='Wire radius (mm; log scale)', ylabel='Inductive EMF (GV per loop)',
           xlim=(wire_mm[0]/1.28, wire_mm[-1]*1.28), ylim=(0, emf_GV.max()*1.24))
    ax.text(.98, .09, f'{width_us:g} μs pulse; {peak_MA:g} MA peak per loop\n'
            'Loop EMF estimate; feed voltage is unmodeled', transform=ax.transAxes,
            ha='right', fontsize=8.8, color='#475569', linespacing=1.5)
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    fig.text(.055, .13, 'Paths: leading flat-exterior approximation. Mouth matching, angular gravity response, field leakage and state preparation remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Conductors: isolated control-pair estimates; original storage-field cross terms, supports, material response and the driver remain unmodeled.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a working device remains unbuilt.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'finite_paths_conductors_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def optical_localized_progress():
    """Render final optical-cutoff and localized-control conductor comparisons."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'optical_matching' not in gravity or 'localized_control' not in future:
        raise ValueError('Final optical and localized-control results are required; no figure written')
    optical = gravity['optical_matching']
    cutoff = optical['outer_cutoff_comparison']
    raw_span = cutoff['raw_span_over_re']
    corrected_span = cutoff['fixed_path_weak_field_span_over_re']
    assert np.isfinite(raw_span) and np.isfinite(corrected_span) and raw_span > corrected_span > 0
    assert np.isclose(raw_span/corrected_span, cutoff['span_reduction_factor'], rtol=1e-13)
    matching_rows = [row for row in optical['matching_cases']
                     if np.isclose(row['inner_throat_rho'], cutoff['inner_rho'])]
    assert len(matching_rows) == len(optical['parameters']['outer_areal_radii_over_re'])
    assert np.isclose(np.ptp([row['leading_flat_length_residual_over_re'] for row in matching_rows]),
                      raw_span, rtol=1e-12)
    assert np.isclose(np.ptp([row['fixed_path_weak_field_length_residual_over_re']
                             for row in matching_rows]), corrected_span, rtol=1e-10)
    localized = future['localized_control']
    old_conductor = future['finite_conductor']
    new_conductor = localized['selected_conductor']
    fixed_conductor = localized['fixed_radius_conductor']

    def unique_row(model, key, value):
        matches = [row for row in model['rows'] if np.isclose(row[key], value, rtol=1e-10, atol=0)]
        assert len(matches) == 1
        return matches[0]

    old_row = unique_row(old_conductor, 'wire_radius_over_loop_radius', .02)
    new_row = unique_row(new_conductor, 'wire_radius_over_loop_radius', .02)
    fixed_row = unique_row(fixed_conductor, 'wire_radius_metres', old_row['wire_radius_metres'])
    models = (old_conductor, new_conductor, fixed_conductor)
    rows = (old_row, new_row, fixed_row)
    trajectories = (future['lower_radiation_control']['selected_control'],
                    localized['selected_control'], localized['selected_control'])
    for model, row, trajectory in zip(models, rows, trajectories):
        parameters = model['parameters']; control_parameters = trajectory['parameters']
        current = control_parameters['peak_impressed_current_per_loop_amperes']
        assert parameters['signed_peak_current_per_loop_amperes'] == [current, current]
        assert parameters['source_pulse_duration_seconds'] == control_parameters['source_pulse_duration_seconds']
        assert parameters['loop_major_radius_metres'] == control_parameters['control_loop_radius_metres']
        assert parameters['loop_centers_metres'] == control_parameters['control_loop_centers_metres']
        assert np.isclose(row['wire_radius_metres']/parameters['loop_major_radius_metres'],
                          row['wire_radius_over_loop_radius'], rtol=1e-12)
        inductance = row['equal_current_quasistatic_inductance_per_loop_henries']
        assert np.isclose(row['control_pair_peak_quasistatic_magnetic_energy_joules'],
                          inductance*current**2, rtol=1e-12)
        assert np.isclose(row['quasistatic_peak_inductive_emf_per_loop_volts'],
                          inductance*parameters['maximum_abs_dI_dt_amperes_per_second'], rtol=1e-12)
    stored_MJ = np.asarray([row['control_pair_peak_quasistatic_magnetic_energy_joules'] for row in rows])/1e6
    emf = np.asarray([row['quasistatic_peak_inductive_emf_per_loop_volts'] for row in rows])
    hoop = np.asarray([row['pair_only_equivalent_hoop_stress_pascals'] for row in rows])
    assert all(np.all(np.isfinite(a)) and np.all(a > 0) for a in (stored_MJ, emf, hoop))
    assert stored_MJ[0] > max(stored_MJ[1:])
    assert old_row['wire_radius_metres'] > new_row['wire_radius_metres']
    emf_ratio = emf[1:]/emf[0]; hoop_ratio = hoop[1:]/hoop[0]
    old_a = old_conductor['parameters']['loop_major_radius_metres']
    new_a = new_conductor['parameters']['loop_major_radius_metres']
    old_I = old_conductor['parameters']['signed_peak_current_per_loop_amperes'][0]/1e6
    new_I = new_conductor['parameters']['signed_peak_current_per_loop_amperes'][0]/1e6
    b_mm = np.asarray([row['wire_radius_metres'] for row in rows])*1000

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Optical matching and smaller control coils: what changed',
                 fontsize=19, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Two independent model checks • approximate geometry and conductor estimates • no built medium',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    segments = [(np.linspace(-1.7, -.45, 70), '#007e87'),
                (np.linspace(-.45, .50, 70), '#b7791f'),
                (np.linspace(.50, 1.7, 70), '#326bb3')]
    for position, color in segments:
        width = .20+.055*(position+1.7)**2
        ax.fill_between(position, -width, width, color=color, alpha=.10)
        ax.plot(position, width, position, -width, color=color, lw=2)
    for position in (-.45, .50):
        ax.plot([position, position], [-.93, .95], color='#64748b', ls=':', lw=1.3)
    ax.text(-1.10, 0, 'Throat', ha='center', va='center', fontsize=10.5, color='#007e87')
    ax.text(.025, 0, 'RN mouth', ha='center', va='center', fontsize=10.5, color='#92400e')
    ax.text(1.14, 0, 'Weak-field\nexterior', ha='center', va='center', fontsize=10.5,
            color='#326bb3', linespacing=1.4)
    ax.text(-.45, 1.05, 'Inner match', ha='center', fontsize=8.8, color='#475569')
    ax.text(.50, 1.05, 'Outer match', ha='center', fontsize=8.8, color='#475569')
    ax.text(0, -1.18, 'Each region uses its own approximation\n'
            'Concept diagram; no completed physical junction', ha='center', va='top',
            fontsize=9, color='#475569', linespacing=1.5)
    ax.set(title='A. Matching three approximate regions', xlim=(-1.85, 1.85), ylim=(-1.78, 1.55))
    ax.axis('off')

    ax = axes[0, 1]
    spans = np.asarray([raw_span, corrected_span])
    ax.bar([0, 1], spans, width=.50, color=['#94a3b8', '#007e87'])
    ax.set_yscale('log'); ax.minorticks_off()
    for i, value in enumerate(spans):
        ax.text(i, value*1.3, f'{value:.5f}', ha='center', fontsize=10, fontweight='bold',
                color='#475569' if i == 0 else '#007e87')
    ax.set_xticks([0, 1], labels=['Flat exterior', 'Exterior optical\ncorrection included'])
    ax.set(title='B. Dependence on outer matching radius',
           ylabel='Cutoff-induced length span / rₑ', ylim=(spans.min()/5, spans.max()*6))
    ax.text(.97, .96, f'Outer radius: {min(optical["parameters"]["outer_areal_radii_over_re"]):g}–'
            f'{max(optical["parameters"]["outer_areal_radii_over_re"]):,.0f} rₑ\n'
            'Fixed inner cutoff; one magnetic path', transform=ax.transAxes,
            ha='right', va='top', fontsize=8.2, color='#475569', linespacing=1.4)
    ax.text(.5, .45, f'{raw_span/corrected_span:.1f}× smaller variation\n'
            'NOT an accuracy or stability proof', transform=ax.transAxes,
            ha='center', va='center', fontsize=9, color='#92400e', fontweight='bold', linespacing=1.5,
            bbox=dict(facecolor='white', alpha=.95, edgecolor='none', pad=3))

    ax = axes[1, 0]
    ax.bar([0, 1, 2], stored_MJ, width=.53, color=['#94a3b8', '#7951b0', '#007e87'])
    for i, value in enumerate(stored_MJ):
        ax.text(i, value+.026*stored_MJ.max(), f'{value:.2f} MJ', ha='center',
                fontsize=9.2, fontweight='bold', color=('#475569', '#7951b0', '#007e87')[i])
    ax.set_xticks([0, 1, 2], labels=[f'Previous pair\nb = {b_mm[0]:g} mm',
                                    f'Same b/a = 0.02\nb = {b_mm[1]:g} mm',
                                    f'Fixed b = {b_mm[2]:g} mm\nb/a = {fixed_row["wire_radius_over_loop_radius"]:g}; less slender'])
    ax.set(title='C. Peak stored energy: control pair only',
           ylabel='Quasistatic magnetic energy (MJ)', ylim=(0, stored_MJ.max()*1.32))
    ax.text(.98, .96, f'Coil radius a: {old_a:g} → {new_a:g} m\n'
            f'Peak current: {old_I:g} → {new_I:g} MA', transform=ax.transAxes,
            ha='right', va='top', fontsize=8.4, color='#475569', linespacing=1.5)

    ax = axes[1, 1]; positions = np.arange(2)
    ax.bar(positions-.16, emf_ratio, width=.30, color='#326bb3', label='Inductive EMF per loop')
    ax.bar(positions+.16, hoop_ratio, width=.30, color='#b7791f', label='Equivalent hoop stress')
    for offset, ratios, color in ((-.16, emf_ratio, '#326bb3'), (.16, hoop_ratio, '#92400e')):
        for x, value in zip(positions+offset, ratios):
            ratio_label = f'{value:.2f}' if value >= 1 else f'{value:.3g}'
            ax.text(x, value*1.20, f'{ratio_label}×', ha='center', fontsize=9, color=color)
    ax.axhline(1, color='#64748b', ls=':', lw=1.2)
    ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks(positions, labels=['Same wire / coil ratio', 'Fixed wire radius'])
    ax.set(title='D. Smaller coils change the conductor tradeoff',
           ylabel='Ratio to previous pair (log scale)',
           ylim=(min(emf_ratio.min(), hoop_ratio.min())/2.5, max(emf_ratio.max(), hoop_ratio.max())*7))
    ax.legend(framealpha=.95, edgecolor='none', fontsize=8, loc='upper right')
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    fig.text(.055, .13, 'Geometry: reduced cutoff sensitivity does not update the quantum source or establish a joined spacetime, stable medium or time offset.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Conductors: stored energy is not dissipation; EMF is not feed voltage. Pair-only hoop estimates assume a solid wire; survival is unproved.',
             fontsize=9, color='#475569')
    fig.text(.055, .059, 'Both wire choices remain cost estimates; their actual electromagnetic fields have not been used to recheck the proton trajectory.',
             fontsize=8.5, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a working device remains unbuilt.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'optical_localized_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def finite_wire_matching_progress():
    """Show final prescribed torus sources, unretuned motion, and one matching cancellation."""
    from matplotlib.patches import Circle

    root = Path(__file__).resolve().parent
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    if 'finite_torus_control' not in future or 'weak_magnetic_response' not in gravity:
        raise ValueError('Final finite-torus and weak-metric results are required; no figure written')
    surface = future['finite_torus_control']
    assert surface['retuning_performed'] is False
    cases = sorted(surface['cases'], key=lambda row: row['minor_radius_metres'])
    radii = np.asarray([row['minor_radius_metres'] for row in cases], dtype=float)
    assert len(cases) == 2 and np.allclose(radii, [.025, .1], rtol=1e-12, atol=0)
    assert np.allclose(radii, sorted(surface['tested_minor_radii_metres']), rtol=1e-12, atol=0)
    baseline = future['localized_control']['selected_control']
    parameters = baseline['parameters']
    major = parameters['control_loop_radius_metres']
    centers = np.asarray(parameters['control_loop_centers_metres'], dtype=float)
    assert centers.shape == (2, 3)
    assert np.allclose(centers[:, :2], centers[0, :2], rtol=0, atol=1e-12)
    half_separation = abs(centers[1, 2]-centers[0, 2])/2
    assert major > radii.max() > 0 and half_separation > radii.max()
    controls = [baseline]+[row['control'] for row in cases]
    labels = ['Filament reference']+[f'Surface: b = {b*1000:g} mm' for b in radii]
    lab_us = []; proper_us = []
    for control in controls:
        current_parameters = control['parameters']
        for key in ('control_loop_radius_metres', 'peak_impressed_current_per_loop_amperes',
                    'source_pulse_duration_seconds', 'starting_azimuth_radians', 'initial_gamma'):
            assert np.isclose(current_parameters[key], parameters[key], rtol=1e-12, atol=0)
        assert current_parameters['control_loop_centers_metres'] == parameters['control_loop_centers_metres']
        endpoint = control['endpoint']; route = control['plot']['route']
        assert np.isclose(endpoint['radius_metres'], parameters['outgoing_checkpoint_metres'], rtol=1e-10)
        time = np.asarray(route['lab_microseconds'], dtype=float)
        radius = np.asarray(route['radius_metres'], dtype=float)
        proper = np.asarray(route['proton_proper_microseconds'], dtype=float)
        assert len(time) == len(radius) == len(proper) and len(time) > 2
        assert all(np.all(np.isfinite(a)) for a in (time, radius, proper))
        assert np.all(np.diff(time) > 0) and np.all(np.diff(proper) > 0)
        assert np.isclose(time[-1]*1e-6, endpoint['lab_seconds'], rtol=1e-10, atol=1e-15)
        assert np.isclose(proper[-1]*1e-6, endpoint['proton_proper_seconds'], rtol=1e-10, atol=1e-15)
        assert np.isclose(radius[-1], endpoint['radius_metres'], rtol=1e-10)
        assert 0 < control['pulse_clear']['lab_seconds'] < endpoint['lab_seconds']
        lab_us.append(endpoint['lab_seconds']*1e6)
        proper_us.append(endpoint['proton_proper_seconds']*1e6)
    shifts = [0.]
    for case in cases:
        shift = case['comparison_to_filament']['endpoint_position_difference_metres']
        difference = (np.asarray(case['control']['endpoint']['position_metres'])-
                      np.asarray(baseline['endpoint']['position_metres']))
        assert np.isclose(np.linalg.norm(difference), shift, rtol=1e-11)
        assert case['control']['parameters']['control_torus_minor_radius_metres'] == case['minor_radius_metres']
        shifts.append(shift)
    annulus = gravity['weak_magnetic_response']['companion_annulus']
    negative = annulus['previously_missing_companion_span_over_re']
    positive = annulus['universal_two_mouth_length_change_over_re']
    residual = annulus['universal_cancellation_residual_over_re']
    assert np.isfinite(negative) and np.isfinite(positive) and negative < 0 < positive
    assert np.isclose(negative+positive, residual, rtol=0, atol=1e-15*positive)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 9))
    grid = fig.add_gridspec(2, 2)
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Finite magnetic source surfaces tested; one matching term cancels',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Forward magnetic controls and a separate static gravity diagnostic • no demonstrated time-travel medium',
             ha='center', fontsize=10.5, color='#475569')

    ax3d = fig.add_subplot(grid[0, 0], projection='3d')
    phi, theta = np.meshgrid(np.linspace(0, 2*np.pi, 73), np.linspace(0, 2*np.pi, 29))
    b_show = radii.max()
    tube_x = (major+b_show*np.cos(theta))*np.cos(phi)
    tube_y = (major+b_show*np.cos(theta))*np.sin(phi)
    centerline_phi = np.linspace(0, 2*np.pi, 160)
    for center_z in (-half_separation, half_separation):
        tube_z = center_z+b_show*np.sin(theta)
        ax3d.plot_surface(tube_x, tube_y, tube_z, color='#7951b0', alpha=.42,
                          edgecolor='none', rstride=1, cstride=1, shade=True, antialiased=True)
        ax3d.plot(major*np.cos(centerline_phi), major*np.sin(centerline_phi),
                  np.full_like(centerline_phi, center_z), color='#475569', ls='--', lw=1)
    limit = max(major+b_show, half_separation+b_show)*1.06
    ax3d.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit))
    ax3d.set_box_aspect((1, 1, 1)); ax3d.view_init(elev=18, azim=-48)
    ax3d.set_title('A. Forward magnetic control sources', pad=12)
    ax3d.set_axis_off()
    ax3d.text2D(.5, -.04, f'Two torus surfaces: a = {major:g} m; separation = {2*half_separation:g} m\n'
                 f'Shown b = {b_show*1000:g} mm; dashed lines = prior filaments',
                 transform=ax3d.transAxes, ha='center', fontsize=8.8, color='#475569', linespacing=1.5)

    ax = fig.add_subplot(grid[0, 1])
    radius_colors = ['#007e87', '#7951b0']
    for b, color in zip(radii, radius_colors):
        ax.add_patch(Circle((0, 0), b*1000, fill=False, edgecolor=color, lw=2,
                            label=f'b = {b*1000:g} mm'))
    ax.scatter([0], [0], s=25, color='#475569', zorder=4, label='Filament')
    radius_mm = radii.max()*1000
    ax.set(title='B. Cross-section of one current tube',
           xlabel='Radial offset (mm)', ylabel='Axial offset (mm)',
           xlim=(-radius_mm*1.13, radius_mm*1.13), ylim=(-radius_mm*1.13, radius_mm*1.13))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([-100, 0, 100]); ax.set_yticks([-100, 0, 100])
    ax.legend(frameon=False, fontsize=8, loc='center left', bbox_to_anchor=(1.01, .5))
    ax.text(.5, -.29, 'Current is prescribed on each surface\nMaterial response is not modeled',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False)

    ax = fig.add_subplot(grid[1, 0]); ax.axis('off')
    ax.set_title('C. Proton reaches 5 km for all three sources', pad=10)
    table_rows = [[label, f'{lab:.6f}', f'{proper:.6f}', f'{shift:.2f}']
                  for label, lab, proper, shift in zip(labels, lab_us, proper_us, shifts)]
    table = ax.table(cellText=table_rows,
                     colLabels=['Source', 'Lab (μs)', 'Proper (μs)', 'Shift at 5 km (m)'],
                     cellLoc='center', colWidths=[.31, .23, .23, .23], bbox=[0, .36, 1, .56])
    table.auto_set_font_size(False); table.set_fontsize(8.5)
    for (row_index, col_index), cell in table.get_celld().items():
        cell.set_edgecolor('#e2e8f0'); cell.set_linewidth(.7)
        if row_index == 0:
            cell.set_facecolor('#f1f5f9'); cell.set_text_props(color='#334155', fontweight='bold')
        else:
            cell.set_facecolor('white'); cell.set_text_props(color=('#475569', '#007e87', '#7951b0')[row_index-1])
    ax.text(.5, .28, f'Same {parameters["source_pulse_duration_seconds"]*1e6:g} μs pulse, '
            f'{parameters["peak_impressed_current_per_loop_amperes"]/1e6:g} MA, starting phase and orbit\n'
            'Shift: position change from the filament at its 5 km crossing\n'
            'Finite-source differences are physical model changes', transform=ax.transAxes,
            ha='center', va='top', fontsize=9, color='#475569', linespacing=1.5)

    ax = fig.add_subplot(grid[1, 1])
    contributions = [negative, positive, residual]
    ax.bar([0, 1, 2], contributions, width=.5, color=['#7951b0', '#007e87', '#94a3b8'])
    ax.axhline(0, color='#64748b', lw=1)
    for i, value in enumerate(contributions):
        offset = positive*.08 if value >= 0 else -positive*.08
        label = f'{value:+.5e}' if i < 2 else f'{value:g}'
        ax.text(i, value+offset, label, ha='center', va='bottom' if value >= 0 else 'top',
                fontsize=8.7, color=('#7951b0', '#007e87', '#475569')[i])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    ax.set_xticks([0, 1, 2], labels=['Exterior term', 'Reference-clock\n+ radial rescaling', 'Leading sum'])
    ax.set(title='D. One leading static matching term cancels',
           ylabel='Optical-length contribution / rₑ', ylim=(-positive*1.65, positive*1.65))
    ax.text(.5, .96, 'Companion contribution on a fixed annulus', transform=ax.transAxes,
            ha='center', va='top', fontsize=8.4, color='#475569')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    fig.text(.055, .13, 'Forward controls: prescribed surface currents now generate the trajectory fields. Feeds, current redistribution and material survival remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Gravity: one leading static clock/coordinate term cancels. Full geometry, quantum source, orbital motion and stability remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; these rings are magnetic controls, not a spacetime portal.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'finite_wire_matching_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def boundary_launch_progress():
    """Render final fixed-background magnetic response and sampled launch sensitivity."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'magnetic_boundary_response' not in gravity or 'torus_launch_perturbations' not in future:
        raise ValueError('Final boundary-response and launch results are required; no figure written')
    magnetic = gravity['magnetic_boundary_response']
    matching = [row for row in magnetic['leading_cases']
                if row['degree'] == 1 and np.isclose(row['throat_cutoff_rho'], 4)]
    assert len(matching) == 1
    mode = matching[0]; cutoff = np.arctan(mode['throat_cutoff_rho'])
    coordinate = np.linspace(-cutoff, cutoff, 401)
    even_cutoff = 1+cutoff*np.tan(cutoff)
    odd_cutoff = np.tan(cutoff)
    profile = ((1+coordinate*np.tan(coordinate))/(2*even_cutoff)
               + np.tan(coordinate)/(2*odd_cutoff))
    derivative = ((np.tan(coordinate)+coordinate/np.cos(coordinate)**2)/(2*even_cutoff)
                  + 1/(2*odd_cutoff*np.cos(coordinate)**2))
    center = mode['right_unit_left_zero_center_value']
    assert np.all(np.isfinite(profile)) and np.all(np.diff(profile) > 0)
    assert np.allclose(profile[[0, -1]], [0, 1], rtol=0, atol=1e-14)
    assert np.isclose(profile[len(profile)//2], center, rtol=1e-13)
    assert np.isclose(derivative[-1], mode['diagonal'], rtol=1e-13)
    assert np.isclose(-derivative[0], mode['off_diagonal'], rtol=1e-12)
    assert mode['even_eigenvalue'] > 0 and mode['odd_eigenvalue'] > 0
    launch = future['torus_launch_perturbations']
    assert launch['source_retuning_performed'] is False
    reference_endpoint = launch['reference']['endpoint']
    audited = {record['case_id']: record for record in launch['audited_cases']}
    assert len(audited) == len(launch['audited_cases'])
    combined = audited['combined_3d']['fine']
    endpoint = combined['endpoint']
    delta_position = (np.asarray(endpoint['position_metres'], dtype=float)
                      - np.asarray(reference_endpoint['position_metres'], dtype=float))
    assert delta_position.shape == (3,) and np.all(np.isfinite(delta_position))
    assert np.linalg.norm(delta_position) > 0
    assert endpoint['lab_seconds'] > endpoint['proton_proper_seconds'] > 0
    assert np.isclose(endpoint['cylindrical_radius_metres'], 5000, rtol=1e-10)
    assert np.isclose(np.linalg.norm(delta_position),
                      combined['arrival_comparison']['endpoint_position_difference_metres'], rtol=1e-11)
    conditions = combined['launch_conditions']
    for key, value in [('radial_offset_metres', .001), ('vertical_offset_metres', .001),
                       ('radial_direction_radians', 1e-5), ('vertical_direction_radians', 1e-5)]:
        assert np.isclose(conditions[key], value, rtol=1e-12, atol=0)
    timing_offsets = [-1e-10, 1e-10, 3.75e-8, 5e-8]
    timing_cases = []
    for offset in timing_offsets:
        matches = [row for row in launch['screening_cases']
                   if str(row['case_id']).startswith('timing_')
                   and np.isclose(row['phase_equivalent_timing_error_seconds'], offset,
                                  rtol=1e-10, atol=1e-18)]
        assert len(matches) == 1
        selected = matches[0]
        if selected['case_id'] in audited:
            selected = audited[selected['case_id']]['fine']
            assert np.isclose(selected['phase_equivalent_timing_error_seconds'], offset, rtol=1e-12, atol=0)
        timing_cases.append(selected)
    successful = [i for i, row in enumerate(timing_cases) if row.get('arrival_comparison') is not None]
    unsuccessful = [i for i in range(len(timing_cases)) if i not in successful]
    shifts = np.asarray([timing_cases[i]['arrival_comparison']['endpoint_position_difference_metres']
                         for i in successful], dtype=float)
    assert np.all(np.isfinite(shifts)) and np.all(shifts > 0)
    assert len(successful) == 3 and unsuccessful == [3]

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Boundary magnetic response checked; launch changes move arrivals',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'A static probe on assumed geometry • separate forward-control launch tests • no demonstrated time machine',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    position = np.linspace(-1, 1, 140); width = .26+.10*position**2
    ax.fill_between(position, -width, width, color='#007e87', alpha=.11)
    ax.plot(position, width, position, -width, color='#007e87', lw=2)
    for side, color in ((-1, '#64748b'), (1, '#326bb3')):
        ax.plot([side, side], [-.66, .66], color=color, lw=2)
        ax.text(side, .87, f'Boundary input\na = {(side+1)//2}', ha='center',
                fontsize=10, color=color, linespacing=1.4)
    ax.scatter([0], [0], s=21, color='#007e87')
    ax.text(0, .48, 'Prescribed throat geometry', ha='center', fontsize=10, color='#007e87')
    ax.text(0, -.71, 'Finite diagnostic boundaries; geometry sketch only\n'
            'Next: solve the magnetic field and gravity together', ha='center', va='top',
            fontsize=9, color='#475569', linespacing=1.5)
    ax.set(title='A. Static magnetic boundary inputs', xlim=(-1.55, 1.55), ylim=(-1.49, 1.48))
    ax.axis('off')

    ax = axes[0, 1]
    ax.plot(coordinate, profile, color='#007e87', lw=2.2)
    ax.scatter([-cutoff, 0, cutoff], [0, center, 1], s=27, color='#007e87', zorder=4)
    ax.annotate(f'Center amplitude = {center:.5f}', xy=(0, center), xytext=(-.90, .49),
                color='#475569', fontsize=9, arrowprops=dict(arrowstyle='->', color='#64748b', lw=1))
    ax.set_xticks([-cutoff, 0, cutoff], labels=['Left boundary', 'Center', 'Right boundary'])
    ax.get_xticklabels()[0].set_ha('left')
    ax.get_xticklabels()[-1].set_ha('right')
    ax.set(title='B. Computed response on fixed geometry',
           ylabel='Normalized magnetic mode a(x)', xlim=(-cutoff*1.07, cutoff*1.07), ylim=(-.035, 1.20))
    ax.text(.03, .95, 'Leading l = 1 pattern; ρ cutoff = 4\n'
            'Static response, not a material-transport fraction',
            transform=ax.transAxes, va='top', fontsize=8.5, color='#475569', linespacing=1.5)

    ax = axes[1, 0]
    ax.bar([0, 1, 2], delta_position, width=.52, color=['#326bb3', '#7951b0', '#007e87'])
    ax.axhline(0, color='#94a3b8', lw=1)
    span = max(float(np.ptp(np.r_[delta_position, 0.])), np.finfo(float).tiny)
    for i, value in enumerate(delta_position):
        offset = .035*span if value >= 0 else -.035*span
        ax.text(i, value+offset, f'{value:+.3g} m', ha='center',
                va='bottom' if value >= 0 else 'top', fontsize=9, color='#475569')
    ax.set_xticks([0, 1, 2], labels=['Δx', 'Δy', 'Δz'])
    ax.set(title='C. A checked three-dimensional launch change',
           ylabel='Arrival-position change at 5 km (m)',
           ylim=(min(delta_position.min(), 0)-.20*span, max(delta_position.max(), 0)+.42*span))
    ax.text(.98, .96, f'Total position shift: {np.linalg.norm(delta_position):.3f} m\n'
            '1 mm radial + vertical offsets; 10 μrad in both directions',
            transform=ax.transAxes, ha='right', va='top', fontsize=8.3, color='#475569', linespacing=1.5)

    ax = axes[1, 1]
    ax.bar(successful, shifts, width=.52, color=['#326bb3', '#326bb3', '#7951b0'])
    for index, value in zip(successful, shifts):
        ax.text(index, value*1.25, f'{value:,.2f} m', ha='center', fontsize=9,
                color='#326bb3' if index < 2 else '#7951b0')
    for index in unsuccessful:
        stopping_us = timing_cases[index]['endpoint']['lab_seconds']*1e6
        assert np.isclose(stopping_us, 30, rtol=1e-12)
        ax.text(index, .51, f'No 5 km arrival\nby {stopping_us:g} μs',
                transform=ax.get_xaxis_transform(), ha='center', va='center',
                fontsize=8.7, color='#92400e', linespacing=1.5,
                bbox=dict(facecolor='#fff7ed', edgecolor='none', pad=5))
    ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks(range(4), labels=[f'{value*1e9:+g}' for value in timing_offsets])
    ax.set(title='D. Selected changes of launch phase',
           xlabel='Phase-equivalent timing offset (ns)',
           ylabel='Position shift at 5 km\n(m; log scale)',
           xlim=(-.6, 3.6), ylim=(shifts.min()/4, shifts.max()*15))
    ax.text(.5, .96, 'Sampled outcomes; no exact timing threshold established',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569')
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)
    fig.text(.055, .13, 'Magnetic probe: unique finite static response on prescribed geometry. Matter passage, the coupled quantum source and full stability remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .085, 'Forward controls: sampled sensitivity only. Arrival changes are physical model responses, not numerical errors or verified receiver tolerances.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     13 → 14 Sep: time dilation is possible; a complete one-day journey remains untested.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'boundary_launch_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def coupled_receiver_progress():
    """Render the coupled scalar and full-trajectory static receiver diagnostics."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'coupled_dipole' not in gravity or 'electrostatic_receiver' not in future:
        raise ValueError('Final coupled-scalar and receiver results are required; no figure written')
    coupled = gravity['coupled_dipole']
    static = coupled['static_profile']
    chi = np.asarray(static['chi'], dtype=float)
    unit_psi = np.asarray(static['unit_right_zero_left_psi'], dtype=float)
    amplitude = float(static['illustrative_amplitude'])
    cutoff = np.arctan(static['cutoff_rho'])
    assert chi.ndim == 1 and chi.shape == unit_psi.shape and len(chi) > 20
    assert np.all(np.isfinite(chi)) and np.all(np.isfinite(unit_psi))
    assert np.all(np.diff(chi) > 0) and np.all(np.diff(unit_psi) > 0)
    assert np.allclose(chi[[0, -1]], [-cutoff, cutoff], rtol=0, atol=1e-13)
    assert np.allclose(unit_psi[[0, -1]], [0, 1], rtol=0, atol=1e-12)
    center_index = int(np.argmin(np.abs(chi)))
    assert abs(chi[center_index]) < 1e-13
    assert np.isclose(unit_psi[center_index], static['center_psi_per_unit_right_boundary'], rtol=1e-12)
    assert np.isclose(amplitude, .001, rtol=1e-12)
    modes = coupled['clamped_scalar_modes']
    assert len(modes) == 4
    assert all(mode['frequency_in_tau_units'] > 0 and
               mode['positive_wave_energy_at_unit_displacement_phase'] > 0 and
               abs(mode['dirichlet_boundary_residual']) < 1e-9 for mode in modes)
    # North pole, first order: the magnetic invariant is gauge invariant;
    # the lapse uses the saved sphere-trace-zero gauge and common time.
    magnetic_ppm = -amplitude*unit_psi*1e6
    lapse_ppm = magnetic_ppm/2

    receiver = future['electrostatic_receiver']
    design = receiver['design']
    fine = receiver['trajectories']['fine']
    assert receiver['capture_demonstrated'] is False
    assert fine['verification']['all_receiver_fields_present_from_initial_time'] is True
    route = fine['plot']
    lab_us = np.asarray(route['lab_microseconds'], dtype=float)
    proper_us = np.asarray(route['proton_proper_microseconds'], dtype=float)
    axial = np.asarray(route['receiver_axial_metres'], dtype=float)
    gamma = np.asarray(route['gamma'], dtype=float)
    positions = np.column_stack([np.asarray(route[f'{axis}_metres'], dtype=float) for axis in 'xyz'])
    origin = np.asarray(design['center_metres'], dtype=float)
    axis_direction = np.asarray(design['axis_unit_direction'], dtype=float)
    frame = np.asarray(design['frame_columns'], dtype=float)
    assert positions.shape == (len(lab_us), 3) and origin.shape == axis_direction.shape == (3,)
    assert frame.shape == (3, 3) and np.allclose(frame.T@frame, np.eye(3), rtol=0, atol=1e-12)
    assert np.isclose(np.linalg.norm(axis_direction), 1, rtol=0, atol=1e-12)
    assert all(values.shape == lab_us.shape and np.all(np.isfinite(values))
               for values in (lab_us, proper_us, axial, gamma))
    assert np.all(np.isfinite(positions)) and np.all(np.diff(lab_us) > 0)
    assert np.all(gamma >= 1) and lab_us[-1] > proper_us[-1] > 0
    assert abs(lab_us[0]) < 1e-12 and abs(proper_us[0]) < 1e-12
    endpoint = fine['endpoint']
    plane = fine['receiver_plane_crossing']
    assert np.isclose(lab_us[-1], endpoint['lab_seconds']*1e6, rtol=1e-12)
    assert np.isclose(proper_us[-1], endpoint['proton_proper_seconds']*1e6, rtol=1e-12)
    assert np.isclose(gamma[-1], endpoint['gamma'], rtol=1e-12)
    assert np.allclose(positions[-1], endpoint['position_metres'], rtol=0, atol=1e-8)
    assert np.min(np.abs(lab_us-plane['lab_seconds']*1e6)) < 1e-10
    assert np.allclose((positions-origin)@axis_direction, axial, rtol=1e-10, atol=1e-8)
    alignment = frame.T@axis_direction
    axial_column = int(np.argmax(np.abs(alignment)))
    assert np.isclose(abs(alignment[axial_column]), 1, rtol=0, atol=1e-12)
    transverse_columns = [column for column in range(3) if column != axial_column]
    local_transverse = ((positions-origin)@frame)[:, transverse_columns]
    inner_radius = float(design['inner_major_radius_metres'])
    outer_radius = float(design['outer_major_radius_metres'])
    tube_radius = float(design['minor_radius_metres'])
    inner_charge = float(design['inner_total_charge_coulombs'])
    outer_charge = float(design['outer_total_charge_coulombs'])
    assert inner_charge > 0 and np.isclose(inner_charge, -outer_charge, rtol=1e-12)
    assert 0 < tube_radius < inner_radius < outer_radius
    near = np.abs(axial) <= 2*outer_radius
    energy_window = (axial >= -200) & (axial <= 200+1e-7)
    assert np.count_nonzero(near) >= 5 and np.count_nonzero(energy_window) >= 20
    assert np.all(np.diff(axial[energy_window]) > 0)
    minimum_index = int(np.argmin(gamma))
    assert near[minimum_index] and gamma[-1] > gamma[minimum_index]+1
    assert np.isclose(gamma[minimum_index], fine['minimum_gamma_state']['gamma'], rtol=1e-12)
    assert np.isclose(axial[-1], 200, rtol=0, atol=1e-7)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 9))
    grid = fig.add_gridspec(2, 2)
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.64, wspace=.30)
    fig.suptitle('Coupled throat response; a receiver that slows but does not capture',
                 fontsize=17.3, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate theoretical tests: leading classical gravity and prescribed static electric charges',
             ha='center', fontsize=10.5, color='#475569')

    ax = fig.add_subplot(grid[0, 0])
    ax.fill_between([-1.08, 1.08], [-.34, -.34], [.34, .34], color='#007e87', alpha=.10)
    ax.plot([-1.08, 1.08], [.34, .34], [-1.08, 1.08], [-.34, -.34], color='#007e87', lw=1.7)
    for side in (-1, 1):
        ax.plot([side, side], [-.60, .60], ls='--', color='#64748b', lw=1.5)
        ax.text(side, .72, 'Fixed cut', ha='center', fontsize=9, color='#64748b')
    ax.text(-.47, 0, 'Magnetic\nfield', ha='center', va='center', fontsize=11,
            color='#007e87', linespacing=1.5)
    ax.text(.47, 0, 'Spacetime\nmetric', ha='center', va='center', fontsize=11,
            color='#7951b0', linespacing=1.5)
    ax.annotate('', xy=(.17, 0), xytext=(-.17, 0),
                arrowprops=dict(arrowstyle='<->', color='#475569', lw=1.7))
    ax.text(0, -.54, 'One physical dipole scalar Ψ couples both', ha='center',
            va='top', fontsize=9, color='#475569')
    ax.text(0, -1.02, 'No growing finite-energy modes for clamped cuts\n'
            'in this linear scalar problem', ha='center', va='center', fontsize=9,
            color='#007e87', linespacing=1.5,
            bbox=dict(facecolor='#eff9f8', edgecolor='none', pad=6))
    ax.text(0, 1.02, 'Assumed leading throat; geometry sketch only', ha='center',
            fontsize=8.7, color='#475569')
    ax.set(title='A. Field and geometry are now coupled', xlim=(-1.5, 1.5), ylim=(-1.38, 1.24))
    ax.axis('off')

    ax = fig.add_subplot(grid[0, 1])
    ax.plot(chi, magnetic_ppm, color='#007e87', lw=2.2, label='Magnetic invariant F²')
    ax.plot(chi, lapse_ppm, color='#7951b0', lw=2.2, ls='--', label='Lapse N: fixed coordinates')
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set_xticks([-cutoff, 0, cutoff], labels=['Left cut', 'Center', 'Right cut'])
    ax.get_xticklabels()[0].set_ha('left'); ax.get_xticklabels()[-1].set_ha('right')
    ax.set(title='B. An illustrative static dipole response', ylabel='Relative change at north pole (ppm)',
           xlim=(-1.04*cutoff, 1.04*cutoff), ylim=(-amplitude*1e6*1.14, amplitude*1e6*.15))
    ax.legend(frameon=False, loc='lower left', fontsize=8.5)
    ax.text(.03, .80, f'Boundary Ψ: 0 → {amplitude:g}; ρ cutoff = {static["cutoff_rho"]:g}\n'
            'First order; sphere-trace-zero gauge and common time', transform=ax.transAxes,
            va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)

    ax = fig.add_subplot(grid[1, 0], projection='3d')
    phi, theta = np.meshgrid(np.linspace(0, 2*np.pi, 89), np.linspace(0, 2*np.pi, 17))
    for radius, color in ((inner_radius, '#b7791f'), (outer_radius, '#326bb3')):
        ax.plot_surface(tube_radius*np.sin(theta), (radius+tube_radius*np.cos(theta))*np.cos(phi),
                        (radius+tube_radius*np.cos(theta))*np.sin(phi), color=color,
                        alpha=.85, linewidth=0, shade=True, antialiased=True)
        angle = np.linspace(0, 2*np.pi, 181)
        ax.plot(np.zeros_like(angle), radius*np.cos(angle), radius*np.sin(angle), color=color, lw=1.7)
    ax.plot(axial[near], local_transverse[near, 0], local_transverse[near, 1],
            color='#475569', lw=1.9, zorder=5)
    center_path_index = int(np.argmin(np.abs(axial)))
    ax.scatter([axial[center_path_index]], [local_transverse[center_path_index, 0]],
               [local_transverse[center_path_index, 1]], color='#475569', s=14, depthshade=False)
    ax.text2D(.01, .96, f'Inner: {inner_radius:g} m, +{inner_charge:.2f} C\n'
              f'Outer: {outer_radius:g} m, {outer_charge:.2f} C; tube radius {tube_radius:g} m', transform=ax.transAxes,
              fontsize=8.5, color='#475569', va='top', linespacing=1.4)
    ax.text2D(.5, -.045, 'Prescribed charged surfaces; field present from t = 0\n'
              'Charge response of a metal electrode is not solved', transform=ax.transAxes,
              ha='center', va='top', fontsize=8.3, color='#475569', linespacing=1.5)
    ax.set(title='C. Static receiver rings and the local proton path',
           xlabel='Axial s (m)', ylabel='Across (m)', zlabel='',
           xlim=(-2*outer_radius, 2*outer_radius),
           ylim=(-1.13*outer_radius, 1.13*outer_radius), zlim=(-1.13*outer_radius, 1.13*outer_radius))
    ax.set_xticks([-40, 0, 40]); ax.set_yticks([-20, 0, 20]); ax.set_zticks([-20, 0, 20])
    ax.tick_params(labelsize=8, pad=0)
    ax.xaxis.label.set_size(9); ax.yaxis.label.set_size(9)
    ax.xaxis.labelpad=0; ax.yaxis.labelpad=0
    ax.view_init(elev=19, azim=-65); ax.set_box_aspect((1.65, 1, 1))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill=False
        axis.pane.set_edgecolor('#e2e8f0')
        axis._axinfo['grid']['color']=(.80, .84, .88, .25)

    ax = fig.add_subplot(grid[1, 1])
    ax.plot(axial[energy_window], gamma[energy_window], color='#7951b0', lw=2.2)
    ax.axvline(0, ls=':', color='#94a3b8', lw=1)
    ax.axhline(1, ls=':', color='#94a3b8', lw=1)
    ax.scatter([axial[minimum_index], axial[-1]], [gamma[minimum_index], gamma[-1]],
               color='#7951b0', s=24, zorder=4)
    ax.annotate(f'At center: γ ≈ {gamma[minimum_index]:.3f}\n'
                f'Still moving at {plane["beta"]:.3f} c',
                xy=(axial[minimum_index], gamma[minimum_index]), xytext=(-180, 2.9),
                fontsize=8.8, color='#475569', linespacing=1.5,
                arrowprops=dict(arrowstyle='->', color='#64748b', lw=1))
    ax.annotate(f'At +{axial[-1]:g} m: γ ≈ {gamma[-1]:.4f}',
                xy=(axial[-1], gamma[-1]), xytext=(53, 6.6),
                fontsize=8.8, color='#475569',
                arrowprops=dict(arrowstyle='->', color='#64748b', lw=1))
    ax.set(title='D. The proton speeds up again after the center',
           xlabel='Along receiver axis, s (m)', ylabel='Energy / rest energy, γ',
           xlim=(-210, 210), ylim=(.5, max(gamma[energy_window])*1.10))
    ax.set_yticks([1, 4, 7, 10])
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)

    fig.text(.055, .125, 'Gravity: a leading classical scalar result. Actual mouths, the quantum source and nonlinear stability remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'Receiver: fixed electric charges act along the full trajectory. Material response, stopping and a clock reunion are not demonstrated.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a complete one-day journey remains untested.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'coupled_receiver_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def exterior_discharge_progress():
    """Render exterior boundary response and the final causal receiver-discharge test."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'coupled_exterior' not in gravity or 'switched_receiver' not in future:
        raise ValueError('Final exterior and switched-receiver results are required; no figure written')
    exterior = gravity['coupled_exterior']
    curve = exterior['plot']
    radius = np.asarray(curve['radius_over_M'], dtype=float)
    static_response = np.asarray(curve['static_impedance_M'], dtype=float)
    potential = np.asarray(curve['potential_U'], dtype=float)
    assert radius.ndim == 1 and len(radius) > 20
    assert radius.shape == static_response.shape == potential.shape
    assert all(np.all(np.isfinite(values)) for values in (radius, static_response, potential))
    assert np.all(np.diff(radius) > 0) and np.all(radius > 1)
    assert np.all(static_response > 0) and np.all(potential > 0)
    static_cases = exterior['static_cases']
    case_radius = np.asarray([row['cut_radius_over_M'] for row in static_cases], dtype=float)
    case_static = np.asarray([row['decaying_master_impedance_M'] for row in static_cases], dtype=float)
    assert np.all(np.diff(case_radius) > 0)
    assert all(np.isclose(row['decaying_master_impedance_M'], row['normalized_exterior_static_energy'],
                          rtol=1e-10, atol=1e-12) for row in static_cases)
    laplace = exterior['laplace_cases']
    assert len(laplace) == 3
    responses = []
    for record in laplace:
        cuts = np.asarray([row['cut_radius_over_M'] for row in record['cases']], dtype=float)
        values = np.asarray([row['impedance_M'] for row in record['cases']], dtype=float)
        excess = np.asarray([row['excess_over_static'] for row in record['cases']], dtype=float)
        assert np.allclose(cuts, case_radius, rtol=0, atol=1e-12)
        assert record['laplace_rate_times_M'] > 0 and np.all(values > case_static)
        assert np.allclose(values-case_static, excess, rtol=1e-10, atol=1e-12)
        assert all(row['combined_outer_boundary_and_step_refinement_difference'] < 1e-9
                   for row in record['cases'])
        responses.append(values)

    static_receiver = future['electrostatic_receiver']
    switched = future['switched_receiver']
    design = static_receiver['design']
    static_fine = static_receiver['trajectories']['fine']
    switched_fine = switched['trajectories']['fine']
    assert switched['capture_demonstrated'] is False
    assert switched['downstream_transport_achieved'] is False
    assert switched['retained_kinetic_energy_reduction_demonstrated'] is True
    assert switched_fine['verification']['solver_success'] is True
    assert switched['definition']['receiver'] == design
    static_plot = switched['plot']['static_reference']
    assert static_plot == static_fine['plot']
    switched_plot = switched_fine['plot']
    reference_us = static_fine['receiver_plane_crossing']['lab_seconds']*1e6
    static_time = np.asarray(static_plot['lab_microseconds'], dtype=float)-reference_us
    switched_time = np.asarray(switched_plot['lab_microseconds'], dtype=float)-reference_us
    static_gamma = np.asarray(static_plot['gamma'], dtype=float)
    switched_gamma = np.asarray(switched_plot['gamma'], dtype=float)
    static_axial = np.asarray(static_plot['receiver_axial_metres'], dtype=float)
    switched_axial = np.asarray(switched_plot['receiver_axial_metres'], dtype=float)
    for times, energy in ((static_time, static_gamma), (switched_time, switched_gamma)):
        assert times.ndim == 1 and times.shape == energy.shape and len(times) > 20
        assert np.all(np.isfinite(times)) and np.all(np.isfinite(energy))
        assert np.all(np.diff(times) > 0) and np.all(energy >= 1)
    assert static_axial.shape == static_time.shape and switched_axial.shape == switched_time.shape
    assert np.all(np.isfinite(static_axial)) and np.all(np.isfinite(switched_axial))
    assert np.isclose(static_gamma[-1], static_fine['endpoint']['gamma'], rtol=1e-12)
    assert np.isclose(switched_gamma[-1], switched_fine['endpoint']['gamma'], rtol=1e-12)
    assert np.isclose(static_axial[-1], 200, rtol=0, atol=1e-8)
    assert np.isclose(switched_axial[-1], -200, rtol=0, atol=1e-8)
    assert abs(switched_time[0]) < 1e-10
    initial = switched_fine['initial']
    reference_state = static_fine['receiver_plane_crossing']
    for key in ('position_metres', 'momentum_over_mc', 'lab_seconds', 'proton_proper_seconds', 'gamma'):
        assert np.allclose(initial[key], reference_state[key], rtol=0, atol=1e-12)
    assert switched['causal_reuse']['exact_saved_state_key'] == 'electrostatic_receiver.trajectories.fine.receiver_plane_crossing'
    assert abs(switched['causal_reuse']['onset_condition_residual_seconds']) < 1e-18
    assert len(switched_fine['axial_turns']) == 1
    turn = switched_fine['axial_turns'][0]
    clear = switched_fine['receiver_pulse_clear']
    turn_time = turn['lab_seconds']*1e6-reference_us
    clear_time = clear['lab_seconds']*1e6-reference_us
    turn_axial = turn['receiver_coordinates_metres'][2]
    assert 0 < turn_time < clear_time < switched_time[-1] and turn_axial > 0
    assert abs(turn['receiver_axial_momentum_over_mc']) < 1e-10
    for event in (turn, clear):
        event_index = int(np.argmin(abs(switched_time-(event['lab_seconds']*1e6-reference_us))))
        assert abs(switched_time[event_index]-(event['lab_seconds']*1e6-reference_us)) < 1e-10
        assert np.isclose(switched_gamma[event_index], event['gamma'], rtol=1e-11)
        assert np.isclose(switched_axial[event_index], event['receiver_coordinates_metres'][2], rtol=0, atol=1e-8)
    after_clear = switched_time >= clear_time-1e-10
    assert np.allclose(switched_gamma[after_clear], clear['gamma'], rtol=1e-11, atol=1e-11)
    assert clear['receiver_electric_field_volts_per_metre'] == 0
    assert clear['receiver_magnetic_field_tesla'] == 0
    inner_radius = float(design['inner_major_radius_metres'])
    outer_radius = float(design['outer_major_radius_metres'])
    tube_radius = float(design['minor_radius_metres'])
    inner_charge = float(design['inner_total_charge_coulombs'])
    outer_charge = float(design['outer_total_charge_coulombs'])
    assert 0 < tube_radius < inner_radius < outer_radius and inner_charge > 0
    assert np.isclose(outer_charge, -inner_charge, rtol=1e-12)
    source = switched['plot']['source']
    source_time = np.asarray(source['lab_microseconds'], dtype=float)
    charge = np.asarray(source['remaining_inner_charge_coulombs'], dtype=float)
    transfer_current = np.asarray(source['outward_transfer_current_amperes'], dtype=float)
    assert source_time.shape == charge.shape == transfer_current.shape
    assert np.all(np.diff(source_time) > 0) and np.all(np.isfinite(transfer_current))
    assert np.all(np.diff(charge) <= 1e-12) and np.all(transfer_current >= 0)
    assert np.isclose(charge[0], inner_charge, rtol=1e-12) and charge[-1] == 0
    peak_current = switched['source_definition']['peak_outward_transfer_current_amperes']
    assert np.isclose(np.max(transfer_current), peak_current, rtol=1e-12)
    duration_ns = switched['definition']['duration_seconds']*1e9
    bound = switched['source_energy_lower_bound']
    initial_energy_bound_GJ = bound['initial_field_energy_outside_sphere_lower_bound_joules']/1e9
    enclosing_radius = bound['enclosing_radius_metres']
    assert enclosing_radius > outer_radius+tube_radius and initial_energy_bound_GJ > 0
    final_kinetic_MeV = switched_fine['endpoint']['kinetic_energy_eV']/1e6
    assert np.isclose(final_kinetic_MeV, switched['energy_accounting']['final_kinetic_energy_eV']/1e6,
                      rtol=1e-12)
    static_window = static_time >= -.20
    switched_window = switched_time >= -.20
    assert np.count_nonzero(static_window) >= 20 and np.count_nonzero(switched_window) >= 20
    # The saved causal-reuse check establishes the identical earlier gamma=10 launch.
    # Show those actual saved samples before the new continuation, without retuning.
    common_prefix = static_window & (static_time < 0)
    display_time = np.r_[static_time[common_prefix], switched_time]
    display_gamma = np.r_[static_gamma[common_prefix], switched_gamma]
    display_axial = np.r_[static_axial[common_prefix], switched_axial]

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 9))
    grid = fig.add_gridspec(2, 2)
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Exterior response tested; receiver discharge reflects the proton',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate gravity and receiver tests • prescribed sources • no built medium or capture',
             ha='center', fontsize=10.5, color='#475569')

    ax = fig.add_subplot(grid[0, 0])
    ax.fill_between([-1.35, -.28], [-.32, -.32], [.32, .32], color='#007e87', alpha=.13)
    ax.plot([-1.35, -.28], [.32, .32], [-1.35, -.28], [-.32, -.32], color='#007e87', lw=1.8)
    ax.fill_between([.28, 1.35], [-.32, -.32], [.32, .32], color='#326bb3', alpha=.10)
    ax.plot([.28, 1.35], [.32, .32], [.28, 1.35], [-.32, -.32], color='#326bb3', lw=1.8)
    for cut in (-.28, .28):
        ax.plot([cut, cut], [-.62, .62], ls='--', color='#64748b', lw=1.5)
    ax.text(-.815, 0, 'Interior\nmodel', ha='center', va='center', fontsize=11,
            color='#007e87', linespacing=1.5)
    ax.text(.815, 0, 'Exterior\nmodel', ha='center', va='center', fontsize=11,
            color='#326bb3', linespacing=1.5)
    ax.text(0, 0, '?', ha='center', va='center', fontsize=22, color='#b7791f')
    ax.text(0, .88, 'Field + metric response in each model', ha='center', fontsize=9, color='#475569')
    ax.text(0, -.69, 'The physical mouth junction is not solved', ha='center', va='top',
            fontsize=9.5, color='#92400e')
    ax.text(0, -1.13, 'Common fields, geometry, clocks and quantum response\n'
            'must be matched; this is a domain sketch', ha='center', va='center',
            fontsize=8.5, color='#475569', linespacing=1.5)
    ax.set(title='A. A new exterior input; no joined medium yet',
           xlim=(-1.53, 1.53), ylim=(-1.49, 1.27))
    ax.axis('off')

    ax = fig.add_subplot(grid[0, 1])
    ax.plot(radius-1, static_response, color='#007e87', lw=2.2, label='Static: exact')
    for record, values, color in zip(laplace, responses, ['#326bb3', '#7951b0', '#b7791f']):
        ax.plot(case_radius-1, values, marker='o', ms=4, lw=1.3, color=color,
                label=f'Rate sM = {record["laplace_rate_times_M"]:g}')
    ax.set_xscale('log'); ax.minorticks_off()
    ax.set_xticks([.01, .1, 1, 9], labels=['1.01', '1.1', '2', '10'])
    ax.set(title='B. Response changes with the sampled test rate',
           xlabel='Exterior cut radius r/M (spacing expanded near 1)',
           ylabel='Exterior boundary response, KM', xlim=(.009, 10),
           ylim=(0, max(float(np.max(values)) for values in responses)*1.30))
    ax.legend(frameon=False, ncol=2, loc='upper right', fontsize=8.3)
    ax.text(.5, -.31, 'Positive real Laplace test rates sM\nNot a complete causal boundary response',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16); ax.set_axisbelow(True)

    ax = fig.add_subplot(grid[1, 0])
    angle = np.linspace(0, 2*np.pi, 241)
    ax.plot(inner_radius*np.cos(angle), inner_radius*np.sin(angle), color='#b7791f', lw=2.0)
    ax.plot(outer_radius*np.cos(angle), outer_radius*np.sin(angle), color='#326bb3', lw=2.0)
    for phi in np.linspace(0, 2*np.pi, 8, endpoint=False):
        start = inner_radius+.23*(outer_radius-inner_radius)
        end = inner_radius+.78*(outer_radius-inner_radius)
        ax.annotate('', xy=(end*np.cos(phi), end*np.sin(phi)),
                    xytext=(start*np.cos(phi), start*np.sin(phi)),
                    arrowprops=dict(arrowstyle='->', color='#64748b', lw=1.1))
    ax.scatter([0], [0], color='#475569', s=15)
    ax.text(0, -4, 'Receiver axis', ha='center', va='top', fontsize=8.5, color='#475569')
    ax.text(0, 1.51*outer_radius, f'Initially +{inner_charge:.2f} C inside and {outer_charge:.2f} C outside\n'
            f'Top view: radii {inner_radius:g}/{outer_radius:g} m; tube radius {tube_radius:g} m',
            ha='center', va='top', fontsize=8.4, color='#475569', linespacing=1.5)
    ax.text(0, -1.28*outer_radius, f'Source duration: {duration_ns:g} ns; peak current {peak_current/1e6:.1f} MA\n'
            'Edge current density diverges; no conductor design', ha='center', va='top',
            fontsize=8.4, color='#475569', linespacing=1.5)
    ax.set(title='C. Discharge has an explicit current path',
           xlim=(-1.65*outer_radius, 1.65*outer_radius),
           ylim=(-1.68*outer_radius, 1.62*outer_radius))
    ax.set_aspect('equal', adjustable='box'); ax.axis('off')

    detail_grid = grid[1, 1].subgridspec(2, 1, height_ratios=[1, 1.15], hspace=.18)
    energy_ax = fig.add_subplot(detail_grid[0, 0])
    position_ax = fig.add_subplot(detail_grid[1, 0], sharex=energy_ax)
    energy_ax.plot(static_time[static_window], static_gamma[static_window], color='#94a3b8',
                   ls='--', lw=1.9, label='Static')
    energy_ax.plot(display_time, display_gamma, color='#7951b0', lw=2.0, label='Discharge')
    energy_ax.scatter([static_time[-1]], [static_gamma[-1]], color='#94a3b8', s=18, zorder=4)
    energy_ax.scatter([switched_time[-1]], [switched_gamma[-1]], color='#7951b0', s=18, zorder=4)
    maximum_time = max(static_time[-1], switched_time[-1])
    maximum_gamma = max(np.max(static_gamma[static_window]), np.max(switched_gamma[switched_window]))
    energy_ax.set(title='D. Lower energy remains; the path reverses',
                  ylabel='Energy / rest\nenergy, γ',
                  xlim=(-.20, maximum_time+.06*(maximum_time+.20)), ylim=(.5, maximum_gamma*1.34))
    energy_ax.set_yticks([1, 5, 10])
    energy_ax.legend(frameon=False, loc='upper right', ncol=2, fontsize=8.1)
    energy_ax.text(.40, .36, f'Receiver fields clear at +{clear_time*1000:.0f} ns: {final_kinetic_MeV:.1f} MeV\n'
                   f'γ = {switched_gamma[-1]:.3f}; still moving', transform=energy_ax.transAxes,
                   fontsize=8.0, color='#7951b0', va='bottom', linespacing=1.4)
    energy_ax.tick_params(axis='x', labelbottom=False)
    position_ax.plot(static_time[static_window], static_axial[static_window], color='#94a3b8', ls='--', lw=1.9)
    position_ax.plot(display_time, display_axial, color='#7951b0', lw=2.0)
    position_ax.scatter([turn_time, switched_time[-1]], [turn_axial, switched_axial[-1]],
                        color='#7951b0', s=18, zorder=4)
    position_ax.scatter([static_time[-1]], [static_axial[-1]], color='#94a3b8', s=18, zorder=4)
    position_ax.axhline(0, color='#94a3b8', lw=.8)
    position_ax.annotate(f'Turn at\n+{turn_axial:.2f} m', xy=(turn_time, turn_axial), xytext=(-.16, 95),
                         fontsize=8.0, color='#475569', linespacing=1.4,
                         arrowprops=dict(arrowstyle='->', color='#64748b', lw=.9))
    position_ax.text(.22, 170, 'Static: +200 m', fontsize=8.1, color='#64748b')
    position_ax.text(.43, -58, 'Discharge: upstream −200 m', fontsize=8.1, color='#7951b0')
    position_ax.set(xlabel='Time after the common center state (μs)', ylabel='Position s (m)', ylim=(-245, 240))
    position_ax.set_yticks([-200, 0, 200])
    for ax in (energy_ax, position_ax):
        ax.axvline(0, ls=':', color='#94a3b8', lw=.8)
        ax.axvline(clear_time, ls=':', color='#c4b5d8', lw=.8)
        ax.tick_params(labelsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)

    fig.text(.055, .125, 'Exterior: static positivity and sampled Laplace response. A physical mouth junction, quantum response and full stability remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, f'Initial receiver field: ≥ {initial_energy_bound_GJ:.3f} GJ outside {enclosing_radius:g} m. '
             'This is a field-energy bound, not dissipation or driver demand.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a complete one-day journey remains untested.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'exterior_discharge_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def finite_boundary_progress():
    """Render finite-radius boundary coefficients and the regularized receiver test."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'coupled_boundary' not in gravity or 'regularized_receiver' not in future:
        raise ValueError('Final finite-boundary and regularized-receiver results are required; no figure written')
    boundary = gravity['coupled_boundary']
    regularized = future['regularized_receiver']
    if 'incoming' not in regularized or 'discharge' not in regularized:
        raise ValueError('The regularized receiver aggregation is incomplete; no figure written')
    assert regularized['current_density_divergence_removed'] is True
    assert regularized['downstream_transport_achieved'] is False
    assert regularized['capture_demonstrated'] is False
    direct = boundary['direct_four_dimensional_check']
    assert len(direct['einstein_components']) == 10 and set(direct['einstein_components']) == {'0'}
    assert len(direct['maxwell_components']) == 4 and set(direct['maxwell_components']) == {'0'}
    assert direct['moncrief_reextraction_residual'] == '0'
    assert set(boundary['residual_gauge']['all_six_boundary_residuals']) == {'0'}
    cases = boundary['static_boundary_cases']
    cut_radius = np.asarray([row['cut_radius_over_M'] for row in cases], dtype=float)
    lapse = np.asarray([row['lapse_log'] for row in cases], dtype=float)
    curvature = np.asarray([row['mixed_tt_extrinsic_curvature'] for row in cases], dtype=float)
    assert len(cases) == 6 and np.all(np.diff(cut_radius) > 0) and np.all(cut_radius > 1)
    assert all(np.all(np.isfinite(values)) for values in (cut_radius, lapse, curvature))
    assert np.all(lapse > 0) and np.all(curvature < 0)
    assert all(np.isclose(row['cap_flux_area_contrast'], .25, rtol=0, atol=1e-12) for row in cases)

    design = regularized['design']
    previous = future['switched_receiver']
    previous_design = previous['definition']['receiver']
    assert design['poloidal_weighting'] == 'cosine_squared'
    for key in ('center_metres', 'frame_columns', 'axis_unit_direction',
                'inner_major_radius_metres', 'outer_major_radius_metres', 'minor_radius_metres',
                'inner_total_charge_coulombs', 'outer_total_charge_coulombs'):
        assert np.allclose(design[key], previous_design[key], rtol=0, atol=1e-12)
    assert design['pulse_parameters'] == previous_design['pulse_parameters']
    new_incoming = regularized['incoming']['trajectories']['fine']
    new_discharge = regularized['discharge']['trajectories']['fine']
    old_discharge = previous['trajectories']['fine']
    assert new_discharge['verification']['solver_success'] is True
    old_plot = old_discharge['plot']
    new_plot = new_discharge['plot']
    reference_us = old_discharge['initial']['lab_seconds']*1e6
    old_time = np.asarray(old_plot['lab_microseconds'], dtype=float)-reference_us
    new_time = np.asarray(new_plot['lab_microseconds'], dtype=float)-reference_us
    old_axial = np.asarray(old_plot['receiver_axial_metres'], dtype=float)
    new_axial = np.asarray(new_plot['receiver_axial_metres'], dtype=float)
    for times, position, record in ((old_time, old_axial, old_discharge), (new_time, new_axial, new_discharge)):
        assert times.ndim == 1 and times.shape == position.shape and len(times) > 20
        assert np.all(np.isfinite(times)) and np.all(np.isfinite(position)) and np.all(np.diff(times) > 0)
        assert np.isclose(position[-1], record['endpoint']['receiver_coordinates_metres'][2], rtol=0, atol=1e-8)
        assert np.isclose(times[-1]+reference_us, record['endpoint']['lab_seconds']*1e6, rtol=1e-12)
    assert np.isclose(new_incoming['initial']['gamma'], 10, rtol=0, atol=1e-12)
    new_center = new_incoming['receiver_plane_crossing']
    for key in ('position_metres', 'momentum_over_mc', 'lab_seconds', 'proton_proper_seconds', 'gamma'):
        assert np.allclose(new_discharge['initial'][key], new_center[key], rtol=0, atol=1e-12)
    assert abs(regularized['causal_connection']['onset_residual_seconds']) < 1e-18
    for record in (old_discharge, new_discharge):
        assert np.isclose(record['endpoint']['receiver_coordinates_metres'][2], -200, rtol=0, atol=1e-8)
        assert record['endpoint']['receiver_axial_momentum_over_mc'] < 0
        assert len(record['axial_turns']) == 1 and record['axial_turns'][0]['receiver_coordinates_metres'][2] > 0
        assert np.isclose(record['endpoint']['gamma'], record['receiver_pulse_clear']['gamma'], rtol=1e-10)

    # Total radial current through a cylinder in the central overlap is fixed.
    # With xi=z/b, each displayed profile integrates to one; they are not fields.
    profiles = regularized['plot']['current_flux_profiles']
    xi = np.asarray(profiles['z_over_b'], dtype=float)
    finite_profile = np.asarray(profiles['regularized_normalized_flux_density'], dtype=float)
    assert xi.shape == finite_profile.shape and len(xi) > 20 and np.all(np.diff(xi) > 0)
    assert np.allclose(xi[[0, -1]], [-1, 1], rtol=0, atol=1e-12)
    assert np.allclose(finite_profile, 2*np.sqrt(np.maximum(0, 1-xi**2))/np.pi, rtol=1e-12, atol=1e-12)
    old_profile = np.full_like(xi, np.nan)
    interior = np.abs(xi) < 1
    old_profile[interior] = 1/(np.pi*np.sqrt(1-xi[interior]**2))
    theta = np.linspace(0, 2*np.pi, 1024, endpoint=False)
    assert np.isclose(np.mean(2*np.cos(theta)**2), 1, rtol=0, atol=1e-14)
    assert np.all(np.isfinite(finite_profile)) and np.all(finite_profile >= 0)
    assert np.all(finite_profile[[0, -1]] == 0)
    assert np.all(np.isfinite(old_profile[interior]))
    assert np.isclose(finite_profile[len(xi)//2], 2*old_profile[len(xi)//2], rtol=1e-12)
    source_audit = regularized['source_audit']
    assert abs(source_audit['normalized_charge_error']) < 1e-12
    assert np.isclose(source_audit['radial_flux_over_negative_qdot'], 1, rtol=1e-12)
    assert source_audit['current_squared_geometry_per_metre'] > 0
    assert np.isclose(source_audit['current_squared_geometry_per_metre'],
                      source_audit['independent_volume_quadrature_geometry_per_metre'], rtol=1e-9)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.20, hspace=.65, wspace=.30)
    fig.suptitle('Finite boundary data and a receiver current without edge divergence',
                 fontsize=17.5, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Static exterior reconstruction and prescribed-source tests • no joined quantum throat or built medium',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0, 0]
    ax.plot([.42, .42], [-.92, .92], ls='--', color='#64748b', lw=1.7)
    for height, text_label, color in ((.65, 'Lapse + induced metric', '#7951b0'),
                                      (0, 'Normal metric derivative', '#007e87'),
                                      (-.65, 'Magnetic flux + normal field', '#326bb3')):
        ax.text(-.55, height, text_label, ha='center', va='center', fontsize=10, color=color,
                bbox=dict(facecolor='#f6f8fa', edgecolor='none', pad=5))
        ax.annotate('', xy=(.38, height), xytext=(.12, height),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.2))
    ax.text(.93, .35, 'Quantum throat\nmodel', ha='center', va='center', fontsize=10,
            color='#475569', linespacing=1.5)
    ax.text(.93, -.42, 'Join\nunresolved', ha='center', va='center', fontsize=10,
            color='#92400e', linespacing=1.5)
    ax.text(0, -1.18, 'Boundary data, not a created tunnel or material', ha='center', va='top',
            fontsize=8.7, color='#475569')
    ax.set(title='A. The finite exterior boundary is reconstructed',
           xlim=(-1.65, 1.52), ylim=(-1.45, 1.24))
    ax.axis('off')

    ax = axes[0, 1]
    ax.plot(cut_radius-1, lapse, 'o-', color='#7951b0', ms=4.5, lw=1.8, label='Lapse change δ log N')
    ax.plot(cut_radius-1, curvature, 'o-', color='#007e87', ms=4.5, lw=1.8,
            label='Normal derivative M δκₜ')
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set_xscale('log'); ax.minorticks_off()
    ax.set_xticks([.01, .1, 1, 9], labels=['1.01', '1.1', '2', '10'])
    ax.set(title='B. Linear response at six finite cuts',
           xlabel='Cut radius r/M (spacing expanded near 1)', ylabel='Response per unit master H(c)',
           xlim=(.009, 10), ylim=(min(curvature)*1.13, max(lapse)*1.9))
    ax.legend(frameon=False, loc='lower right', fontsize=8.6)
    ax.text(.5, -.31, 'Each cut: H(c) = 1, I(c) = 0; no surface displacement or relabeling\n'
            'Normal points toward +r; these are linear coefficients', transform=ax.transAxes,
            ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16); ax.set_axisbelow(True)

    ax = axes[1, 0]
    ax.plot(xi, old_profile, color='#94a3b8', ls='--', lw=1.8, label='Earlier: divergent edges')
    ax.plot(xi, finite_profile, color='#007e87', lw=2.2, label='Reweighted: finite')
    for side in (-1, 1):
        ax.annotate('', xy=(side*.992, 1.63), xytext=(side*.965, 1.30),
                    arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=1.2))
    ax.scatter([-1, 1], [0, 0], color='#007e87', s=20, zorder=4)
    ax.set(title='C. The current profile is finite at its edges',
           xlabel='Height through current cross-section, z / b', ylabel='Unit-integral radial current profile',
           xlim=(-1.06, 1.06), ylim=(-.035, 1.78))
    ax.legend(frameon=False, loc='upper center', fontsize=8.4)
    ax.text(.5, .62, 'Equal total current\nCentral radial overlap of the annulus', transform=ax.transAxes,
            ha='center', va='center', fontsize=8.6, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)

    ax = axes[1, 1]
    ax.plot(old_time, old_axial, color='#94a3b8', ls='--', lw=2.0, label='Earlier profile')
    ax.plot(new_time, new_axial, color='#007e87', lw=2.0, label='Reweighted')
    ax.scatter([old_time[-1]], [old_axial[-1]], color='#94a3b8', s=20, zorder=4)
    ax.scatter([new_time[-1]], [new_axial[-1]], color='#007e87', s=20, zorder=4)
    ax.axhline(0, color='#94a3b8', lw=.8)
    maximum_time = max(old_time[-1], new_time[-1])
    low = min(old_axial.min(), new_axial.min()); high = max(old_axial.max(), new_axial.max())
    span = max(high-low, 1)
    ax.set(title='D. Both cases still return upstream',
           xlabel='Time after the earlier center crossing (μs)', ylabel='Receiver-axis position s (m)',
           xlim=(-.025, maximum_time*1.06), ylim=(low-.15*span, high+.45*span))
    ax.legend(frameon=False, loc='upper right', fontsize=8.6)
    old_energy = old_discharge['endpoint']['kinetic_energy_eV']/1e6
    new_energy = new_discharge['endpoint']['kinetic_energy_eV']/1e6
    ax.text(.05, .30, f'Endpoint kinetic energy\nEarlier: {old_energy:.3f} MeV\n'
            f'Reweighted: {new_energy:.3f} MeV', transform=ax.transAxes,
            ha='left', va='center', fontsize=8.6, color='#475569', linespacing=1.5)
    ax.text(.57, .47, 'Reflected; no capture', transform=ax.transAxes,
            fontsize=8.8, color='#92400e', ha='center')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.16, axis='y'); ax.set_axisbelow(True)

    fig.text(.055, .125, 'Boundary: direct linear Einstein–Maxwell and surface-gauge checks pass. The quantum-supported background and its physical junction remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'Receiver: bounded prescribed current density is a mathematical improvement; conductor response, source preparation and permanent capture remain unproved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a complete one-day journey remains untested.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'finite_boundary_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def taper_timing_progress():
    """Render imposed-taper pressure requirements and the verified receiver-duration scan."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'quantum_taper' not in gravity or 'receiver_duration_scan' not in future:
        raise ValueError('Final quantum-taper and receiver-duration results are required; no figure written')
    taper = gravity['quantum_taper']
    assert str(taper['direct_four_dimensional_divergence_residual']) == '0'
    cases = sorted(taper['cases'], key=lambda row: row['end_x']-row['start_x'])
    assert len(cases) == 3
    for case in cases:
        arrays = case['plot']
        coordinate = np.asarray(arrays['x'], dtype=float)
        window = np.asarray(arrays['window'], dtype=float)
        state_pressure = np.asarray(arrays['state_required_angular_pressure_hat'], dtype=float)
        whole_pressure = np.asarray(arrays['whole_required_angular_pressure_hat'], dtype=float)
        radial_pressure = np.asarray(arrays['original_radial_pressure_hat'], dtype=float)
        assert coordinate.ndim == 1 and len(coordinate) > 20
        assert all(values.shape == coordinate.shape and np.all(np.isfinite(values))
                   for values in (coordinate, window, state_pressure, whole_pressure, radial_pressure))
        assert np.all(np.diff(coordinate) > 0) and np.all(np.diff(window) <= 1e-12)
        assert np.allclose(coordinate[[0, -1]], [case['start_x'], case['end_x']], rtol=0, atol=1e-12)
        assert np.allclose(window[[0, -1]], [1, 0], rtol=0, atol=1e-12)
        u = (coordinate-case['start_x'])/(case['end_x']-case['start_x'])
        assert np.allclose(window, 1-10*u**3+15*u**4-6*u**5, rtol=1e-10, atol=1e-12)
        assert np.all(state_pressure >= -1e-12) and np.all(whole_pressure >= state_pressure-1e-12)
        assert np.all(radial_pressure < 0)
        assert case['any_monotone_whole_taper_max_pressure_ratio_lower_bound'] > 0
        assert case['whole_tensor_taper_peak']['required_angular_pressure_hat'] >= np.max(whole_pressure)-1e-9
    widest = cases[-1]
    assert np.isclose(widest['start_x'], .5) and np.isclose(widest['end_x'], 1.2)
    wide_plot = widest['plot']
    x = np.asarray(wide_plot['x'], dtype=float)
    state_pressure = np.asarray(wide_plot['state_required_angular_pressure_hat'], dtype=float)
    whole_pressure = np.asarray(wide_plot['whole_required_angular_pressure_hat'], dtype=float)
    peak = widest['whole_tensor_taper_peak']
    ratio_bound = widest['any_monotone_whole_taper_max_pressure_ratio_lower_bound']
    assert np.isclose(ratio_bound, 1/(2*widest['delta_log_areal_radius']), rtol=1e-12)

    scan = future['receiver_duration_scan']
    if 'selected' not in scan or 'screening' not in scan:
        raise ValueError('The selected receiver-duration refinement is incomplete; no figure written')
    assert scan['capture_demonstrated'] is False and scan['controlled_rest_demonstrated'] is False
    selected = scan['selected']['trajectory']
    records = [dict(record) for record in scan['screening']['records']]
    selected_duration = selected['definition']['duration_seconds']
    assert np.isclose(selected_duration, scan['selected']['duration_seconds'], rtol=1e-12, atol=0)
    selected_indices = [i for i, record in enumerate(records)
                        if np.isclose(record['duration_seconds'], selected_duration, rtol=1e-12, atol=0)]
    assert len(selected_indices) == 1
    records[selected_indices[0]]['endpoint'] = selected['endpoint']
    records[selected_indices[0]]['receiver_pulse_clear'] = selected['receiver_pulse_clear']
    records.sort(key=lambda record: record['duration_seconds'])
    durations_ns = np.asarray([record['duration_seconds']*1e9 for record in records], dtype=float)
    endpoint_energy_MeV = np.asarray([record['endpoint']['kinetic_energy_eV']/1e6 for record in records], dtype=float)
    endpoint_axial = np.asarray([record['endpoint']['receiver_coordinates_metres'][2] for record in records], dtype=float)
    assert len(records) == 8 and np.all(np.diff(durations_ns) > 0)
    assert np.all(endpoint_energy_MeV > 0) and np.all(np.isfinite(endpoint_energy_MeV))
    assert np.all(np.isclose(np.abs(endpoint_axial), 200, rtol=0, atol=1e-7))
    downstream = endpoint_axial > 0
    upstream = endpoint_axial < 0
    assert np.any(downstream) and np.any(upstream)
    selected_index = int(np.argmin(abs(durations_ns-selected_duration*1e9)))
    assert downstream[selected_index]
    assert endpoint_energy_MeV[selected_index] <= np.min(endpoint_energy_MeV[downstream])+1e-9
    baseline = future['regularized_receiver']
    assert selected['definition']['receiver'] == baseline['design']
    assert np.isclose(selected['definition']['onset_seconds'], baseline['discharge']['definition']['onset_seconds'],
                      rtol=0, atol=1e-18)
    baseline_center = baseline['incoming']['trajectories']['fine']['receiver_plane_crossing']
    assert selected['verification']['solver_success'] is True
    for key in ('position_metres', 'momentum_over_mc', 'lab_seconds', 'proton_proper_seconds', 'gamma'):
        assert np.allclose(selected['initial'][key], baseline_center[key], rtol=0, atol=1e-12)
    # The fixed charge waveform only changes duration; its peak current scales as 1/T.
    baseline_duration = baseline['discharge']['definition']['duration_seconds']
    peak_current_GA = scan['selected']['source_cost']['peak_transfer_current_amperes']/1e9
    assert np.isclose(peak_current_GA, baseline['source_audit']['peak_transfer_current_amperes']
                      * baseline_duration/selected_duration/1e9, rtol=1e-12, atol=0)
    path = selected['plot']
    axial = np.asarray(path['receiver_axial_metres'], dtype=float)
    transverse = np.asarray(path['receiver_transverse_metres'], dtype=float)
    lab_us = np.asarray(path['lab_microseconds'], dtype=float)
    proper_us = np.asarray(path['proton_proper_microseconds'], dtype=float)
    gamma = np.asarray(path['gamma'], dtype=float)
    assert axial.shape == transverse.shape == lab_us.shape == proper_us.shape == gamma.shape
    assert np.all(np.diff(lab_us) > 0) and np.all(transverse >= 0) and np.all(gamma >= 1)
    assert all(np.all(np.isfinite(values)) for values in (axial, transverse, lab_us, proper_us, gamma))
    endpoint = selected['endpoint']
    clear = selected['receiver_pulse_clear']
    assert np.isclose(axial[-1], 200, rtol=0, atol=1e-7)
    assert np.isclose(transverse[-1], endpoint['receiver_transverse_offset_metres'], rtol=1e-10)
    assert np.isclose(lab_us[-1], endpoint['lab_seconds']*1e6, rtol=1e-12)
    assert np.isclose(proper_us[-1], endpoint['proton_proper_seconds']*1e6, rtol=1e-12)
    assert endpoint['receiver_axial_momentum_over_mc'] > 0 and endpoint['beta'] > 0
    assert clear['receiver_electric_field_volts_per_metre'] == 0
    assert clear['receiver_magnetic_field_tesla'] == 0
    after_clear = lab_us >= clear['lab_seconds']*1e6-1e-10
    assert np.allclose(gamma[after_clear], endpoint['gamma'], rtol=1e-10, atol=1e-12)
    clear_after_center_ns = (clear['lab_seconds']-selected['initial']['lab_seconds'])*1e9
    selected_keV = endpoint['kinetic_energy_eV']/1e3
    turns = selected['axial_turns']
    assert len(turns) == 2
    assert [turn['axial_momentum_direction_change'] for turn in turns] == ['decreasing', 'increasing']
    turn_axial = np.asarray([turn['receiver_coordinates_metres'][2] for turn in turns], dtype=float)
    turn_transverse = np.asarray([turn['receiver_transverse_offset_metres'] for turn in turns], dtype=float)
    assert np.all(np.isfinite(turn_axial)) and np.all(np.isfinite(turn_transverse))
    assert selected['initial']['lab_seconds'] < turns[0]['lab_seconds'] < turns[1]['lab_seconds'] < endpoint['lab_seconds']

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.22, hspace=.65, wspace=.30)
    fig.suptitle('Source tapering requires pressure; receiver timing changes the outcome',
                 fontsize=17.3, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate theoretical tests • held-fixed throat geometry and prescribed receiver currents • no built time machine',
             ha='center', fontsize=10.3, color='#475569')

    ax = axes[0, 0]
    final_x = max(case['end_x'] for case in cases)
    for case, color in zip(cases, ['#94a3b8', '#7951b0', '#007e87']):
        xx = np.asarray(case['plot']['x'], dtype=float)
        ww = np.asarray(case['plot']['window'], dtype=float)
        ax.plot(np.r_[.45, xx, final_x], np.r_[1, ww, 0], color=color, lw=2,
                label=f'{case["start_x"]:g} → {case["end_x"]:g}')
    ax.set(title='A. Three imposed fading windows', xlabel='Conformal coordinate x', ylabel='Remaining source factor w',
           xlim=(.45, 1.23), ylim=(-.035, 1.16))
    ax.legend(title='Window interval', frameon=False, fontsize=8.5, title_fontsize=8.5, loc='upper right')
    ax.text(.5, -.31, 'The background geometry is kept fixed\nA consistent quantum source remains unproved',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[0, 1]
    ax.plot(x, state_pressure, color='#7951b0', ls='--', lw=2, label='State term only')
    ax.plot(x, whole_pressure, color='#007e87', lw=2.2, label='Whole stress tensor')
    ax.scatter([peak['x']], [peak['required_angular_pressure_hat']], color='#007e87', s=25, zorder=4)
    ax.annotate(f'Peak = {peak["required_angular_pressure_hat"]:.3f}',
                xy=(peak['x'], peak['required_angular_pressure_hat']), xytext=(.56, 1.07),
                fontsize=8.8, color='#007e87',
                arrowprops=dict(arrowstyle='->', color='#007e87', lw=1))
    ax.set(title='B. Required side pressure in the widest window',
           xlabel='Conformal coordinate x', ylabel='Normalized side pressure\n' + r'$8\pi G r_e^2 p_\perp/c^4$',
           xlim=(.48, 1.22), ylim=(-.03, peak['required_angular_pressure_hat']*1.40))
    ax.legend(frameon=False, loc='upper right', fontsize=8.5)
    ax.text(.5, -.31, f'Any monotone whole-stress taper on x = 0.5–1.2 needs\n'
            f'max(side pressure / |original radial pressure|) ≥ {ratio_bound/1e8:.2f} × 10⁸',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[1, 0]
    ax.scatter(durations_ns[downstream], endpoint_energy_MeV[downstream], marker='o', s=45,
               color='#007e87', label='Downstream +200 m', zorder=3)
    ax.scatter(durations_ns[upstream], endpoint_energy_MeV[upstream], marker='v', s=47,
               color='#7951b0', label='Upstream −200 m', zorder=3)
    ax.scatter([durations_ns[selected_index]], [endpoint_energy_MeV[selected_index]], marker='o', s=145,
               facecolors='none', edgecolors='#b7791f', lw=1.6, zorder=4)
    ax.annotate(f'{selected_duration*1e9:g} ns: {selected_keV:.1f} keV\nLowest tested downstream energy',
                xy=(durations_ns[selected_index], endpoint_energy_MeV[selected_index]),
                xytext=(6.5, .09), fontsize=8.3, color='#92400e', linespacing=1.5,
                arrowprops=dict(arrowstyle='->', color='#b7791f', lw=1))
    ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100], labels=['1', '2', '5', '10', '20', '50', '100'])
    ax.set(title='C. Eight tested durations; direction also changes',
           xlabel='Prescribed source duration (ns; log scale)', ylabel='Endpoint kinetic energy\n(MeV; log scale)',
           xlim=(.8, 130), ylim=(endpoint_energy_MeV.min()/4, endpoint_energy_MeV.max()*5))
    ax.legend(frameon=False, loc='upper left', fontsize=8.3)

    ax = axes[1, 1]
    ax.plot(axial, transverse, color='#007e87', lw=2.2)
    ax.axhline(0, color='#94a3b8', ls=':', lw=.9)
    ax.axvline(200, color='#94a3b8', ls='--', lw=.9)
    ax.scatter([axial[-1]], [transverse[-1]], color='#007e87', s=28, zorder=4)
    ax.scatter(turn_axial, turn_transverse, color='#b7791f', marker='x', s=34, zorder=4)
    ax.annotate('Two axial reversals near center\n(overlap at this scale)',
                xy=(float(np.mean(turn_axial)), float(np.mean(turn_transverse))),
                xytext=(.32, .58), textcoords='axes fraction', fontsize=8.0, color='#92400e', linespacing=1.5,
                arrowprops=dict(arrowstyle='->', color='#b7791f', lw=1))
    y_max = max(float(np.max(transverse)), 1e-6)
    ax.annotate(f'{transverse[-1]:.3f} m from the axis', xy=(axial[-1], transverse[-1]),
                xytext=(112, y_max*.28), fontsize=8.8, color='#007e87',
                arrowprops=dict(arrowstyle='->', color='#007e87', lw=1))
    ax.text(.03, .94, f'Source: {selected_duration*1e9:g} ns; {peak_current_GA:.2f} GA peak\n'
            f'Receiver fields clear +{clear_after_center_ns:.1f} ns after center\n'
            f'{selected_keV:.1f} keV at the +200 m plane; still moving', transform=ax.transAxes,
            ha='left', va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    ax.set(title='D. A plane crossing with a transverse offset',
           xlabel='Receiver-axis position s (m)', ylabel='Transverse offset (m)',
           xlim=(min(float(np.min(axial)), 0)-5, 210), ylim=(-.04*y_max, 1.29*y_max))
    ax.text(.5, -.25, 'Different axis scales; no point-target acceptance or capture test',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16); ax.set_axisbelow(True)
    fig.text(.055, .125, f'Tapers: α = {taper["alpha"]:.0e}; geometry held fixed. Pressure alone supplies no consistent quantum theory, gravitational response or stable join.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'Receiver: a classical test proton with prescribed currents. Capture, target acceptance, practical current generation and reunion remain unproved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a complete one-day journey remains untested.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'taper_timing_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def loop_source_radiation_progress():
    """Render the channel-source identity tests and full receiver-discharge radiation."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'magnetic_loop_source' not in gravity or 'receiver_radiation_spectrum' not in future:
        raise ValueError('Final magnetic-source and receiver-radiation blocks are required; no figure written')
    source = gravity['magnetic_loop_source']
    radiation_result = future['receiver_radiation_spectrum']
    if 'radiation' not in radiation_result or 'plot' not in radiation_result:
        raise ValueError('The receiver-radiation result is incomplete; no figure written')
    paths = gravity['magnetic_return_paths']['field_line_samples']
    labels = (.5, 1., 1.5)
    colors = ('#64748b', '#7951b0', '#007e87')
    path_rows = []
    for nu in labels:
        matches = [row for row in paths if np.isclose(row['nu'], nu, rtol=0, atol=1e-12)]
        assert len(matches) == 1
        row = matches[0]
        rho = np.asarray(row['rho_over_d'], dtype=float)
        z = np.asarray(row['z_over_d'], dtype=float)
        assert rho.ndim == 1 and rho.shape == z.shape and len(rho) > 20
        assert np.all(np.isfinite(rho)) and np.all(np.isfinite(z)) and np.all(rho >= 0)
        assert np.allclose(rho[[0, -1]], 0, rtol=0, atol=1e-12)
        assert np.allclose(z[[0, -1]], [-.5, .5], rtol=0, atol=1e-12)
        path_rows.append((rho, z))
    assert source['parameters']['flux_interval'] == [.5, 1.5]
    assert np.isclose(source['parameters']['fraction_of_total_flux'], .5)
    middle_profile = next(row for row in source['profile'] if np.isclose(row['nu'], 1.))
    middle = min(middle_profile['samples'], key=lambda row: abs(row['s']))
    assert abs(middle['s']) < 1e-12
    force = np.asarray(middle['force_times_d5_over_q'], dtype=float)
    assert force[0] > 0 and np.allclose(force[1:], 0, rtol=0, atol=1e-15)
    assert middle['azimuthal_current_times_d3'] < 0

    action_cases = source['shape_variation_checks']
    assert len(action_cases) == 9
    convergence_rows = []
    for nu in labels:
        matches = [row for row in action_cases if np.isclose(row['nu'], nu, rtol=0, atol=1e-12)]
        assert len(matches) == 3
        case_steps, case_errors = [], []
        for case in matches:
            predicted = case['action_derivative_predicted']
            assert predicted != 0 and np.isfinite(predicted)
            finite = sorted(case['finite_difference_derivatives'], key=lambda row: row['step'])
            steps = np.asarray([row['step'] for row in finite], dtype=float)
            errors = np.asarray([abs(row['full']-predicted)/abs(predicted) for row in finite], dtype=float)
            assert len(steps) == 3 and np.all(np.diff(steps) > 0)
            assert np.all(np.isfinite(errors)) and np.all(errors > 0)
            assert np.isclose(errors[0], case['relative_action_error'], rtol=1e-8, atol=1e-16)
            case_steps.append(steps)
            case_errors.append(errors)
        assert np.allclose(case_steps, case_steps[0], rtol=0, atol=1e-15)
        convergence_rows.append((case_steps[0], np.max(case_errors, axis=0)))
    ward = source['ward_checks']
    assert len(ward) == 9
    ward_errors = np.asarray([error for row in ward for error in row['relative_errors']], dtype=float)
    assert ward_errors.size == 18 and np.all(np.isfinite(ward_errors)) and np.all(ward_errors >= 0)
    ward_max = float(np.max(ward_errors))

    radiation = radiation_result['radiation']
    spectral = radiation_result['plot']
    frequency_x = np.asarray(spectral['frequency_x'], dtype=float)
    frequency_hz = np.asarray(spectral['frequency_hertz'], dtype=float)
    per_x = np.asarray(spectral['spectral_energy_joules_per_x'], dtype=float)
    per_hz = np.asarray(spectral['spectral_energy_joules_per_hertz'], dtype=float)
    cumulative_x = np.asarray(spectral['cumulative_frequency_x'], dtype=float)
    cumulative_energy = np.asarray(spectral['cumulative_energy_joules'], dtype=float)
    assert frequency_x.ndim == 1 and len(frequency_x) > 20
    assert frequency_x.shape == frequency_hz.shape == per_x.shape == per_hz.shape
    assert cumulative_x.shape == cumulative_energy.shape and len(cumulative_x) > 20
    assert all(np.all(np.isfinite(v)) for v in
               (frequency_x, frequency_hz, per_x, per_hz, cumulative_x, cumulative_energy))
    assert np.all(np.diff(frequency_x) > 0) and np.all(np.diff(cumulative_x) > 0)
    assert np.all(per_x >= 0) and np.all(per_hz >= 0)
    duration = future['receiver_duration_scan']['selected']['duration_seconds']
    assert np.isclose(duration, 3.5e-9, rtol=1e-12, atol=0)
    jacobian = 2*np.pi*duration
    assert np.allclose(frequency_hz*jacobian, frequency_x, rtol=1e-12, atol=1e-14)
    assert np.allclose(per_hz, per_x*jacobian, rtol=1e-12, atol=0)
    evaluated = radiation['evaluated_spectrum_joules']
    tail_bound = radiation['omitted_positive_tail_upper_bound_joules']
    cutoff_hz = radiation['cutoff_frequency_hertz']
    low_band = radiation['evaluated_low_band_joules']
    assert evaluated > 0 and 0 <= low_band <= evaluated and tail_bound > 0
    # Spectrum samples are interior quadrature nodes; cumulative values use bin edges.
    assert 0 < frequency_hz[0] < frequency_hz[-1] < cutoff_hz
    assert np.isclose(cumulative_x[-1]/jacobian, cutoff_hz, rtol=1e-12, atol=0)
    assert np.isclose(cumulative_energy[-1], evaluated, rtol=1e-10, atol=1e-6)
    assert np.all(cumulative_energy >= 0) and np.all(np.diff(cumulative_energy) >= -evaluated*1e-12)
    energy_units = ((1e12, 'TJ'), (1e9, 'GJ'), (1e6, 'MJ'), (1e3, 'kJ'), (1., 'J'))
    energy_scale, energy_unit = next((scale, unit) for scale, unit in energy_units if evaluated >= scale)

    def energy_text(value, upper_bound=False):
        scale, unit = next(((s, u) for s, u in energy_units if value >= s), (1., 'J'))
        scaled = value/scale
        if upper_bound:
            precision = 10.**(3-int(np.floor(np.log10(scaled))))
            scaled = np.ceil(scaled*precision)/precision
        return f'{scaled:.4g} {unit}'

    frequency_GHz = frequency_hz/1e9
    cumulative_GHz = cumulative_x/(jacobian*1e9)
    positive_frequency = frequency_GHz > 0
    positive_cumulative = cumulative_GHz > 0
    spectral_scaled = per_hz*1e9/energy_scale

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.22, hspace=.65, wspace=.30)
    fig.suptitle('Testing a magnetic source and the receiver’s radiation cost',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate theoretical calculations • exterior magnetic paths and prescribed discharge currents • no built time machine',
             ha='center', fontsize=10.1, color='#475569')

    ax = axes[0, 0]
    for nu, color, (rho, z) in zip(labels, colors, path_rows):
        ax.plot(rho, z, color=color, lw=2, label=f'ν = {nu:g}')
    ax.scatter([0, 0], [-.5, .5], color='#334155', s=30, zorder=4)
    ax.annotate('', xy=(middle['rho_over_d']+.38, middle['z_over_d']),
                xytext=(middle['rho_over_d'], middle['z_over_d']),
                arrowprops=dict(arrowstyle='->', color='#b7791f', lw=1.7))
    ax.text(middle['rho_over_d']+.10, .14, 'Force direction', color='#92400e', fontsize=8.5)
    ax.set(title='A. Exterior paths only; throat not solved',
           xlabel='Cylindrical radius ρ / d', ylabel='Axis position z / d',
           xlim=(-.18, 2.60), ylim=(-.92, .92))
    ax.set_aspect('equal', adjustable='box')
    ax.legend(frameon=False, loc='upper right', fontsize=8.5)
    ax.text(.5, -.31, 'Dots: assumed mouth locations. The arrow shows direction only.\n'
            'The derived current still needs a new magnetic field solution.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[0, 1]
    for nu, color, (steps, errors) in zip(labels, colors, convergence_rows):
        ax.plot(steps*1e4, errors, 'o-', color=color, lw=1.8, ms=4, label=f'ν = {nu:g}')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks([.5, 1, 2], labels=['0.5', '1', '2'])
    ax.set(title='B. Smaller test steps improve the agreement',
           xlabel='Deformation size / d (×10⁻⁴)', ylabel='Relative action-response error')
    ax.legend(frameon=False, loc='upper left', fontsize=8.5)
    ax.text(.5, -.31, 'Worst of three test metrics; one is deliberately amplified.\n'
            f'Independent sampled stress differentiation: max error {ward_max:.2g}.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[1, 0]
    ax.plot(frequency_GHz[positive_frequency],
            np.where(spectral_scaled[positive_frequency] > 0, spectral_scaled[positive_frequency], np.nan),
            color='#7951b0', lw=1.5)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
    ax.set(title='C. Radiation from the receiver source',
           xlabel='Frequency (GHz; log scale)',
           ylabel=f'Spectral energy ({energy_unit} / GHz; log scale)')
    ax.text(.5, -.25, f'{duration*1e9:g} ns charge fall; full discharge, all directions at infinity.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')

    ax = axes[1, 1]
    ax.plot(cumulative_GHz[positive_cumulative], cumulative_energy[positive_cumulative]/energy_scale,
            color='#007e87', lw=2.2)
    ax.scatter([cutoff_hz/1e9], [evaluated/energy_scale], color='#007e87', s=28, zorder=4)
    ax.set_xscale('log'); ax.minorticks_off()
    ax.set(title='D. Accumulated radiation up to the cutoff',
           xlabel='Upper frequency (GHz; log scale)', ylabel=f'Evaluated radiation ({energy_unit})',
           ylim=(0, 1.29*evaluated/energy_scale))
    ax.text(.03, .95, f'Evaluated energy ≈ {energy_text(evaluated)}\n'
            f'Omitted high-frequency tail ≤ {energy_text(tail_bound, upper_bound=True)}',
            transform=ax.transAxes, va='top', fontsize=9, color='#007e87', linespacing=1.6)
    ax.text(.5, -.25, 'The tail bound is separate from numerical integration error.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.16); ax.set_axisbelow(True)
    fig.text(.055, .125, 'Source: selected flux channels on prescribed geometry. Full field and gravity responses, a global quantum state and a stable throat remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'Radiation: the isolated receiver discharge. Initial full field energy, driver work and a material or circuit response remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; a complete journey and receiver capture remain unproved.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'loop_source_radiation_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def patch_rf_progress():
    """Render the conditional magnetic patch response and separate synchronous RF ramp."""
    from matplotlib.patches import Rectangle

    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'magnetic_patch_response' not in gravity or 'staged_rf_deceleration' not in future:
        raise ValueError('Final patch-response and RF-ramp blocks are required; no figure written')
    patch = gravity['magnetic_patch_response']
    rf = future['staged_rf_deceleration']
    if 'plot' not in rf:
        raise ValueError('The RF-ramp result is incomplete; no figure written')
    base = max(patch['grid_results'], key=lambda row: row['grid_count'])
    larger = max(patch['enlarged_grid_results'], key=lambda row: row['grid_count'])
    assert np.allclose(base['spacing'], larger['spacing'], rtol=1e-12, atol=1e-14)
    assert base['radial_interval'] == [.7, 1.1] and base['axial_interval'] == [-.2, .2]
    assert larger['radial_interval'] == [.6, 1.2] and larger['axial_interval'] == [-.3, .3]
    field_plot = patch['plot']
    rho = np.asarray(field_plot['radial'], dtype=float)
    z = np.asarray(field_plot['axial'], dtype=float)
    w = np.asarray(field_plot['w'], dtype=float)
    delta_b = np.asarray(field_plot['delta_b_over_q'], dtype=float)
    background_b = np.asarray(field_plot['background_b_over_q'], dtype=float)
    assert w.shape == (len(rho), len(z))
    assert delta_b.shape == background_b.shape == (len(rho), len(z), 2)
    assert all(np.all(np.isfinite(a)) for a in (rho, z, w, delta_b, background_b))
    assert np.all(np.diff(rho) > 0) and np.all(np.diff(z) > 0)
    assert np.allclose(w[[0, -1], :], 0, rtol=0, atol=1e-20)
    assert np.allclose(w[:, [0, -1]], 0, rtol=0, atol=1e-20)
    assert np.allclose(w, w[:, ::-1], rtol=1e-10, atol=1e-18)
    center_index = (len(rho)//2, len(z)//2)
    assert np.isclose(rho[center_index[0]], .9) and abs(z[center_index[1]]) < 1e-12
    assert np.allclose(delta_b[center_index], base['centre_delta_b_over_q'], rtol=1e-12, atol=1e-25)
    center_fields = np.asarray([base['centre_delta_b_over_q'][1], larger['centre_delta_b_over_q'][1]])
    assert np.all(center_fields < 0)
    boundary_ratio = abs(center_fields[1]/center_fields[0])
    assert boundary_ratio > 1
    path_rows = []
    for nu in (.5, 1., 1.5):
        row = next(row for row in gravity['magnetic_return_paths']['field_line_samples']
                   if np.isclose(row['nu'], nu, rtol=0, atol=1e-12))
        path_rho = np.asarray(row['rho_over_d'], dtype=float)
        path_z = np.asarray(row['z_over_d'], dtype=float)
        assert path_rho.shape == path_z.shape and np.all(np.isfinite(path_rho)) and np.all(np.isfinite(path_z))
        path_rows.append((nu, path_rho, path_z))

    ramp = rf['plot']
    lab = np.asarray(ramp['lab_seconds'], dtype=float)
    proper = np.asarray(ramp['proper_seconds'], dtype=float)
    energy = np.asarray(ramp['kinetic_energy_MeV'], dtype=float)
    gamma = np.asarray(ramp['gamma'], dtype=float)
    magnet = np.asarray(ramp['field_tesla'], dtype=float)
    frequency = np.asarray(ramp['rf_frequency_hz'], dtype=float)
    assert lab.ndim == 1 and len(lab) > 20
    assert all(a.shape == lab.shape and np.all(np.isfinite(a)) for a in (lab, proper, energy, gamma, magnet, frequency))
    assert np.all(np.diff(lab) > 0) and np.all(np.diff(proper) > 0)
    assert np.all(np.diff(proper) < np.diff(lab))
    assert np.all(energy > 0) and np.all(np.diff(energy) < 0) and np.all(gamma > 1)
    assert np.all(magnet > 0) and np.all(frequency > 0)
    assert np.isclose(gamma[0], 10, rtol=1e-12)
    elapsed_lab = lab-lab[0]
    elapsed_proper = proper-proper[0]
    time_scale, time_unit = (60., 'min') if elapsed_lab[-1] >= 180 else (1., 's')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.22, hspace=.65, wspace=.30)
    fig.suptitle('A conditional field correction and a separate RF slowdown',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Two idealized calculations • a fixed exterior source and an assumed synchronized orbit • no constructed medium',
             ha='center', fontsize=10.2, color='#475569')

    ax = axes[0, 0]
    for nu, path_rho, path_z in path_rows:
        ax.plot(path_rho, path_z, color='#7951b0' if nu == 1 else '#94a3b8', lw=1.7 if nu == 1 else 1.3)
    ax.scatter([0, 0], [-.5, .5], color='#334155', s=28, zorder=4)
    for row, color, linestyle, label in ((base, '#007e87', '-', 'Base patch'),
                                         (larger, '#b7791f', '--', 'Enlarged patch')):
        left, right = row['radial_interval']; low, high = row['axial_interval']
        ax.add_patch(Rectangle((left, low), right-left, high-low,
                               fill=False, edgecolor=color, ls=linestyle, lw=1.8, label=label, zorder=4))
    ax.set(title='A. Where the exterior calculation is made',
           xlabel='Cylindrical radius ρ / d', ylabel='Axis position z / d',
           xlim=(-.18, 2.60), ylim=(-.92, .92))
    ax.set_aspect('equal', adjustable='box')
    ax.legend(frameon=False, loc='upper right', fontsize=8.2)
    ax.text(.5, -.31, 'Curves: ideal exterior paths. Dots: assumed mouth locations.\n'
            'Rectangles are artificial calculation boundaries, not objects.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[0, 1]
    magnification = 1e14
    scaled_bz = delta_b[:, :, 1]*magnification
    color_limit = float(np.max(abs(scaled_bz)))
    mesh = ax.pcolormesh(rho, z, scaled_bz.T, cmap='RdBu_r', shading='auto',
                         vmin=-color_limit, vmax=color_limit, rasterized=True)
    ax.scatter([rho[center_index[0]]], [z[center_index[1]]], color='#334155', s=20, zorder=4)
    ax.set(title='B. Axial correction, magnified ×10¹⁴',
           xlabel='Cylindrical radius ρ / d', ylabel='Axis position z / d')
    ax.set_aspect('equal', adjustable='box')
    panel_box = ax.get_position(original=True)
    square_width = panel_box.height*fig.get_figheight()/fig.get_figwidth()
    heatmap_left = panel_box.x0+(panel_box.width-square_width-.078)/2
    ax.set_position([heatmap_left, panel_box.y0, square_width, panel_box.height])
    colorbar_ax = fig.add_axes([heatmap_left+square_width+.018, panel_box.y0, .010, panel_box.height])
    colorbar = fig.colorbar(mesh, cax=colorbar_ax)
    colorbar.set_label('10¹⁴ δBz d² / q')
    ax.text(.5, -.31, 'This shows the field correction, not the total field.\n'
            'Zero flux correction was imposed at every patch edge.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[1, 0]
    positions = np.arange(2)
    scaled_centers = center_fields*magnification
    ax.bar(positions, scaled_centers, width=.48, color=['#007e87', '#b7791f'])
    ax.axhline(0, color='#64748b', lw=.8)
    ax.set_xticks(positions, labels=['Base patch', 'Enlarged patch'])
    ax.set(title='C. The artificial boundary changes the answer',
           ylabel='10¹⁴ δBz d² / q at center',
           ylim=(float(np.min(scaled_centers))*1.32, .25))
    for position, value in zip(positions, scaled_centers):
        ax.text(position, value-.25, f'{value:.3f}', ha='center', va='top', fontsize=9, color='#334155')
    ax.text(.03, .08, f'{boundary_ratio:.2f}× larger center correction',
            transform=ax.transAxes, ha='left', va='bottom', fontsize=9, color='#92400e')
    ax.text(.5, -.25, 'Same grid spacing; different artificial boundaries. This is not a stability test.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')

    ax = axes[1, 1]
    ax.plot(elapsed_lab/time_scale, energy, color='#7951b0', lw=2.2)
    ax.scatter([elapsed_lab[-1]/time_scale], [energy[-1]], color='#7951b0', s=28, zorder=4)
    ax.set_yscale('log'); ax.minorticks_off()
    ax.set(title='D. A separate ideal RF energy ramp',
           xlabel=f'Elapsed laboratory time ({time_unit})', ylabel='Proton kinetic energy (MeV; log scale)',
           ylim=(energy[-1]/2, energy[0]*1.7))
    ax.text(.035, .12, f'Final kinetic energy: {energy[-1]:g} MeV\n'
            f'Laboratory time: {elapsed_lab[-1]:.2f} s\n'
            f'Proton proper time: {elapsed_proper[-1]:.2f} s',
            transform=ax.transAxes, va='bottom', fontsize=8.6, color='#475569', linespacing=1.5)
    ax.text(.5, -.25, 'Synchronization is assumed; capture and beam stability are not demonstrated.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .125, 'Magnetic response: current and geometry held fixed. Physical boundaries, a consistent global source, gravity response and a stable throat remain unsolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'RF ramp: a separate prescribed orbit, not the earlier receiver trajectory. Complete transport, source hardware and permanent capture remain unproved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; no complete human journey has been demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'patch_rf_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def sourced_mouth_rf_progress():
    """Render a sourced RN transfer and the accepted reduced RF-spread diagnostics."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'sourced_rn_transfer' not in gravity or 'rf_ramp_spread' not in future:
        raise ValueError('Final sourced-transfer and RF-spread results are required; no figure written')
    rn = gravity['sourced_rn_transfer']
    spread = future['rf_ramp_spread']
    if 'verification' not in spread or str(spread.get('status', '')).startswith('partial'):
        raise ValueError('The RF-spread aggregation is incomplete; no figure written')
    assert all(value == 0 for value in rn['verification']['direct_four_dimensional_Einstein_and_Maxwell_residuals'])
    assert rn['verification']['radial_source_Ward_identity'] is True
    assert rn['verification']['residual_coordinate_gauge_identity'] is True
    source = rn['source']
    amplitude = source['retained_amplitude_at_common_clock']
    assert amplitude > 0 and np.isclose(source['retained_fraction_of_flux'], .5)
    assert np.isclose(amplitude, source['amplitude_per_total_moment_at_common_clock']
                      * source['retained_dipole_moment'], rtol=1e-12, atol=0)
    radial = np.asarray(rn['plot']['radius_over_re'], dtype=float)
    source_h = np.asarray(rn['plot']['source_response_H_per_A'], dtype=float)
    retained_h = np.asarray(rn['plot']['retained_band_H_contribution'], dtype=float)
    assert radial.shape == source_h.shape == retained_h.shape and len(radial) > 20
    assert all(np.all(np.isfinite(a)) for a in (radial, source_h, retained_h))
    assert np.all(np.diff(radial) > 0) and np.all(source_h >= 0)
    assert np.allclose(radial[[0, -1]], [1.01, 4], rtol=0, atol=1e-12)
    assert np.allclose(retained_h, amplitude*source_h, rtol=1e-12, atol=0)
    assert source_h[-1] == 0
    inner_case = min(rn['transfer_cases'], key=lambda row: row['inner_radius_over_re'])
    assert np.isclose(source_h[0], inner_case['source_response_per_A'][0], rtol=1e-10, atol=0)

    off_crest = spread['fine']
    crest = spread['crest_fine']
    assert off_crest['crest_control'] is False and crest['crest_control'] is True
    assert len(off_crest['cases']) == 4 and len(crest['cases']) == 1
    assert off_crest['outcome'] == 'four_sample_envelopes_within_aperture'
    assert crest['outcome'] == 'guaranteed_outside_in_reduced_model'
    cp = crest['plot']
    crest_time = np.asarray(cp['lab_seconds'], dtype=float)
    crest_offset = np.asarray(cp['guiding_offset_metres'], dtype=float)
    crest_radial = np.asarray(cp['radial_betatron_envelope_metres'], dtype=float)
    assert crest_offset.shape == (1, len(crest_time)) and crest_radial.shape == crest_time.shape
    assert np.all(np.diff(crest_time) > 0) and np.all(crest_radial >= 0)
    assert all(np.all(np.isfinite(a)) for a in (crest_time, crest_offset, crest_radial))
    crest_lower = np.maximum(abs(crest_offset[0])-crest_radial, 0)*1e3
    crest_upper = (abs(crest_offset[0])+crest_radial)*1e3
    departure_time = crest['guaranteed_outside_times_seconds'][0]
    assert departure_time is not None and crest_time[0] < departure_time <= crest_time[-1]+1e-10
    departure_ms = (departure_time-crest_time[0])*1e3
    # Use adaptive-solver extrema, not the under-resolved 24-second display trace.
    corner_cases = sorted(off_crest['cases'], key=lambda row:
                          (row['initial_phase_error_radians'], row['initial_relative_momentum_error']))
    margins_mm = np.asarray([row['minimum_sampled_radial_envelope_margin_metres']*1e3
                             for row in corner_cases], dtype=float)
    phases = np.asarray([row['initial_phase_error_radians'] for row in corner_cases], dtype=float)
    momenta = np.asarray([row['initial_relative_momentum_error'] for row in corner_cases], dtype=float)
    assert np.all(np.isfinite(margins_mm)) and np.all(margins_mm > 0)
    assert np.allclose(abs(phases), abs(phases[0]), rtol=1e-12, atol=0) and phases[0] != 0
    assert np.allclose(abs(momenta), abs(momenta[0]), rtol=1e-12, atol=0) and momenta[0] != 0
    assert len(set(zip(np.sign(phases), np.sign(momenta)))) == 4
    corner_labels = [f'{"+" if phase > 0 else "−"}φ, {"+" if momentum > 0 else "−"}p'
                     for phase, momentum in zip(phases, momenta)]
    ramp_time = off_crest['final_lab_seconds']
    assert ramp_time > crest_time[-1]
    phase_degrees = abs(phases[0])*180/np.pi
    late_check = max(spread['verification']['lorentz_window_comparisons'],
                     key=lambda row: row['nominal_start_turn'])
    late_window = [row for row in spread['windows']
                   if row['fine'] and row['nominal_start_turn'] == late_check['nominal_start_turn']]
    assert len(late_window) == 1 and late_window[0]['turns'] == 16
    late_model_difference_mm = late_check['guiding_model_difference_metres']*1e3
    assert late_model_difference_mm > 0
    assert np.isclose(late_model_difference_mm,
                      late_window[0]['maximum_projected_guiding_difference_from_average_metres']*1e3,
                      rtol=1e-12, atol=0)
    amp_mantissa, amp_exponent = f'{amplitude:.2e}'.split('e')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.22, hspace=.65, wspace=.30)
    fig.suptitle('A sourced exterior response and a reduced RF spread test',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Physical matching remains open • RF results use a turn-averaged model and a separate linear transverse envelope',
             ha='center', fontsize=10.1, color='#475569')

    ax = axes[0, 0]
    inner, outer = radial[0], radial[-1]
    ax.plot([.80, inner], [.56, .56], color='#94a3b8', ls='--', lw=1.5)
    ax.plot([outer, 4.25], [.56, .56], color='#94a3b8', ls='--', lw=1.5)
    ax.plot([inner, outer], [.56, .56], color='#b7dfdf', lw=18, solid_capstyle='butt')
    for boundary in (inner, outer):
        ax.plot([boundary, boundary], [.43, .70], color='#007e87', lw=1.8)
    ax.text(inner, .80, 'Inner matching open', ha='left', color='#475569', fontsize=9)
    ax.text(outer, .80, 'Outer data open', ha='right', color='#475569', fontsize=9)
    ax.text((inner+outer)/2, .56, 'Finite RN exterior', ha='center', va='center', fontsize=10, color='#005d63')
    ax.text(.5, .25, rf'Known band only: $A_{{known}}\approx {amp_mantissa}\times10^{{{int(amp_exponent)}}}$',
            transform=ax.transAxes, ha='center', fontsize=10, color='#007e87')
    ax.text(.5, .09, 'Total amplitude also contains an unknown channel contribution.',
            transform=ax.transAxes, ha='center', fontsize=8.8, color='#475569')
    ax.set(title='A. The finite exterior region in this test',
           xlabel='Areal radius r / rₑ', xlim=(.8, 4.25), ylim=(0, 1))
    ax.set_xticks([1.01, 2, 3, 4], labels=['1.01', '2', '3', '4'])
    ax.set_yticks([])
    ax.spines[['left', 'right', 'top']].set_visible(False)
    ax.text(.5, -.31, 'The retained band covers half the flux range.\n'
            'Missing channels and physical end data have not been set to zero.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[0, 1]
    ax.plot(radial-1, source_h, color='#007e87', lw=2.2)
    ax.scatter([radial[0]-1, radial[-1]-1], [source_h[0], source_h[-1]], color='#007e87', s=25, zorder=4)
    ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=.01, linscale=.5); ax.minorticks_off()
    ax.set_xticks([.01, .1, .5, 1, 3], labels=['1.01', '1.1', '1.5', '2', '4'])
    ax.set(title='B. The additive response per unit source amplitude',
           xlabel='Areal radius r / rₑ\n(spacing logarithmic in r/rₑ − 1)',
           ylabel=r'Gauge-invariant $H_{\mathrm{source}}/A$',
           xlim=(.009, 3.2), ylim=(-.003, source_h.max()*2))
    ax.text(.5, -.36, 'Source column only; the physical end-data response is not included.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')

    ax = axes[1, 0]
    time_ms = (crest_time-crest_time[0])*1e3
    ax.fill_between(time_ms, crest_lower, crest_upper, color='#7951b0', alpha=.12)
    ax.plot(time_ms, crest_upper, color='#7951b0', ls='--', lw=1.7, label='|Offset| + linear envelope')
    ax.plot(time_ms, crest_lower, color='#7951b0', lw=2, label='|Offset| − linear envelope')
    ax.axhline(10, color='#64748b', ls=':', lw=1.4, label='1 cm aperture')
    ax.axvline(departure_ms, color='#b7791f', ls=':', lw=1.1)
    ax.text(departure_ms-.6, .97, f'{departure_ms:.2f} ms\nreduced-model departure',
            transform=ax.get_xaxis_transform(), ha='right', va='top', fontsize=8.1, color='#92400e', linespacing=1.4)
    ax.set(title='C. Crest control: the lower bound reaches 1 cm',
           xlabel='Elapsed laboratory time (ms)', ylabel='Absolute radial excursion bounds (mm)',
           ylim=(0, max(float(crest_upper.max()), 10)*1.28))
    ax.legend(frameon=True, facecolor='white', edgecolor='none', framealpha=.95, loc='upper left', fontsize=8.0)
    ax.text(.5, -.25, 'The lower-bound crossing implies departure only within this reduced model.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569')

    ax = axes[1, 1]
    positions = np.arange(4)
    ax.bar(positions, margins_mm, color='#007e87', width=.56)
    margin_max = float(margins_mm.max())
    for position, value in zip(positions, margins_mm):
        ax.text(position, value+margin_max*.025, f'{value:.2f}', ha='center', va='bottom', fontsize=8.8, color='#007e87')
    ax.set_xticks(positions, labels=corner_labels)
    ax.set(title='D. Four off-crest samples retain a radial margin',
           xlabel='Signs of initial phase and momentum offsets',
           ylabel='Smallest sampled radial margin (mm)', ylim=(0, margin_max*1.67))
    ax.text(.03, .95, f'{ramp_time:.2f} s ramp; 1 cm aperture\n'
            f'Initial offsets: ±{phase_degrees:.3g}°; Δp/p = ±{abs(momenta[0]):.1g}\n'
            f'Late 16-turn check: {late_model_difference_mm:.3f} mm kick-versus-average difference',
            transform=ax.transAxes, va='top', fontsize=8.2, color='#475569', linespacing=1.5)
    ax.text(.5, -.25, 'Sampled minima include a linear transverse envelope; this is not whole-beam acceptance.',
            transform=ax.transAxes, ha='center', va='top', fontsize=7.9, color='#475569')
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .125, 'Exterior: coupled equations include the source within a finite RN region. A full quantum state, physical matching and a stable global geometry remain open.',
             fontsize=9, color='#475569')
    fig.text(.055, .081, 'RF: four reduced-model samples and a crest control. Complete electromagnetic transport, whole-beam acceptance and permanent capture remain unproved.',
             fontsize=9, color='#475569')
    fig.text(.055, .033, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; no complete human journey has been demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'sourced_mouth_rf_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def _incident_orbit_progress_axes(incident):
    """Create the frozen classical-incident panels; the RF helper fills row two."""
    assert all(incident['verification'].values())
    orbit = incident['orbit']
    selected_orbit = [row for row in orbit['cases']
                      if row['model'] == 'EM_plus_standard_leading_GW']
    assert len(selected_orbit) == 1
    selected_orbit = selected_orbit[0]
    phase = selected_orbit['orbital_phase_radians']
    phase_degrees = selected_orbit['orbital_phase_degrees']
    radial_fraction = selected_orbit['fractional_separation_decrease']
    assert phase > 0 and 0 < radial_fraction < 1
    assert np.isclose(phase_degrees, np.degrees(phase), rtol=1e-12, atol=0)
    assert np.isclose(orbit['interval_over_re'], np.pi*orbit['ell_over_re'], rtol=1e-12, atol=0)
    incident_amplitude = incident['incident']['leading_master_growing_amplitude']
    assert incident_amplitude < 0
    assert np.isclose(incident_amplitude, -2/orbit['separation_over_re']**2, rtol=1e-12, atol=0)
    cuts = incident['finite_cut_relation']['cuts']
    cut_radii = np.asarray([row['radius_over_re'] for row in cuts], dtype=float)
    classical_terms = np.asarray([row['leading_classical_incident_forcing'] for row in cuts], dtype=float)
    assert np.allclose(cut_radii, [1.01, 1.1, 4], rtol=0, atol=1e-12)
    assert np.all(np.isfinite(classical_terms)) and np.all(classical_terms < 0)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9))
    fig.subplots_adjust(left=.085, right=.975, top=.835, bottom=.23, hspace=.66, wspace=.30)
    fig.suptitle('An incoming field, orbital motion and an RF event test',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .915, 'Separate conditional models • No joined rotating throat or built transport system',
             ha='center', fontsize=10.3, color='#475569')

    ax = axes[0, 0]
    centre_x, orbit_radius, angle_magnification = -.60, .63, 30
    theta = np.linspace(0, 2*np.pi, 240)
    ax.plot(centre_x+orbit_radius*np.cos(theta), orbit_radius*np.sin(theta),
            color='#94a3b8', ls=':', lw=1.2)
    enlarged_angle = phase*angle_magnification
    initial = np.array([[centre_x, orbit_radius], [centre_x, -orbit_radius]])
    final = np.array([[centre_x-orbit_radius*np.sin(enlarged_angle), orbit_radius*np.cos(enlarged_angle)],
                      [centre_x+orbit_radius*np.sin(enlarged_angle), -orbit_radius*np.cos(enlarged_angle)]])
    ax.scatter(initial[:, 0], initial[:, 1], s=55, color='#007e87', zorder=4)
    ax.scatter(final[:, 0], final[:, 1], s=55, facecolors='white', edgecolors='#7951b0', lw=1.7, zorder=5)
    ax.scatter([centre_x], [0], s=10, color='#64748b')
    ax.text(centre_x+.06, -.045, 'orbit centre', fontsize=7.4, color='#64748b', va='top')
    for sign in (-1, 1):
        ax.annotate('', xy=(centre_x, sign*.20), xytext=(centre_x, sign*.43),
                    arrowprops={'arrowstyle': '->', 'color': '#b7791f', 'lw': 1.5})
    ax.text(-1.30, 0, 'Attraction', color='#92400e', fontsize=8.0, va='center')
    arc_theta = np.linspace(np.pi/2, np.pi/2+enlarged_angle, 45)
    arc_x = centre_x+.78*np.cos(arc_theta)
    arc_y = .78*np.sin(arc_theta)
    ax.plot(arc_x, arc_y, color='#7951b0', lw=1.6)
    ax.annotate('', xy=(arc_x[-1], arc_y[-1]), xytext=(arc_x[-6], arc_y[-6]),
                arrowprops={'arrowstyle': '->', 'color': '#7951b0', 'lw': 1.5})
    ax.text(centre_x+.13, .61, '+Q', fontsize=8.7, color='#007e87')
    ax.text(centre_x-.29, -.68, '−Q', fontsize=8.7, color='#007e87')
    ax.plot(centre_x+.105*np.cos(theta), orbit_radius+.105*np.sin(theta),
            color='#007e87', lw=1, alpha=.5)
    ax.annotate('Local RN exterior\n1.01 ≤ r/rₑ ≤ 4', xy=(centre_x+.10, orbit_radius+.02),
                xytext=(.56, .93), textcoords='axes fraction', fontsize=8.7, color='#007e87', va='top',
                arrowprops={'arrowstyle': '-', 'color': '#007e87', 'lw': .9})
    ax.text(.56, .62, f'Orbital orientation change\n{phase_degrees:.4f}°', transform=ax.transAxes, va='top', color='#7951b0', fontsize=10, linespacing=1.5)
    ax.text(.56, .29, f'Fractional distance decrease\n{radial_fraction:.3g}',
            transform=ax.transAxes, va='top', color='#475569', fontsize=9.5, linespacing=1.5)
    ax.text(centre_x, -.92, f'Schematic; angle shown {angle_magnification}×',
            ha='center', fontsize=8, color='#475569')
    ax.set(title='A. Nearly fixed distance does not fix orientation', xlim=(-1.40, 1.60), ylim=(-1.08, 1.05))
    ax.set_aspect('equal', adjustable='datalim'); ax.axis('off')
    ax.text(.5, -.22, 'Assumed mouths with opposite magnetic charges; not to scale.\n'
            'Values use the initial throat light-crossing scale πℓ/c, not a payload transit.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[0, 1]
    cut_positions = np.arange(len(cuts))
    ax.bar(cut_positions, abs(classical_terms), width=.52, color='#007e87')
    ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks(cut_positions, labels=[f'{r:g}' for r in cut_radii])
    for position, value in zip(cut_positions, abs(classical_terms)):
        ax.text(position, value*1.45, f'{value:.3g}', ha='center', va='bottom', fontsize=8.5, color='#007e87')
    ax.set(title='B. One classical incoming term is now fixed',
           xlabel='Finite exterior cut: areal radius r / rₑ',
           ylabel='Magnitude of known boundary term', ylim=(1e-22, 1.5e-14))
    ax.text(.03, .95, f'Leading companion amplitude Cg = {incident_amplitude:.3g}',
            transform=ax.transAxes, va='top', fontsize=8.5, color='#475569')
    ax.text(.5, -.22, 'All three signed terms are negative. Quantum incoming data, the full\n'
            'source amplitude and physical throat matching remain unresolved.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    return fig, axes


def incident_orbit_rf_progress():
    """Render frozen companion-orbit data and the event-resolved RF comparison."""
    root = Path(__file__).resolve().parent
    gravity = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    future = json.loads((root/'future_time_results.json').read_text(encoding='utf-8'))
    if 'companion_incident' not in gravity or 'rf_event_ramp' not in future:
        raise ValueError('Final companion and RF-event results are required; no figure written')
    event = future['rf_event_ramp']
    if str(event.get('status', '')).startswith('partial'):
        raise ValueError('The RF-event result is not frozen; no figure written')
    native = event['extended']
    averaged = event['matched_average']
    plot = averaged['plot']
    time = np.asarray(plot['lab_seconds'], dtype=float)
    raw = np.asarray(plot['raw_guiding_offset_metres'], dtype=float)
    midpoint = np.asarray(plot['midpoint_guiding_offset_metres'], dtype=float)
    matched = np.asarray(plot['averaged_matched_guiding_offset_metres'], dtype=float)
    difference = np.asarray(plot['midpoint_minus_matched_average_metres'], dtype=float)
    assert all(array.shape == time.shape for array in (raw, midpoint, matched, difference))
    assert all(np.all(np.isfinite(array)) for array in (time, raw, midpoint, matched, difference))
    assert len(time) > 100 and np.all(np.diff(time) > 0)
    assert np.allclose(difference, midpoint-matched, rtol=1e-10, atol=1e-16)
    maximum_difference = averaged['maximum_sampled_midpoint_guiding_difference_metres']
    assert np.isclose(np.max(abs(difference)), maximum_difference, rtol=1e-12, atol=0)
    summary = native['summary']
    sampled_pre_post = summary['maximum_pre_post_guiding_offset_metres']
    assert sampled_pre_post > 0 and maximum_difference > 0
    assert native['parameters']['extended_precision'] is True
    assert native['parameters']['extrema_stride_turns'] == 16
    assert np.isclose(native['parameters']['initial_phase_error_radians'], .02, rtol=0, atol=1e-14)
    assert np.isclose(native['parameters']['initial_relative_momentum_error'], 1e-7, rtol=1e-12, atol=0)
    final_time = native['endpoint']['lab_seconds']
    assert 24 < final_time < 25 and time[0] > 0 and time[-1] < final_time
    assert summary['gaps'] > 30_000_000

    fig, axes = _incident_orbit_progress_axes(gravity['companion_incident'])
    ax = axes[1, 0]
    # These are sparse events over the full ramp, with a dense final window.
    # Unconnected markers do not imply that every synchrotron cycle is resolved.
    ax.scatter(time, raw*1e3, s=3, alpha=.50, color='#007e87', linewidths=0)
    ax.axhline(0, color='#94a3b8', lw=.7, zorder=0)
    raw_limit_mm = 1.35*max(float(np.max(abs(raw)))*1e3, sampled_pre_post*1e3)
    ax.set(title=f'C. One event map, {summary["gaps"]/1e6:.2f} million ideal RF gaps',
           xlabel='Elapsed laboratory time (s)', ylabel='Raw post-kick guiding offset (mm)',
           xlim=(0, final_time*1.015), ylim=(-raw_limit_mm, raw_limit_mm))
    ax.text(.03, .96, f'Sampled pre/post maximum: {sampled_pre_post*1e3:.2f} mm\n'
            '(native extrema sampled every 16 turns)',
            transform=ax.transAxes, va='top', fontsize=8.2, color='#475569', linespacing=1.5)
    ax.text(.5, -.25, 'Points are saved post-kick samples, not a resolved whole-ramp trace.\n'
            'Prescribed zero-length gaps; transverse motion is outside this model.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[1, 1]
    ax.scatter(time, difference*1e6, s=3, alpha=.55, color='#7951b0', linewidths=0)
    ax.axhline(0, color='#94a3b8', lw=.7, zorder=0)
    difference_limit_um = maximum_difference*1e6*1.35
    ax.set(title='D. Compare the models at the same RF sections',
           xlabel='Elapsed laboratory time (s)',
           ylabel='Midpoint map − matched average (µm)',
           xlim=(0, final_time*1.015), ylim=(-difference_limit_um, difference_limit_um))
    ax.text(.03, .96, f'Largest sampled difference: {maximum_difference*1e6:.2f} µm\n'
            'Same source times and matched initialization',
            transform=ax.transAxes, va='top', fontsize=8.2, color='#475569', linespacing=1.5)
    ax.text(.5, -.25, 'Use the energy midpoint of each kick and a consistently initialized average.\n'
            'This sampled comparison is not a continuous error bound.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .115, 'Orbit: equal mouths share the leading clock rescaling; no relative offset is created. The rotating quantum state and throat matching remain open.',
             fontsize=8.9, color='#475569')
    fig.text(.055, .077, 'RF: one initial condition in a longitudinal reduction. Finite cavities, transverse coupling, a full beam distribution and permanent capture remain unproved.',
             fontsize=8.9, color='#475569')
    fig.text(.055, .030, '13 → 12 Sep 2027: unproved.     Forward time dilation is supported; no complete human journey has been demonstrated.',
             fontsize=10, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'incident_orbit_rf_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def rotating_screening_progress():
    """Render the rotating-source balance test and finite-screening diagnostics."""
    from matplotlib.patches import FancyBboxPatch

    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['rotating_source']
    verified = record['verification']
    for key in ('direct_four_dimensional_naive_and_flux_divergence', 'dipole_angular_laplacian',
                'trace_free', 'particular_and_exponents_exact'):
        assert verified[key] is True
    assert max(verified['independent_collocation_max_basis_errors']) < 1e-8
    parameters = record['parameters']
    kappa = parameters['kappa']
    assert np.isclose(kappa, 10, rtol=0, atol=1e-12)
    assert np.isclose(kappa, parameters['alpha']*parameters['q']**2, rtol=1e-12, atol=0)
    screening = record['screening']
    cases = screening['comparison_cases']
    unscreened = [row for row in cases if row['kappa'] == 0]
    actual = [row for row in cases if row['kappa'] == kappa]
    assert len(unscreened) == len(actual) == 1
    unscreened, actual = unscreened[0], actual[0]
    finite_ratio = screening['screened_to_unscreened_particular_ratio']
    relative_error = screening['large_kappa_relative_particular_error']
    approximate_ratio = finite_ratio*(1+relative_error)
    ratios = np.array([unscreened['screened_to_unscreened_particular_ratio'],
                       finite_ratio, approximate_ratio], dtype=float)
    assert np.allclose(ratios, [1, 1/6, 1/5], rtol=1e-12, atol=0)
    assert np.isclose(approximate_ratio, 2/kappa, rtol=1e-12, atol=0)
    assert np.isclose(finite_ratio, abs(screening['actual_particular_over_q_omega']
                                      /unscreened['particular_over_q_omega']), rtol=1e-12, atol=0)
    radius = np.asarray(record['plot']['radius_over_re'], dtype=float)
    inner = np.asarray(record['plot']['inner_boundary_weight'], dtype=float)
    outer = np.asarray(record['plot']['outer_boundary_weight'], dtype=float)
    assert radius.shape == inner.shape == outer.shape and len(radius) > 20
    assert all(np.all(np.isfinite(array)) for array in (radius, inner, outer))
    assert np.all(np.diff(radius) > 0)
    assert np.all((inner >= 0) & (inner <= 1)) and np.all((outer >= 0) & (outer <= 1))
    assert np.all(np.diff(inner) <= 0) and np.all(np.diff(outer) >= 0)
    assert np.allclose(radius[[0, -1]], parameters['diagnostic_buffer_radii_over_re'], rtol=1e-12, atol=0)
    assert np.allclose(inner[[0, -1]], [1, 0], rtol=0, atol=1e-12)
    assert np.allclose(outer[[0, -1]], [0, 1], rtol=0, atol=1e-12)
    midpoint = int(np.argmin(abs(np.log(radius)-np.log(np.sqrt(radius[0]*radius[-1])))))
    assert np.isclose(inner[midpoint], actual['midpoint_inner_boundary_weight'], rtol=1e-10, atol=0)
    assert np.isclose(outer[midpoint], actual['midpoint_outer_boundary_weight'], rtol=1e-10, atol=0)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2))
    fig.subplots_adjust(left=.055, right=.969, top=.76, bottom=.32, wspace=.39)
    fig.suptitle('Rotation needs energy flow; screening needs boundary data',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .88, 'A fixed-background conservation test and a separate monopole screening calculation',
             ha='center', fontsize=10.5, color='#475569')

    ax = axes[0]
    ax.set(title='A. Rotating density\nneeds energy flow', xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    cards = [(.70, .23, '#fff0ee', '#b45346', 'Rotate this density alone',
              'Energy balance fails at first order.'),
             (.37, .23, '#e9f5f3', '#007e87', 'Add angular energy flow',
              'First-order balance is restored.'),
             (.03, .24, '#fff5e4', '#a56c15', 'Second-order forces remain',
              'Required stress and charge current\nare still unresolved.')]
    for y, height, face, color, title, detail in cards:
        ax.add_patch(FancyBboxPatch((.02, y), .96, height,
                                   boxstyle='round,pad=.012,rounding_size=.035',
                                   facecolor=face, edgecolor=color, linewidth=1.0))
        ax.text(.5, y+height-.045, title, ha='center', va='top', color=color,
                fontsize=9.6, fontweight='bold')
        detail_offset = .11 if '\n' in detail else .13
        ax.text(.5, y+height-detail_offset, detail, ha='center', va='top', color='#475569',
                fontsize=8.4, linespacing=1.45)
    for top, bottom in ((.685, .615), (.355, .285)):
        ax.annotate('', xy=(.5, bottom), xytext=(.5, top),
                    arrowprops={'arrowstyle': '->', 'color': '#64748b', 'lw': 1.3})
    ax.text(.5, -.20, 'A necessary balance example on fixed RN geometry.\n'
            'The full charged quantum source has not been constructed.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569', linespacing=1.5)

    ax = axes[1]
    bars = ax.bar(np.arange(3), ratios, width=.56, color=['#94a3b8', '#007e87', '#c89135'])
    bars[2].set_hatch('//')
    for position, value, label in zip(range(3), ratios, ('1', '1/6', '1/5')):
        ax.text(position, value+.035, label, ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set(title='B. Screening of the\nparticular solution',
           ylabel='Potential particular amplitude /\nunscreened particular amplitude', ylim=(0, 1.2))
    ax.set_xticks(range(3), labels=['Unscreened\nκ = 0', 'Finite κ\nκ = 10', 'Large-κ formula\nat κ = 10'])
    ax.tick_params(axis='x', labelsize=8.5)
    ax.text(.67, .74, f'Large-κ estimate:\n{relative_error*100:.0f}% too large',
            transform=ax.transAxes, ha='center', va='top', color='#92400e', fontsize=9.4, linespacing=1.5)
    ax.annotate('', xy=(2, ratios[2]+.10), xytext=(1.72, .58),
                arrowprops={'arrowstyle': '->', 'color': '#b7791f', 'lw': 1.2})
    ax.text(.5, -.20, 'Ratios of particular solutions, not the full potential.\n'
            'The 20% error is relative to the finite-κ particular.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569', linespacing=1.5)

    ax = axes[2]
    ax.plot(radius, inner, color='#007e87', lw=2.1, label='Inner-boundary weight')
    ax.plot(radius, outer, color='#7951b0', lw=2.1, label='Outer-boundary weight')
    ax.set_xscale('log'); ax.minorticks_off()
    ax.set_xticks([100, 1000, 10000], labels=['100', '1,000', '10,000'])
    ax.set(title='C. Physical boundary\nvalues are still needed',
           xlabel='Radius r / rₑ (log scale)', ylabel='Unit boundary-response weight',
           xlim=(radius[0], radius[-1]), ylim=(-.025, 1.12))
    ax.legend(frameon=False, loc='upper center', fontsize=8.1)
    ax.text(.5, .52, f'At r = {radius[midpoint]:,.0f} rₑ:\n'
            f'inner ≈ {inner[midpoint]:.0e}\nouter ≈ {outer[midpoint]:.0e}',
            transform=ax.transAxes, ha='center', va='center', fontsize=8.9, color='#475569', linespacing=1.5)
    ax.text(.5, -.20, 'These are weights multiplying unknown end potentials.\n'
            'The physical end potentials have not been set to zero.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569', linespacing=1.5)

    for ax in axes[1:]:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .150, 'Source equation: Maldacena, Milekhin & Popov, Appendix G (G.38).',
             fontsize=9, color='#475569', url='https://arxiv.org/html/1807.04726v3')
    fig.text(.055, .112, 'Finite-κ solution, boundary-response weights and conservation checks: our derivations.',
             fontsize=9, color='#475569')
    fig.text(.055, .052, 'No stable medium or 13 → 12 Sep 2027 passage established.     Forward branch unchanged.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'rotating_screening_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def loop_charge_progress():
    """Render complete-loop charge projection, exterior parity and optical weights."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['loop_charge_response']
    cases = sorted(record['cases'], key=lambda row: row['nu'])
    nu = np.asarray([row['nu'] for row in cases], dtype=float)
    throat = np.asarray([row['throat_mean_weight'] for row in cases], dtype=float)
    exterior = np.asarray([row['exterior_mean_weight'] for row in cases], dtype=float)
    assert np.allclose(nu, [.5, 1, 1.5], rtol=0, atol=1e-12)
    assert np.all((throat > 0) & (throat < 1)) and np.all((exterior > 0) & (exterior < 1))
    assert np.allclose(throat+exterior, 1, rtol=0, atol=1e-12)
    for row in cases:
        total = row['total_optical_length_over_d']
        assert np.isclose(total, row['throat_optical_length_over_d']
                          + row['exterior_optical_length_over_d'], rtol=1e-12, atol=0)
        assert np.isclose(row['throat_mean_weight'], row['throat_optical_length_over_d']/total,
                          rtol=1e-12, atol=0)
        assert max(abs(value) for value in row['exterior_trial_potential_integrals']) < 1e-14
        assert row['charge_balance_residual'] < 1e-12
        assert row['constrained_solve_max_error'] < 1e-12
    trial = [row for row in record['plot'] if row['nu'] == 1]
    assert len(trial) == 1
    coordinate = np.asarray(trial[0]['prolate_coordinate'], dtype=float)
    potential = np.asarray(trial[0]['bare_potential_over_qOmega'], dtype=float)
    assert coordinate.shape == potential.shape and len(coordinate) > 20
    assert np.all(np.isfinite(coordinate)) and np.all(np.isfinite(potential))
    assert np.all(np.diff(coordinate) > 0)
    assert np.allclose(coordinate, -coordinate[::-1], rtol=0, atol=1e-12)
    assert np.allclose(potential, -potential[::-1], rtol=1e-12, atol=1e-14)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2))
    fig.subplots_adjust(left=.055, right=.966, top=.76, bottom=.32, wspace=.41)
    fig.suptitle('Charge response depends on the whole loop',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .88, 'Neutral equilibrium for each independent channel • Leading rotation order • Prescribed optical lengths',
             ha='center', fontsize=10.1, color='#475569')

    ax = axes[0]
    ax.set(title='A. The complete optical\naverage enters the response', xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    theta = np.linspace(0, np.pi, 150)
    ax.plot(.5+.43*np.cos(theta), .65+.23*np.sin(theta), color='#007e87', lw=2.3)
    ax.plot(.5+.43*np.cos(theta), .65-.23*np.sin(theta), color='#7951b0', lw=2.3, ls='--')
    ax.scatter([.07, .93], [.65, .65], color='#64748b', s=30, zorder=4)
    ax.text(.5, .94, 'Exterior trial', ha='center', va='top', color='#007e87', fontsize=10)
    ax.text(.5, .65, 'Complete loop', ha='center', va='center', fontsize=8.7, color='#475569')
    ax.text(.5, .36, 'Throat: mean potential unknown', ha='center', va='top', color='#7951b0', fontsize=9.4)
    ax.text(.5, .21, r'$j_{\rm neutral}=(A_t-\langle A_t\rangle_{\rm loop})/\pi$',
            ha='center', va='center', fontsize=11.2, color='#005d63')
    ax.text(.5, .055, 'Subtract the whole-loop optical mean.',
            ha='center', va='center', fontsize=8.6, color='#475569')
    ax.text(.5, -.20, 'j is channel charge density in optical coordinates.\n'
            'Throat attachment and state unresolved; diagram not to scale.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[1]
    ax.plot(coordinate, potential, color='#007e87', lw=2.2)
    ax.fill_between(coordinate, 0, potential, where=(potential >= 0), color='#007e87', alpha=.12)
    ax.fill_between(coordinate, 0, potential, where=(potential <= 0), color='#007e87', alpha=.12)
    ax.axhline(0, color='#64748b', lw=.9)
    ax.axvline(0, color='#94a3b8', lw=.8, ls=':')
    ax.set(title='B. An odd exterior trial\ncancels by reflection',
           xlabel='Exterior coordinate s (prolate)',
           ylabel='Trial potential Aₜ / (qΩ)', xlim=(-1, 1), ylim=(-.63, .63))
    ax.text(.04, .96, 'Saved bare trial, channel ν = 1', transform=ax.transAxes,
            va='top', fontsize=8.3, color='#475569',
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 1.5})
    ax.text(.5, -.20, 'Odd potential × even optical weight: exterior integral = 0.\n'
            'Exterior parity alone does not fix the throat mean.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[2]
    positions = np.arange(len(cases))
    ax.barh(positions, throat*100, color='#7951b0', height=.55, label='Throat')
    ax.barh(positions, exterior*100, left=throat*100, color='#007e87', height=.55, label='Exterior')
    for position, weight in zip(positions, throat):
        ax.text(weight*50, position, f'{weight*100:.2f}%', ha='center', va='center',
                fontsize=10.2, fontweight='bold', color='white')
    ax.set_yticks(positions, labels=[f'ν = {value:g}' for value in nu])
    ax.set(title='C. Most prescribed optical\nlength lies in the throat',
           xlabel='Fraction of total optical length (%)',
           xlim=(0, 100), ylim=(-.95, 2.6))
    ax.invert_yaxis()
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.legend(loc='upper center', ncol=2, frameon=False, fontsize=8.7)
    ax.text(.5, -.20, 'Fractions of prescribed optical length, not research\n'
            'completion or solved field strengths. The throat mean is unknown.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    for ax in axes[1:]:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .150, 'Bosonized action: MMP, Appendix G. Closed neutral-loop projection and parity checks: our constrained application.',
             fontsize=9, color='#475569', url='https://arxiv.org/html/1807.04726v3')
    fig.text(.055, .112, 'The cancellation uses a symmetric zeroth-order optical measure and an odd first-order potential; no full rotating state is established.',
             fontsize=9, color='#475569')
    fig.text(.055, .052, 'No stable medium or 13 → 12 Sep 2027 passage established.     Forward branch unchanged.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'loop_charge_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def loop_retarded_progress():
    """Render causal charge profiles and the scope of a matched finite ramp."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['loop_retarded_response']
    plot = record['plot']
    y = np.asarray(plot['y'], dtype=float)
    tau = np.asarray(plot['tau'], dtype=float)
    step = np.asarray(plot['step_first_mode'], dtype=float)
    matched = np.asarray(plot['matched_first_mode'], dtype=float)
    drive = np.asarray(plot['matched_drive'], dtype=float)
    assert len(y) > 20 and len(tau) > 20
    assert np.all(np.diff(y) > 0) and np.all(np.diff(tau) > 0)
    assert np.allclose(y[[0, -1]], [0, 1], rtol=0, atol=1e-12)
    assert all(a.shape == tau.shape and np.all(np.isfinite(a)) for a in (step, matched, drive))
    profiles = []
    for time in (.12, .24):
        rows = [row for row in plot['step_profiles'] if np.isclose(row['tau'], time, rtol=0, atol=1e-12)]
        assert len(rows) == 1
        density = np.asarray(rows[0]['density'], dtype=float)
        assert density.shape == y.shape and np.all(np.isfinite(density))
        assert np.isclose(density[0], density[-1], rtol=0, atol=1e-12)
        profiles.append((time, density))
    assert record['propagation']['max_charge_error'] < 1e-12
    control = record['matched_control']
    duration = control['total_duration_over_reference_length']
    assert np.isclose(duration, 1.3, rtol=0, atol=1e-12)
    assert np.allclose(matched[tau >= duration], 1, rtol=0, atol=1e-8)
    assert np.allclose(drive[tau >= duration], 1, rtol=0, atol=1e-12)
    assert max(control['matched_residual_amplitudes']) < 1e-9
    cases = sorted(control['actual_loop_cases'], key=lambda row: row['nu'])
    nu = np.asarray([row['nu'] for row in cases], dtype=float)
    residual_percent = np.asarray([row['residual_amplitudes'][0]*100 for row in cases], dtype=float)
    assert np.allclose(nu, [.5, 1, 1.5], rtol=0, atol=1e-12)
    assert np.all(np.isfinite(residual_percent)) and np.all(residual_percent >= 0)
    assert residual_percent[1] < 1e-9 and residual_percent[0] > residual_percent[2] > 1
    for row in cases:
        assert abs(row['residual_amplitudes'][0]-row['predicted_amplitudes'][0]) < 1e-8
    feedback = record['maxwell_comparison']['cases']
    assert feedback[0]['mass_times_Lref'] == 0
    assert all(row['residual_relative_to_screened_static'] > .01 for row in feedback[1:])

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2))
    fig.subplots_adjust(left=.070, right=.970, top=.76, bottom=.32, wspace=.42)
    fig.suptitle('A timed drive cancels mean charge oscillations in one ideal loop',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .88, 'Prescribed static optical loop • External potential drive • No damping or physical source apparatus',
             ha='center', fontsize=10.3, color='#475569')

    ax = axes[0]
    for (time, density), color, style in zip(profiles, ['#007e87', '#7951b0'], ['-', '--']):
        ax.plot(y, density, color=color, ls=style, lw=2, label=f'τ = {time:.2f}')
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set(title='A. Charge waves travel\naround the fixed loop',
           xlabel='Optical position x / L (periodic)', ylabel='Charge density (normalized)',
           xlim=(0, 1), ylim=(-.63, 1.19))
    ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=8.4)
    ax.text(.5, -.20, 'τ = ct/L; the plot ends are the same point on the loop.\n'
            'Saved step-response profiles show causal propagation.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[1]
    ax.plot(tau, step, color='#c18b32', ls='--', lw=1.7, label='Sudden step response')
    ax.plot(tau, matched, color='#007e87', lw=2.1, label='Matched-ramp response')
    ax.plot(tau, drive, color='#7951b0', ls=':', lw=1.7, label='Matched external drive')
    ax.axhline(1, color='#64748b', ls='-.', lw=.9, label='Final equilibrium = 1')
    ax.axvline(duration, color='#94a3b8', lw=.7, ls=':')
    ax.set(title='B. A specially timed ramp\ncancels the free transients',
           xlabel=r'Time $ct/L_{\rm ref}$', ylabel='First-mode response or drive\n(normalized to final value)',
           xlim=(tau[0], tau[-1]), ylim=(-.06, 2.55))
    ax.legend(loc='upper center', ncol=2, frameon=True, facecolor='white', edgecolor='none',
              fontsize=7.0, columnspacing=.7, handlelength=1.7)
    ax.text(.5, -.20, r'Ramp ends at $1.3 L_{\rm ref}/c$; first charge mode shown.' '\n'
            'Cancellation is analytic for the free reference harmonics.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)

    ax = axes[2]
    positions = np.arange(len(cases))
    ax.bar(positions, residual_percent, color=['#7951b0', '#007e87', '#7951b0'], width=.58)
    ax.scatter([1], [residual_percent[1]], color='#007e87', s=28, zorder=4)
    for position, value in zip(positions, residual_percent):
        label = '≈0%' if position == 1 else f'{value:.2f}%'
        ax.text(position, value+.15, label, ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set_xticks(positions, labels=['ν = 0.5', 'ν = 1\nreference', 'ν = 1.5'])
    ax.set(title='C. The same ramp does not\ncancel unequal loops',
           ylabel='First-mode residual\n(% of static amplitude)',
           ylim=(0, float(residual_percent.max())*1.25))
    ax.text(.5, -.20, 'One common ramp on three prescribed loop lengths.\n'
            'Residuals are relative to each static first mode.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .150, 'Action: MMP, Appendix G. Retarded response and timed-ramp calculations: our application of that model.',
             fontsize=9, color='#475569', url='https://arxiv.org/html/1807.04726v3')
    fig.text(.055, .112, 'A separate electromagnetic feedback test shifts the frequencies and leaves residual waves; actual coupled fields and throat matching remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .052, 'No stable medium or 13 → 12 Sep 2027 passage established.     Forward branch unchanged.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'loop_retarded_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def throat_charge_progress():
    """Render the fixed-throat electric-displacement boundary diagnostics."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['throat_charge_response']
    plot = record['plot']
    chi = np.asarray(plot['chi'], dtype=float)
    restoring = np.asarray(plot['potential'], dtype=float)
    same = np.asarray(plot['same_end_displacement'], dtype=float)
    opposite = np.asarray(plot['opposite_end_displacement'], dtype=float)
    assert len(chi) > 20 and np.all(np.diff(chi) > 0)
    assert all(a.shape == chi.shape and np.all(np.isfinite(a)) for a in (restoring, same, opposite))
    assert np.allclose(chi, -chi[::-1], rtol=0, atol=1e-12)
    assert np.allclose(restoring, record['parameters']['kappa']/np.cos(chi)**2, rtol=1e-12, atol=0)
    assert np.allclose(same, same[::-1], rtol=1e-12, atol=1e-12)
    assert np.allclose(opposite, -opposite[::-1], rtol=1e-12, atol=1e-12)
    assert np.allclose(same[[0, -1]], [1, 1], rtol=0, atol=1e-12)
    assert np.allclose(opposite[[0, -1]], [-1, 1], rtol=0, atol=1e-12)
    cases = sorted(record['static_cases'], key=lambda row: row['rho_cut'])
    cuts = np.asarray([row['rho_cut'] for row in cases], dtype=float)
    ratios = np.asarray([row['common_endpoint_displacement_to_center'] for row in cases], dtype=float)
    assert np.allclose(cuts, [4, 10, 100, 1000], rtol=0, atol=1e-12)
    assert np.all(np.isfinite(ratios)) and np.all(ratios > 0) and np.all(np.diff(ratios) < 0)
    center = int(np.argmin(abs(chi)))
    assert abs(chi[center]) < 1e-12
    assert np.isclose(same[center], ratios[0], rtol=1e-10, atol=0)
    assert np.allclose(chi[[0, -1]], [-cases[0]['chi_cut'], cases[0]['chi_cut']], rtol=0, atol=1e-12)
    assert record['independent_bvp']['profile_max_error'] < 1e-8

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2))
    fig.subplots_adjust(left=.070, right=.969, top=.76, bottom=.32, wspace=.45)
    fig.suptitle('Electric displacement in a prescribed curved throat',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .88, 'D = electric displacement, not voltage • Spherical charged mode ℓ = 0 • Fixed leading geometry',
             ha='center', fontsize=10.3, color='#475569')

    ax = axes[0]
    ax.plot(chi, restoring, color='#007e87', lw=2.2)
    ax.fill_between(chi, 0, restoring, color='#007e87', alpha=.07)
    for boundary in chi[[0, -1]]:
        ax.axvline(boundary, color='#94a3b8', ls=':', lw=1)
    ax.set(title='A. Curved restoring term;\nphysical end data unknown',
           xlabel='Optical coordinate χ', ylabel='Wave coefficient κ sec²χ',
           xlim=(chi[0]-.04, chi[-1]+.04), ylim=(0, float(restoring.max())*1.42))
    ax.text(.03, .96, 'Left D:\nunknown', transform=ax.transAxes,
            ha='left', va='top', fontsize=9, color='#475569', linespacing=1.4)
    ax.text(.97, .96, 'Right D:\nunknown', transform=ax.transAxes,
            ha='right', va='top', fontsize=9, color='#475569', linespacing=1.4)
    ax.text(.5, -.20, 'The curve is a wave-equation coefficient, not a voltage.\n'
            'These finite cuts are not established physical mouth joins.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[1]
    ax.plot(chi, same, color='#007e87', lw=2.2, label='Left D = +1, right D = +1')
    ax.plot(chi, opposite, color='#7951b0', ls='--', lw=2.0, label='Left D = −1, right D = +1')
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set(title='B. Arbitrary unit boundary\nvalues test the static response',
           xlabel='Optical coordinate χ', ylabel='Normalized electric displacement D',
           xlim=(chi[0], chi[-1]), ylim=(-1.13, 1.43))
    ax.legend(loc='upper center', frameon=True, facecolor='white', edgecolor='none', fontsize=8.1)
    ax.text(.5, -.20, 'The imposed ±1 values are arbitrary test units.\n'
            'Actual end displacement and its source remain unknown.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[2]
    ax.plot(cuts, ratios, marker='o', color='#007e87', lw=1.8, markersize=5)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
    ax.set_xticks(cuts, labels=['4', '10', '100', '1,000'])
    for cut, ratio in zip(cuts, ratios):
        ax.annotate(f'{ratio:.2e}', xy=(cut, ratio), xytext=(0, 9), textcoords='offset points',
                    ha='left' if cut == cuts[0] else 'center', va='bottom', fontsize=8.3, color='#007e87')
    ax.set(title='C. Equal-end tests at\ndifferent cutoff positions',
           xlabel=r'Cutoff $\rho_{\rm cut}$ (log scale)',
           ylabel=r'$D(0)/D_{\rm end}$ (equal end values)',
           xlim=(3, 1500), ylim=(ratios[-1]/5, ratios[0]*5))
    ax.text(.5, -.20, 'A static centre-to-boundary ratio in four saved tests.\n'
            'This is not a transmission probability or stability result.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .150, 'Charged collective action: MMP, Appendix A. Finite-cut electric response and checks: our application.',
             fontsize=9, color='#475569', url='https://arxiv.org/html/1807.04726v3')
    fig.text(.055, .112, 'This spherical ℓ = 0 test does not solve the rotating dipole, physical exterior source or nonlinear metric feedback.',
             fontsize=9, color='#475569')
    fig.text(.055, .052, 'No stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'throat_charge_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def charged_mouth_progress():
    """Show the saved leading static throat/RN match and its unit tests."""
    from matplotlib.patches import FancyBboxPatch

    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['charged_mouth_matching']
    plot = record['plot']
    radial_offset = np.asarray(plot['radial_offset_over_re'], dtype=float)
    displacement = np.asarray(plot['unit_outer_static_displacement'], dtype=float)
    cases = sorted(record['cases'], key=lambda row: row['rho_cut'])
    cuts = np.asarray([row['rho_cut'] for row in cases], dtype=float)
    centers = np.asarray([row['corrected_center_per_unit_outer_displacement'] for row in cases], dtype=float)
    bias = np.asarray([row['extremal_omission_fraction'] for row in cases], dtype=float)
    extremal_centers = np.asarray([row['extremal_center_per_unit_outer_displacement'] for row in cases], dtype=float)
    profile_case = next(row for row in cases if row['rho_cut'] == plot['rho_cut'])
    spread = float(record['corrected_center_relative_cut_spread'])
    assert radial_offset.shape == displacement.shape and len(radial_offset) > 20
    assert np.all(np.isfinite(radial_offset)) and np.all(np.diff(radial_offset) > 0)
    assert np.all(np.isfinite(displacement)) and np.all(displacement > 0)
    assert np.allclose(cuts, [10, 100, 1000], rtol=0, atol=0)
    assert np.all(centers > 0) and np.all(bias < 0)
    assert np.allclose(bias, extremal_centers/centers-1, rtol=1e-9, atol=1e-15)
    assert np.isclose(radial_offset[0], profile_case['inner_radius_minus_re_over_re'], rtol=1e-12, atol=0)
    assert np.isclose(radial_offset[-1], record['parameters']['outer_areal_radius_over_re']-1, rtol=1e-12, atol=0)
    assert np.isclose(displacement[0], profile_case['inner_displacement_per_unit_outer_displacement'], rtol=1e-12, atol=0)
    assert np.isclose(displacement[-1], 1, rtol=0, atol=1e-12)
    assert np.isclose(plot['center_displacement'], profile_case['corrected_center_per_unit_outer_displacement'], rtol=1e-12, atol=0)
    assert 0 < spread < 1e-9 and record['parameters']['kappa'] == 10

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 6.2))
    fig.subplots_adjust(left=.048, right=.970, top=.76, bottom=.32, wspace=.48)
    fig.suptitle('Static electrical model — not a physical medium',
                 fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .88, 'Spherical mode ℓ = 0 • Fixed approximate geometry • D = electric displacement, not voltage',
             ha='center', fontsize=10.3, color='#475569')

    ax = axes[0]
    ax.set(xlim=(0, 1), ylim=(0, 1), title='A. Join finite regions;\nouter source still unknown')
    ax.axis('off')
    for x, width, label, color in ((.01, .25, 'Finite RN\nregion', '#007e87'),
                                    (.355, .29, 'Leading\nthroat', '#7951b0'),
                                    (.74, .25, 'Finite RN\nregion', '#007e87')):
        ax.add_patch(FancyBboxPatch((x, .64), width, .17,
                                   boxstyle='round,pad=0.008,rounding_size=0.012',
                                   facecolor=color, alpha=.10, edgecolor=color, linewidth=1.5))
        ax.text(x+width/2, .725, label, ha='center', va='center', fontsize=9, color=color)
    for left, right in ((.26, .355), (.645, .74)):
        ax.plot([left, right], [.725, .725], color='#64748b', ls='--', lw=1.2)
    for x in (.135, .865):
        ax.text(x, .96, 'Outer D:\nUNKNOWN', ha='center', va='top', fontsize=8.8,
                color='#475569', linespacing=1.4)
    ax.text(.5, .545, 'Response matrices • Leading overlap\nSchematic regions, not to scale',
            ha='center', va='center', fontsize=8.1, color='#475569', linespacing=1.5)
    ax.add_patch(FancyBboxPatch((.035, .055), .93, .365,
                               boxstyle='round,pad=0.015,rounding_size=0.018',
                               facecolor='#f1f5f9', edgecolor='#cbd5e1', linewidth=.8))
    ax.text(.5, .35, 'Arbitrary outer test: (1, 1)', ha='center', va='center', fontsize=9.2)
    ax.text(.5, .235, r'$D_{\rm center}/D_{\rm outer}\approx1.399\times10^{-33}$',
            ha='center', va='center', fontsize=11, color='#007e87')
    ax.text(.5, .115, 'Static screening only', ha='center', va='center', fontsize=9,
            fontweight='bold', color='#007e87')
    ax.text(.5, -.20, 'All D values use one direction through the model.\n'
            'Equal test values mean opposite outward radial flux.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[1]
    ax.loglog(radial_offset, displacement, color='#007e87', lw=2.2)
    ax.scatter(radial_offset[[0, -1]], displacement[[0, -1]], color='#007e87', s=23, zorder=3)
    ax.minorticks_off()
    ax.set_xticks([1e-7, 1e-4, 1e-1, 1e2])
    ax.set_yticks([1e-27, 1e-18, 1e-9, 1])
    ax.set(title='B. One radial RN segment:\nunit outer displacement',
           xlabel=r'Radial offset $(r-r_e)/r_e$', ylabel=r'$D/D_{\rm outer}$ (unit test)',
           xlim=(radial_offset[0]/2, radial_offset[-1]*2), ylim=(1e-29, 30))
    ax.annotate('inner match', xy=(radial_offset[0], displacement[0]), xytext=(6, 12),
                textcoords='offset points', ha='left', va='bottom', fontsize=8.4, color='#007e87',
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .90, 'pad': 1})
    ax.annotate('outer D = 1', xy=(radial_offset[-1], displacement[-1]), xytext=(-5, -11),
                textcoords='offset points', ha='right', va='top', fontsize=8.4, color='#007e87',
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .90, 'pad': 1})
    ax.text(.5, -.20, 'ρcut = 100; the curve stops at the throat match.\n'
            'The throat center is not on this radial curve.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)

    ax = axes[2]
    percent = -100*bias
    ax.loglog(cuts, percent, marker='o', color='#7951b0', lw=1.8, markersize=5)
    ax.minorticks_off()
    ax.set_xticks(cuts, labels=['10', '100', '1,000'])
    for index, (cut, value) in enumerate(zip(cuts, percent)):
        label = np.format_float_positional(value, precision=4, unique=False, fractional=False, trim='-')
        ax.annotate(f'−{label}%', xy=(cut, value), xytext=(0, 10), textcoords='offset points',
                    ha=('left' if index == 0 else 'right' if index == 2 else 'center'),
                    va='bottom', fontsize=8.7, color='#7951b0',
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .90, 'pad': 1})
    ax.set(title='C. Omitting the mass deficit\nunderestimates the center response',
           xlabel=r'Overlap cutoff $\rho_{\rm cut}$', ylabel='Magnitude of omitted-mass bias (%)',
           xlim=(7, 1500), ylim=(1e-5, 4))
    ax.text(.5, -.20, f'Corrected cut spread: {spread:.2e} (relative).\n'
            'Internal consistency, not physical accuracy.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569', linespacing=1.5)
    for ax in axes[1:]:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.055, .150, 'Charged action and overlap: MMP, Appendix A, Eqs. 5.39–5.42. Finite-region response matrices: our application.',
             fontsize=9, color='#475569', url='https://arxiv.org/html/1807.04726v3')
    fig.text(.055, .112, 'Sphere-profile, lapse and radial corrections remain omitted; no actual outer source or self-consistent geometry is fixed.',
             fontsize=9, color='#475569')
    fig.text(.055, .052, 'No stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'charged_mouth_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def nonstationary_support_progress():
    """Render the saved local plateau and fixed-throat support-duration test."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['nonstationary_support']
    x = np.asarray(record['plot']['chiral_coordinate'], dtype=float)
    stress = np.asarray(record['plot']['stress_over_required_magnitude'], dtype=float)
    plateau = record['plateau_candidate']
    width = float(plateau['guaranteed_local_window'])
    case = next(row for row in record['cases'] if row['rho_cut'] == 100)
    durations = np.array([case['maximum_unchanged_support_duration'], case['null_crossing_time']])
    assert x.shape == stress.shape and len(x) > 20 and np.all(np.diff(x) > 0)
    assert np.all(np.isfinite(stress)) and np.isclose(width, np.pi, rtol=1e-12)
    assert np.allclose(stress[np.abs(x) < width/2-1e-8], -1, rtol=1e-11, atol=1e-11)
    assert plateau['local_window_pass'] and not plateau['whole_throat_for_crossing_pass']
    assert not case['unchanged_support_for_full_crossing'] and 0 < durations[0] < durations[1]
    assert np.isclose(durations.sum(), record['maximum_chiral_support_width'], rtol=1e-12)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.subplots_adjust(left=.082, right=.966, top=.70, bottom=.29, wspace=.39)
    fig.suptitle('Local negative-energy window passes; unchanged-throat strategy fails',
                 fontsize=16, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Fixed free quantum-field model • Dynamic geometry untested',
             ha='center', fontsize=11, color='#475569')
    ax = axes[0]
    ax.axvspan(-width/2, width/2, color='#007e87', alpha=.07)
    ax.plot(x, stress, color='#007e87', lw=1.6, label='Saved plateau state')
    ax.axhline(-1, color='#b7791f', ls='--', lw=1.1, label='Required stress = −1')
    ax.set_yscale('symlog', linthresh=2)
    ax.set(title='A. A local interval has the required stress', xlabel='Chiral optical coordinate',
           ylabel='Stress / |required stress|', xlim=(x[0], x[-1]))
    ax.annotate('', xy=(width/2, .89), xytext=(-width/2, .89),
                xycoords=ax.get_xaxis_transform(), arrowprops={'arrowstyle': '<->', 'color': '#475569'})
    ax.text(.5, .92, 'Guaranteed local window W = π', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=8.5, color='#475569')
    ax.legend(loc='lower center', framealpha=.95, fontsize=8, ncol=2, edgecolor='none')
    ax.text(.5, -.25, 'Symmetric log scale shows large positive and negative transition spikes.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569')
    ax = axes[1]
    ax.barh([1, 0], durations, height=.28, color=['#7951b0', '#007e87'])
    for y, duration, label in zip([1, 0], durations, ['Whole-throat support limit', 'One light crossing']):
        ax.text(.02, y+.23, label, ha='left', va='bottom', fontsize=10)
        ax.text(duration+.06, y, f'{duration:.6f}', ha='left', va='center', fontsize=10)
    ax.set(title='B. The whole-throat window is too short (ρcut = 100)',
           xlabel='Duration (same optical units)', xlim=(0, 3.95), ylim=(-.42, 1.60), yticks=[])
    ax.text(.5, -.25, 'Full light crossing fails only for the unchanged-support strategy.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .118, 'Fixed-circle mean-stress test (Fewster–Hollands; MMP). No preparation protocol or seconds/calendar conversion.',
             fontsize=9, color='#475569')
    fig.text(.05, .050, 'No physical medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'nonstationary_support_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def plateau_gravity_progress():
    """Render the saved local leading-JT response without inferring collapse."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['plateau_gravity_response']
    rows = record['plot']
    time = np.array([row['global_time'] for row in rows])
    right = np.array([row['right_formal_coefficient_ratio'] for row in rows])
    left = np.array([row['left_formal_coefficient_ratio'] for row in rows])
    phi = np.array([row['right_fixed_cut_phi_over_alpha'] for row in rows])
    alpha = record['parameters']['alpha']
    minimum_area = record['minimum_sampled_total_area_over_extremal_area']
    assert len(time) > 20 and np.all(np.diff(time) > 0)
    assert all(np.all(np.isfinite(values)) for values in (time, right, left, phi))
    assert np.allclose(left, 1, rtol=0, atol=1e-12) and np.isclose(right[0], 1)
    assert abs(right[-1]) < 1e-12 and np.isclose(time[-1], record['formal_positive_right_coefficient_gate'])
    assert minimum_area > 0 and np.isclose(1+alpha*phi[-1], minimum_area, rtol=0, atol=3e-16)
    assert record['verification']['local_JT_max_scaled_residual'] < 1e-10
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.subplots_adjust(left=.077, right=.965, top=.71, bottom=.30, wspace=.33)
    fig.suptitle('Prescribed quantum state changes the mouth conditions',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Local leading-JT model • Incoming state prescribed • Boundary and state feedback unsolved',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]
    ax.plot(time, left, color='#7951b0', ls='--', lw=1.8, label='Left coefficient')
    ax.plot(time, right, color='#007e87', lw=2.1, label='Right coefficient')
    ax.scatter([time[-1]], [right[-1]], color='#007e87', s=25, zorder=3)
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set(title='A. Formal mouth response', ylabel='Coefficient / initial coefficient', ylim=(-.08, 1.26))
    ax.legend(loc='upper center', ncol=2, fontsize=8.2, frameon=False)
    ax.text(.97, .16, f'Formal zero at τ = {time[-1]:.6f}', transform=ax.transAxes,
            ha='right', fontsize=9, color='#007e87')
    ax.text(.5, -.25, 'A zero coefficient is not a horizon or vanished area.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.3, color='#475569')
    ax = axes[1]
    ax.plot(time, phi, color='#007e87', lw=2.1)
    ax.scatter(time[[0, -1]], phi[[0, -1]], color='#007e87', s=24, zorder=3)
    ax.annotate(f'Initial: {phi[0]:.3f}', (time[0], phi[0]), xytext=(7, -18),
                textcoords='offset points', fontsize=8.5, color='#007e87')
    ax.annotate(f'End: {phi[-1]:.5f}', (time[-1], phi[-1]), xytext=(-7, 11),
                textcoords='offset points', ha='right', fontsize=8.5, color='#007e87')
    ax.text(.56, .89, f'Smallest sampled total-area ratio:\n{minimum_area:.14f}  (> 0)',
            transform=ax.transAxes, ha='center', va='top', fontsize=9.1, color='#475569', linespacing=1.5,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .92, 'pad': 2})
    ax.set(title='B. Right fixed-cut area correction (ρcut = 100)', ylabel='Leading correction Φ / α', ylim=(-10, 190))
    ax.text(.5, -.25, 'Total area / reference area = 1 + Φ; Φ is only a correction.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.3, color='#475569')
    for ax in axes:
        ax.set(xlabel='Formal global time τ (optical units)', xlim=(0, time[-1]*1.025))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.14); ax.set_axisbelow(True)
    fig.text(.05, .12, 'No physical collapse time follows: exterior initial data are incomplete. Optical units do not specify clock days.',
             fontsize=9, color='#475569')
    fig.text(.05, .055, 'No physical medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'plateau_gravity_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def quantum_seed_feedback_progress():
    """Plot saved reduced-feedback clock shifts and launched-probe returns."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['quantum_seed_feedback']
    cases = record['cases']
    keys = sorted(cases, key=float)
    window = record['parameters']['final_exterior_time']
    delay = record['parameters']['exterior_delay']
    assert np.allclose([float(key) for key in keys], [0, .001, .01, .05, .1])
    assert np.isclose(delay, .1) and np.isclose(window, 2*record['parameters']['original_loop_time'])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    fig.subplots_adjust(left=.075, right=.966, top=.72, bottom=.30, wspace=.36)
    fig.suptitle('Clock offsets appear; no tested return reaches the past',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Initial quantum pulse • Coupled reduced feedback • Fixed exterior delay D = 0.1',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]
    for key, color in zip(['0.01', '0.05', '0.1'], ['#007e87', '#7951b0', '#b7791f']):
        plot = cases[key]['plot']
        time = np.asarray(plot['exterior_time'], dtype=float)
        offset = np.asarray(plot['right_clock_shift'])-np.asarray(plot['left_clock_shift'])
        assert time.shape == offset.shape and len(time) > 20 and np.all(np.diff(time) > 0)
        assert np.all(np.isfinite(offset)) and abs(offset[0]) < 1e-12 and np.isclose(time[-1], window)
        assert np.isclose(offset[-1], cases[key]['final_clock_offset'], rtol=1e-10, atol=1e-12)
        ax.plot(time, offset, color=color, lw=2, label=f'Seed amplitude {key}')
    ax.axhline(0, color='#94a3b8', lw=.8)
    ax.set(title='A. The quantum pulse creates a clock offset', xlabel='Exterior time (optical units)',
           ylabel='Right − left clock shift', xlim=(0, window))
    ax.legend(loc='lower left', fontsize=8.3, framealpha=.95, edgecolor='none')
    ax.text(.5, -.25, 'A coordinate clock difference alone does not test past travel.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.3, color='#475569')
    ax = axes[1]
    for index, key in enumerate(keys):
        probes = sorted((p for p in cases[key]['launched_probes'] if p['launch_exterior_time'] == 0),
                        key=lambda p: p['launch_side'])
        assert [p['launch_side'] for p in probes] == [0, 1]
        for side, (probe, color) in enumerate(zip(probes, ['#007e87', '#7951b0'])):
            y = index+(-.16 if side == 0 else .16)
            if probe['exit_reached_in_interval']:
                value = probe['complete_return_delay']
                assert value > 0 and np.isclose(value, probe['exit_exterior_time']+delay)
                ax.barh(y, value, height=.27, color=color, label=f'Launch side {side}' if index == 0 else None)
                ax.text(value+.07, y, f'{value:.3f}', ha='left', va='center', fontsize=8.1, color=color)
            else:
                assert probe['complete_return_delay'] is None and probe['exit_exterior_time'] is None
                ax.text(.12, y, f'Side {side}: no exit in window', ha='left', va='center', fontsize=8.2, color=color)
    ax.set_yticks(range(len(keys)), labels=keys)
    ax.set(title='B. Complete return delays for launch at t = 0', xlabel='Return delay (optical units)',
           ylabel='Seed amplitude', xlim=(0, 7.0), ylim=(len(keys)-.45, -.65))
    ax.legend(loc='upper right', fontsize=8.2, frameon=False)
    ax.text(.5, -.25, f'Window ends at {window:.6f}; missing exits do not prove trapping.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.3, color='#475569')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .12, 'Exterior gravity and physical preparation remain unresolved. Optical units are not seconds; the pulse supplies initial energy.',
             fontsize=9, color='#475569')
    fig.text(.05, .055, 'No stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'quantum_seed_feedback_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def causal_foliation_progress():
    """Contrast one-way ordering with complete return in saved affine tests."""
    from matplotlib.patches import Circle, FancyArrowPatch
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['causal_foliation']
    delay = float(record['parameters']['exterior_delay'])
    examples = record['affine_examples']
    throat = float(examples[0]['one_way_exterior_elapsed'])
    assert np.isclose(throat, np.pi) and np.isclose(delay, .1)
    for row in examples:
        assert np.isclose(row['one_way_exterior_elapsed'], throat-row['imposed_clock_offset'], atol=1e-12)
        assert np.isclose(row['complete_return_elapsed'], throat+delay-row['imposed_clock_offset'], atol=1e-12)
    example = next(row for row in examples if np.isclose(row['imposed_clock_offset'], throat+.05))
    assert example['one_way_exterior_elapsed'] < 0 < example['complete_return_elapsed']
    assert example['throat_lambda_elapsed'] > 0 and example['exterior_lambda_elapsed'] > 0
    offset = np.linspace(throat-.1, throat+delay+.05, 301)
    one_way, complete = throat-offset, throat+delay-offset
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12, 5.6))
    ax = fig.add_axes([.080, .30, .575, .43])
    inset = fig.add_axes([.725, .33, .25, .38])
    fig.suptitle('Solver time-ordering limit is not physical chronology protection',
                 fontsize=16, fontweight='bold', y=.965)
    fig.text(.5, .875, 'Prescribed toy clocks • No self-consistent source • Optical units, not seconds',
             ha='center', fontsize=10.5, color='#475569')
    ax.axvspan(throat, throat+delay, color='#007e87', alpha=.09)
    ax.plot(offset, one_way, color='#7951b0', lw=2, label='Throat leg: π − S (exterior clock)')
    ax.plot(offset, complete, color='#007e87', lw=2, label='Complete return: π + D − S')
    ax.axhline(0, color='#64748b', lw=.9)
    for boundary, label in ((throat, 'One-way zero\nS = π'), (throat+delay, 'Closed null loop\nS = π + D')):
        ax.axvline(boundary, color='#64748b', ls=':', lw=1)
        ax.text(boundary, 1.025, label, transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', fontsize=9, color='#475569')
    ax.text(throat+delay/2, .205, 'Alternative λ time\nremains forward',
            ha='center', va='center', fontsize=9.2, color='#007e87', linespacing=1.5)
    ax.scatter([example['imposed_clock_offset']]*2,
               [example['one_way_exterior_elapsed'], example['complete_return_elapsed']],
               color=['#7951b0', '#007e87'], s=28, zorder=3)
    ax.set(xlabel='Imposed clock offset S', ylabel='Elapsed exterior-clock time',
           xlim=(offset[0], offset[-1]), ylim=(-.18, .27))
    ax.legend(loc='lower left', fontsize=8.2, framealpha=.95, edgecolor='none')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.13); ax.set_axisbelow(True)
    inset.set(xlim=(0, 1), ylim=(0, 1), title='Example: S = π + 0.05')
    inset.axis('off')
    for x, label in ((.15, 'L'), (.85, 'R')):
        inset.add_patch(Circle((x, .49), .075, facecolor='#eef4f6', edgecolor='#64748b', lw=1))
        inset.text(x, .49, label, ha='center', va='center', fontsize=11)
    for start, end, color in (((.22, .53), (.78, .53), '#7951b0'), ((.78, .45), (.22, .45), '#007e87')):
        inset.add_patch(FancyArrowPatch(start, end, connectionstyle='arc3,rad=-.55',
                                       arrowstyle='-|>', mutation_scale=13, lw=1.6, color=color))
    inset.text(.5, .83, f"Throat: {example['one_way_exterior_elapsed']:+.2f}",
               ha='center', fontsize=11, color='#7951b0')
    inset.text(.5, .18, f'Exterior: +{delay:.2f}', ha='center', fontsize=11, color='#007e87')
    inset.text(.5, -.03, f"Full return: +{example['complete_return_elapsed']:.2f}",
               ha='center', fontsize=11, fontweight='bold', color='#007e87')
    fig.text(.05, .185, 'Current solver additionally requires transit > D/2; π is not its exact failure point. Coupled λ-feedback is not yet implemented.',
             fontsize=9, color='#475569')
    fig.text(.05, .120, 'The shaded interval has a negative one-way exterior-clock elapsed time but a positive complete return in this prescribed graph.',
             fontsize=9, color='#475569')
    fig.text(.05, .055, 'No 4D stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'causal_foliation_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def lambda_feedback_progress():
    """Compare saved lambda/exterior-time solutions at common physical times."""
    root = Path(__file__).resolve().parent
    data = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))
    record = data['lambda_feedback']
    new = record['cases']['0.05']
    old = data['quantum_seed_feedback']['cases']['0.05']
    time = np.asarray(new['plot']['exterior_time'], dtype=float)
    assert np.array_equal(time, np.asarray(old['plot']['exterior_time']))
    assert len(time) > 20 and np.all(np.diff(time) > 0)
    old_clocks = np.array([old['plot'][side+'_clock_shift'] for side in ('right', 'left')])
    new_clocks = np.array([new['plot'][side+'_clock_shift'] for side in ('right', 'left')])
    assert old_clocks.shape == new_clocks.shape == (2, len(time))
    assert np.all(np.isfinite(new_clocks)) and np.all(np.isfinite(old_clocks))
    assert np.allclose(new_clocks[:, 0], 0, rtol=0, atol=1e-12)
    assert np.allclose(old_clocks[:, 0], 0, rtol=0, atol=1e-12)
    difference = new_clocks-old_clocks
    displayed_max = float(np.max(np.abs(difference)))
    overlap = next(row for row in record['verification']['overlap'] if row['amplitude'] == .05)
    assert displayed_max < 1e-8 and overlap['maximum_clock_state_difference'] < 1e-8
    assert record['verification']['default_quantum_records_exact']
    for case in record['cases'].values():
        for probe in case['launched_probes']:
            assert probe['complete_return_delay'] > 0 if probe['exit_reached_in_interval'] else probe['complete_return_delay'] is None
    foliation = new['lambda_foliation']
    right_time, left_time = foliation['final_boundary_exterior_times']
    ledger = foliation['maximum_slice_energy_ledger_error']
    assert ledger < 1e-8 and foliation['minimum_evaluated_retarded_lambda_gap'] > 0
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    fig.subplots_adjust(left=.078, right=.964, top=.73, bottom=.33, wspace=.35)
    fig.suptitle('New coupled solver agrees; no backward return found',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Same initial pulse: amplitude 0.05 • Clocks synchronized initially • Fixed-delay exterior',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]
    ax.plot(time, old_clocks[0]-old_clocks[1], color='#94a3b8', lw=2.6, label='Exterior-time solver')
    ax.plot(time, new_clocks[0]-new_clocks[1], color='#007e87', lw=1.9, ls='--', label='λ-time solver')
    ax.set(title='A. The same physical clock offset', ylabel='Right − left clock shift')
    ax.legend(loc='lower left', fontsize=8.5, framealpha=.95, edgecolor='none')
    ax.text(.5, -.25, 'Curves overlap at shared exterior times.', transform=ax.transAxes,
            ha='center', va='top', fontsize=8.4, color='#475569')
    ax = axes[1]
    for values, side, color in zip(difference, ['Right clock', 'Left clock'], ['#007e87', '#7951b0']):
        ax.plot(time, values, color=color, lw=1.5, label=side)
    ax.axhline(0, color='#94a3b8', lw=.7)
    ax.set(title='B. New − old clocks at matched times', ylabel='Clock shift difference (optical units)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    ax.legend(loc='best', fontsize=8.3, framealpha=.95, edgecolor='none')
    ax.text(.5, -.25, f'Max displayed |difference|: {displayed_max:.2e} (optical units).',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569')
    for ax in axes:
        ax.set(xlabel='Exterior time (optical units)', xlim=(time[0], time[-1]))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .17, f'Final equal-λ slice: sR = {right_time:.6f}, sL = {left_time:.6f}; normalized energy-ledger residual {ledger:.2e}.',
             fontsize=9, color='#475569')
    fig.text(.05, .115, 'These pulses test only the shared regime; the extra affine domain remains untested. Exterior gravity remains unresolved.',
             fontsize=9, color='#475569')
    fig.text(.05, .050, 'No 4D stable medium or 13 → 12 Sep 2027 passage established. Optical units are not seconds.',
             fontsize=10.2, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'lambda_feedback_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def wide_packet_progress():
    """Compare matched launches, without mixing censored return subsets."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['wide_packet']
    keys = ['narrow', 'wide', 'wide_reversed']
    cases = [record['cases'][key] for key in keys]
    common = record['common_complete_probes']
    assert len(common) == 1 and common[0]['launch_phase'] == 0 and common[0]['launch_side'] == 0
    returns = np.asarray(common[0]['normalized_return_delays'], dtype=float)
    energy = record['verification']['matched_initial_energy_over_b0']
    counts, exits = [], []
    for index, case in enumerate(cases):
        probes = case['causal_diagnostics']['phase_matched_probes']
        summary = case['probe_summary']
        assert len(probes) == summary['total'] == 18
        assert sorted((p['launch_phase'], p['launch_side']) for p in probes) == [(i/8, side) for i in range(9) for side in (0, 1)]
        probe = next(p for p in probes if p['launch_phase'] == 0 and p['launch_side'] == 0)
        assert probe['complete_return_within_interval'] and probe['exit_reached_in_interval']
        assert np.isclose(returns[index], probe['complete_return_delay']/case['packet_parameters']['initial_loop_time'], rtol=1e-12)
        assert np.isclose(case['packet_parameters']['initial_energy_over_b0'], energy, rtol=1e-12)
        assert np.isclose(case['causal_diagnostics']['normalized_exterior_time'][-1], 2)
        counts.append(sum(p['complete_return_within_interval'] for p in probes))
        exits.append(sum(p['exit_reached_in_interval'] for p in probes))
        assert counts[-1] == summary['full_returns_within_interval'] and exits[-1] == summary['exits_reached']
    assert np.all(np.isfinite(returns)) and np.all(returns > 0)
    assert counts == [1, 13, 14] and exits == [1, 17, 17]
    positions = np.arange(3)
    colors = ['#7951b0', '#007e87', '#c18b32']
    labels = ['Narrow\nD = 0.1', 'Wide\nD = 1.2', 'Wide reversed\nD = 1.2']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    fig.subplots_adjust(left=.075, right=.967, top=.73, bottom=.34, wspace=.33)
    fig.suptitle('Improved passage in the model; no past return found',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .88, f'Matched initial energy I/b₀ = {energy:.6f} • Same 18 launch phases and sides • Two initial loops',
             ha='center', fontsize=10.4, color='#475569')
    ax = axes[0]
    ax.bar(positions, returns, color=colors, width=.6)
    ax.axhline(1, color='#64748b', ls='--', lw=1.2, label='Unperturbed loop = 1')
    ax.axhline(0, color='#64748b', lw=.9)
    for position, value in zip(positions, returns):
        ax.text(position, value+.06, f'{value:.6f}', ha='center', va='bottom', fontsize=9.3)
    ax.set(title='A. Same completed launch: right side, s/L = 0',
           ylabel='Complete return delay C / L', ylim=(-.12, 2.3))
    ax.legend(loc='upper right', fontsize=8.3, frameon=False)
    ax.text(.5, -.28, 'C = 0 is closed-null; a past return would require C < 0.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569')
    ax = axes[1]
    ax.bar(positions, [18]*3, color='#f1f5f9', width=.6)
    ax.bar(positions, counts, color=colors, width=.6)
    for position, count in zip(positions, counts):
        ax.text(position, count+.6, f'{count}/18', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set(title='B. Full returns within two initial loops',
           ylabel='Complete returns among 18 tests', ylim=(0, 20), yticks=[0, 6, 12, 18])
    ax.text(.5, -.28, 'Throat exits found: 1/18, 17/18, 17/18, respectively.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569')
    for ax in axes:
        ax.set_xticks(positions, labels=labels)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .16, 'Changing D also changes the assumed loop geometry and background support; this is not just pulse shaping.',
             fontsize=9, color='#475569')
    fig.text(.05, .105, 'L is each initial loop time. Missing returns concern only the finite interval; all times use optical model units.',
             fontsize=9, color='#475569')
    fig.text(.05, .045, 'No 4D source, stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'wide_packet_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def long_packet_progress():
    """Plot matched late returns and same-event local source diagnostics."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['long_packet']
    parameters = record['parameters']
    length, delay, alpha = (parameters[key] for key in ('initial_loop_time', 'exterior_delay', 'source_alpha'))
    shortcut = 2*delay/length
    assert np.isclose(length, np.pi+delay) and np.isclose(alpha, 1e-9)
    assert np.isclose(parameters['final_exterior_time'], 4*length)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
    fig.subplots_adjust(left=.075, right=.963, top=.73, bottom=.33, wspace=.33)
    fig.suptitle('Later returns improve; no past return or exterior shortcut found',
                 fontsize=16, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Four initial loops • Identical launch phases and sides • Declared source scale α = 10⁻⁹',
             ha='center', fontsize=10.5, color='#475569')
    completed, greatest = [], 1.2
    for key, label, color in (('wide', 'Wide', '#007e87'), ('wide_reversed', 'Reversed', '#7951b0')):
        case = record['cases'][key]
        probes = case['causal_diagnostics']['phase_matched_probes']
        summary = case['long_run_summary']
        assert len(probes) == summary['total_launches'] == 50
        completed.append(sum(p['complete_return_within_interval'] for p in probes))
        assert completed[-1] == summary['complete_returns_within_interval']
        for side, linestyle in ((0, '-'), (1, '--')):
            selected = sorted((p for p in probes if p['launch_side'] == side), key=lambda p: p['launch_phase'])
            phase = np.array([p['launch_phase'] for p in selected])
            values = np.array([p['complete_return_delay']/length if p['exit_reached_in_interval'] else np.nan for p in selected])
            assert np.allclose(phase, np.arange(25)/8)
            assert all(p['complete_return_delay'] is None for p in selected if not p['exit_reached_in_interval'])
            assert np.all(values[np.isfinite(values)] > shortcut)
            greatest = max(greatest, float(np.nanmax(values)))
            axes[0].plot(phase, values, color=color, ls=linestyle, marker='o', ms=2.8, lw=1.6,
                         label=f'{label}, '+('right' if side == 0 else 'left'))
        diagnostic = case['source_diagnostics']
        time = np.asarray(diagnostic['center_global_time'], dtype=float)
        strength = np.asarray(diagnostic['center_chiral_strength'], dtype=float)
        source = np.asarray(diagnostic['center_chiral_sources_over_alpha'], dtype=float)
        area = np.asarray(diagnostic['center_area_ratio'], dtype=float)
        assert np.isclose(diagnostic['alpha'], alpha) and source.shape == (len(time), 2)
        assert time.shape == strength.shape == area.shape and np.all(np.diff(time) > 0)
        assert np.all(np.isfinite(strength)) and np.all(strength > 0) and np.all(area > 0)
        assert np.allclose(strength, alpha*np.sum(np.abs(source), axis=1)/area, rtol=1e-12, atol=0)
        axes[1].semilogy(time, strength, color=color, lw=1.7, label=label)
    ax = axes[0]
    ax.axhline(1, color='#64748b', ls='--', lw=1, label='Unperturbed loop = 1')
    ax.axhline(shortcut, color='#b7791f', ls=':', lw=1.3, label=f'Exterior-only round trip = {shortcut:.3f}')
    ax.axhline(0, color='#64748b', lw=.8)
    ax.set(title='A. Return delays for matched launches', xlabel='Launch phase s / L',
           ylabel='Complete return delay C / L', xlim=(0, 3), ylim=(-.10, greatest*1.12))
    ax.legend(loc='lower left', ncol=2, fontsize=7.6, framealpha=.95, edgecolor='none', columnspacing=.7)
    ax.text(.5, -.24, 'Gaps: no throat exit in the window. C < 0 would mean a past return.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569')
    ax = axes[1]
    ax.set(title='B. Local source at the throat center', xlabel='Center global time (optical units)',
           ylabel='Local chiral-strength diagnostic (log scale)')
    ax.legend(loc='best', fontsize=8.3, framealpha=.95, edgecolor='none')
    ax.text(.5, -.24, 'Both chiral components evaluated at the same center event.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.2, color='#475569')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .17, f'Full returns within four initial loops: {completed[0]}/50 wide; {completed[1]}/50 reversed. Found exits determine the plotted delays.',
             fontsize=9, color='#475569')
    fig.text(.05, .115, 'Source strength of order 1 invalidates the weak-source approximation. Small center values do not prove global stability.',
             fontsize=9, color='#475569')
    fig.text(.05, .050, 'Optical model units only. No 4D source, stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.1, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'long_packet_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def throat_strip_progress():
    """Display saved throat maps; peak annotations use the densest-grid summaries."""
    from matplotlib.colors import LogNorm
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['throat_strip']
    prepared = []
    for key, label in (('wide', 'Wide packet'), ('wide_reversed', 'Reversed wide packet')):
        diagnostic = record['cases'][key]['strip_diagnostics']
        plot, summary = diagnostic['plot'], diagnostic['summary']
        time = np.asarray(plot['global_time'], dtype=float)
        chi = np.asarray(plot['chi'], dtype=float)
        strength = np.asarray(plot['chiral_strength'], dtype=float)
        phi = np.asarray(plot['phi_over_alpha'], dtype=float)
        alpha = float(diagnostic['alpha'])
        assert np.isclose(alpha, 1e-9) and np.all(np.diff(time) > 0) and np.all(np.diff(chi) > 0)
        assert strength.shape == phi.shape == (len(time), len(chi))
        assert np.all(np.isfinite(strength)) and np.all(strength > 0) and np.all(np.isfinite(phi))
        assert np.all(1+alpha*phi > 0) and summary['minimum_sampled_area_ratio'] > 0
        assert np.allclose(chi[[0, -1]], [-np.arctan(4), np.arctan(4)], rtol=0, atol=1e-12)
        center = int(np.argmin(np.abs(chi)))
        assert abs(chi[center]) < 1e-12
        center_peak = float(summary['maximum_sampled_center_strength'])
        strip_peak = float(summary['maximum_sampled_strip_strength'])
        ratio = float(summary['strip_to_center_peak_ratio'])
        assert strip_peak >= strength.max()*(1-1e-12) and center_peak >= strength[:, center].max()*(1-1e-12)
        assert center_peak > 0 and np.isclose(ratio, strip_peak/center_peak, rtol=1e-12)
        if not diagnostic['sampling_converged']:
            label += ' (peak unresolved)'
        prepared.append((label, time, chi, strength, center_peak, strip_peak, ratio))
    assert np.array_equal(prepared[0][1], prepared[1][1]) and np.array_equal(prepared[0][2], prepared[1][2])
    reversed_diagnostic = record['cases']['wide_reversed']['strip_diagnostics']
    peak_change_percent = 100*reversed_diagnostic['summary']['peak_relative_change']
    envelope = reversed_diagnostic['source_envelope']['resolution_checks'][-1]['chiral_strength_envelope']
    assert not reversed_diagnostic['sampling_converged'] and peak_change_percent > .5
    assert reversed_diagnostic['source_envelope']['numerically_converged'] and envelope > prepared[1][5]
    norm = LogNorm(vmin=min(row[3].min() for row in prepared), vmax=max(row[5] for row in prepared))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.1))
    fig.subplots_adjust(left=.072, right=.882, top=.76, bottom=.36, wspace=.24)
    fig.suptitle('Source strength across a finite throat region', fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Fixed reduced geometry • Static observers • α = 10⁻⁹ • |χ| ≤ arctan 4',
             ha='center', fontsize=10.5, color='#475569')
    for ax, (label, time, chi, strength, center_peak, strip_peak, ratio) in zip(axes, prepared):
        mesh = ax.pcolormesh(time, chi, strength.T, shading='nearest', cmap='viridis', norm=norm, rasterized=True)
        ax.set(title=label, xlabel='Global time τ (optical units)', ylabel='Throat coordinate χ',
               xlim=(time[0], time[-1]), ylim=(chi[0], chi[-1]))
        ax.axhline(0, color='white', lw=.6, ls=':', alpha=.5)
        ax.text(.5, -.22, f'Sampled maxima (densest grid):\n'
                f'center {center_peak:.2e}; strip {strip_peak:.2e}; ratio {ratio:.2f}', transform=ax.transAxes,
                ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.5)
    cax = fig.add_axes([.914, .36, .016, .40])
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label('Chiral-strength diagnostic (log scale)', fontsize=9.5)
    colorbar.ax.tick_params(labelsize=8.5)
    fig.text(.05, .195, 'Downsampled heatmaps; sampled maxima use the densest saved grids, not display pixels.',
             fontsize=9, color='#475569')
    fig.text(.05, .150, f'Reversed peak unresolved: {peak_change_percent:.4f}% change versus 0.5% target. A separate source envelope was checked.',
             fontsize=9, color='#9a6700')
    fig.text(.05, .105, f'Conditional reversed envelope: {envelope:.3e}; a numerical estimate, not a certified continuum bound or stability proof.',
             fontsize=9, color='#475569')
    fig.text(.05, .050, 'No 4D source, stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'throat_strip_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def energy_sweep_progress():
    """Compare saved fixed-geometry energy trials and flag uncontrolled late output."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['energy_sweep']
    groups = [sorted((row for row in record['rows'] if row['orientation'] == orientation),
                     key=lambda row: row['energy_ratio']) for orientation in ('wide', 'reversed')]
    assert all(len(group) == 4 for group in groups)
    energy = np.array([row['energy_ratio'] for row in groups[0]])
    assert np.all(np.diff(energy) > 0) and np.allclose(energy, [row['energy_ratio'] for row in groups[1]])
    for group in groups:
        for row in group:
            assert row['total_launches'] == 50 and 0 <= row['complete_returns_within_interval'] <= 50
            assert row['common_launch_return_over_loop'] > 0 and row['maximum_sampled_center_strength'] > 0
            onset = row['first_sampled_order_one_exterior_phase']
            assert onset is None or row['common_launch_return_over_loop'] < onset
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.7))
    fig.subplots_adjust(left=.061, right=.968, top=.74, bottom=.34, wspace=.43)
    fig.suptitle('Energy sweep: return counts fall as source strength grows',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Fixed exterior delay D = 1.2 • 50 matched launches per state • Declared α = 10⁻⁹',
             ha='center', fontsize=10.5, color='#475569')
    fields = ['complete_returns_within_interval', 'common_launch_return_over_loop', 'maximum_sampled_center_strength']
    for group, label, color, marker, size, linestyle in zip(groups, ['Wide', 'Reversed'],
            ['#007e87', '#7951b0'], ['o', 's'], [45, 24], ['-', '--']):
        for column, (ax, field) in enumerate(zip(axes, fields)):
            values = np.array([row[field] for row in group], dtype=float)
            assert np.all(np.isfinite(values))
            ax.plot(energy, values, color=color, ls=linestyle, lw=1.4, label=label)
            for x, value, row in zip(energy, values, group):
                formal = column != 1 and row['first_sampled_order_one_exterior_phase'] is not None
                ax.scatter(x, value, marker=marker, s=size, facecolor='white' if formal else color,
                           edgecolor=color, linewidth=1.2, zorder=4)
    ax = axes[0]
    for index, x in enumerate(energy):
        counts = [group[index]['complete_returns_within_interval'] for group in groups]
        ax.text(x, max(counts)+2.0, f'{counts[0]}, {counts[1]}', ha='center', fontsize=8.3)
    ax.set(title='A. Full returns in the test window', ylabel='Complete returns out of 50', ylim=(-2, 56))
    ax.legend(loc='upper right', fontsize=8.3, frameon=False)
    ax.text(.5, -.24, 'Labels: wide, reversed. Missing exits are censored.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.0, color='#475569')
    ax = axes[1]
    ax.axhline(1, color='#64748b', ls=':', lw=1.1)
    values = [row[fields[1]] for group in groups for row in group]
    ax.set(title='B. Same early right-side launch', ylabel='Return delay C / initial loop L',
           ylim=(min(1, min(values))-.015, max(values)+.025))
    ax.text(.5, -.24, 'Launch phase 0; before sampled order-one onset.\nDotted line: unperturbed loop = 1; axis is enlarged.',
            transform=ax.transAxes, ha='center', va='top', fontsize=7.8, color='#475569', linespacing=1.4)
    ax = axes[2]
    ax.set_yscale('log')
    ax.axhline(1, color='#b7791f', ls=':', lw=1.3)
    ax.set(title='C. Local source-control diagnostic', ylabel='Maximum sampled center strength')
    ax.text(.03, 1, 'Screen: S = 1', transform=ax.get_yaxis_transform(),
            va='bottom', fontsize=8.4, color='#9a6700')
    ax.text(.5, -.24, 'Observer-dependent screen, not a 4D residual test.\nSmall center values alone do not prove stability.',
            transform=ax.transAxes, ha='center', va='top', fontsize=7.8, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.set(xlabel='Initial energy ratio I / b₀', xlim=(energy[0]-.18, energy[-1]+.18))
        ax.set_xticks(energy, labels=[f'{value:.3f}' if index == 0 else f'{value:g}' for index, value in enumerate(energy)])
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.05, .16, 'Open late-run points fail the assumed weak-source screen; those late counts remain formal model outputs.',
             fontsize=9, color='#9a6700')
    fig.text(.05, .105, 'All compared early return delays remain positive. Counts concern a finite window; times use optical model units.',
             fontsize=9, color='#475569')
    fig.text(.05, .045, 'No 4D source, stable medium or 13 → 12 Sep 2027 passage established.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'energy_sweep_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def center_residual_progress():
    """Contrast center stress with specified equation residuals, not metric error."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['center_residual']
    rows = {row['case']: row for row in record['rows']}
    closure = record['conditional_polyakov']
    assert set(rows) == {'3_wide', '3_reversed', '4_wide', '4_reversed'}
    assert np.isclose(record['parameters']['alpha'], 1e-9)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 6.1))
    fig.subplots_adjust(left=.055, right=.968, top=.74, bottom=.32, wspace=.46)
    fig.suptitle('Large stress does not imply a large equation residual',
                 fontsize=18, fontweight='bold', y=.965)
    fig.text(.5, .88, 'Conditional reconstruction • Sampled center events only • Four initial loops • α = 10⁻⁹',
             ha='center', fontsize=10.5, color='#475569')
    for key, row in rows.items():
        case = record['cases'][key]
        phase = np.asarray(case['exterior_phase'], dtype=float)
        strength = np.asarray(case['center_chiral_strength'], dtype=float)
        residual = np.asarray(case['absolute_observer_residual'], dtype=float)
        phi = np.asarray(case['center_phi'], dtype=float)
        gradient = np.asarray(case['center_phi_u_v'], dtype=float)
        assert np.all(np.diff(phase) > 0) and np.allclose(phase[[0, -1]], [0, 4])
        assert strength.shape == phi.shape == phase.shape and residual.shape == gradient.shape == (len(phase), 2)
        assert np.all(np.isfinite(residual)) and np.all(residual >= 0) and np.all(strength > 0)
        assert np.all(1+phi > 0) and np.allclose(residual, 2*gradient**2/(1+phi[:, None])**2, rtol=1e-12, atol=0)
        null_max = np.max(residual, axis=1)
        assert np.isclose(null_max.max(), row['maximum_absolute_radial_null_residual'], rtol=1e-12, atol=0)
        assert np.isclose(strength.max(), row['maximum_sampled_center_strength'], rtol=1e-12)
        for field, summary_field in (('mixed_observer_residual', 'maximum_absolute_mixed_residual'),
                ('angular_residual', 'maximum_absolute_angular_residual'), ('angular_pulse_minus_static', 'maximum_absolute_angular_pulse_change')):
            values = np.asarray(case['conditional_polyakov'][field])
            assert values.shape == phase.shape and np.isclose(np.max(np.abs(values)), row[summary_field], rtol=1e-12, atol=0)
        color = '#007e87' if key.startswith('3_') else '#7951b0'
        style = '--' if key.endswith('reversed') else '-'
        label = f"η = {row['energy_ratio']:g}, "+('reversed' if key.endswith('reversed') else 'wide')
        axes[1].semilogy(phase, strength, color=color, ls=style, lw=1.6, label=label)
        axes[2].semilogy(phase, np.ma.masked_equal(null_max, 0), color=color, ls=style, lw=1.6)
    ax = axes[0]
    ax.set(xlim=(0, 1), ylim=(0, 1), title='A. Center probe in an assumed throat'); ax.axis('off')
    x = np.linspace(.03, .97, 101); radius = .07+.11*(2*x-1)**2
    ax.fill_between(x, .77-radius, .77+radius, color='#e2e8f0')
    ax.plot(x, .77-radius, color='#94a3b8', lw=1); ax.plot(x, .77+radius, color='#94a3b8', lw=1)
    ax.scatter([.5], [.77], color='#007e87', s=34, zorder=3)
    ax.text(.5, .97, 'Schematic — no physical medium', ha='center', va='top', fontsize=8.5, color='#475569')
    ax.annotate('Same event at χ = 0', xy=(.5, .755), xytext=(.5, .55), ha='center', fontsize=9,
                color='#007e87', arrowprops={'arrowstyle': '-', 'color': '#007e87', 'lw': .8})
    ax.text(.5, .43, 'Specified Polyakov closure\nSampled absolute maxima', ha='center', fontsize=9,
            fontweight='bold', color='#475569', linespacing=1.4)
    table = [('Mixed', max(row['maximum_absolute_mixed_residual'] for row in rows.values())),
             ('Static mixed calibration', abs(closure['static_mixed_residual'])),
             ('Angular', max(row['maximum_absolute_angular_residual'] for row in rows.values())),
             ('Angular pulse − static', max(row['maximum_absolute_angular_pulse_change'] for row in rows.values()))]
    for y, (label, value) in zip([.31, .22, .13, .04], table):
        ax.text(.01, y, label, ha='left', fontsize=8.3, color='#475569')
        ax.text(.99, y, f'{value:.2e}', ha='right', fontsize=9, color='#007e87')
    ax.text(.5, -.22, 'Mixed defect is already present in the static calibration.',
            transform=ax.transAxes, ha='center', va='top', fontsize=7.9, color='#475569')
    axes[1].set(title='B. Center stress diagnostic S', ylabel='S (log scale)')
    axes[1].legend(loc='upper left', fontsize=7.6, framealpha=.95, edgecolor='none')
    axes[2].set(title='C. Exact radial-null residual', ylabel='Max of two absolute observer components', ylim=(1e-22, 1e-10))
    for ax in axes[1:]:
        ax.set(xlabel='Exterior phase s / L', xlim=(0, 4))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    peak_stress = max(row['maximum_sampled_center_strength'] for row in rows.values())
    peak_null = max(row['maximum_absolute_radial_null_residual'] for row in rows.values())
    axes[1].text(.5, -.22, f'Largest sampled S: {peak_stress:.2f}; no validity cutoff imposed.',
                 transform=axes[1].transAxes, ha='center', va='top', fontsize=7.9, color='#475569')
    axes[2].text(.5, -.22, f'Largest sampled residual: {peak_null:.3e}.\nNear-zero start below plotted range; zero omitted.',
                 transform=axes[2].transAxes, ha='center', va='top', fontsize=7.9, color='#475569', linespacing=1.4)
    fig.text(.045, .145, r'Residuals are dimensionless $r_e^2(G-8\pi G T)$ with specified observer contractions, not metric-error percentages.',
             fontsize=9, color='#475569')
    fig.text(.045, .095, 'Reconstructed metric and specified matter closure only; residuals are evaluated, not solved away. No global stability proof.',
             fontsize=9, color='#475569')
    fig.text(.045, .040, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'center_residual_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def spatial_residual_progress():
    """Plot finite model coverage and residual envelopes conditional on numeric jets."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['spatial_residual']
    assert set(record['cases']) == {'4_wide', '4_reversed'}
    assert np.isclose(record['parameters']['alpha'], 1e-9)
    assert record['verification']['unresolved_sampled_peak_checks'] == ['4_wide']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.08, right=.975, top=.76, bottom=.34, wspace=.30)
    fig.suptitle('Extending the equation check beyond the center', fontsize=18, fontweight='bold', y=.955)
    fig.text(.5, .875, 'Finite interior only: |\u03c7| \u2264 arctan 4  \u2022  \u03b7 = 4  \u2022  \u03b1 = 10\u207b\u2079  \u2022  Specified Polyakov closure',
             ha='center', fontsize=10.5, color='#475569')
    labels, envelope_peaks, time_notes = [], [], []
    for index, (key, label, color, style) in enumerate((('4_wide', 'Wide', '#007e87', '-'),
                                                       ('4_reversed', 'Reversed', '#7951b0', '--'))):
        case = record['cases'][key]; summary = case['summary']
        chi = np.asarray(case['domain_plot']['chi']); maximum = np.asarray(case['domain_plot']['maximum_global_time'])
        rmax, lmax = case['entry_clock_upper_limits']
        rectangle = case['full_width_rectangle_final_time']; center = case['latest_center_global_time']
        envelope = np.asarray(summary['absolute_residual_envelope']); sampled = np.asarray(summary['sampled_absolute_residual_maxima'])
        assert case['components'] == ['null_plus', 'null_minus', 'mixed', 'angular']
        assert chi.shape == maximum.shape and np.all(np.diff(chi) > 0)
        assert np.allclose(chi[[0, -1]], [-np.arctan(4), np.arctan(4)], rtol=0, atol=1e-14)
        assert np.allclose(maximum, np.pi/2+np.minimum(rmax-chi, lmax+chi), rtol=0, atol=1e-13)
        assert np.isclose(rectangle, maximum.min()) and np.isclose(center, np.pi/2+min(rmax, lmax))
        assert center > rectangle > 0 and envelope.shape == sampled.shape == (4,)
        assert np.all(np.isfinite(envelope)) and np.all(envelope >= sampled) and np.all(sampled > 0)
        assert summary == case['resolution_checks'][-1] and case['numerical_envelope_converged']
        assert case['sampled_peaks_converged'] == (key == '4_reversed')
        assert np.isclose(summary['minimum_area_ratio'], 1+summary['phi_interval'][0])
        axes[0].fill_between(chi, 0, maximum, color=color, alpha=.04)
        axes[0].fill_between(chi, rectangle, maximum, color=color, alpha=.15)
        axes[0].plot(chi, maximum, color=color, ls=style, lw=2, label=label)
        axes[0].hlines(rectangle, chi[0], chi[-1], color=color, ls=':', lw=1.5)
        axes[0].scatter([0], [center], color=color, s=35, zorder=4)
        position = np.arange(4)+(index-.5)*.27
        axes[1].bar(position, envelope-1e-16, bottom=1e-16, width=.23, facecolor='none', edgecolor=color, lw=1.5)
        axes[1].scatter(position, sampled, color=color, s=34, zorder=4)
        labels.append(label); envelope_peaks.append(envelope.max())
        time_notes.append(f'{label}: {rectangle:.3f} \u2192 {center:.3f}')
    axes[0].axvline(0, color='#94a3b8', ls=':', lw=.8, zorder=0)
    axes[0].plot([], [], color='#64748b', ls=':', label='Rectangle end')
    axes[0].plot([], [], color='#64748b', marker='o', ls='none', label='Latest center event')
    axes[0].set(title='A. Available model domain', xlabel='Interior position \u03c7', ylabel='Global time \u03c4 (optical units)',
                xlim=(-np.arctan(4), np.arctan(4)), ylim=(0, 4.65))
    axes[0].legend(loc='upper right', fontsize=8, framealpha=.95, edgecolor='none')
    axes[0].text(.5, -.19, 'A rectangle misses the darker late region.\nAt the center: '+ '; '.join(time_notes),
                 transform=axes[0].transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    axes[1].set(title='B. Sampled maxima and conditional envelopes', yscale='log', ylim=(1e-16, 1e-6),
                ylabel='Absolute observer residual', xticks=np.arange(4), xticklabels=['Null +', 'Null \u2212', 'Mixed', 'Angular'])
    axes[1].text(.035, .975, 'Open bars: conditional envelopes\nDots: sampled maxima*',
                 transform=axes[1].transAxes, va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    axes[1].text(.5, -.19, 'Largest component envelopes\n'+ '  /  '.join(f'{label}: {value:.3e}' for label, value in zip(labels, envelope_peaks)),
                 transform=axes[1].transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13); ax.set_axisbelow(True)
    change = record['cases']['4_wide']['summary']['sampled_peak_relative_change']
    fig.text(.065, .18, f'* Wide peak scan unresolved: last change {100*change:.3f}% exceeds the 0.5% target. Dots are not certified true peaks.',
             fontsize=9, color='#475569')
    fig.text(.065, .135, 'Analytic inequalities use numerical maxima of six field quantities; these are not certified continuum bounds.',
             fontsize=9, color='#475569')
    fig.text(.065, .09, r'Residuals: $r_e^2(G-8\pi G T)$, specified contractions. Hypothetical geometry; near-mouth/exterior regions and stability unresolved.',
             fontsize=9, color='#475569')
    fig.text(.065, .035, 'No constructed medium or 13 \u2192 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'spatial_residual_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def coupled_response_progress():
    """Keep the local coupled benchmark distinct from conditional exit kinematics."""
    root = Path(__file__).resolve().parent
    response = json.loads((root/'polyakov_response_results.json').read_text(encoding='utf-8'))
    timing = json.loads((root/'mmp_pulse_results.json').read_text(encoding='utf-8'))['timing_response']
    data = response['plot']; parameters = response['parameters']
    time = np.asarray(data['center_time']); radius = np.asarray(data['center_radius_response'])
    conformal = np.asarray(data['center_conformal_response']); scale = parameters['alpha']*parameters['infinitesimal_amplitude_epsilon']
    assert response['acceptance']['all_pass'] and np.isclose(scale, 1e-12, rtol=1e-12, atol=0)
    assert radius.shape == conformal.shape == time.shape and np.all(np.diff(time) > 0)
    assert np.allclose(time[[0, -1]], parameters['null_coordinate_interval'])
    assert np.allclose(radius, np.diag(data['radius_response_over_alpha_epsilon']), rtol=1e-12, atol=0)
    assert np.allclose(conformal, np.diag(data['conformal_response_over_alpha_epsilon']), rtol=1e-12, atol=0)
    assert np.all(np.isfinite(radius)) and np.all(np.isfinite(conformal)) and radius[0] == conformal[0] == 0
    assert np.all((np.asarray(response['observed_refinement_ratios']) > 3.5) & (np.asarray(response['observed_refinement_ratios']) < 4.5))
    assert np.isclose(response['response']['maximum_conformal_log_change'], scale*response['refinement'][-1]['maximum_normalized_conformal_response'], rtol=1e-12, atol=0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4), gridspec_kw={'width_ratios': [1.05, 1]})
    fig.subplots_adjust(left=.085, right=.975, top=.755, bottom=.345, wspace=.26)
    fig.suptitle('A first local coupled geometry response', fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Two independent diagnostics \u2022 No matched physical mouth clocks', ha='center', fontsize=11, color='#475569')
    ax = axes[0]
    ax.axvspan(*parameters['pulse_u_interval'], color='#cbd5e1', alpha=.45, label='Prescribed incoming pulse')
    ax.plot(time, radius, color='#007e87', lw=2, label=r'Radius: $\delta R/(\alpha\epsilon)$')
    ax.plot(time, conformal, color='#7951b0', ls='--', lw=2, label=r'Metric: $\delta\sigma/(\alpha\epsilon)$')
    ax.axhline(0, color='#64748b', lw=.65)
    ax.set(title='A. Local linear benchmark', xlabel='Center time t = u = v (model units)', ylabel='Normalized field response', xlim=(0, 1.5))
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.13); ax.set_axisbelow(True)
    ax.legend(loc='best', fontsize=8, framealpha=.95, edgecolor='none')
    ax.text(.5, -.205, '\u03b1\u03b5 = 10\u207b\u00b9\u00b2; curves show fields, not area.\nLocal linear response; no nonlinear finite-amplitude test.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.45)
    ax = axes[1]; ax.set(xlim=(0, 1), ylim=(0, 1), title='B. Conditional exit sensitivity'); ax.axis('off')
    ax.text(.5, .955, 'Different state: \u03b7 = 4; completed side-0, phase-0 exit', ha='center', fontsize=8.7, color='#475569')
    ax.text(.015, .705, 'Receiving clock rate', fontsize=9.3, color='#475569')
    ax.text(.015, .555, 'Timing gain = 1 / rate', fontsize=9.3, color='#475569')
    ax.text(.5, .375, 'Assumed combined phase error: \u00b10.001', ha='center', fontsize=9.3, fontweight='bold', color='#475569')
    ax.text(.015, .22, 'Exit shift lower bound', fontsize=9.3, color='#475569')
    ax.text(.015, .075, 'Exit shift upper bound', fontsize=9.3, color='#475569')
    for key, label, color, xpos in [('4_wide', 'Wide', '#007e87', .64), ('4_reversed', 'Reversed', '#7951b0', .91)]:
        completed = [probe for probe in timing['cases'][key]['probes'] if probe['exit_reached_in_interval']]
        assert len(completed) == 1
        probe = completed[0]
        assert probe['launch_side'] == 0 and probe['launch_phase'] == 0 and probe['complete_return_delay'] > 0
        gain = probe['linear_time_gain_per_phase']; rate = probe['receiving_clock_rate']
        assert np.isclose(gain*rate, 1, rtol=1e-13)
        bounds = [row for row in probe['conditional_response_intervals'] if row['assumed_combined_phase_error'] == .001]
        assert len(bounds) == 1
        interval = bounds[0]; low, high = interval['exit_shift_interval']
        assert low < 0 < high and np.allclose(np.asarray(interval['exit_interval'])-probe['exit_exterior_time'], [low, high], rtol=1e-12, atol=1e-14)
        ax.text(xpos, .83, label, ha='center', fontsize=10, fontweight='bold', color=color)
        for ypos, value in [(.705, f'{rate:.5f}'), (.555, f'{gain:.5f}'), (.22, f'{low:+.8f}'), (.075, f'{high:+.8f}')]:
            ax.text(xpos, ypos, value, ha='center', fontsize=9.5, color=color)
    ax.hlines([.785, .46, .315], .015, .99, color='#cbd5e1', lw=.7)
    ax.text(.5, -.205, 'Exit shifts are in exterior optical units, not seconds.\nAn earlier exit shift is not a past return.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.45)
    fig.text(.065, .175, 'No mouth-matching calculation connects A to B. The phase bound in B is assumed, not predicted by A.', fontsize=9.2, color='#475569')
    fig.text(.065, .12, 'Finite local region and specified Polyakov matter; exterior completion, physical sources and stability remain unresolved.', fontsize=9.2, color='#475569')
    fig.text(.065, .045, 'No constructed medium or 13 \u2192 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'coupled_response_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def finite_energy_response_progress():
    """Show a finite-state local benchmark without implying throat formation."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_finite_results.json').read_text(encoding='utf-8'))
    parameters = record['parameters']; case = record['physical_alpha']; data = case['plot']
    fine = case['refinement'][-1]; clock = fine['local_clock']; alpha = parameters['alpha']
    u = np.asarray(record['plot']['source_u']); source = np.asarray(record['plot']['delta_t_u'])
    time = np.asarray(data['center_time']); radius = np.asarray(data['center_radius_response'])
    conformal = np.asarray(data['center_conformal_response'])
    linear_radius = np.asarray(data['linear_geometry_center_radius_response'])
    linear_conformal = np.asarray(data['linear_geometry_center_conformal_response'])
    assert record['acceptance']['all_pass'] and np.isclose(alpha, 1e-9, rtol=1e-12, atol=0)
    assert u.shape == source.shape and np.all(np.diff(u) > 0) and np.all(np.isfinite(source))
    assert radius.shape == conformal.shape == linear_radius.shape == linear_conformal.shape == time.shape
    assert np.all(np.diff(time) > 0) and np.allclose(time[[0, -1]], parameters['null_coordinate_interval'])
    assert np.allclose(radius, np.diag(data['radius_response_over_alpha']), rtol=1e-12, atol=0)
    assert np.allclose(conformal, np.diag(data['conformal_response_over_alpha']), rtol=1e-12, atol=0)
    assert parameters['minimum_finite_map_derivative'] > 0 and source.min() < 0 < source.max()
    integral = parameters['integrated_delta_t']
    assert integral > 0 and np.isclose(integral, 1.5*parameters['initial_energy_helper_normalization'], rtol=1e-12)
    assert np.isclose(integral, parameters['independent_integrated_delta_t'], rtol=1e-12)
    assert np.isclose(clock['proper_duration_change'], alpha*clock['proper_duration_change_over_alpha'], rtol=1e-12, atol=0)
    assert clock['static_proper_duration']+clock['proper_duration_change'] > 0
    assert fine['normalized_geometry_nonlinear_remainder'] < case['common_grid_state_differences'][-1]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(13.7, 6.4), gridspec_kw={'width_ratios': [1, 1, .70]})
    fig.subplots_adjust(left=.06, right=.97, top=.75, bottom=.35, wspace=.36)
    fig.suptitle('A finite-energy pulse with coupled geometry response', fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Conditional Einstein–Maxwell–Polyakov model • α = 10⁻⁹ • A pre-existing throat is assumed',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]
    ax.fill_between(u, 0, source, where=source >= 0, color='#007e87', alpha=.16)
    ax.fill_between(u, 0, source, where=source < 0, color='#7951b0', alpha=.15)
    ax.plot(u, source, color='#007e87', lw=1.7); ax.axhline(0, color='#64748b', lw=.7)
    ax.set(title='A. Exact incoming state source', xlabel='Incoming null coordinate u', ylabel=r'State-source excess $\delta t_u$', xlim=(0, 1.5))
    ax.text(.5, -.23, f'Integrated excess = {integral:.3f} > 0.\nModel normalization; local lobes may be negative.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.plot(time, radius, color='#007e87', lw=2, label=r'Radius $\delta R/\alpha$')
    ax.plot(time, conformal, color='#7951b0', lw=2, label=r'Metric $\delta\sigma/\alpha$')
    ax.plot(time, linear_radius, color='#334155', ls='--', lw=.8, alpha=.8, label='Same-source linear geometry')
    ax.plot(time, linear_conformal, color='#334155', ls='--', lw=.8, alpha=.8)
    ax.axhline(0, color='#64748b', lw=.6)
    ax.set(title='B. Center geometry response', xlabel='Center time t = u = v', ylabel='Response / α (radius is not area)', xlim=(0, 1.5))
    ax.legend(loc='lower left', fontsize=7.6, framealpha=.95, edgecolor='none')
    ax.text(.5, -.23, 'Full and linear curves overlap at this scale.\nTheir small difference is below the mesh difference.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color='#475569', linespacing=1.4)
    ax = axes[2]
    ax.fill([0, .75, 0, -.75], [0, .75, 1.5, .75], color='#e2e8f0', alpha=.65)
    ax.plot([.2, .2], [.2, 1.3], color='#7951b0', ls=':', lw=1.3)
    ax.plot([-.2, -.2], [.2, 1.3], color='#94a3b8', ls=':', lw=1)
    ax.plot([.2, -.2, .2], [.4, .8, 1.2], color='#007e87', lw=1.9, marker='o', markersize=4)
    for start, end in [((.12, .48), (-.06, .66)), ((-.12, .88), (.06, 1.06))]:
        ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle': '->', 'color': '#007e87', 'lw': 1.4})
    for x, t, label, align in [(.25, .4, 'Launch 0.4', 'left'), (-.26, .8, 'Reflect\n0.8', 'right'), (.25, 1.2, 'Return 1.2', 'left')]:
        ax.text(x, t, label, ha=align, va='center', fontsize=7.8, color='#475569')
    ax.set(title='C. Prescribed local light trip', xlabel='Position x', ylabel='Time t', xlim=(-.8, .9), ylim=(0, 1.5))
    ax.text(.5, -.23, f"Static duration: {clock['static_proper_duration']:.6f} rₑ/c\nChange: {clock['proper_duration_change']:+.3e} rₑ/c; total stays positive.",
            transform=ax.transAxes, ha='center', va='top', fontsize=8.1, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.12); ax.set_axisbelow(True)
    fig.text(.05, .17, 'This evolves an assumed local throat; it does not form one from ordinary space or supply a physical preparation method.', fontsize=9, color='#475569')
    fig.text(.05, .115, 'C uses one prescribed emitter clock and an ideal reflection. The trip is causal by construction; no physical mouth matching or global stability.', fontsize=8.8, color='#475569')
    fig.text(.05, .04, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'finite_energy_response_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def radius_boundary_progress():
    """Track local constant-area levels without presenting them as material mouths."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_boundary_results.json').read_text(encoding='utf-8'))
    case = record['refinement'][-1]; parameters = record['parameters']; probe = case['probe']
    assert record['verification']['all_declared_checks_pass'] and case['intervals'] == 1024
    assert set(case['sides']) == {'-1', '1'} and np.isclose(parameters['alpha'], 1e-9, rtol=1e-12, atol=0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.085, right=.97, top=.76, bottom=.345, wspace=.27)
    fig.suptitle('Tracking constant-area levels through a quantum pulse', fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Local nonlinear model • A pre-existing throat is assumed • These levels are not material mouths',
             ha='center', fontsize=10.5, color='#475569')
    final_times = []; all_positions = []
    for side, label, color in [('1', 'Right level', '#007e87'), ('-1', 'Left level', '#7951b0')]:
        boundary = case['sides'][side]; data = boundary['plot']
        time = np.asarray(data['time']); position = np.asarray(data['x']); speed = np.asarray(data['velocity'])
        assert time.shape == position.shape == speed.shape and np.all(np.diff(time) > 0)
        assert np.all(np.isfinite(position)) and np.all(np.isfinite(speed))
        assert np.isclose(time[0], parameters['initial_time'])
        assert np.isclose(position[0], int(side)*parameters['cutoff_x'], rtol=0, atol=1e-8)
        assert np.isclose(time[-1], boundary['final_time']) and np.isclose(position[-1], boundary['final_x'])
        assert np.all(np.abs(speed[:-1]) < 1+1e-8) and np.isclose(abs(speed[-1]), 1, rtol=0, atol=1e-6)
        axes[0].plot(position, time, color=color, lw=2, label=label)
        axes[0].scatter(position[[0, -1]], time[[0, -1]], s=[20, 55], color=color, zorder=4)
        axes[1].plot(time, speed, color=color, lw=2, label=label)
        axes[1].scatter([time[-1]], [speed[-1]], s=55, color=color, zorder=4)
        axes[1].annotate(f'Null at t={time[-1]:.4f}', xy=(time[-1], speed[-1]),
                         xytext=(-6, -19 if speed[-1] > 0 else 45), textcoords='offset points',
                         ha='right', fontsize=8.2, color=color)
        final_times.append(time[-1]); all_positions.extend(position)
    if probe['return_time'] is not None:
        ray_time = np.array([probe['launch_time'], probe['reflection_time'], probe['return_time']])
        ray_x = np.array([probe['launch_x'], probe['reflection_x'], probe['return_x']])
        assert np.all(np.diff(ray_time) > 0) and probe['proper_round_trip_duration'] > 0
        assert probe['return_time'] < case['sides']['1']['final_time']
        assert np.allclose(np.abs(np.diff(ray_x)), np.diff(ray_time), rtol=0, atol=1e-10)
        axes[0].plot(ray_x, ray_time, color='#475569', ls='--', lw=1.3, marker='o', markersize=3, label='Ideal test ray')
        for start, end in zip(np.c_[ray_x, ray_time][:-1], np.c_[ray_x, ray_time][1:]):
            axes[0].annotate('', xy=start+.65*(end-start), xytext=start+.35*(end-start),
                             arrowprops={'arrowstyle': '->', 'color': '#475569', 'lw': 1.2})
        clock_note = f"Completed local return: proper duration {probe['proper_round_trip_duration']:.6f} rₑ/c > 0."
    else:
        clock_note = 'No return before the tracked right level reaches its null limit.'
    span = max(abs(np.asarray(all_positions)))*1.22
    axes[0].set(title='A. Continuously tracked local levels', xlabel='Position x', ylabel='Coordinate time t (model units)',
                xlim=(-span, span), ylim=(parameters['initial_time']-.015, max(final_times)+.025))
    axes[0].legend(loc='best', fontsize=8, framealpha=.95, edgecolor='none')
    axes[0].text(.5, -.205, f"Initialized at x = ±{parameters['cutoff_x']:g}, t = {parameters['initial_time']:g}.\n"+clock_note,
                 transform=axes[0].transAxes, ha='center', va='top', fontsize=8.4, color='#475569', linespacing=1.4)
    axes[1].axhspan(-1, 1, color='#e2e8f0', alpha=.24)
    for limit in [-1, 1]:
        axes[1].axhline(limit, color='#64748b', ls=':', lw=1.1)
    axes[1].set(title='B. Timelike levels approach the null limit', xlabel='Coordinate time t (model units)', ylabel='Level velocity dx / dt',
                ylim=(-1.14, 1.14), xlim=(parameters['initial_time']-.008, max(final_times)+.015))
    axes[1].legend(loc='center left', fontsize=8, framealpha=.95, edgecolor='none')
    axes[1].text(.5, -.205, 'Dotted limits: dx/dt = ±1. Timelike levels lie between.\nThis tracks area levels, not trajectories of material objects.',
                 transform=axes[1].transAxes, ha='center', va='top', fontsize=8.4, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.06, .17, 'The plotted branches stop at their first null limit; this does not establish an event horizon or a physical entrance.', fontsize=9.1, color='#475569')
    fig.text(.06, .115, 'Local prescribed source and ideal reflection only. Exterior matching, preparation and global stability remain unresolved.', fontsize=9.1, color='#475569')
    fig.text(.06, .045, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'radius_boundary_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def source_retaining_junction_progress():
    """Display source-retaining boundary derivatives and their numerical checks."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_junction_results.json').read_text(encoding='utf-8'))
    parameters = record['parameters']; fine = record['refinement'][-1]
    assert record['verification']['all_declared_checks_pass'] and fine['intervals'] == 1024
    assert np.isclose(parameters['alpha'], 1e-9, rtol=1e-12, atol=0)
    assert set(fine['sides']) == {'-1', '1'}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.085, right=.97, top=.755, bottom=.345, wspace=.28)
    fig.suptitle('The boundary source requires a radial lapse gradient', fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Conditional local model • Pre-existing throat assumed • Boundary values and first derivatives only',
             ha='center', fontsize=10.5, color='#475569')
    for side, label, color, style in [('1', 'Right: x = +0.05', '#007e87', '-'),
                                      ('-1', 'Left: x = −0.05', '#7951b0', '--')]:
        rows = fine['sides'][side]
        time = np.array([row['time'] for row in rows]); values = np.array([row['lapse_radial_derivative_times_alpha'] for row in rows])
        assert np.all(np.diff(time) > 0) and np.allclose(time[[0, -1]], parameters['time_interval'])
        assert np.all(np.isfinite(values)) and np.all(values < 0)
        assert np.allclose([row['x'] for row in rows], int(side)*parameters['x_cutoff'], rtol=0, atol=1e-14)
        axes[0].plot(time, values, color=color, ls=style, lw=2, label=label)
    control = record['zero_source_control']['sides']['1']
    control_time = np.array([row['time'] for row in control]); control_values = np.array([row['lapse_radial_derivative_times_alpha'] for row in control])
    assert np.allclose(control_values, control_values[0], rtol=0, atol=1e-10) and control_values[0] < 0
    axes[0].plot(control_time, control_values, color='#475569', ls=':', lw=1.3, label='Static quantum-state control')
    axes[0].axhline(0, color='#64748b', ls='-.', lw=1.2)
    axes[0].text(.5, .955, r'Restricted Vaidya: $\psi_R=0$', transform=axes[0].transAxes,
                 ha='center', va='top', fontsize=9, color='#475569')
    axes[0].set(title='A. Required radial lapse derivative', xlabel='Boundary coordinate time t',
                ylabel=r'$\alpha\,\partial\psi/\partial R$  ($r_e=1$)', xlim=parameters['time_interval'], ylim=(-4200, 550))
    axes[0].legend(loc='lower left', fontsize=8, framealpha=.95, edgecolor='none')
    axes[0].text(.5, -.205, 'A scaled coordinate derivative, not a measured field.\nStatic control overlaps the left-boundary trace.',
                 transform=axes[0].transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    grids = np.array([row['intervals'] for row in record['refinement']])
    for key, label, color, marker in [('quantum_einstein_scaled', 'Quantum–Einstein stress check', '#007e87', 'o'),
                                      ('lapse_source_scaled', 'Radial lapse/source check', '#7951b0', 's')]:
        residual = np.array([row['maximum_checks'][key] for row in record['refinement']])
        assert np.all(residual > 0) and np.all(np.diff(residual) < 0)
        assert np.isclose(residual[-1], record['verification']['source_maximum'][key], rtol=1e-13)
        assert residual[-1] < record['verification']['thresholds']['final_source_scaled']
        axes[1].plot(grids, residual, color=color, marker=marker, lw=1.8, ms=5, label=label)
    axes[1].set_xscale('log', base=2); axes[1].set_yscale('log')
    axes[1].set(title='B. Independent source checks converge', xlabel='Grid intervals per null direction', ylabel='Maximum normalized check residual',
                xticks=grids, xticklabels=[str(n) for n in grids], ylim=(1e-4, .2))
    axes[1].legend(loc='upper right', fontsize=8, framealpha=.95, edgecolor='none')
    errors = fine['maximum_checks']
    axes[1].text(.5, -.205, f"At 1024: stress {errors['quantum_einstein_scaled']:.3e}; lapse {errors['lapse_source_scaled']:.3e}.\nRefinement is not a certified continuum error bound.",
                 transform=axes[1].transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.06, .17, 'R is measured in units of rₑ; α = 10⁻⁹. The plotted gradient is required within this specified source and coordinate chart.', fontsize=9, color='#475569')
    fig.text(.06, .115, 'Matching first derivatives supplies boundary data. The full exterior, its quantum state, preparation and stability remain unsolved.', fontsize=9, color='#475569')
    fig.text(.06, .045, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'source_retaining_junction_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def quantum_collar_progress():
    """Show a finite static continuation while keeping the returning path open."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_collar_results.json').read_text(encoding='utf-8'))
    parameters = record['parameters']; case = record['tighter_integration']; rows = case['plot']
    alpha = parameters['alpha']; radius = np.array([row['radius'] for row in rows])
    correction = np.array([row['mass_excess_over_alpha'] for row in rows])
    leading = -5/(18*radius**3)+1/(3*radius**2)-1/18
    last = case['checkpoints'][-1]; budget = last['remaining_optical_path']
    assert record['verification']['all_declared_checks_pass'] and case['integration_refinement'] == 2
    assert np.isclose(alpha, 1e-9, rtol=1e-12, atol=0) and len(rows) == 161
    assert np.all(np.diff(radius) > 0) and np.allclose(radius[[0, -1]], [1.01, parameters['outer_radius']])
    assert np.all(np.isfinite(correction)) and np.allclose(correction, np.array([r['mass_excess'] for r in rows])/alpha, rtol=1e-12, atol=0)
    assert np.max(np.abs(correction-leading)) < 1e-9
    assert np.isclose(last['radius'], parameters['outer_radius']) and budget > 0
    assert np.isclose(budget, parameters['total_assumed_optical_period']-2*last['optical_x'], rtol=0, atol=1e-14)
    assert case['maximum_constraint_scaled'] < record['verification']['thresholds']['constraint']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.08, right=.97, top=.755, bottom=.345, wspace=.27)
    fig.suptitle('A static exterior with the quantum source retained', fontsize=18, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Existing throat assumed • Same reduced action • Angular equation retained • α = 10⁻⁹',
             ha='center', fontsize=10.5, color='#475569')
    ax = axes[0]
    ax.semilogx(radius, correction, color='#007e87', lw=2.3, label='Static coupled solution')
    ax.semilogx(radius, leading, color='#7951b0', ls='--', lw=1.2, label='Leading small-α formula')
    ax.axhline(0, color='#94a3b8', lw=.7)
    ax.set(title='A. Mass correction along the exterior', xlabel='Areal radius / rₑ (log scale)',
           ylabel=r'$(M-1)/\alpha$', xlim=(1.01, 100), xticks=[1.01, 2, 5, 10, 100], xticklabels=['1.01', '2', '5', '10', '100'])
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.13); ax.set_axisbelow(True)
    ax.legend(loc='upper right', fontsize=8.3, framealpha=.95, edgecolor='none')
    ax.text(.5, -.205, r'$m_1(R)=-5/(18R^3)+1/(3R^2)-1/18$'+'\nCurves overlap at this scale; fixed finite R expansion.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.6, color='#475569', linespacing=1.4)
    ax = axes[1]; ax.set(xlim=(0, 1), ylim=(0, 1), title='B. Finite continuation; return path still open'); ax.axis('off')
    ax.text(.5, .965, 'Region map only — schematic, not to scale', ha='center', va='top', fontsize=8.7, color='#475569')
    for xpos, label, face, edge, style in [(.13, 'Assumed\nstatic throat', '#f1f5f9', '#94a3b8', '-'),
                                           (.50, 'Same-source\nfinite exterior', '#e0f2f1', '#007e87', '-'),
                                           (.87, 'Returning path\nnot solved', 'white', '#7951b0', '--')]:
        ax.text(xpos, .76, label, ha='center', va='center', fontsize=9, color='#334155',
                bbox={'boxstyle': 'round,pad=.5', 'facecolor': face, 'edgecolor': edge, 'linestyle': style, 'linewidth': 1.2})
    ax.annotate('', xy=(.35, .76), xytext=(.27, .76), arrowprops={'arrowstyle': '->', 'color': '#007e87', 'lw': 1.3})
    ax.annotate('', xy=(.72, .76), xytext=(.65, .76), arrowprops={'arrowstyle': '->', 'color': '#7951b0', 'ls': '--', 'lw': 1.2})
    ax.text(.50, .535, f"Solved from x = {parameters['start_x']:g}\nto radius = {parameters['outer_radius']:g} rₑ", ha='center', fontsize=8.5, color='#007e87', linespacing=1.3)
    ax.text(.5, .40, 'Remaining optical budget at the outer boundary', ha='center', fontsize=9.5, color='#475569')
    ax.text(.5, .26, f'{budget:.11f}', ha='center', fontsize=20, fontweight='bold', color='#007e87')
    ax.text(.5, .15, 'Center-normalized model units', ha='center', fontsize=8.5, color='#475569')
    ax.text(.5, -.205, 'Positive budget; the actual return path is unsolved.\nPositivity alone does not close the quantum state.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.6, color='#475569', linespacing=1.4)
    fig.text(.06, .17, r'$M=Gm/(c^2r_e)$. This is a static finite-region solution, not evolution of the incoming pulse or a dynamic stability test.', fontsize=9, color='#475569')
    fig.text(.06, .115, 'The full returning magnetic path, four-dimensional quantum source, physical preparation and global matching remain unresolved.', fontsize=9, color='#475569')
    fig.text(.06, .045, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'quantum_collar_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def magnetic_state_closure_progress():
    """Compare the saved angular source with its leading fixed-gauge particular."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_closure_results.json').read_text(encoding='utf-8'))
    fine = next(row for row in record['quadrature_refinement'] if row['order'] == 128)
    case = next(row for row in fine['roots'] if row['cutoff'] == 100)
    profile = case['profile']; mu = np.asarray(profile['nu'])-1
    weights = np.asarray(profile['weights']); omega2 = np.asarray(profile['channel_omega_squared'])
    source = (omega2-1)/3; mean = float(weights@source)
    rms = float(np.sqrt(weights@((source-mean)**2)))
    assert record['verification']['all_declared_checks_pass'] and len(mu) == 128
    assert mu.shape == weights.shape == source.shape and np.all(np.diff(mu) > 0)
    assert np.all(np.abs(mu) < 1) and np.all(weights > 0) and np.isclose(weights.sum(), 1, atol=1e-14, rtol=0)
    assert np.all(np.asarray(profile['periods']) > 0) and np.all(np.isfinite(source))
    assert np.isclose(mean, case['channel_source_mean_beta'], rtol=1e-12)
    assert np.isclose(rms, case['channel_source_rms_beta'], rtol=1e-12)
    assert np.isclose(rms/mean, case['source_rms_over_mean'], rtol=1e-12)
    modes = {row['degree']: row for row in case['source_angular_modes']}
    for degree, mode in modes.items():
        polynomial = np.polynomial.legendre.legval(mu, [0]*degree+[1])
        coefficient = (2*degree+1)*float(weights@(source*polynomial))
        assert np.isclose(coefficient, mode['source_coefficient'], rtol=1e-11, atol=1e-14)
        assert np.isclose(coefficient**2/((2*degree+1)*rms**2), mode['variance_fraction'], rtol=1e-10, atol=1e-14)
    shares = 100*np.array([modes[1]['variance_fraction'], modes[2]['variance_fraction'],
                           1-modes[1]['variance_fraction']-modes[2]['variance_fraction']])
    assert np.all(shares > 0) and np.isclose(shares.sum(), 100)
    old = (fine['old_assumed_omega_squared']-1)/3
    angular = json.loads((root/'polyakov_angular_results.json').read_text(encoding='utf-8'))
    geometry = next(row for row in angular['quadrature_refinement'] if row['order'] == 128)
    checks = angular['verification']; alpha = geometry['alpha']
    assert checks['all_declared_checks_pass'] and checks['ten_einstein_residuals'] == ['0']*10
    assert checks['four_maxwell_residuals'] == ['0']*4 and checks['delta_f_squared'] == '0'
    assert alpha == 1e-9 and np.isclose(geometry['angular_source_rms'], rms, rtol=1e-13)
    coefficients = np.asarray(geometry['metric_coefficients_over_alpha'])
    assert len(coefficients) == 7 and coefficients[0] == 0
    for degree in range(1, 7):
        assert np.isclose(coefficients[degree], -2*modes[degree]['source_coefficient']/(degree*(degree+1)), rtol=1e-13)
    angle = np.asarray(geometry['profile']['mu']); metric = np.asarray(geometry['profile']['metric_over_alpha'])
    assert angle.shape == metric.shape == (401,) and np.all(np.diff(angle) > 0)
    assert np.allclose(metric, np.polynomial.legendre.legval(angle, coefficients), atol=1e-14, rtol=0)
    bound = geometry['uniform_all_modes_metric_bound']
    assert np.isclose(bound, 2*alpha*rms, rtol=1e-13)
    assert np.isclose(alpha*np.max(np.abs(metric)), geometry['maximum_six_mode_metric_magnitude'], rtol=1e-13)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.08, right=.97, top=.755, bottom=.345, wspace=.29)
    fig.suptitle('Angular quantum source and leading metric response', fontsize=17.5, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Assumed throat geometry and prescribed magnetic paths • Cutoff radius 100 rₑ • 128-point angular quadrature',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]
    ax.plot(mu, source, color='#007e87', lw=2, label='Loop-dependent source')
    ax.axhline(mean, color='#7951b0', ls='--', lw=1.5, label=f'New angular mean: {mean:.6f}')
    ax.axhline(old, color='#64748b', ls=':', lw=1.3, label='Earlier single-delay assumption')
    ax.axhline(0, color='#94a3b8', lw=.6)
    ax.set(title='A. Source varies with direction', xlabel='Angular label μ = ν − 1',
           ylabel=r'$S=(\omega_{\rm channel}^2-1)/3$', xlim=(-1, 1), ylim=(-.38, .63))
    ax.legend(loc='lower right', fontsize=8, framealpha=.95, edgecolor='none')
    ax.text(.5, -.205, f'Mean {mean:.6f}; angular RMS deviation {rms:.6f}.\nOnly the angularly averaged state is closed here.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.plot(angle, metric, color='#007e87', lw=2, label='Six-mode reconstruction')
    ax.axhline(bound/alpha, color='#7951b0', ls='--', lw=1.2, label='Conditional all-mode envelope')
    ax.axhline(-bound/alpha, color='#7951b0', ls='--', lw=1.2)
    ax.axhline(0, color='#94a3b8', lw=.6)
    ax.set(title='B. Leading throat metric particular', xlabel='Angular coordinate μ',
           ylabel=r'$H/\alpha$', xlim=(-1, 1), ylim=(-.64, .64))
    ax.legend(loc='upper right', fontsize=8, framealpha=.95, edgecolor='none')
    ax.text(.5, -.205, r'Fixed gauge: $H=h_{tt}=h_{xx}$; this is not a throat shape.'+'\n'+
            f'All-mode |H| ≤ 5.556×10⁻¹⁰; source dipole + quadrupole: {shares[:2].sum():.3f}%.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.13); ax.set_axisbelow(True)
    fig.text(.06, .17, f'RMS/mean = {rms/mean:.3f} does not prove instability: the absolute stress remains proportional to α = 10⁻⁹.', fontsize=9.2, color='#475569')
    fig.text(.06, .115, 'Leading throat particular only; mouths and stability unresolved. Bounds use numerical moments, not certified continuum enclosures.', fontsize=9, color='#475569')
    fig.text(.06, .045, 'No constructed medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'magnetic_state_closure_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig)
    print(outputs[0])
    return outputs


def flux_channel_current_progress():
    """Show the saved local channel current and magnetic principal coefficient."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_channel_results.json').read_text(encoding='utf-8'))
    checks = record['verification']; params = record['parameters']; alpha = params['alpha']
    assert checks['all_declared_checks_pass'] and checks['negative_pressure_and_positive_principal']
    assert checks['action_variation_residual'] == checks['maxwell_reduction_residual'] == '0'
    assert checks['direct_radial_and_angular_ward_residuals'] == ['0', '0'] and alpha == 1e-9
    assert checks['ten_einstein_equation_reduction_residuals'] == ['0']*10
    assert checks['four_maxwell_equation_reduction_residuals'] == ['0']*4 and checks['background_quantum_trace_terms_retained']
    rows = [record['throat']]+record['collar_samples']
    radius = np.array([row['radius'] for row in rows]); pressure = np.array([row['quantum_pressure'] for row in rows])
    correction = np.array([row['magnetic_principal_correction'] for row in rows])
    scaled = np.array([row['principal_correction_over_alpha'] for row in rows])
    lapse = np.array([row['lapse'] for row in rows]); gradient = np.array([row['lapse_log_derivative'] for row in rows])
    assert np.all(np.diff(radius) > 0) and np.all(pressure < 0) and np.all(correction > 0)
    assert np.allclose(correction, -radius**4*pressure/(2*params['Q']**2), rtol=1e-12, atol=0)
    assert np.allclose(scaled, correction/alpha, rtol=1e-12, atol=0)
    expected = alpha*radius**2*(gradient**2+params['omega_squared']/lapse**2)/(6*params['Q']**2)
    assert np.allclose(correction, expected, rtol=1e-12, atol=0)
    assert np.isclose(radius[-1], params['outer_radius'], rtol=1e-13)
    assert np.isclose(correction.max(), record['maximum_sampled_principal_correction'], rtol=1e-13)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.055, right=.97, top=.755, bottom=.31, wspace=.25)
    fig.suptitle('Quantum channel stress and current are now linked', fontsize=17.5, fontweight='bold', y=.96)
    fig.text(.5, .87, 'Derived local channel equations • Assumed throat and finite static exterior • α = 10⁻⁹',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Field and source must respond together', pad=12)
    x = np.linspace(.05, .95, 200); y = .84+.045*np.sin(2*np.pi*x)
    for shift in (-.022, 0, .022):
        ax.plot(x, y+shift, color='#007e87', lw=1.8 if shift == 0 else .8, alpha=1 if shift == 0 else .3)
    ax.annotate('', xy=(.53, .84+.045*np.sin(2*np.pi*.53)), xytext=(.47, .84+.045*np.sin(2*np.pi*.47)),
                arrowprops={'arrowstyle': '->', 'color': '#007e87', 'lw': 1.8})
    ax.text(.5, .97, 'SCHEMATIC: bent magnetic channel', ha='center', va='top', fontsize=9, color='#475569')
    boxes = [( .22, .59, 'Metric change', '#edf7f8'), (.76, .59, 'Quantum stress\nDensity + tilt terms', '#edf7f8'),
             (.76, .20, 'Induced current', '#f3eff9'), (.22, .20, 'Physical mouth\nmatching\nunresolved', '#f3eff9')]
    for xx, yy, label, color in boxes:
        ax.text(xx, yy, label, ha='center', va='center', fontsize=10, linespacing=1.5,
                bbox={'boxstyle': 'round,pad=.6', 'facecolor': color, 'edgecolor': '#cbd5e1'})
    for start, end in [((.405, .59), (.55, .59)), ((.76, .44), (.76, .32)), ((.57, .20), (.41, .20))]:
        ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle': '->', 'color': '#64748b', 'lw': 1.3})
    ax.text(.5, -.215, 'Action variation and both conservation identities pass.\nThis diagram represents equations, not a material recipe.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.loglog(radius[1:], scaled[1:], color='#007e87', lw=2, label='Saved exterior samples')
    ax.scatter(radius[0], scaled[0], color='#7951b0', s=24, zorder=3, label='Throat sample')
    ax.set(title='B. Positive correction to the field equation', xlabel='Areal radius / rₑ (log scale)',
           ylabel=r'$(C-1)/\alpha$', xlim=(.9, 115), ylim=(9e-6, .5))
    ax.annotate('Throat: C − 1 = 2.699×10⁻¹⁰', xy=(radius[0], scaled[0]), xytext=(2.4, .26),
                fontsize=8.5, color='#475569', arrowprops={'arrowstyle': '-', 'lw': .7, 'color': '#94a3b8'})
    ax.annotate('At 100 rₑ: C − 1 = 1.667×10⁻¹⁴', xy=(radius[-1], scaled[-1]), xytext=(2.0, 2.4e-5),
                fontsize=8.5, color='#475569', arrowprops={'arrowstyle': '-', 'lw': .7, 'color': '#94a3b8'})
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.15, which='major'); ax.set_axisbelow(True)
    ax.text(.5, -.215, 'C > 1 avoids this static magnetic degeneracy only.\nIt does not establish dynamic or global stability.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.055, .12, 'Seed test uses α = 10⁻³. Full geometry, complete-loop state, physical matching and stability remain unresolved.',
             fontsize=9.2, color='#475569')
    fig.text(.055, .045, 'No stable medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'flux_channel_current_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs:
        fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs




def invariant_collar_progress():
    """Plot the saved invariant seed and boundary-basis tests on one collar."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_invariant_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass'] and record['algebra_checks']['all_declared_checks_pass']
    seed = record['target_seed_refinement'][-1]; basis = record['boundary_map_refinement'][-1]
    assert seed['alpha'] == record['parameters']['alpha'] == 1e-9
    traces = []
    for case in (seed, basis):
        rows = case['plot']; distance = np.array([row['proper_distance'] for row in rows])
        radius = np.array([row['radius'] for row in rows]); field = np.array([row['invariant_master'] for row in rows])
        assert np.all(np.isfinite(field)) and np.all(np.diff(distance) > 0) and np.all(np.diff(radius) > 0)
        assert distance[0] == 0 and radius[0] > 1
        traces.append((distance, radius, field))
    sl, sr, sp = traces[0]; bl, br, bp = traces[1]
    assert np.isclose(sr[-1], seed['gate_radius'], rtol=1e-13) and np.isclose(abs(sp[-1]), .01, rtol=1e-12)
    assert np.isclose(br[-1], record['parameters']['outer_radius'], rtol=1e-13)
    assert np.isclose(np.max(abs(bp)), basis['maximum_sampled_forced_master'], rtol=1e-13)
    assert abs(bp[0])+abs(bp[-1]) < 1e-40
    assert np.allclose(bp, basis['source_multiplier_for_saved_dipole']*np.array([row['unit_source_master'] for row in basis['plot']]), rtol=1e-13, atol=0)
    assert np.asarray(basis['endpoint_matrix']).shape == (2, 2)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.5))
    fig.subplots_adjust(left=.055, right=.975, top=.75, bottom=.31, wspace=.25)
    fig.suptitle('Collar response now uses coordinate-independent variables', fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Static linear model • Same α = 10⁻⁹ and prescribed source • Different endpoint assumptions',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. One side of an assumed throat', pad=12)
    x = np.linspace(.04, .94, 180); widening = np.clip((x-.30)/.64, 0, 1)
    halfwidth = .055+.245*widening**1.4
    ax.fill_between(x, .53-halfwidth, .53+halfwidth, color='#e5f3f4')
    ax.plot(x, .53-halfwidth, color='#007e87', lw=1.5); ax.plot(x, .53+halfwidth, color='#007e87', lw=1.5)
    for cut in (.30, .94): ax.plot([cut, cut], [.16, .85], ls='--', color='#7951b0', lw=1)
    ax.text(.30, .92, 'Inner cut\nx = 0.5', ha='center', va='center', fontsize=9.5)
    ax.text(.94, .95, 'Outer cut\nradius = 100 rₑ', ha='right', va='center', fontsize=9.5)
    ax.annotate('', xy=(.04, .53), xytext=(.20, .53), arrowprops={'arrowstyle': '->', 'color': '#007e87'})
    ax.text(.135, .39, 'Toward\nthroat', ha='center', va='top', fontsize=9, color='#475569')
    ax.text(.66, .53, 'Finite exterior collar', ha='center', va='center', fontsize=10, color='#007e87')
    ax.text(.53, .075, 'P(inner) and P(outer): physically undetermined', ha='center', fontsize=9.2, color='#475569')
    ax.text(.5, -.215, 'Schematic, not to scale; the other mouth is not shown.\nThe response map still needs physical end conditions.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.plot(sl, abs(sp), color='#007e87', lw=2, label='Arbitrary local seed; stopped at 1%')
    ax.plot(bl, abs(bp), color='#7951b0', lw=1.8, label='Zero P at both ends; diagnostic basis')
    ax.set_xscale('symlog', linthresh=.5, linscale=.7); ax.set_yscale('log')
    ax.set(title='B. End assumptions determine the response', xlabel='Outward proper distance / rₑ',
           ylabel=r'Invariant amplitude $|P|$', xlim=(0, 135), ylim=(1e-31, .05))
    ax.set_xticks([0, .5, 2, 10, 50, 126.3]); ax.set_xticklabels(['0', '0.5', '2', '10', '50', '126'])
    ax.scatter([sl[-1]], [abs(sp[-1])], color='#007e87', s=21, zorder=3)
    ax.annotate('Seed: |P| = 0.01\nradius = 1.206 rₑ', xy=(sl[-1], abs(sp[-1])), xytext=(.55, 2e-10),
                fontsize=8.3, color='#475569', arrowprops={'arrowstyle': '-', 'color': '#94a3b8', 'lw': .7})
    peak = int(np.argmax(abs(bp)))
    ax.annotate('Forced basis sampled maximum\n|P| = 3.334×10⁻²¹', xy=(bl[peak], abs(bp[peak])), xytext=(.065, 1e-15),
                fontsize=8.3, color='#475569', arrowprops={'arrowstyle': '-', 'color': '#94a3b8', 'lw': .7})
    ax.legend(loc='upper left', fontsize=7.8, framealpha=.95, edgecolor='none')
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.15); ax.set_axisbelow(True)
    ax.text(.5, -.215, 'Inner distances expanded (symlog); near-zero values below range.\nDifferent end assumptions: neither curve establishes stability.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.055, .12, 'The 1% invariant stop differs from the earlier coordinate stop. Global state feedback and physical mouth matching remain unresolved.',
             fontsize=9, color='#475569')
    fig.text(.055, .045, 'No constructed stable medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'invariant_collar_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def embedded_collar_progress():
    """Render one collar's embedded-boundary map and saved optical sensitivity."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_embedded_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass']
    geometry = record['geometry_checks']
    assert geometry['geometry_residuals'] == geometry['gauge_residuals'] == ['0']*6
    assert geometry['optical_gauge_residual'] == geometry['optical_P_factorization_residual'] == '0'
    assert record['input_order'] == ['P_inner', 'P_outer', 'Zeta_inner', 'Zeta_outer', 'beta_inner', 'beta_outer', 'delta_w']
    assert len(record['observable_order']) == 6 and record['parameters']['alpha'] == 1e-9
    fine = record['refinement'][-1]; alpha = record['parameters']['alpha']; rows = fine['plot']
    distance = np.array([row['proper_distance'] for row in rows]); radius = np.array([row['radius'] for row in rows])
    response = np.array([row['cumulative_optical_change_per_alpha_delta_w'] for row in rows])
    raw = np.array([row['cumulative_optical_change_per_delta_w'] for row in rows])
    optical = np.asarray(fine['optical_row']); matrix = np.asarray(fine['assembled_response_matrix'])
    assert matrix.shape == (13, 7) and np.array_equal(matrix[-1], optical)
    assert optical[4] == optical[5] == 0 and np.all(np.isfinite(response))
    assert len(rows) == 501 and distance[0] == response[0] == 0
    assert np.all(np.diff(distance) > 0) and np.all(np.diff(radius) > 0)
    assert np.allclose(raw, alpha*response, rtol=1e-13, atol=0)
    assert np.isclose(raw[-1], optical[-1], rtol=1e-12) and np.isclose(radius[-1], 100, rtol=1e-13)
    assert np.isclose(optical[2], -fine['cuts'][0]['inverse_lapse'], rtol=1e-13)
    assert np.isclose(optical[3], fine['cuts'][1]['inverse_lapse'], rtol=1e-13)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.7))
    fig.subplots_adjust(left=.055, right=.97, top=.76, bottom=.365, wspace=.25)
    fig.suptitle('Boundary response derived; actual matching unresolved', fontsize=17.5, fontweight='bold', y=.96)
    fig.text(.5, .88, 'Assumed static throat and channel model • One collar only • Common center-time normalization; rₑ = 1',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Geometry and fields at movable cuts', pad=12)
    x = np.linspace(.04, .94, 150); width = .045+.13*np.clip((x-.26)/.68, 0, 1)**1.3
    ax.fill_between(x, .67-width, .67+width, color='#e5f3f4')
    ax.plot(x, .67-width, color='#007e87', lw=1.5); ax.plot(x, .67+width, color='#007e87', lw=1.5)
    for cut, shifted, label in ((.26, .31, r'$Z_{\rm in}$'), (.86, .91, r'$Z_{\rm out}$')):
        ax.plot([cut, cut], [.46, .88], color='#7951b0', lw=1)
        ax.plot([shifted, shifted], [.46, .88], color='#7951b0', lw=1, ls='--')
        ax.annotate('', xy=(shifted, .93), xytext=(cut, .93), arrowprops={'arrowstyle': '<->', 'color': '#7951b0', 'lw': 1})
        ax.text((cut+shifted)/2, .98, label, ha='center', va='bottom', fontsize=10, color='#7951b0')
    ax.text(.56, .67, 'Finite collar', ha='center', color='#007e87')
    ax.text(.5, .405, 'Six boundary quantities (schematic, not to scale)', ha='center', fontsize=9, color='#475569')
    labels = ['Clock rate (lapse)', 'Angular metric', 'Time curvature', 'Angular curvature', 'Flux / area contrast', 'Magnetic boundary response']
    for i, label in enumerate(labels):
        xx = .035 if i < 3 else .53; yy = .29-.09*(i % 3)
        ax.text(xx, yy, '\u2022 '+label, va='center', fontsize=9.2, color='#334155')
    ax.text(.5, -.215, 'Inner cut x = 0.5; outer cut radius = 100 rₑ.\nZ denotes invariant cut displacement; end values remain inputs.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.plot(distance, response, color='#007e87', lw=2)
    ax.axhline(0, color='#94a3b8', lw=.7); ax.axhline(response[-1], color='#7951b0', ls='--', lw=1)
    ax.set_xscale('symlog', linthresh=.5, linscale=.7)
    ax.set(title='B. Source sensitivity of this collar', xlabel='Outward proper distance / rₑ (symlog)',
           ylabel=r'$\delta L/(\alpha\,\delta w)$', xlim=(0, 135), ylim=(-.12, .015))
    ax.set_xticks([0, .5, 2, 10, 50, 126.3]); ax.set_xticklabels(['0', '0.5', '2', '10', '50', '126'])
    ax.text(.50, .66, f'Final normalized sensitivity\n{response[-1]:.9f}', transform=ax.transAxes, ha='center', fontsize=10, color='#475569')
    ax.text(.50, .43, r'$q_w=-1.0834\times10^{-10}$', transform=ax.transAxes, ha='center', fontsize=11, color='#7951b0')
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.14); ax.set_axisbelow(True)
    ax.text(.5, -.215, 'P and cut displacements are held fixed as the source varies.\nThis is a static sensitivity, not a backward-time trajectory.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.5, .185, r'$\delta L=q_{\rm in}P_{\rm in}+q_{\rm out}P_{\rm out}+q_w\delta w+Z_{\rm out}/N_{\rm out}-Z_{\rm in}/N_{\rm in}$',
             ha='center', fontsize=13, color='#334155')
    fig.text(.055, .11, 'Physical end conditions, outside path, other mouth, complete-loop state and stability remain unresolved.', fontsize=9.2, color='#475569')
    fig.text(.055, .043, 'No constructed stable medium or 13 → 12 Sep 2027 passage established.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'embedded_collar_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def center_exterior_progress():
    """Show the reflection-even center and saved local RN omission diagnostic."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_center_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass']
    assert record['refinement'][-1]['parity_maximum_residual'] == 0
    assert record['rn_algebra']['residuals'] == ['0']*8
    diagnostic = record['rn_omission_diagnostic']; assert diagnostic['all_declared_checks_pass']
    assert diagnostic['parameters']['alpha'] == record['parameters']['alpha'] == 1e-9
    rows = diagnostic['refinement'][-1]['rows']; radius = np.array([row['radius'] for row in rows])
    signed = np.array([row['pulled_signed_ratios_Pprime_and_R_Xprime'] for row in rows])
    assert signed.shape == (3, 2) and np.allclose(radius, [10, 30, 100], rtol=1e-12, atol=0)
    for row in rows:
        omitted = np.asarray(row['omitted_rhs']); source = np.asarray(row['retained_state_rhs'])
        assert source[0] != 0 and source[1] == 0
        assert np.allclose(np.asarray(row['matrix_difference'])@row['classical_invariants'], omitted, rtol=1e-12, atol=0)
        legacy = [omitted[0]/source[0], row['radius']*omitted[1]/source[0]]
        assert np.allclose(legacy, row['signed_ratios_Pprime_and_R_Xprime'], rtol=1e-12, atol=0)
        assert np.isclose(max(abs(np.array(legacy))), row['weighted_infinity_ratio'], rtol=1e-12)
        pulled = np.asarray(row['pulled_classical_operator_mismatch'])
        assert np.allclose(pulled, omitted+row['radial_correction'], rtol=1e-12, atol=0)
        expected = [pulled[0]/source[0], row['radius']*pulled[1]/source[0]]
        assert np.allclose(expected, row['pulled_signed_ratios_Pprime_and_R_Xprime'], rtol=1e-12, atol=0)
        assert row['pulled_direct_equality_relative_error'] < 1e-40
    ratios = abs(signed); assert np.all(np.isfinite(ratios)) and np.all(ratios > 0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.055, right=.97, top=.755, bottom=.34, wspace=.25)
    fig.suptitle('Smooth center fixed; the exterior still needs quantum terms', fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Assumed static throat • Reflection-even sector • Target α = 10⁻⁹ • No material or machine established',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Center symmetry supplies an inner condition', pad=12)
    x = np.linspace(.04, .77, 150); width = .055+.21*((x-.18)/.59)**2
    ax.fill_between(x, .54-width, .54+width, color='#e5f3f4')
    ax.plot(x, .54-width, color='#007e87', lw=1.5); ax.plot(x, .54+width, color='#007e87', lw=1.5)
    ax.plot([.18, .18], [.35, .75], color='#7951b0', ls='--', lw=1)
    ax.text(.18, .84, 'Smooth center\nX = 0', ha='center', va='center', fontsize=10, color='#7951b0')
    ax.text(.47, .54, 'Computed collar', ha='center', color='#007e87', fontsize=10)
    ax.plot([.77, .77], [.23, .85], color='#7951b0', lw=1)
    ax.text(.77, .93, 'Outer cut\n100 rₑ', ha='center', va='center', fontsize=9.5)
    ax.plot([.77, .97], [.805, .855], color='#94a3b8', ls='--', lw=1.3)
    ax.plot([.77, .97], [.275, .225], color='#94a3b8', ls='--', lw=1.3)
    ax.text(.89, .54, 'Exterior\nunresolved', ha='center', va='center', fontsize=9.2, color='#64748b')
    ax.text(.48, .12, 'Theoretical geometry; schematic, not to scale', ha='center', fontsize=9, color='#475569')
    ax.text(.5, -.21, 'X = 0 follows from the chosen reflection symmetry.\nThe physical exterior response remains an unknown input.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]; locations = np.arange(3); width = .30
    for j, (color, label) in enumerate((('#007e87', r'$|O^{\rm pull}_P/S_P|$'), ('#7951b0', r'$|R\,O^{\rm pull}_X/S_P|$'))):
        centers = locations+(j-.5)*width
        ax.bar(centers, ratios[:, j], width=width, color=color, label=label)
        for xx, value in zip(centers, ratios[:, j]):
            ax.text(xx, value+.08, f'{value:.3g}', ha='center', va='bottom', fontsize=9, color='#334155')
    ax.axhline(1, color='#94a3b8', ls='--', lw=1)
    ax.text(.01, 1.05, 'Same magnitude as state source', transform=ax.get_yaxis_transform(), fontsize=8, color='#64748b', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': .3})
    ax.set(title='B. Pulled operator mismatch vs state forcing', xlabel='Areal radius / rₑ',
           ylabel='Pulled mismatch / state source', xticks=locations, xticklabels=['10', '30', '100'], ylim=(0, 2.8))
    ax.legend(loc='upper left', fontsize=10, framealpha=.95, edgecolor='none', ncol=2)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    ax.text(.5, -.21, 'O includes the radial-coordinate pullback.\nBoth ratios use nonzero Sₚ; the X mismatch is weighted by R.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.055, .16, 'Growing-only classical input: the operator mismatch is not a full forced residual, global solution error or instability test.',
             fontsize=9.2, color='#475569')
    fig.text(.055, .108, 'Physical matching, exterior quantum response, complete-loop state and dynamic stability remain unresolved.', fontsize=9.2, color='#475569')
    fig.text(.055, .043, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'center_exterior_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def overlap_transport_progress():
    """Plot retained finite-collar transport and its asymptotic cancellation."""
    from fractions import Fraction
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_overlap_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass'] and record['asymptotic_checks']['all_declared_checks_pass']
    checks = record['asymptotic_checks']; large = checks['large_radius']
    assert checks['background_equation_residuals'] == ['0']*3 and checks['growing_projection_residual'] == '0'
    assert checks['inverse_residual'] == [['0', '0'], ['0', '0']] and large['leading_cancellation'] == '0'
    leading = [Fraction(x) for x in large['growing_adjoint_terms_times_r_squared']]
    assert leading == [Fraction(-4, 27), Fraction(4, 27)] and sum(leading) == 0
    assert Fraction(large['growing_projection_times_r_cubed']) == Fraction(-31, 27)
    fine = record['refinement'][-1]; rows = fine['plot']; amplitude = record['parameters']['source_amplitude']
    radius = np.array([row['radius'] for row in rows]); correction = np.array([row['growing_shift_from_incident_over_A'] for row in rows])
    raw = np.array([row['growing_shift_from_incident'] for row in rows])
    cg = record['incident_growing_amplitude']; separation = record['parameters']['separation_over_re']
    assert record['parameters']['alpha'] == 1e-9 and amplitude > 0 and np.isclose(cg, -2/separation**2, rtol=1e-13)
    assert np.allclose(radius[[0, -1]], [10, 100], rtol=1e-13) and np.all(np.diff(radius) > 0)
    assert np.all(np.isfinite(correction)) and correction[0] == 0 and np.all(np.diff(correction) > 0)
    assert np.allclose(raw, amplitude*correction, rtol=1e-13, atol=0)
    assert np.isclose(raw[-1], fine['transfer_correction'][0][0]*cg, rtol=1e-13)
    assert np.isclose(correction[-1], fine['source_order_incident_growing_shift'], rtol=1e-13)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.055, right=.97, top=.755, bottom=.34, wspace=.25)
    fig.suptitle('Local terms partly cancel; finite boundary correction remains', fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Static linear collar model • Computed region: radii 10–100 rₑ • α = 10⁻⁹',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Leading terms cancel after projection', pad=12)
    x = np.linspace(.10, .74, 90); halfwidth = .06+.075*(x-.10)/.64
    ax.fill_between(x, .70-halfwidth, .70+halfwidth, color='#e5f3f4')
    ax.plot(x, .70-halfwidth, color='#007e87', lw=1.5); ax.plot(x, .70+halfwidth, color='#007e87', lw=1.5)
    for cut, label in ((.10, '10 rₑ'), (.74, '100 rₑ')):
        ax.plot([cut, cut], [.51, .91], color='#7951b0', lw=1)
        ax.text(cut, .95, label, ha='center', va='bottom', fontsize=10)
    ax.text(.41, .70, 'Computed overlap', ha='center', va='center', color='#007e87')
    ax.plot([.74, .97], [.835, .87], color='#94a3b8', ls='--', lw=1.3)
    ax.plot([.74, .97], [.565, .53], color='#94a3b8', ls='--', lw=1.3)
    ax.text(.875, .70, 'Exterior\nunresolved', ha='center', va='center', fontsize=9, color='#64748b')
    ax.text(.5, .37, r'$-\frac{4}{27r^2}+\frac{4}{27r^2}=0$', ha='center', va='center', fontsize=19, color='#7951b0')
    ax.text(.5, .16, r'Next term: $-31/(27r^3)$', ha='center', fontsize=12, color='#334155')
    ax.text(.5, -.21, 'Schematic; first-order large-radius coefficients per k c_g, k = α/3.\nThese expansion terms do not determine the far exterior.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.plot(radius, correction, color='#007e87', lw=2)
    ax.axhline(correction[-1], color='#7951b0', ls='--', lw=1)
    ax.scatter(radius[-1], correction[-1], color='#7951b0', s=24, zorder=3)
    ax.set(title='B. Cumulative correction for one test input', xlabel='Areal radius / rₑ',
           ylabel=r'$\Delta c_g/A_{\rm dip}$', xlim=(10, 103), ylim=(0, .0205), xticks=[10, 30, 50, 70, 100])
    ax.text(.57, .48, f'At 100 rₑ\nDelta c_g / A_dip = {correction[-1]:.5f}', transform=ax.transAxes,
            ha='center', va='center', fontsize=11, color='#475569', linespacing=1.5)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.14); ax.set_axisbelow(True)
    ax.text(.5, -.21, 'Curve: incident-induced correction only, normalized by A_dip.\nSource forcing and an unknown decaying input are separate terms.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.055, .16, r'Declared growing-only input: $c_g=-2/d^2$. RN functions serve as a basis for the retained collar equations.',
             fontsize=9.2, color='#475569')
    fig.text(.055, .108, 'Physical endpoint data, the exterior return region beyond 100 rₑ, complete-loop state and stability remain unresolved.', fontsize=9.2, color='#475569')
    fig.text(.055, .043, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'overlap_transport_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def return_current_progress():
    """Render prescribed return paths and discrete saved current quadrature."""
    from mmp_pulse import _magnetic_return_path
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_return_current_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass'] and record['verification']['default_source_compatibility']
    fine = record['refinement'][-1]; nu = np.asarray(fine['nu']); weights = np.asarray(fine['weights'])
    integrand = np.asarray(fine['flux_integrand']); contributions = np.asarray(fine['weighted_contribution'])
    cumulative = np.asarray(fine['cumulative_integral']); integral = record['integral']
    assert fine['flux_order'] == len(nu) == 128 and fine['path_order'] == 256
    assert np.all(np.diff(nu) > 0) and np.all((nu > 0) & (nu < 2)) and np.all(weights > 0)
    assert np.isclose(weights.sum(), 1, rtol=0, atol=1e-14) and np.all(integrand < 0)
    assert np.allclose(contributions, weights*integrand, rtol=1e-13, atol=0)
    assert np.allclose(cumulative, np.cumsum(contributions), rtol=1e-13, atol=0)
    assert np.isclose(cumulative[-1], integral, rtol=1e-13, atol=0) and np.isclose(fine['integral'], integral, rtol=1e-13, atol=0)
    alpha = record['parameters']['alpha']; separation = record['parameters']['separation_over_re']
    assert alpha == 1e-9 and np.isclose(record['magnetic_growing_contribution'], 16*np.pi**2*alpha*integral/separation**2, rtol=1e-13, atol=0)
    assert np.isclose(record['magnetic_contribution_over_A'], record['magnetic_growing_contribution']/record['source_amplitude'], rtol=1e-13, atol=0)
    paths = []
    for label in (.5, 1., 1.5):
        values = np.array([_magnetic_return_path(label, float(s))[:2] for s in np.linspace(-.999999, .999999, 501)])
        rho, z = values.T
        flux = (z+.5)/np.hypot(rho, z+.5)-(z-.5)/np.hypot(rho, z-.5)
        assert np.all(rho > 0) and np.all(np.isfinite(values)) and np.allclose(flux, label, atol=1e-8, rtol=0)
        paths.append((label, np.r_[0, rho, 0], np.r_[-.5, z, .5]))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.6))
    fig.subplots_adjust(left=.075, right=.97, top=.755, bottom=.34, wspace=.28)
    fig.suptitle('Returning channels give a negative magnetic contribution', fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'Prescribed two-mouth magnetic geometry • Ideal flat-region model • Full saved channel weights • α = 10⁻⁹',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]
    for (label, rho, z), color in zip(paths, ('#007e87', '#7951b0', '#b07d24')):
        ax.plot(rho, z, color=color, lw=1.5, label=f'ν = {label:g}')
        ax.plot(-rho, z, color=color, lw=1.5)
    ax.scatter([0, 0], [-.5, .5], s=28, color='#334155', zorder=4)
    ax.annotate('Ideal mouth A', xy=(0, .5), xytext=(.2, .32), fontsize=8, color='#334155',
                arrowprops={'arrowstyle': '-', 'lw': .7, 'color': '#94a3b8'})
    ax.annotate('Ideal mouth B', xy=(0, -.5), xytext=(.2, -.33), fontsize=8, color='#334155',
                arrowprops={'arrowstyle': '-', 'lw': .7, 'color': '#94a3b8'})
    ax.set(title='A. Reference magnetic return paths', xlabel='Transverse position / separation d', ylabel='Axial position / d')
    ax.set_aspect('equal', adjustable='datalim'); ax.legend(loc='upper right', fontsize=8, framealpha=.95, edgecolor='none')
    ax.text(.5, -.21, 'Meridional slice of the theoretical field; mouth centers are idealized.\nSmaller ν gives larger loops; ν → 2 approaches the direct path.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]
    ax.step(np.r_[0, nu], np.r_[0, -cumulative], where='post', color='#007e87', lw=1.6)
    ax.scatter(nu[-1], -cumulative[-1], color='#7951b0', s=25, zorder=3)
    ax.set(title='B. Accumulated magnitude of the contribution', xlabel='Field-line label (ν)',
           ylabel=r'$-\sum_i w_iG(\nu_i)$', xlim=(0, 2.03), ylim=(0, .0035))
    ax.text(.67, .25, f'Negative total contribution\nI = {integral:.8f}', transform=ax.transAxes,
            ha='center', fontsize=11, color='#475569', linespacing=1.5)
    ax.text(.5, -.21, 'Discrete accumulated Gauss weights; not a continuous partial integral.\nThe positive curve displays the magnitude of a negative contribution.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.12); ax.set_axisbelow(True)
    fig.text(.055, .16, r'$\delta c_{g,\rm magnetic}=16\pi^2\alpha I/d^2$ : a source-generated field contribution, not a complete coupled boundary solution.',
             fontsize=9.2, color='#475569')
    fig.text(.055, .108, 'Finite-channel endpoint physics, curved-region/interface fields, the full metric and dynamical stability remain unresolved.', fontsize=9.2, color='#475569')
    fig.text(.055, .043, 'No material recipe, constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'return_current_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def return_feedback_progress():
    """Compare the saved fixed-collar forcing and conditional mean response."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_feedback_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass'] and record['parameters']['alpha'] == 1e-9
    fine = record['refinement'][-1]; mean = record['mean_closure']
    weights = np.asarray(fine['weights']); magnetic = np.asarray(fine['channel_omega_squared_change'])
    forcing = mean['fixed_collar_forcing']; adjusted = mean['mean_omega_squared_shift']
    derivative = mean['derivative_refinement'][-1]; slope = derivative['mean_map_derivative']
    assert len(weights) == fine['flux_order'] == len(magnetic) == 128
    assert np.all(weights > 0) and np.isclose(weights.sum(), 1, atol=1e-14, rtol=0)
    assert np.all(np.asarray(fine['length_response_per_coupling']) < 0)
    assert np.all(np.asarray(fine['fractional_period_change']) < 0) and np.all(magnetic > 0)
    assert np.isclose(weights@magnetic, forcing, rtol=1e-12, atol=0)
    assert np.isclose(fine['fixed_collar_mean_omega_squared_change'], forcing, rtol=1e-13, atol=0)
    channel_derivative = np.asarray(derivative['channel_derivative'])
    assert np.isclose(weights@channel_derivative, slope, rtol=1e-12)
    assert np.isclose(mean['implicit_denominator'], 1-slope, rtol=1e-13)
    assert np.isclose(adjusted, forcing/(1-slope), rtol=1e-13, atol=0) and 0 < adjusted < forcing
    assert np.isclose(adjusted/forcing, mean['retained_fraction_of_fixed_collar_forcing'], rtol=1e-13)
    total = np.asarray(mean['channel_omega_squared_total_change'])
    assert np.allclose(total, magnetic+channel_derivative*adjusted, rtol=1e-12, atol=1e-24)
    assert np.isclose(weights@total, adjusted, rtol=1e-12, atol=0)
    assert np.isclose(mean['clock_fractional_shift'], derivative['clock_log_derivative']*adjusted, rtol=1e-13, atol=0)
    reduction = 100*(1-adjusted/forcing)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.055, right=.97, top=.755, bottom=.335, wspace=.25)
    fig.suptitle('Mean throat and clock response partly offsets magnetic feedback', fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .875, 'First static response only • Prescribed flat-region magnetic model • α = 10⁻⁹',
             ha='center', fontsize=10, color='#475569')
    ax = axes[0]; ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Competing responses in the model', pad=12)
    boxes = [(.29, .83, 'Quantum current\nchanges the field', '#edf7f8'),
             (.29, .51, 'Returning paths\nbecome shorter', '#edf7f8'),
             (.29, .18, 'Stronger negative\nCasimir contribution', '#edf7f8'),
             (.80, .51, 'Allow mean throat\n+ clock response', '#f3eff9')]
    for xx, yy, label, color in boxes:
        ax.text(xx, yy, label, ha='center', va='center', fontsize=10, linespacing=1.4,
                bbox={'boxstyle': 'round,pad=.55', 'facecolor': color, 'edgecolor': '#cbd5e1'})
    for start, end in [((.29, .70), (.29, .64)), ((.29, .38), (.29, .31))]:
        ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle': '->', 'color': '#007e87', 'lw': 1.4})
    ax.annotate('', xy=(.53, .19), xytext=(.73, .37), arrowprops={'arrowstyle': '->', 'color': '#7951b0', 'lw': 1.4})
    ax.text(.78, .25, 'Partly offsets\nthe mean change', ha='center', va='top', fontsize=9, color='#7951b0', linespacing=1.4)
    ax.text(.5, -.21, 'Fixed-flux paths shorten for the prescribed magnetic response.\nThe additional mean throat/clock adjustment opposes part of it.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    ax = axes[1]; values = np.array([forcing, adjusted])/1e-10
    ax.bar([0, 1], values, width=.52, color=['#007e87', '#7951b0'])
    for xx, value in enumerate(values):
        ax.text(xx, value+.065, f'{value:.3f}', ha='center', va='bottom', fontsize=11, color='#334155')
    ax.set(title='B. A smaller conditional mean-state correction', ylabel=r'Change in mean $\omega^2$ ($\times10^{-10}$)',
           xticks=[0, 1], xticklabels=['Collar fixed', 'Mean collar allowed\nto respond'], xlim=(-.55, 1.65), ylim=(0, 2.6))
    ax.text(.66, .81, f'About {reduction:.0f}% smaller\nin this restricted calculation', transform=ax.transAxes,
            ha='center', va='center', fontsize=10.5, color='#475569', linespacing=1.4)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    ax.text(.5, -.21, f'Mean correction = fixed-collar forcing / {mean["implicit_denominator"]:.3f}.\nThis static comparison does not establish dynamical stability.',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color='#475569', linespacing=1.4)
    fig.text(.055, .15, 'Longest channels exceed the perturbative/localized-channel regime; nonspherical metric and interface feedback remain unknown.',
             fontsize=9.2, color='#475569')
    fig.text(.055, .043, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'return_feedback_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def gravity_feedback_progress():
    """Show selected static path corrections and a separate matching-clock test."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'polyakov_gravity_feedback_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass']
    assert record['finite_cut_checks']['verification']['all_declared_checks_pass']
    mean = record['conditional_mean']; match = record['matching_sensitivity']; par = record['parameters']
    assert par['alpha'] == 1e-9 and mean['implicit_denominator'] > 1
    keys = ('magnetic', 'gravitational', 'combined')
    shifts = np.array([mean[f'{key}_mean_shift'] for key in keys])
    forcings = np.array([mean[f'{key}_forcing'] for key in keys])
    assert np.all(np.isfinite(shifts)) and np.all(shifts > 0)
    assert np.allclose(shifts, forcings/mean['implicit_denominator'], rtol=1e-13, atol=0)
    assert np.isclose(shifts[0]+shifts[1], shifts[2], rtol=1e-13, atol=0)
    assert np.isclose(forcings[0]+forcings[1], forcings[2], rtol=1e-13, atol=0)
    epsilon = par['epsilon_m_over_d']; clock_test = match['mean_shift_for_test_clock_offset_epsilon']
    assert np.isclose(epsilon, par['gravitational_mass_over_re']/par['separation_over_re'], rtol=1e-13, atol=0)
    assert np.isclose(clock_test, match['mean_shift_per_independent_log_clock_offset']*epsilon, rtol=1e-13, atol=0)
    assert record['causality_scope']['throat_optical_time_over_separation'] > 0
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 6.4))
    fig.suptitle('Two path corrections are comparable; clock matching still matters',
                 fontsize=16, fontweight='bold', y=.955)
    fig.text(.5, .865, 'Conditional static model  |  Same radial mean response applied to all three bars  |  α = 10⁻⁹',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.085, .345, .49, .43]); values = shifts/1e-10
    ax.bar(np.arange(3), values, width=.55, color=['#007e87', '#7951b0', '#334155'])
    for xx, value in enumerate(values):
        ax.text(xx, value+.055, f'{value:.3f}', ha='center', va='bottom', fontsize=12, color='#334155')
    ax.set(xticks=np.arange(3), xticklabels=['Quantum-current\npath', 'Gravity-induced\npath', 'Selected sum'],
           ylabel=r'Change in mean $\omega^2$ ($\times10^{-10}$)', ylim=(0, 1.9), xlim=(-.55, 2.55))
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    box = fig.add_axes([.64, .345, .305, .43]); box.set(xlim=(0, 1), ylim=(0, 1)); box.axis('off')
    box.add_patch(plt.Rectangle((.015, .015), .97, .97, facecolor='#fff8df', edgecolor='#e7cf83', lw=1.1))
    box.text(.5, .85, 'Separate clock matching test', ha='center', va='center', fontsize=11.5, fontweight='bold', color='#7c5d13')
    box.text(.5, .665, r'Impose test offset $+\eta=m/d$', ha='center', fontsize=11, color='#7c5d13')
    box.text(.5, .49, f'Additional mean change\n{clock_test/1e-10:.3f} × 10⁻¹⁰', ha='center', va='center',
             fontsize=13, linespacing=1.45, color='#7c5d13')
    box.text(.5, .205, 'Physical offset remains undetermined.\nThis is not an uncertainty bound\nand is not included in the bars.',
             ha='center', va='center', fontsize=9.5, linespacing=1.45, color='#7c5d13')
    fig.text(.085, .218, 'Shapiro delay is already in the baseline; the new gravity term adds path bending only.', fontsize=10, color='#475569')
    fig.text(.085, .161, 'Physical mouth matching and nonspherical response remain unresolved; these bars are not a stability test.',
             fontsize=9.5, color='#475569')
    fig.text(.085, .096, 'A positive static lapse and shared time exclude backward travel within this assumed model.',
             fontsize=10, fontweight='bold', color='#007e87')
    fig.text(.085, .041, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10, color='#475569')
    outputs = tuple(root/f'gravity_feedback_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def clock_response_progress():
    """Separate stationary-offset tests from the zero-offset linear spectrum."""
    root = Path(__file__).resolve().parent
    static = json.loads((root/'stationary_offset_gate_results.json').read_text(encoding='utf-8'))
    spectral = json.loads((root/'distributed_clock_spectrum_results.json').read_text(encoding='utf-8'))
    assert static['verification']['all_declared_checks_pass']
    assert spectral['verification']['all_declared_checks_pass']
    assert spectral['verification']['all_identity_errors_within_analytic_truncation_bounds']
    assert spectral['proof']['antisymmetric_sign'] and spectral['proof']['symmetric_sign']
    assert all(c['b0'] > 0 and np.all(np.asarray(c['weights']) > 0) for c in spectral['cases'])
    for case in spectral['cases']:
        for sector in ('antisymmetric', 'symmetric'):
            assert case[f'{sector}_error'] <= case[f'{sector}_tail_bound']
    fine = sorted((c for c in static['cases'] if c['order'] == 128), key=lambda c: c['offset_over_d'])
    coarse = sorted((c for c in static['cases'] if c['order'] == 64), key=lambda c: c['offset_over_d'])
    offset = np.array([c['offset_over_d'] for c in fine])
    crossing = np.array([c['throat_time_over_d'] for c in fine])
    delay = np.array([c['earlier_direction_exit_delay_over_d'] for c in fine])
    assert len(fine) == len(coarse) == 4 and np.array_equal(offset, [0, 1, 3, 10])
    assert np.all(np.isfinite(crossing)) and np.all(delay > 0)
    assert np.allclose(crossing-offset, delay, rtol=1e-13, atol=0)
    assert np.max(np.abs(crossing-np.array([c['throat_time_over_d'] for c in coarse]))) < 1e-10
    assert all(c['minimum_short_period'] > 0 and c['fastest_complete_return_over_d'] > 0 for c in fine)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 6.8))
    fig.suptitle('Tested static clock offsets still give later exits', fontsize=18, fontweight='bold', y=.955)
    fig.text(.5, .875, 'Separate assumed states, not a preparation sequence  |  Times in d/c units; no calendar conversion',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.05, .375, .345, .405]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Assumed geometry — schematic only', fontsize=11.5, pad=9)
    t = np.linspace(0, 1, 121); x = .13+.74*t; width = .09+.15*np.abs(2*t-1)**1.7
    ax.fill_between(x, .55-width, .55+width, color='#edf7f8')
    for sign in (-1, 1): ax.plot(x, .55+sign*width, color='#007e87', lw=1.6)
    angle = np.linspace(0, 2*np.pi, 100)
    for xx, label in ((.13, 'Mouth A'), (.87, 'Mouth B')):
        ax.plot(xx+.035*np.cos(angle), .55+.24*np.sin(angle), color='#007e87', lw=1.6)
        ax.text(xx, .91, label, ha='center', color='#334155')
    ax.annotate('', xy=(.82, .55), xytext=(.18, .55), arrowprops={'arrowstyle': '->', 'lw': 1.5, 'color': '#7951b0'})
    ax.text(.5, .66, 'Light crossing', ha='center', fontsize=9.5, color='#7951b0')
    ax.text(.5, .34, 'Imposed relative clock offset S', ha='center', fontsize=10, color='#475569')
    ax.text(.5, .095, 'Arrival delay = crossing time − offset\nPositive means arrival is still later.',
            ha='center', va='center', fontsize=10, color='#334155', linespacing=1.5)
    ax = fig.add_axes([.475, .375, .485, .405]); xx = np.arange(4); w = .235
    ax.bar(xx-w, crossing, width=w, color='#007e87', label='Crossing time')
    ax.bar(xx, offset, width=w, color='#c2cbd5', label='Imposed offset')
    ax.bar(xx+w, delay, width=w, color='#7951b0', label='Arrival delay')
    for xval, value in zip(xx+w, delay): ax.text(xval, value+.27, f'{value:.2f}', ha='center', fontsize=9, color='#7951b0')
    ax.set(title='B. Static offset: no past exit in the four tests', ylabel='Time / (d/c)',
           xticks=xx, xticklabels=['0', '1', '3', '10'], xlabel='Imposed time offset / (d/c)', ylim=(0, 20.5))
    ax.legend(loc='upper left', ncol=1, frameon=False, fontsize=9)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    fig.text(.5, .264, 'The static mean-throat crossing time increases as the imposed offset increases.',
             ha='center', fontsize=10, color='#475569')
    fig.text(.5, .196, 'Small clock perturbations: no exponentially growing mode in either linear sector.',
             ha='center', fontsize=11, fontweight='bold', color='#007e87')
    fig.text(.5, .151, 'Zero-offset NS vacuum only • fixed delays • b₀ > 0. Neutral modes remain; this does not establish damping or nonlinear stability.',
             ha='center', fontsize=9.2, color='#475569')
    fig.text(.05, .080, 'Moving exterior and physical preparation remain unresolved; the pictured geometry is assumed, not constructed.',
             fontsize=10, color='#475569')
    fig.text(.05, .030, 'No stable material or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, color='#007e87', fontweight='bold')
    outputs = tuple(root/f'clock_response_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def distributed_sampling_progress():
    """Compare verified saved quadratures without rerunning the delay dynamics."""
    root = Path(__file__).resolve().parent
    names = ('distributed_loop_response_results.json', 'distributed_loop_panel_results.json', 'distributed_loop_panel_fine_results.json')
    records = [json.loads((root/name).read_text(encoding='utf-8')) for name in names]
    final = records[-1]; verification = final['verification']
    assert verification['all_declared_checks_pass'] and verification['flux_convergence_pass']
    assert verification['timestep_convergence_pass']
    assert verification['one_channel_D01_old_source_maximum_state_difference'] is not None
    assert not records[0]['verification']['all_declared_checks_pass']
    assert not records[1]['verification']['all_declared_checks_pass']
    comparisons = ('flux_64_128', 'flux_8_16', 'flux_16_32')
    errors, scales, tolerances = [], [], []
    for record, key in zip(records, comparisons):
        v = record['verification']; observable = v['observable_maximum_absolute_differences'][key]
        error = max(observable['mean_j_plus_response_over_alpha'], observable['mean_j_minus_response_over_alpha'])
        assert np.isclose(error, v[f'{key}_maximum_absolute_difference'], rtol=1e-12, atol=0)
        scale = record['refinement'][1]['maximum_absolute_mean_stress_response']
        assert scale > 0 and error > 0
        assert np.isclose(error/scale, v['flux_error_over_peak_mean_stress'], rtol=1e-9, atol=0)
        errors.append(error); scales.append(scale); tolerances.append(v['acceptance_thresholds'][f'{key}_absolute'])
    assert np.array_equal(tolerances, [2e-7]*3)
    assert errors[0] > tolerances[0] and errors[1] > tolerances[1] and errors[2] < tolerances[2]
    assert all(r['parameters']['pulse_epsilon'] == 1e-4 and r['parameters']['shortest_loops'] == 4 for r in records)
    percentages = 100*np.array(errors)/np.array(scales); tolerance_percent = 100*tolerances[-1]/scales[-1]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 6.6))
    fig.suptitle('Same equations; finer arrival-time sampling meets the check', fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .875, 'Finite linear pulse-response test  |  The acceptance threshold was not relaxed',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.05, .355, .36, .415]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Unequal return paths — assumed geometry', fontsize=11, pad=10)
    t = np.linspace(0, 1, 161); x = .13+.74*t
    colors = ['#007e87', '#7951b0', '#aa7900']
    for height, color, pulse_t in zip((.20, .36, .53), colors, (.76, .61, .49)):
        ax.plot(x, .35+height*np.sin(np.pi*t), color=color, lw=1.7)
        ax.plot(.13+.74*pulse_t, .35+height*np.sin(np.pi*pulse_t), 'o', color=color, ms=6)
    ax.plot(x, .35-.10*np.sin(np.pi*t), '--', color='#64748b', lw=1.2)
    ax.plot([.13, .87], [.35, .35], 'o', color='#334155', ms=8)
    ax.text(.5, .30, 'Assumed throat', ha='center', fontsize=9.5, color='#475569')
    ax.text(.06, .35, 'A', ha='center', va='center', fontweight='bold', color='#334155')
    ax.text(.94, .35, 'B', ha='center', va='center', fontweight='bold', color='#334155')
    ax.annotate('', xy=(.93, .065), xytext=(.10, .065), arrowprops={'arrowstyle': '->', 'color': '#64748b'})
    for center, color in zip((.35, .57, .80), colors):
        pulse_x = np.linspace(center-.06, center+.06, 41)
        ax.plot(pulse_x, .065+.085*np.exp(-((pulse_x-center)/.022)**2), color=color, lw=1.5)
    ax.text(.5, -.08, 'Different path lengths spread the pulse arrival times.\nSampling is refined within each arrival-time interval.',
            transform=ax.transAxes, ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.5)
    ax = fig.add_axes([.50, .355, .455, .415])
    ax.bar(np.arange(3), percentages, width=.52, color=['#94a3b8', '#aa7900', '#007e87'])
    for xx, value in enumerate(percentages):
        ax.text(xx, value*1.25, f'{value:.3g}%', ha='center', fontsize=11, color='#334155')
    ax.set(yscale='log', ylim=(.01, 10), xlim=(-.55, 2.55), yticks=[.01, .1, 1, 10], yticklabels=['0.01', '0.1', '1', '10'],
           xticks=np.arange(3), xticklabels=['Global\n64 → 128', 'Panels\n8 → 16', 'Panels\n16 → 32'],
           ylabel='Sampling difference / peak mean source (%)', title='B. Successive sampling comparisons')
    ax.axhline(tolerance_percent, color='#334155', ls='--', lw=1.2)
    ax.text(2.48, tolerance_percent*1.22, f'Unchanged tolerance: {tolerance_percent:.3f}%', ha='right', fontsize=9, color='#334155')
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    fig.text(.728, .239, 'Panel numbers are points per arrival-time interval.\nAbsolute source tolerance: 2 × 10⁻⁷ (model units).',
             ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.5)
    fig.text(.05, .130, 'Time-step refinement and the independent one-channel control also pass the saved checks.', fontsize=10, color='#007e87')
    fig.text(.05, .080, 'Finite sampled convergence is not a stability proof; exterior motion, physical matching and preparation remain unresolved.',
             fontsize=9.5, color='#475569')
    fig.text(.05, .030, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'distributed_sampling_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def moving_exterior_progress():
    """Plot the saved prescribed-motion test, with first-order timing scope."""
    root = Path(__file__).resolve().parent
    maps = json.loads((root/'moving_exterior_map_results.json').read_text(encoding='utf-8'))
    wrapper = json.loads((root/'moving_exterior_response_results.json').read_text(encoding='utf-8'))
    assert wrapper['status'] == 'VERIFIED' and maps['verification']['all_declared_checks_pass']
    record = wrapper['result']; verification = record['verification']; par = record['parameters']
    assert verification['all_declared_checks_pass']
    assert verification['timestep_halving_maximum_absolute_difference'] < verification['timestep_threshold']
    assert wrapper['extraction_regression']['state_maximum_difference'] == 0
    assert par['external_pulse'] == 0 and par['displacement_fraction'] == .001
    assert np.isclose(par['period'], np.pi+par['exterior_delay'], rtol=1e-14, atol=0)
    fine = record['refinement'][-1]; data = {k: np.asarray(v) for k, v in fine['plot'].items()}
    time = data['exterior_time']; launch = data['return_launch_time']; period = par['period']
    assert np.all(np.diff(time) > 0) and np.all(np.diff(launch) > 0)
    assert np.isclose(par['peak_displacement'], par['displacement_fraction']*par['exterior_delay'], rtol=1e-14, atol=0)
    z = time/par['duration']; expected = np.zeros_like(z); active = (z > 0)&(z < 1)
    expected[active] = par['peak_displacement']*256*z[active]**4*(1-z[active])**4
    assert np.max(np.abs(expected-data['displacement'])) < 64*np.finfo(float).eps*par['peak_displacement']
    for side in ('A', 'B'):
        assert np.all(data[f'clock_rate_{side}'] > 0) and np.all(data[f'full_return_{side}'] > 0)
        total = data[f'return_change_{side}']; summed = data[f'clock_return_change_{side}']+data[f'direct_return_displacement_{side}']
        assert np.allclose(total, summed, rtol=1e-12, atol=0)
        assert np.max(np.abs(data[f'full_return_{side}']-period-total)) < 8*np.finfo(float).eps*period
    minimum_return = min(np.min(data[f'full_return_{side}']) for side in ('A', 'B'))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.5})
    fig = plt.figure(figsize=(14, 6.9)); colors = ['#007e87', '#7951b0']
    fig.suptitle('Prescribed motion changes the clocks; full returns remain positive', fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .88, 'A fixed; B moves by 0.1% of initial separation  |  Assumed throat  |  Linear response in model time units',
             ha='center', fontsize=10, color='#475569')
    geometry = fig.add_axes([.055, .625, .235, .155]); geometry.set(xlim=(0, 1), ylim=(0, 1)); geometry.axis('off')
    geometry.set_title('A. Prescribed B motion', fontsize=11, pad=11)
    geometry.plot([.12, .85], [.65, .65], color='#94a3b8', lw=1.2)
    u = np.linspace(0, 1, 101); geometry.plot(.12+.73*u, .65-.36*np.sin(np.pi*u), '--', color='#94a3b8', lw=1)
    geometry.plot([.12, .85], [.65, .65], 'o', color=colors[0], ms=7)
    geometry.text(.12, .96, 'A fixed', ha='center', fontsize=9)
    geometry.text(.85, .96, 'B moves', ha='center', fontsize=9)
    geometry.annotate('', xy=(.98, .80), xytext=(.73, .80), arrowprops={'arrowstyle': '<->', 'color': colors[1], 'lw': 1.2})
    geometry.text(.50, .20, 'Assumed throat', ha='center', va='center', fontsize=8.5, color='#64748b')
    geometry.text(.50, -.06, 'Schematic, not to scale', ha='center', fontsize=8.5, color='#64748b')
    motion = fig.add_axes([.055, .33, .235, .225]); mask = (time >= 0)&(time <= 1.1*par['duration'])
    motion.plot(time[mask], 100*data['displacement'][mask]/par['exterior_delay'], color=colors[1], lw=1.7)
    motion.set(xlabel='Exterior time', ylabel='B displacement / D (%)', xlim=(0, 1.1*par['duration']), ylim=(-.004, .115))
    source = fig.add_axes([.385, .585, .255, .195]); source.set_title('B. Quantum source and clock response', fontsize=11, pad=11)
    clock = fig.add_axes([.385, .33, .255, .185]); physical = time >= 0
    for side, color, key in zip(('A', 'B'), colors, ('induced_plus_source', 'induced_minus_source')):
        source.plot(time[physical], data[key][physical]*1e3, color=color, lw=1.25, label=side)
        clock.plot(time[physical], (data[f'clock_rate_{side}'][physical]-1)*1e3, color=color, lw=1.25)
    source.set(ylabel='Source change\n(×10⁻³)', xlim=(0, time[-1])); source.tick_params(labelbottom=False)
    source.legend(loc='upper right', ncol=2, frameon=False, fontsize=8)
    clock.set(ylabel='Clock-rate change\n(×10⁻³)', xlabel='Exterior time', xlim=(0, time[-1]))
    returns = fig.add_axes([.735, .33, .24, .45]); returns.set_title('C. Full-return time change', fontsize=11, pad=11)
    for side, color in zip(('A', 'B'), colors):
        returns.plot(launch, data[f'return_change_{side}']*1e3, color=color, lw=1.4, label=f'Launch at {side}')
    returns.set(xlabel='Launch time', ylabel='Return change\n(×10⁻³ model units)', xlim=(0, launch[-1]))
    returns.legend(loc='upper right', frameon=False, fontsize=8)
    for ax in (motion, source, clock, returns):
        ax.axhline(0, color='#94a3b8', lw=.7); ax.spines[['top', 'right']].set_visible(False); ax.grid(alpha=.12)
    fig.text(.172, .245, f'Prescribed duration: {par["duration"]:g}.\nMotion generation is not solved.', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.512, .245, 'Source-chart changes alone\ndo not establish emitted radiation.', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.855, .245, f'Baseline return: {period:.4f}\nShortest sampled return: {minimum_return:.4f} > 0', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.055, .130, 'Return curves include direct exterior-path change and induced clock response; all transit times are first-order predictions.',
             fontsize=10, color='#475569')
    fig.text(.055, .079, 'Physical mouth matching, exterior gravity and preparation remain unresolved; finite sampled behavior does not prove stability.',
             fontsize=9.5, color='#475569')
    fig.text(.055, .028, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'moving_exterior_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def finite_motion_progress():
    """Display the saved finite-motion energy balance and complete return times."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'moving_exterior_nonlinear_results.json').read_text(encoding='utf-8'))
    verification = record['verification']
    assert verification['all_declared_checks_pass']
    assert record['exact_map_controls']['all_declared_checks_pass']
    assert record['energy_identity_controls']['all_declared_checks_pass']
    assert all(v == '0' for v in record['energy_identity_controls']['residuals'].values())
    assert max(verification['moderate_timestep_differences'].values()) < verification['maximum_step_threshold']
    fine = record['cases'][-1]; par = fine['parameters']; data = {k: np.asarray(v) for k, v in fine['plot'].items()}
    assert par['amplitude'] == .05 and par['maximum_step'] == .025 and fine['censored_return_probes'] == 0
    assert np.isclose(par['length'], np.pi+par['delay'], rtol=1e-14, atol=0)
    time = data['lab_time']; launch = data['return_launch_lab_time']
    assert np.all(np.diff(time) > 0) and np.all(np.diff(launch) > 0)
    assert fine['minimum_sampled_proper_clock_rate'] > 0 and fine['minimum_sampled_lab_clock_rate'] > 0
    bulk_change = data['reduced_bulk_energy']-2*par['delay']*par['q0']
    total = data['improved_total_energy_change']; drive = data['reduced_drive']
    sum_scale = max(np.max(np.abs(bulk_change)), np.max(np.abs(data['improved_boundary_energy'])), 1)
    assert np.max(np.abs(total-bulk_change-data['improved_boundary_energy'])) < 8*np.finfo(float).eps*sum_scale
    assert np.array_equal(total-drive, data['energy_balance_residual'])
    residual = np.max(np.abs(total-drive))
    assert np.isclose(residual, fine['maximum_reduced_balance_residual'], rtol=1e-13, atol=0)
    assert np.isclose(drive[-1], fine['final_reduced_drive'], rtol=1e-14, atol=0)
    returns = np.vstack([data['return_A'], data['return_B']])
    assert np.all(np.isfinite(returns)) and np.all(returns > 0)
    assert np.allclose(returns.min(axis=1), fine['minimum_sampled_complete_return'], rtol=1e-13, atol=0)
    minimum = returns.min(); maximum = returns.max()
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.5})
    fig = plt.figure(figsize=(14, 6.7)); colors = ['#007e87', '#7951b0']
    fig.suptitle('Finite prescribed motion: energy balance closes; no past return in this test', fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .878, 'Composed rays and proper clocks  |  Prescribed 5% excursion  |  Assumed throat and flat exterior',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.045, .355, .27, .415]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Moving exterior energy domain', fontsize=11, pad=11)
    ax.add_patch(plt.Rectangle((.14, .49), .72, .20, facecolor='#edf7f8', edgecolor='none'))
    ax.plot([.14, .86], [.59, .59], color='#94a3b8', lw=1.1)
    ax.plot([.14, .86], [.59, .59], 'o', color=colors[0], ms=8)
    ax.plot([.80, .80], [.43, .75], '--', color='#94a3b8', lw=1)
    ax.text(.14, .90, 'A fixed', ha='center', fontsize=10)
    ax.text(.84, .90, 'B prescribed', ha='center', fontsize=10)
    ax.annotate('', xy=(.98, .80), xytext=(.74, .80), arrowprops={'arrowstyle': '<->', 'color': colors[1], 'lw': 1.3})
    ax.text(.49, .71, 'Exterior field energy', ha='center', fontsize=9, color='#475569')
    u = np.linspace(0, 1, 101); ax.plot(.14+.72*u, .59-.24*np.sin(np.pi*u), '--', color='#94a3b8', lw=1)
    ax.text(.5, .29, 'Assumed throat; schematic only', ha='center', fontsize=9, color='#64748b')
    ax.text(.5, .085, 'Exterior energy change\n+ boosted boundary charge', ha='center', va='center', fontsize=10,
            linespacing=1.4, color='#334155', bbox={'boxstyle': 'round,pad=.4', 'facecolor': '#f3f6f9', 'edgecolor': '#d4dce5'})
    energy = fig.add_axes([.395, .355, .255, .415]); energy.set_title('B. The reduced balance closes', fontsize=11, pad=11)
    energy.plot(time, total, color=colors[0], lw=2.3, label='Total energy change')
    energy.plot(time, drive, color='#aa7900', ls='--', lw=1.5, label='Integrated drive')
    energy.axvline(par['duration'], color='#94a3b8', ls=':', lw=1)
    energy.set(xlabel='Lab time (model units)', ylabel='Dimensionless reduced energy', xlim=(0, time[-1]))
    energy.legend(loc='upper right', frameon=False, fontsize=8.5)
    ret = fig.add_axes([.755, .355, .22, .415]); ret.set_title('C. Complete elapsed return time', fontsize=11, pad=11)
    for side, values, color in zip(('A', 'B'), returns, colors):
        ret.plot(launch, values, color=color, lw=1.5, label=f'Launch {side}')
    ret.axhline(par['length'], color='#94a3b8', ls='--', lw=1)
    ret.set(xlabel='Launch time (model units)', ylabel='Elapsed return time (model units)', xlim=(0, launch[-1]),
            ylim=(minimum-.008, maximum+.015))
    ret.ticklabel_format(axis='y', useOffset=False); ret.legend(loc='upper right', frameon=False, fontsize=8)
    for chart in (energy, ret):
        chart.spines[['top', 'right']].set_visible(False); chart.grid(alpha=.12); chart.set_axisbelow(True)
    fig.text(.18, .255, f'Motion lasts {par["duration"]:g} model units, then\nreturns to the initial position.', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.5225, .255, f'Final reduced transfer: {drive[-1]:.4g}\nSampled balance residual: {residual:.3g}', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.865, .255, f'Zoomed vertical scale; all returns positive.\nMinimum {minimum:.4f}; baseline {par["length"]:.4f}.', ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.4)
    fig.text(.045, .133, 'Energy uses the M₀ = 0 boundary-charge reference; this dimensionless balance is not an apparatus energy budget in joules.',
             fontsize=10, color='#475569')
    fig.text(.045, .080, 'Physical mouth matching, motion generation, material preparation and nonlinear stability remain unresolved.',
             fontsize=10, color='#475569')
    fig.text(.045, .028, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'finite_motion_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def magnetic_polarization_progress():
    """Separate the local polarization solution from the classical-throat trace audit."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'magnetic_polarization_results.json').read_text(encoding='utf-8'))
    verification = record['verification']; par = record['parameters']; product = record['product_solution']
    assert verification['all_declared_checks_pass']
    assert all(v == '0' for v in verification['symbolic_residuals'].values())
    assert all(row['radial_nec_times_8piG_re2'] == 0 for row in record['finite_collar_table'])
    epsilon = par['epsilon']; charge = par['magnetic_charge_q']; coupling = par['gauge_coupling_at_Xstar']
    assert charge == 100000 and np.isclose(epsilon, coupling**2/(24*np.pi**2), rtol=1e-13, atol=0)
    radius_ratio = product['ads_radius_over_re']/product['sphere_radius_over_re']
    assert np.isclose(radius_ratio**2, 1/(1-2*epsilon), rtol=1e-13, atol=0)
    trace = next(row for row in record['classical_rn_trace_audit']['table'] if row['radius_over_re'] == 1)
    gauge = trace['gauge_trace_times_re4']; gravity = trace['gravitational_trace_times_re4']
    lll = trace['lowest_landau_trace_times_re4']; missing = trace['required_additional_trace_times_re4']
    assert np.isclose(gauge, -charge**2/(48*np.pi**2), rtol=1e-13, atol=0)
    assert np.isclose(gravity, 11/(720*np.pi**2), rtol=1e-13, atol=0)
    assert np.isclose(lll, -charge/(48*np.pi**2), rtol=1e-13, atol=0)
    assert np.isclose(missing, gravity-lll, rtol=1e-13, atol=0)
    assert np.isclose(trace['full_required_trace_times_re4'], gauge+gravity, rtol=1e-13, atol=0)
    assert np.isclose(trace['local_log_plus_lll_trace_times_re4'], gauge+lll, rtol=1e-13, atol=0)
    values = np.array([gravity, lll, missing])/abs(lll)
    assert 0 < values[0] < 1e-4 and values[1] == -1 and 1 < values[2] < 1.0001
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 6.8))
    fig.suptitle('Local magnetic polarization alone cannot open a regular static throat', fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .875, f'Two separate diagnostics  |  One massless Dirac species  |  q = {charge:,}; g = {coupling:.5f}',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.05, .345, .39, .43]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('A. Background change versus opening stress', fontsize=11.5, pad=14)
    ax.add_patch(plt.Rectangle((.025, .58), .95, .35, facecolor='#edf7f8', edgecolor='#c4dfe2', lw=1))
    ax.text(.5, .855, 'Local magnetic vacuum polarization', ha='center', va='center', fontsize=10.5, color='#334155')
    ax.text(.5, .742, r'$\rho+p_r=0$', ha='center', va='center', fontsize=17, color='#007e87')
    ax.text(.5, .625, 'Changes the background; supplies no opening stress.', ha='center', va='center', fontsize=9.5, color='#475569')
    ax.add_patch(plt.Rectangle((.025, .12), .95, .35, facecolor='#f3eff9', edgecolor='#d8cce9', lw=1))
    ax.text(.5, .398, 'Loop Casimir contribution', ha='center', va='center', fontsize=10.5, color='#334155')
    ax.text(.5, .282, r'Candidate $\rho+p_r<0$', ha='center', va='center', fontsize=15, color='#7951b0')
    ax.text(.5, .165, 'A complete four-dimensional source is still needed.', ha='center', va='center', fontsize=9.5, color='#475569')
    ax = fig.add_axes([.565, .345, .39, .43]); xx = np.arange(3)
    ax.bar(xx, values, width=.53, color=['#94a3b8', '#7951b0', '#b38516'])
    ax.plot(0, values[0], 'o', color='#64748b', ms=5)
    ax.axhline(0, color='#475569', lw=1)
    ax.text(0, .11, f'+{values[0]:.2e}', ha='center', fontsize=10, color='#475569')
    ax.text(1, -1.12, '−1', ha='center', va='top', fontsize=11, color='#7951b0')
    ax.text(2, 1.09, f'+{values[2]:.6f}', ha='center', va='bottom', fontsize=10, color='#8a650f')
    ax.set(title='B. A missing trace of comparable size', ylabel='Trace / |lowest-mode trace|',
           xticks=xx, xticklabels=['Full target\nminus gauge', 'Lowest magnetic\nmodes (LLL)', 'Required\nadditional trace'],
           ylim=(-1.4, 1.45), xlim=(-.55, 2.55))
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.14); ax.set_axisbelow(True)
    fig.text(.245, .25, f'Local-model radius ratio: Lₐ / Rₛ = {radius_ratio:.9f}\nThis background shift does not demonstrate an open throat.',
             ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.45)
    fig.text(.76, .25, 'Classical equal-radius RN throat, not the corrected metric.\nCommon gauge trace removed; vertical units are |LLL trace|.',
             ha='center', va='top', fontsize=9.5, color='#475569', linespacing=1.45)
    fig.text(.05, .13, 'The additional trace is a requirement, not an identified source; matching a trace does not determine a conserved stress tensor.',
             fontsize=9.5, color='#475569')
    fig.text(.05, .08, 'Renormalization, double counting, full four-dimensional matching and stability remain unresolved.', fontsize=10, color='#475569')
    fig.text(.05, .03, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'magnetic_polarization_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def dirac_flux_progress():
    """Compare a finite Dirac flux contrast against two specified local comparators."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'dirac_flux_response_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and all(record['tests'].values())
    cases = sorted(record['cases'], key=lambda c: c['q']); flux = np.array([c['q'] for c in cases])
    assert np.array_equal(flux, [2, 3])
    thresholds = record['thresholds']
    for case in cases:
        for field, limit in [('response_main_to_fine_absolute_change', 'response_refinement_absolute'),
                             ('action_main_to_fine_absolute_change', 'action_refinement_absolute'),
                             ('direct_vs_ibp_absolute_difference', 'direct_vs_ibp_absolute'),
                             ('richardson_vs_ibp_absolute_difference', 'richardson_vs_ibp_absolute')]:
            assert case[field] < thresholds[limit]
    signed = np.array([[c['local_log_response'], c['local_plus_minimal_trace_repair_response'], c['response']] for c in cases])
    assert np.allclose(signed[:, 0], -flux**2*np.log(flux)/(24*np.pi**2), rtol=1e-13, atol=0)
    assert np.allclose(signed[:, 1], signed[:, 0]+(flux-flux**2)/(24*np.pi**2), rtol=1e-13, atol=0)
    assert np.allclose(signed[:, 2], [-.03242773217598199, -.1031821229596531], rtol=1e-9, atol=0)
    assert np.all(signed < 0) and np.all(np.diff(np.abs(signed), axis=1) > 0)
    magnitudes = np.abs(signed); missed = 100*(1-magnitudes[:, 1]/magnitudes[:, 2])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13.2, 7.1))
    fig.suptitle('A minimal trace repair still misses part of the shape response', fontsize=18, fontweight='bold', y=.955)
    fig.text(.5, .88, 'One minimal massless Dirac field  |  Controlled low-flux tests  |  Euclidean product geometry',
             ha='center', fontsize=10, color='#475569')
    shape = fig.add_axes([.045, .38, .39, .405]); shape.set(xlim=(0, 1.45), ylim=(0, 1)); shape.set_aspect('equal'); shape.axis('off')
    shape.set_title('A. Opposite radius changes, fixed volume', fontsize=11.5, pad=12)
    shape.text(.28, .93, r'$H^2$ factor', ha='center', fontsize=11, color='#334155')
    shape.text(1.13, .93, r'$S^2$ sphere', ha='center', fontsize=11, color='#334155')
    u = np.linspace(-.18, .18, 101)
    for offset in (-.14, 0, .14):
        shape.plot(.28+u, .55+offset+.7*u*u, color='#007e87', lw=1.25)
        shape.plot(.28+offset+.7*u*u, .55+u, color='#007e87', lw=1.25)
    shape.annotate('', xy=(.39, .32), xytext=(.17, .32), arrowprops={'arrowstyle': '<->', 'color': '#b8c2cc', 'lw': 1})
    shape.annotate('', xy=(.50, .27), xytext=(.06, .27), arrowprops={'arrowstyle': '<->', 'color': '#007e87', 'lw': 1.5})
    shape.text(.28, .18, 'a grows', ha='center', fontsize=10, color='#007e87')
    shape.text(.28, .06, r'$a=e^{\eta}$', ha='center', fontsize=11)
    angle = np.linspace(0, 2*np.pi, 161)
    shape.plot(1.13+.22*np.cos(angle), .55+.22*np.sin(angle), '--', color='#b8c2cc', lw=1)
    shape.fill(1.13+.15*np.cos(angle), .55+.15*np.sin(angle), color='#f1ecf8')
    shape.plot(1.13+.15*np.cos(angle), .55+.15*np.sin(angle), color='#7951b0', lw=1.6)
    shape.plot(1.13+.15*np.cos(angle), .55+.045*np.sin(angle), color='#7951b0', lw=.8)
    shape.plot(1.13+.045*np.cos(angle), .55+.15*np.sin(angle), color='#7951b0', lw=.8)
    shape.annotate('', xy=(.985, .55), xytext=(.87, .55), arrowprops={'arrowstyle': '->', 'color': '#7951b0'})
    shape.annotate('', xy=(1.275, .55), xytext=(1.39, .55), arrowprops={'arrowstyle': '->', 'color': '#7951b0'})
    shape.text(1.13, .18, 'r shrinks', ha='center', fontsize=10, color='#7951b0')
    shape.text(1.13, .06, r'$r=e^{-\eta}$', ha='center', fontsize=11)
    shape.text(.72, .55, r'$\times$', ha='center', va='center', fontsize=22, color='#64748b')
    bars = fig.add_axes([.555, .375, .405, .41]); bars.set_title('B. Magnitudes of the finite shape response', fontsize=11.5, pad=12)
    xx = np.arange(2); width = .23
    for j, (label, color) in enumerate(zip(('Local', '+ minimal trace repair', 'Full spectrum'), ('#94a3b8', '#b38516', '#007e87'))):
        positions = xx+(j-1)*width
        bars.bar(positions, magnitudes[:, j], width=width*.91, label=label, color=color)
        for xpos, value in zip(positions, magnitudes[:, j]):
            bars.text(xpos, value+.0026, f'{value:.4f}', ha='center', va='bottom', fontsize=9, color='#334155')
    bars.set(xticks=xx, xticklabels=[r'$q=2$', r'$q=3$'], ylabel=r'Dimensionless $|R_q|$', ylim=(0, .128), xlim=(-.52, 1.52))
    bars.legend(loc='upper left', frameon=False, fontsize=9)
    bars.spines[['top', 'right']].set_visible(False); bars.grid(axis='y', alpha=.14); bars.set_axisbelow(True)
    fig.text(.24, .285, r'Fixed volume element: $a^2r^2=1$; reference $a=r=1$.', ha='center', fontsize=10, color='#475569')
    fig.text(.24, .242, 'Schematic mathematical deformation, not a formed tunnel.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.7575, .285, f'Minimal repair misses {missed[0]:.0f}% (q = 2) and {missed[1]:.0f}% (q = 3).', ha='center', fontsize=10, color='#475569')
    fig.text(.7575, .242, 'Signed responses are negative; bars show their magnitudes.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.5, .177, r'$C_q=w(q)-q^2w(1)+(q^2-1)w(0),\qquad R_q=\left.dC_q/d\eta\right|_{\eta=0}$',
             ha='center', fontsize=13, color='#334155')
    fig.text(.045, .116, 'w is action density; q = 1 is the fixed reference. This combination is not a single-q pressure or extra throat-opening stress.',
             fontsize=9.7, color='#475569')
    fig.text(.045, .073, 'Next: a conserved stress tensor and physical throat matching. These low-flux comparisons give no high-flux error bar.',
             fontsize=9.7, color='#475569')
    fig.text(.045, .029, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.', fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'dirac_flux_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def high_flux_progress():
    """Display the saved high-flux contrast and the resolved linear correction."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'dirac_high_flux_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) == 7 and all(record['tests'].values())
    asym = record['asymptotic']; c0 = asym['quadratic_coefficient']
    ch = asym['linear_hyperbolic']; cs = asym['linear_spherical']; c1 = asym['linear_total']
    assert np.isclose(c1, ch+cs, rtol=1e-14, atol=0) and c0 < 0 and ch < 0 and cs < 0
    cases = sorted((c for c in record['cases'] if c['q'] >= 100), key=lambda c: c['q'])
    flux = np.array([c['q'] for c in cases], dtype=float)
    assert np.array_equal(flux, [100, 1000, 10000, 100000])
    coefficient = np.array([c['fine']['curvature_coefficient'] for c in cases])
    calculated = np.array([(c['fine']['response_over_q2']-c0)*c['q'] for c in cases])
    assert np.allclose(coefficient, calculated, rtol=1e-13, atol=0)
    assert all(c['fine']['direct_vs_ibp_over_q2'] < 2e-13 for c in cases)
    difference = np.abs(coefficient-c1)
    refinement = np.array([c['linear_coefficient_refinement_change'] for c in cases])
    assert np.all(difference > 0) and np.all(refinement > 0)
    assert np.isclose(coefficient[-1], c1, rtol=2e-6, atol=0)
    target = flux[-1]; parts = np.abs(np.array([ch, cs]))*target
    assert np.allclose(parts, [844.3431970194815, 925.5221127717903], rtol=1e-13, atol=0)
    fraction = 100*parts[1]/sum(parts); ppm = 1e6*parts[1]/abs(c0*target**2)
    assert 52 < fraction < 53 and 2.9 < ppm < 3
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.1))
    fig.suptitle('An omitted quantum contribution survives at the model’s magnetic flux',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .88, 'Minimal massless Dirac field  |  Homogeneous product vacuum  |  Full angular spectrum retained',
             ha='center', fontsize=10, color='#475569')
    bar = fig.add_axes([.075, .39, .39, .40])
    bar.set_title('A. Linear correction at q = 100,000', fontsize=11.5, pad=15)
    bar.barh(.64, parts[0], height=.29, color='#b38516')
    bar.barh(.64, parts[1], left=parts[0], height=.29, color='#007e87')
    bar.text(parts[0]/2, .64, f'Minimal repair\n{parts[0]:.3f}', ha='center', va='center',
             color='white', fontsize=11, linespacing=1.5)
    bar.text(parts[0]+parts[1]/2, .64, f'Higher angular modes\n{parts[1]:.3f}', ha='center', va='center',
             color='white', fontsize=10.5, linespacing=1.5)
    bar.text(sum(parts)/2, .96, 'Known quadratic part removed', ha='center', fontsize=10, color='#475569')
    bar.text(sum(parts)/2, .29, f'Total linear correction magnitude: {sum(parts):.3f}',
             ha='center', fontsize=10, color='#334155')
    bar.set(xlim=(0, 2000), ylim=(0, 1.15), yticks=[], xticks=[0, 500, 1000, 1500, 2000],
            xlabel=r'Response magnitude [$\hbar c/L^4$]')
    bar.spines[['top', 'right', 'left']].set_visible(False)
    bar.text(.5, -.20, 'Both signed contributions are negative.', transform=bar.transAxes,
             ha='center', fontsize=9.5, color='#475569')
    chart = fig.add_axes([.62, .39, .335, .40])
    chart.set_title('B. Numerical approach to the analytic coefficient', fontsize=11.5, pad=15)
    chart.loglog(flux, difference, 'o-', color='#007e87', lw=1.8, ms=6, label='Distance from analytic limit')
    chart.loglog(flux, refinement, 's--', color='#94a3b8', lw=1.3, ms=4, label='Main-to-fine numerical change')
    chart.scatter([target], [difference[-1]], s=105, facecolors='none', edgecolors='#007e87', zorder=5)
    chart.annotate('Model flux', xy=(target, difference[-1]), xytext=(-8, -21), textcoords='offset points',
                   ha='right', fontsize=9.5, color='#007e87')
    chart.set(xlabel=r'Magnetic flux label $q$ (log scale)', ylabel='Absolute coefficient difference',
              xlim=(70, 150000), ylim=(2e-13, 1e-6))
    chart.legend(loc='upper right', frameon=False, fontsize=8.8)
    chart.spines[['top', 'right']].set_visible(False); chart.grid(alpha=.12); chart.set_axisbelow(True)
    fig.text(.27, .255, f'Omitted: {fraction:.1f}% of this linear correction', ha='center',
             fontsize=12, fontweight='bold', color='#007e87')
    fig.text(.27, .213, f'Only {ppm:.2f} ppm of the leading quadratic response.', ha='center', fontsize=10, color='#475569')
    fig.text(.7875, .265, r'Extracted coefficient: $q(R/q^2-C_0)$', ha='center', fontsize=10, color='#475569')
    fig.text(.7875, .223, f'Analytic limit C₁ = {c1:.12f}', ha='center', fontsize=10, color='#475569')
    fig.text(.7875, .181, 'Observed convergence, not a certified total error bound.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.05, .134, r'$R$ is the fixed-volume shape derivative of $3W(3q)-8W(2q)+5W(q)$; these bars are not a single-flux pressure.',
             fontsize=10, color='#475569')
    fig.text(.05, .085, r'Homogeneous vacuum: $\rho+p_r=0$. This correction supplies no throat-opening stress or material prescription.',
             fontsize=10, color='#475569')
    fig.text(.05, .032, 'No constructed stable medium or 13 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'high_flux_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def matched_source_progress():
    """Show the matched leading source and its unequal-radius spectral check."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'dirac_matched_source_results.json').read_text(encoding='utf-8'))
    high = json.loads((root/'dirac_high_flux_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) == 8 and all(record['tests'].values())
    assert high['passed'] and all(high['tests'].values())
    assert all(value == '0' for value in record['symbolic_residuals'].values())
    ch = high['asymptotic']['linear_hyperbolic']; cs = high['asymptotic']['linear_spherical']
    assert record['d_gap'] > 0 and np.isclose(cs, -8*record['d_gap'], rtol=1e-13, atol=0)
    cases = sorted((row for row in record['cases'] if row['q'] == 100000), key=lambda row: row['eta'])
    assert len(cases) == 5
    ratio = np.array([row['fine']['a_over_r'] for row in cases])
    measured = np.array([row['fine']['linear_coefficient'] for row in cases])
    predicted = np.array([row['fine']['predicted_linear_coefficient'] for row in cases])
    assert np.allclose(ratio, np.exp(2*np.array([row['eta'] for row in cases])), rtol=1e-13, atol=0)
    assert np.allclose(predicted, ch+cs*ratio**2, rtol=1e-13, atol=0)
    assert np.all(measured < 0) and np.all(np.diff(measured) < 0)
    mismatch = float(np.max(np.abs(measured-predicted)))
    assert 0 < mismatch < 6e-8
    assert all(row['fine']['direct_vs_ibp'] < 2e-13 and row['linear_coefficient_refinement'] < 3e-9 for row in cases)
    assert all(row['absolute_error'] < 1e-11 for row in record['action_derivative_controls'])
    stress = np.array([1., -1., 1.]); assert stress[0]+stress[1] == 0
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.1))
    fig.suptitle('A matched quantum source changes geometry without opening the throat',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Leading O(q) bulk action  |  Fixed smooth geometry  |  Model validation, not a material construction',
             ha='center', fontsize=10, color='#475569')
    ax = fig.add_axes([.075, .385, .36, .40]); positions = np.arange(3)
    ax.set_title('A. Added energy and pressures', fontsize=11.5, pad=15)
    ax.bar(positions, stress, width=.58, color=['#b38516', '#7951b0', '#007e87'])
    ax.axhline(0, color='#64748b', lw=1)
    for x, value, label in zip(positions, stress, (r'$+u$', r'$-u$', r'$+u$')):
        ax.text(x, value+(.09 if value > 0 else -.09), label, ha='center',
                va='bottom' if value > 0 else 'top', fontsize=14, color='#334155')
    ax.set(xticks=positions, xticklabels=['Energy density\nρ', 'Radial pressure\npᵣ', 'Angular pressure\npₜ'],
           ylabel='Added stress / u', ylim=(-1.45, 1.45), yticks=[-1, 0, 1], xlim=(-.6, 2.6))
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', alpha=.13); ax.set_axisbelow(True)
    chart = fig.add_axes([.59, .385, .37, .40])
    chart.set_title('B. Unequal-radius test at q = 100,000', fontsize=11.5, pad=15)
    grid = np.linspace(ratio[0], ratio[-1], 201)
    chart.plot(grid, ch+cs*grid**2, color='#007e87', lw=2, label='Matched leading prediction')
    chart.scatter(ratio, measured, s=60, facecolors='white', edgecolors='#7951b0', linewidths=1.6,
                  zorder=4, label='Full spectrum')
    chart.axhline(ch, color='#b38516', lw=1.3, ls='--', label='Minimal repair alone')
    chart.set(xlabel=r'Curvature-radius ratio $a/r$', ylabel='Extracted linear shape coefficient',
              xlim=(.46, 1.91), ylim=(-.044, -.004), xticks=[.5, 1, 1.5, 1.8])
    chart.legend(loc='lower left', frameon=False, fontsize=9)
    chart.spines[['top', 'right']].set_visible(False); chart.grid(alpha=.13); chart.set_axisbelow(True)
    fig.text(.255, .265, r'$u=d_{\rm gap}q/R^4>0$', ha='center', fontsize=12, color='#475569')
    fig.text(.255, .210, r'$\rho+p_r=0$: no added radial opening stress.', ha='center',
             fontsize=11.5, fontweight='bold', color='#007e87')
    fig.text(.775, .265, r'Prediction: $C_{1H}+C_{1S}(a/r)^2$', ha='center', fontsize=12, color='#475569')
    fig.text(.775, .213, f'Largest sampled coefficient mismatch: {mismatch:.2e}',
             ha='center', fontsize=9.7, color='#475569')
    fig.text(.05, .138, r'Right panel: fixed-volume derivative of $3W(3q)-8W(2q)+5W(q)$, with the known quadratic term removed.',
             fontsize=10, color='#475569')
    fig.text(.05, .089, 'The source shifts the product reference. Physical interfaces, the returning-path state and full stability remain unresolved.',
             fontsize=9.8, color='#475569')
    fig.text(.05, .032, 'No constructed stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'matched_source_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def return_closure_progress():
    """Illustrate the saved conditional loop closure and path-scale sensitivity."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'dirac_return_closure_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) == 8 and all(record['tests'].values())
    fine = max(record['records'], key=lambda row: (row['order'], row['refinement']))
    state = max(fine['roots'], key=lambda row: row['cutoff'])
    assert fine['order'] == 128 and state['cutoff'] == 100
    assert np.isclose(state['omega_squared'], record['corrected_root'], rtol=1e-14, atol=0)
    weights = np.asarray(state['weights']); periods = np.asarray(state['periods'])
    assert np.all(weights > 0) and np.all(periods > 0) and np.isclose(sum(weights), 1, rtol=1e-13, atol=0)
    assert np.isclose(np.dot(weights, (2*np.pi/periods)**2), state['omega_squared'], rtol=1e-12, atol=0)
    throat = state['throat_time']; exterior = state['shortest_exterior_time']
    ratio = state['throat_to_exterior_ratio']
    assert throat > exterior > 0 and np.isclose(ratio, throat/exterior, rtol=1e-13, atol=0)
    assert 4.62 < ratio < 4.63 and record['root_refinement'] < 2e-10
    rows = fine['path_scale_sensitivity']+[{'path_scale': 1., 'central_flare_over_alpha': state['background']['central_flare_over_alpha']}]
    rows = sorted(rows, key=lambda row: row['path_scale'])
    scale = np.array([row['path_scale'] for row in rows]); flare = np.array([row['central_flare_over_alpha'] for row in rows])
    assert np.array_equal(scale, [.999, 1., 1.001]) and np.all(flare > 0) and np.all(np.diff(flare) < 0)
    assert np.isclose(record['parameters']['alpha'], 1e-9, rtol=1e-14, atol=0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.1))
    fig.suptitle('Conditional averaged-source opening; the throat route remains slower',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Finite spherical collar + prescribed exterior paths  |  Common optical units, not clock days',
             ha='center', fontsize=10, color='#475569')
    diagram = fig.add_axes([.05, .385, .425, .405])
    diagram.set(xlim=(-.05, 1.7), ylim=(-.03, 1.03)); diagram.set_aspect('equal'); diagram.axis('off')
    diagram.set_title('A. One throat pass + one exterior return', fontsize=11.5, pad=15)
    angle = np.linspace(0, np.pi, 181)
    x = .825+.625*np.cos(angle)
    diagram.plot(x, .5+.29*np.sin(angle), color='#007e87', lw=5, solid_capstyle='round')
    diagram.plot(x, .5-.29*np.sin(angle), color='#b38516', lw=2, ls='--')
    diagram.annotate('', xy=(.91, .79), xytext=(.73, .79),
                     arrowprops={'arrowstyle': '-|>', 'color': '#004b50', 'lw': 1.5, 'mutation_scale': 14})
    diagram.annotate('', xy=(.73, .21), xytext=(.91, .21),
                     arrowprops={'arrowstyle': '-|>', 'color': '#8a650f', 'lw': 1.5, 'mutation_scale': 14})
    for xpos, letter in ((.2, 'A'), (1.45, 'B')):
        diagram.add_patch(plt.Circle((xpos, .5), .075, facecolor='white', edgecolor='#475569', lw=1.7, zorder=4))
        diagram.text(xpos, .5, letter, ha='center', va='center', fontsize=11, color='#334155', zorder=5)
        diagram.text(xpos, .63, 'Mouth', ha='center', fontsize=9.5, color='#475569')
    diagram.text(.825, .93, f'Modeled throat time: {throat:.6f}', ha='center', fontsize=11, color='#007e87')
    diagram.text(.825, .035, f'Shortest exterior time: {exterior:.6f}', ha='center', fontsize=10.5, color='#8a650f')
    chart = fig.add_axes([.62, .385, .335, .405])
    chart.set_title('B. Central opening in the averaged model', fontsize=11.5, pad=15)
    chart.plot(scale, flare, 'o-', color='#007e87', lw=1.8, ms=7)
    for xx, yy in zip(scale, flare):
        chart.annotate(f'{yy:.8f}', (xx, yy), xytext=(0, 12), textcoords='offset points', ha='center', fontsize=10, color='#334155', bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.5})
    chart.set(xlabel='Prescribed exterior-path scale', ylabel=r'Central radial flare / $\alpha$',
              xticks=scale, xticklabels=['0.999', '1.000', '1.001'],
              xlim=(.99872, 1.00128), ylim=(.20632, .20670))
    chart.ticklabel_format(axis='y', useOffset=False)
    chart.spines[['top', 'right']].set_visible(False); chart.grid(alpha=.13); chart.set_axisbelow(True)
    fig.text(.2625, .285, 'Schematic assumed geometry, not a physical blueprint.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.2625, .23, f'Throat crossing takes {ratio:.3f}× the exterior time.', ha='center',
             fontsize=11.5, fontweight='bold', color='#007e87')
    fig.text(.7875, .285, 'Zoomed vertical scale; all three values are positive.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.7875, .23, f'Closed mean-state parameter: w = {state["omega_squared"]:.10f}',
             ha='center', fontsize=10, color='#475569')
    fig.text(.05, .151, 'NS spin/Wilson phase is assumed. The ±0.1% path changes are sensitivity probes, not physical uncertainty bounds.',
             fontsize=9.8, color='#475569')
    fig.text(.05, .097, 'The complete nonradial exterior, state preparation and dynamical stability remain unresolved.',
             fontsize=10, color='#475569')
    fig.text(.05, .039, 'No constructed stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'return_closure_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def magnetic_path_progress():
    """Display saved first-order exterior paths, with explicit visual enlargement."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'magnetic_path_polarization_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) == 8 and all(record['tests'].values())
    fine = record['records'][-1]; epsilon = record['epsilon']
    assert fine['configuration']['outer'] == 40000 and fine['configuration']['modes'] == 36
    paths = sorted(fine['representative_paths'], key=lambda row: row['nu'])
    nu = np.array([row['nu'] for row in paths]); coefficient = np.array([row['relative_length_coefficient'] for row in paths])
    assert np.array_equal(nu, [.02, .1, .5, 1., 1.5, 1.9]) and np.all(coefficient < 0)
    assert np.allclose(coefficient, [row['delta_length_per_epsilon']/row['base_length'] for row in paths], rtol=1e-13, atol=0)
    ppm = epsilon*coefficient*1e6
    screen = np.interp(nu, fine['nu'], fine['gradient_turning_screen']); flagged = screen >= .1
    assert np.array_equal(flagged, [True, False, False, False, False, False])
    representative = next(row for row in paths if row['nu'] == 1.)
    profile = representative['profile']; v_half = np.asarray(profile['v']); u_half = np.asarray(profile['u'])
    du_half = np.asarray(profile['delta_u'])
    assert v_half.shape == u_half.shape == du_half.shape and np.all(np.diff(v_half) > 0)
    assert np.all((v_half > 0) & (v_half < 1)) and np.all(u_half > 1) and np.all(du_half < 0)
    v = np.r_[-1., -v_half[::-1], v_half, 1.]
    u = np.r_[1., u_half[::-1], u_half, 1.]; du = np.r_[0., du_half[::-1], du_half, 0.]
    enlargement = 500; shown_u = u+enlargement*epsilon*du
    assert np.all(shown_u >= 1)
    rho = .5*np.sqrt((u*u-1)*(1-v*v)); z = .5*u*v
    shown_rho = .5*np.sqrt((shown_u*shown_u-1)*(1-v*v)); shown_z = .5*shown_u*v
    assert np.array_equal(z[[0, -1]], [-.5, .5]) and np.array_equal(shown_z[[0, -1]], [-.5, .5])
    representative_ppm = epsilon*representative['relative_length_coefficient']*1e6
    assert -17.19 < representative_ppm < -17.17
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.1))
    fig.suptitle('Local magnetic polarization slightly shortens the modeled exterior paths',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Computed first-order correction  |  Flat local-action approximation  |  Fixed point-source endpoints',
             ha='center', fontsize=10, color='#475569')
    path_ax = fig.add_axes([.065, .38, .405, .405])
    path_ax.set_title(r'A. Exterior field line $\nu=1$', fontsize=11.5, pad=15)
    path_ax.plot(z, rho, color='#94a3b8', lw=2.5, label='Baseline')
    path_ax.plot(shown_z, shown_rho, color='#007e87', lw=1.8, label='Correction enlarged 500×')
    path_ax.scatter([-.5, .5], [0, 0], s=45, facecolors='white', edgecolors='#334155', zorder=5)
    path_ax.set(xlabel=r'Axial coordinate $z/d$', ylabel=r'Radial coordinate $\rho/d$',
                xlim=(-.61, .61), ylim=(-.065, 1.0))
    path_ax.set_aspect('equal', adjustable='box')
    path_ax.legend(loc='center', bbox_to_anchor=(.5, .45), frameon=False, fontsize=8.7)
    path_ax.spines[['top', 'right']].set_visible(False); path_ax.grid(alpha=.12); path_ax.set_axisbelow(True)
    chart = fig.add_axes([.62, .38, .335, .405])
    chart.set_title('B. Actual predicted length changes', fontsize=11.5, pad=15)
    chart.semilogx(nu, ppm, '--', color='#b8c2cc', lw=1.2)
    chart.scatter(nu[~flagged], ppm[~flagged], s=50, color='#007e87', zorder=3, label='Within chosen screen')
    chart.scatter(nu[flagged], ppm[flagged], s=75, facecolors='white', edgecolors='#b38516',
                  linewidths=1.8, zorder=4, label='Outside chosen screen')
    chart.axhline(0, color='#64748b', lw=.8)
    chart.annotate(f'{representative_ppm:.2f} ppm at ν = 1', xy=(1, representative_ppm),
                   xytext=(.06, .76), textcoords='axes fraction', ha='left', fontsize=10, color='#007e87',
                   arrowprops={'arrowstyle': '-', 'color': '#007e87', 'lw': .8})
    chart.set(xlabel=r'Field-line label $\nu$ (log scale)', ylabel='Relative length change (ppm)',
              xlim=(.014, 2.5), ylim=(-130, 12), xticks=[.02, .1, .5, 1, 2], xticklabels=['0.02', '0.1', '0.5', '1', '2'])
    chart.tick_params(axis='x', labelsize=9)
    chart.legend(loc='lower right', frameon=False, fontsize=9)
    chart.spines[['top', 'right']].set_visible(False); chart.grid(alpha=.12); chart.set_axisbelow(True)
    fig.text(.2675, .275, 'Displacement magnified 500× for visibility.', ha='center',
             fontsize=11, fontweight='bold', color='#007e87')
    fig.text(.2675, .228, 'Assumed mouth points; no finite entrances or throat shown.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.7875, .275, 'No enlargement in this panel; negative means shorter.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.7875, .228, 'The flagged ν = 0.02 result remains a formal prediction.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.05, .157, 'Chosen screen: full gradient parameter < 0.1 at turning points. This is not a uniform physical error bound.',
             fontsize=10, color='#475569')
    fig.text(.05, .102, 'The nonlocal current, curved exterior geometry, finite-mouth matching and stability remain unresolved.',
             fontsize=10, color='#475569')
    fig.text(.05, .041, 'No constructed stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'magnetic_path_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def magnetic_validity_progress():
    """Show full-path derivative diagnostics and their conditional flux screen."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'magnetic_path_validity_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) >= 7 and all(record['tests'].values())
    screen = record['screen']; threshold = screen['eta1_limit']; cut = screen['joint_nu_cut']
    assert record['q'] == 100000 and np.isclose(threshold, .1, rtol=1e-14, atol=0)
    assert np.isclose(screen['eta2_limit'], threshold**2, rtol=1e-14, atol=0)
    assert cut == max(screen['nu_cut_eta1'], screen['nu_cut_eta2'])
    fraction = screen['continuum_flux_below_cut']
    assert np.isclose(fraction, cut/2, rtol=1e-14, atol=0) and .0232 < fraction < .0233
    extrema = sorted(record['extrema'], key=lambda row: row['nu'])
    nu = np.array([row['nu'] for row in extrema])
    maxima = np.array([max(row['max_eta1'], np.sqrt(row['max_eta2'])) for row in extrema])/threshold
    assert np.all(np.diff(nu) > 0) and np.all(np.isfinite(maxima)) and np.all(maxima > 0)
    profiles = [next(row for row in record['profiles'] if row['nu'] == n) for n in (.02, .05, 1.9)]
    counter = record['eta2_endpoint_counterexample']
    numerator, denominator = map(int, counter['exact_q2_eta2_squared'].split('/'))
    assert np.isclose(counter['eta2'], np.sqrt(numerator/denominator)/record['q'], rtol=1e-13, atol=0)
    assert counter['eta2'] > max(counter['eta2_turning'], counter['eta2_mouth'])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.2))
    fig.suptitle(f'At least {100*fraction:.2f}% of the modeled flux fails the chosen derivative screen',
                 fontsize=16.5, fontweight='bold', y=.955)
    fig.text(.5, .885, 'First and second Cartesian derivatives  |  Flat opposite point charges  |  q = 100,000',
             ha='center', fontsize=10, color='#475569')
    left = fig.add_axes([.075, .38, .37, .405])
    left.set_title('A. Diagnostics along three field lines', fontsize=11.5, pad=15)
    handles = []
    for profile, color in zip(profiles, ('#b38516', '#007e87', '#7951b0')):
        v = np.asarray(profile['v']); eta1 = np.asarray(profile['eta1'])/threshold
        eta2 = np.asarray(profile['sqrt_eta2'])/threshold
        assert v.shape == eta1.shape == eta2.shape and np.all(np.diff(v) > 0)
        assert v[0] == 0 and v[-1] < 1 and np.all(eta1 > 0) and np.all(eta2 > 0)
        line, = left.plot(v, eta2, color=color, lw=1.8, label=f'ν = {profile["nu"]:g}')
        left.plot(v, eta1, color=color, lw=1.4, ls='--'); handles.append(line)
    left.axhline(1, color='#475569', lw=1, ls=':')
    left.text(.98, 1.045, 'Chosen screen', ha='right', fontsize=9.5, color='#475569')
    left.set(xlabel='Path coordinate v (turning point → mouth)', ylabel='Derivative diagnostic / 0.1',
             xlim=(0, 1), ylim=(0, 1.7), xticks=[0, .25, .5, .75, 1])
    left.spines[['top', 'right']].set_visible(False); left.grid(alpha=.12); left.set_axisbelow(True)
    right = fig.add_axes([.60, .38, .355, .405])
    right.set_title('B. Largest diagnostic on each full path', fontsize=11.5, pad=15)
    right.axvspan(nu[0], cut, facecolor='#b38516', alpha=.10)
    right.loglog(nu, maxima, color='#007e87', lw=1.8)
    right.axhline(1, color='#475569', lw=1, ls=':')
    right.axvline(cut, color='#b38516', lw=1, ls='--')
    right.scatter([cut], [1], color='#b38516', s=40, zorder=4)
    right.annotate(f'ν = {cut:.5f}', xy=(cut, 1), xytext=(10, 13), textcoords='offset points',
                   fontsize=9.5, color='#8a650f')
    right.text(4e-5, .17, 'Guaranteed failure\nfrom turning-point test', fontsize=9.2, color='#8a650f', linespacing=1.45)
    right.set(xlabel=r'Field-line label $\nu$ (log scale)', ylabel='max(η₁, √η₂) / 0.1',
              xlim=(1e-5, 2.2), ylim=(.1, 100))
    right.spines[['top', 'right']].set_visible(False); right.grid(alpha=.12); right.set_axisbelow(True)
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(.26, .29), ncol=3, frameon=False, fontsize=9.5)
    fig.text(.26, .235, r'Solid: $\sqrt{\eta_2}/0.1$; dashed: $\eta_1/0.1$.', ha='center', fontsize=10, color='#475569')
    fig.text(.7775, .29, f'At least {100*fraction:.2f}% of continuous flux measure dν/2.', ha='center',
             fontsize=11, fontweight='bold', color='#007e87')
    fig.text(.7775, .235, 'This is a flagged flux fraction, not an error fraction.', ha='center', fontsize=9.7, color='#475569')
    fig.text(.05, .166, f'An exact counterexample at ν = {counter["nu"]:g} has η₂ larger inside the path than at either endpoint.',
             fontsize=10, color='#475569')
    fig.text(.05, .111, 'Higher-ν whole-path maxima are numerical scans. These chosen screens prove neither physical accuracy nor stability.',
             fontsize=9.8, color='#475569')
    fig.text(.05, .047, 'No constructed stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'magnetic_validity_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def spinful_return_progress():
    """Plot exact saved transverse modes and a verified short spinful bulk packet."""
    root = Path(__file__).resolve().parent
    record = json.loads((root/'spinful_return_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and len(record['tests']) == 9 and all(record['tests'].values())
    assert record['operator_square_residuals'] == ['0', '0'] and record['marginal_k'] == 0
    charge = record['q']; assert charge == 100000
    selected = [next(p for p in record['profiles'] if p['k'] == k) for k in (1, 16)]
    evolution = next(e for e in record['evolution'] if e['grid'] == 256 and e['step'] == .001)
    assert evolution['mu'] == 16.5 and evolution['finish'] == .4 and evolution['inner'] == .03
    rho = np.asarray(evolution['rho']); z = np.asarray(evolution['z']); density = np.asarray(evolution['final_density'])
    assert density.shape == (len(rho), len(z)) and np.all(np.isfinite(density)) and np.all(density >= 0)
    dr = rho[1]-rho[0]; dz = z[1]-z[0]; R, Z = np.meshgrid(rho, z, indexing='ij')
    final = evolution['history'][-1]
    assert np.isclose(density.sum()*dr*dz, final['norm'], rtol=1e-13, atol=0)
    assert np.isclose((R*density).sum()*dr*dz, final['mean_rho'], rtol=1e-13, atol=0)
    assert np.isclose((Z*density).sum()*dr*dz, final['mean_z'], rtol=1e-13, atol=0)
    assert evolution['maximum_norm_error'] < 2e-11 and evolution['maximum_edge_probability'] < 2e-6
    mu = evolution['mu']; initial = np.exp(2*(-mu*np.log(R)-mu/R+mu-Z*Z/(4*.12**2)))
    initial /= initial.sum()*dr*dz
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.2))
    fig.suptitle('Spinful transverse modes and a computed two-dimensional quantum packet',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Exact transverse benchmark + short local evolution  |  Neither panel gives a return or escape probability',
             ha='center', fontsize=10, color='#475569')
    left = fig.add_axes([.075, .38, .37, .405])
    left.set_title('A. Exact transverse radial modes', fontsize=11.5, pad=15)
    maximum = 0.
    for profile, color in zip(selected, ('#b38516', '#007e87')):
        radius = np.asarray(profile['r_over_qd']); probability = np.asarray(profile['density'])
        peak = .5*np.sqrt(1/profile['mu']**2-1/charge**2)
        assert np.all(radius > 0) and np.all(probability >= 0) and abs(radius[np.argmax(probability)]/peak-1) < .02
        scaled_radius = radius/peak; scaled_density = probability*peak
        maximum = max(maximum, float(scaled_density.max()))
        left.semilogx(scaled_radius, scaled_density, color=color, lw=2, label=f'k = {profile["k"]}')
    left.axvline(1, color='#94a3b8', lw=1, ls=':')
    left.set(xlabel='Radius / each mode’s peak radius', ylabel='Density in the rescaled radius',
             xlim=(.15, 20), ylim=(0, maximum*1.15), xticks=[.2, .5, 1, 2, 5, 10, 20],
             xticklabels=['0.2', '0.5', '1', '2', '5', '10', '20'])
    left.legend(loc='upper right', frameon=False, fontsize=10)
    left.spines[['top', 'right']].set_visible(False); left.grid(alpha=.12); left.set_axisbelow(True)
    right = fig.add_axes([.585, .38, .325, .405])
    right.set_title('B. Evolved packet at t = 0.4 ℓ/c', fontsize=11.5, pad=15)
    mesh = right.pcolormesh(z, rho, density, shading='auto', cmap='viridis', vmin=0,
                           vmax=float(density.max()), rasterized=True)
    right.contour(z, rho, initial, levels=[float(initial.max())/2], colors='white', linewidths=1.4, linestyles='--')
    right.set(xlabel=r'Axial coordinate $z/\ell$', ylabel=r'Radius $\rho/\ell$',
              xlim=(-.9, .6), ylim=(.45, 1.8))
    right.set_aspect('equal', adjustable='box')
    right.text(.22, 1.47, 'Initial outline', color='white', fontsize=9, ha='center')
    color_ax = fig.add_axes([.926, .38, .013, .405]); colorbar = fig.colorbar(mesh, cax=color_ax)
    colorbar.set_label(r'Packet density $\times\ell^2$', fontsize=9.5)
    fig.text(.26, .285, 'Each full curve is normalized; radii are scaled separately.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.26, .235, 'The k = 1 mode has a long radial tail.', ha='center', fontsize=10.5, color='#8a650f')
    fig.text(.7475, .285, 'Dashed: initial half-maximum density; local view shown.', ha='center', fontsize=9.5, color='#475569')
    fig.text(.7475, .235, r'$\mu=16.5,\quad \ell=qd/(2\mu)$; includes both spinor components.',
             ha='center', fontsize=10, color='#475569')
    fig.text(.05, .168, 'The k = 0 slice mode is not normalizable. This does not determine a three-dimensional escape probability.',
             fontsize=10, color='#475569')
    fig.text(.05, .112, 'The packet evolves briefly in a periodic box; no mouth-to-mouth or outgoing-scattering conditions are imposed.',
             fontsize=9.9, color='#475569')
    fig.text(.05, .047, 'No constructed stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'spinful_return_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def scattering_vacuum_progress():
    """Compare saved local toy stresses while keeping physical scattering unresolved."""
    from matplotlib.patches import FancyBboxPatch
    root = Path(__file__).resolve().parent
    record = json.loads((root/'scattering_vacuum_results.json').read_text(encoding='utf-8'))
    incoming = json.loads((root/'incoming_modes_results.json').read_text(encoding='utf-8'))
    assert record['status'] == 'passed' and incoming['passed'] and all(incoming['tests'].values())
    assert incoming['tests']['monopole_stationary_dirac_equation'] and incoming['tests']['positive_spectrum_norm_and_moments']
    assert record['checks']['energy_integrals'] == 6
    assert max(record['checks']['max_absolute_errors'].values()) < 2e-12
    names = ('partial_NS', 'same_probability_positive_phase', 'same_probability_causal_delay')
    cases = [next(row for row in record['examples'] if row['name'] == name) for name in names]
    for row in cases:
        assert row['physical_loop_length'] == 1 and row['return_probability'] == .25
        assert np.isclose(row['return_modulus']**2, row['return_probability'], rtol=1e-14, atol=0)
        assert np.isclose(row['physical_loop_interaction_energy'],
                          2*row['physical_loop_length']*row['physical_loop_Tuu_equals_Tvv'], rtol=1e-14, atol=0)
    stress = np.array([row['physical_loop_Tuu_equals_Tvv'] for row in cases])
    assert stress[0] < stress[2] < 0 < stress[1]
    assert np.isclose(cases[1]['phase_radians'], 2*np.pi/3, rtol=1e-14, atol=0)
    assert cases[2]['causal_vertex_delay'] == .8
    assert np.isclose(stress[2], stress[0]/1.8**2, rtol=1e-14, atol=0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.2))
    fig.suptitle('Return probability alone does not determine vacuum stress',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Conditional model — no physical medium', ha='center',
             fontsize=12, fontweight='bold', color='#007e87')
    left = fig.add_axes([.065, .29, .38, .49]); left.set_axis_off()
    left.set_title('A. The missing physical connection', fontsize=11.5, pad=16)
    stages = [(.78, '#e6f4f3', '#007e87', 'Exact positive-energy reference',
               'Isolated-monopole incoming mode'),
              (.45, '#fff5df', '#a77612', 'Physical exterior return matrix',
               'Amplitude and phase still unsolved'),
              (.12, '#f0eaf8', '#7951b0', 'Conditional vacuum-stress bridge',
               'Toy requires a state and reference')]
    for y, fill, color, title, detail in stages:
        left.add_patch(FancyBboxPatch((.025, y-.105), .95, .21,
                       boxstyle='round,pad=.012,rounding_size=.025',
                       facecolor=fill, edgecolor=color, linewidth=1.25,
                       linestyle='--' if y == .45 else '-'))
        left.text(.5, y+.025, title, ha='center', va='center', fontsize=11, fontweight='bold', color=color)
        left.text(.5, y-.043, detail, ha='center', va='center', fontsize=10, color='#475569')
    for top, bottom in ((.66, .57), (.33, .24)):
        left.annotate('', xy=(.5, bottom), xytext=(.5, top),
                      arrowprops={'arrowstyle': '->', 'color': '#64748b', 'lw': 1.6})
    right = fig.add_axes([.59, .35, .355, .43])
    right.set_title('B. Same toy return probability: 25%', fontsize=11.5, pad=16)
    right.axhline(0, color='#64748b', lw=1)
    colors = ('#007e87', '#b38516', '#7951b0')
    right.bar(np.arange(3), stress, width=.56, color=colors, zorder=3)
    for x, value, color in zip(range(3), stress, colors):
        right.text(x, value+(.009 if value > 0 else -.01), f'{value:+.4f}', ha='center',
                   va='bottom' if value > 0 else 'top', fontsize=10.5, color=color, fontweight='bold')
    right.set(xticks=range(3), xticklabels=['Phase 0\nNo delay', 'Phase 2π/3\nNo delay', 'Phase 0\nDelay 0.8 L'],
              ylabel=r'Local null stress $T_{uu}L^2$ ($\hbar=c=1$)', ylim=(-.185, .105), xlim=(-.65, 2.65))
    right.spines[['top', 'right']].set_visible(False); right.grid(axis='y', alpha=.14); right.set_axisbelow(True)
    right.tick_params(axis='x', length=0, pad=10)
    fig.text(.7675, .235, 'Phase reverses the sign; delay changes the magnitude.',
             ha='center', fontsize=10.3, color='#475569')
    fig.text(.05, .174, 'The reference input and toy are separate tests; no physical complex return amplitude has been computed.',
             fontsize=10, color='#475569')
    fig.text(.05, .117, 'Toy: a filled reservoir sea, NS loop phase and fixed subtraction. This is not the physical throat stress.',
             fontsize=10, color='#475569')
    fig.text(.05, .052, 'No stable medium or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'scattering_vacuum_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def dipole_collar_progress():
    """Compare saved collar residuals with the analytic first-corrector bound."""
    import hashlib
    from matplotlib.patches import Circle, FancyArrowPatch
    root = Path(__file__).resolve().parent
    angular_path = root/'dipole_collar_modes_results.json'
    angular = json.loads(angular_path.read_text(encoding='utf-8'))
    bound = json.loads((root/'collar_dressing_bound_results.json').read_text(encoding='utf-8'))
    assert angular['passed'] and all(angular['tests'].values())
    assert bound['passed'] and all(bound['tests'].values())
    assert bound['input_sha256'][angular_path.name] == hashlib.sha256(angular_path.read_bytes()).hexdigest()
    parameters = angular['saved_collar_parameters']; assert parameters == bound['parameters']
    assert parameters['q'] == 100000 and 0 < parameters['isotropic_b'] < parameters['separation_d']
    ks = (0, 16, 50000)
    modes = [next(row for row in angular['actual_collar'] if row['k'] == k) for k in ks]
    dressed = [next(row for row in bound['actual_collar'] if row['k'] == k) for k in ks]
    previous = np.array([row['previous_potential_residual'] for row in modes])
    exact = np.array([row['radial_residual_norm'] for row in modes])
    upper = np.array([row['flux_normalized_residual_bound'] for row in dressed])
    assert np.all(np.isfinite([previous, exact, upper])) and np.all(0 < upper)
    assert np.all(upper < exact) and np.all(exact < previous)
    for mode, row in zip(modes, dressed):
        assert row['q'] == parameters['q'] and row['b'] == parameters['isotropic_b']
        assert np.isclose(row['saved_undressed_radial_residual'], mode['radial_residual_norm'], rtol=1e-14, atol=0)
        zero = next(item for item in row['energy_comparisons'] if item['label'] == 'zero')
        assert zero['energy'] == 0 and np.isclose(zero['residual_bound'], row['flux_normalized_residual_bound'], rtol=1e-14, atol=0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13, 7.2))
    fig.suptitle('Correcting the local collar mode reduces its Dirac residual',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Flat point-pair model  |  q = 100,000  |  A local mathematical test',
             ha='center', fontsize=11, color='#475569')
    left = fig.add_axes([.045, .32, .40, .46]); left.set_axis_off()
    left.set_title('A. Local collar only — schematic, not to scale', fontsize=11, pad=17)
    left.set(xlim=(0, 1), ylim=(0, 1)); left.set_aspect('equal', adjustable='box')
    left.add_patch(Circle((.22, .53), .145, facecolor='#e6f4f3', edgecolor='#007e87', linewidth=1.7))
    left.scatter([.22, .82], [.53, .53], s=45, c=['#007e87', '#475569'], zorder=5)
    left.annotate('', xy=(.325, .625), xytext=(.22, .53), arrowprops={'arrowstyle': '->', 'color': '#007e87', 'lw': 1.3})
    left.text(.31, .65, 'b', fontsize=11, color='#007e87')
    left.text(.22, .29, 'Selected focus', ha='center', fontsize=10)
    left.text(.82, .29, 'Partner focus', ha='center', fontsize=10)
    left.annotate('', xy=(.82, .39), xytext=(.22, .39), arrowprops={'arrowstyle': '<->', 'color': '#64748b', 'lw': 1.2})
    left.text(.52, .41, 'Separation d', ha='center', fontsize=9.8, color='#475569')
    left.add_patch(FancyArrowPatch((.23, .70), (.82, .60), connectionstyle='arc3,rad=-.32',
                                  arrowstyle='-', linestyle='--', linewidth=1.5, color='#a77612'))
    left.text(.56, .92, 'Global return / escape unsolved', ha='center', fontsize=10, color='#8a650f')
    left.text(.5, .16, 'Angular solution valid inside 0 < b < d', ha='center', fontsize=10.3, color='#007e87')
    left.text(.5, .065, 'Saved collar: b/d = 2.36 × 10⁻⁸', ha='center', fontsize=10, color='#475569')
    right = fig.add_axes([.575, .355, .375, .425])
    right.set_title('B. Local residuals and an upper bound', fontsize=11.5, pad=17)
    series = [(previous, '#64748b', 'o', 'Isolated mode residual'),
              (exact, '#007e87', 'D', 'Exact-angular residual'),
              (upper, '#7951b0', 'v', 'First-correction bound (E = 0)')]
    for shift, (values, color, marker, label) in zip((-.19, 0, .19), series):
        right.scatter(np.arange(3)+shift, values, s=65, marker=marker, color=color, label=label, zorder=4)
    right.set_yscale('log')
    right.set(xticks=range(3), xticklabels=['k = 0', 'k = 16', 'k = 50,000'],
              ylabel='Local residual norm (1/rₑ)', xlim=(-.55, 2.55),
              ylim=(float(upper.min())/4, float(previous.max())*4))
    right.spines[['top', 'right']].set_visible(False); right.grid(axis='y', which='major', alpha=.15)
    right.set_axisbelow(True); right.tick_params(axis='x', length=0, pad=9)
    fig.legend(loc='center', bbox_to_anchor=(.7625, .255), ncol=1, frameon=False, fontsize=9.6,
               handletextpad=.6, labelspacing=.55)
    fig.text(.05, .176, 'E = 0 isolates the derivative correction; a finite-energy input adds an energy-dependent term to the bound.',
             fontsize=10, color='#475569')
    fig.text(.05, .119, 'Small local residuals do not establish transport accuracy, global return, escape or a stable medium.',
             fontsize=10, color='#475569')
    fig.text(.05, .052, 'Model calculation, not a new physics discovery. 13 Sep 2027 → 12 Sep 2027 remains unproved.',
             fontsize=10.3, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'dipole_collar_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def exterior_boundary_progress():
    """Plot the saved conditional two-mouth exterior field and numerical ledger."""
    from matplotlib.colors import LogNorm
    from matplotlib.tri import Triangulation
    from matplotlib.patches import Circle
    root = Path(__file__).resolve().parent
    record = json.loads((root/'exterior_boundary_map_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and all(record['tests'].values())
    assert record['parameters'] == dict(q=1, mu=.5, energy=.4, d=1., b=.2)
    case = record['representative']; field = case['field']
    vertices = np.asarray(field['vertices']); triangles = np.asarray(field['triangles'], dtype=int)
    density = np.asarray(field['rho_u_squared'])
    assert vertices.ndim == 2 and vertices.shape[1] == 2 and triangles.shape == (len(density), 3)
    assert triangles.min() >= 0 and triangles.max() < len(vertices)
    assert np.all(np.isfinite(density)) and np.all(density >= 0) and density.max() > 0
    assert case['energy'] == .4 and case['R'] == 4 and case['h'] == .035
    omitted = sum(case['omitted_outgoing_norm'].values())
    flux = np.array([case['lower_zero_flux'], case['retained_free_outer_flux'], omitted, case['ledger']['numerical_loss']])
    assert np.all(flux > 0) and abs(case['mouth_nonzero_net_flux']) < 1e-10
    assert abs(flux.sum()+case['mouth_nonzero_net_flux']-1) < 1e-8
    assert abs(case['reduced_physical_plus_numerical_balance_residual']) < 1e-8
    amplitude = complex(case['lower_zero_amplitude']['real'], case['lower_zero_amplitude']['imag'])
    assert np.isclose(abs(amplitude)**2, flux[0], rtol=1e-13, atol=0)
    mesh = sorted(record['mesh_refinement'], key=lambda row: row['h'], reverse=True)
    assert len(mesh) >= 2 and len({row['R'] for row in mesh}) == 1 and len({row['nmax'] for row in mesh}) == 1
    h = np.array([row['h'] for row in mesh])
    numerical = np.array([row['ledger']['numerical_loss'] for row in mesh])
    modal = np.array([sum(row['omitted_outgoing_norm'].values()) for row in mesh])
    assert np.all(np.diff(h) < 0) and np.all(numerical > 0) and np.all(modal > 0)
    radius = sorted(record['radius_refinement'], key=lambda row: row['R'])
    assert [row['R'] for row in radius] == [2., 4., 8.]
    fixed_radius_mesh = len({row['h'] for row in radius}) == len({row['nmax'] for row in radius}) == 1
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13.5, 8.1))
    fig.suptitle('A computed exterior wave, with explicit model boundaries', fontsize=17, fontweight='bold', y=.962)
    fig.text(.5, .906, 'q = 1  |  μ = 1/2  |  E d = 0.4  |  Unit incident current  |  No physical throat matching',
             ha='center', fontsize=10.5, color='#475569')
    image_ax = fig.add_axes([.067, .30, .32, .525])
    image_ax.set_title('A. Computed wave intensity', fontsize=11.5, pad=13)
    triangulation = Triangulation(vertices[:, 0], vertices[:, 1], triangles)
    positive = density[density > 0]; low = max(float(positive.min()), float(density.max())*1e-4)
    mesh_image = image_ax.tripcolor(triangulation, facecolors=density, shading='flat', cmap='viridis',
                                   norm=LogNorm(vmin=low, vmax=float(density.max())), rasterized=True)
    for z in (.5, -.5):
        image_ax.add_patch(Circle((0, z), .2, facecolor='white', edgecolor='#475569', lw=1))
    image_ax.annotate('Input lead', xy=(.18, .57), xytext=(.53, .93), color='white', fontsize=10,
                      arrowprops={'arrowstyle': '->', 'color': 'white', 'lw': 1.2})
    image_ax.annotate('Other-mouth lead', xy=(.18, -.55), xytext=(.50, -.92), color='white', fontsize=9.5,
                      arrowprops={'arrowstyle': '->', 'color': 'white', 'lw': 1.2})
    image_ax.set(xlabel='Cylindrical radius ρ / d', ylabel='Axial position z / d', xlim=(0, 1.4), ylim=(-1.2, 1.2))
    image_ax.set_aspect('equal', adjustable='box')
    color_ax = fig.add_axes([.40, .34, .014, .44]); colorbar = fig.colorbar(mesh_image, cax=color_ax)
    colorbar.set_label(r'Weighted intensity $\rho |u|^2$ (log scale)', fontsize=9.7)
    fig.text(.242, .222, f'Near-mouth crop; whole-domain radius R = {case["R"]:g} d.', ha='center', fontsize=10, color='#475569')
    ledger_ax = fig.add_axes([.675, .595, .285, .23])
    ledger_ax.set_title('B. Fractions of unit model input', fontsize=11.5, pad=13)
    colors = ('#007e87', '#326bb3', '#7951b0', '#b38516')
    bars = ledger_ax.barh(np.arange(4), flux, color=colors, height=.57)
    for item in bars[2:]: item.set_hatch('//')
    for y, value in enumerate(flux):
        ledger_ax.text(value+.018, y, f'{100*value:.2f}%', va='center', fontsize=9.8, color='#334155')
    ledger_ax.set(yticks=range(4), yticklabels=['Other-mouth mode', 'Free outer lead', 'Omitted mode norm', 'Numerical loss (DG)'],
                  xlim=(0, max(.9, float(flux.max())+.20)), xticks=[0, .25, .5, .75])
    ledger_ax.invert_yaxis(); ledger_ax.spines[['top', 'right']].set_visible(False)
    ledger_ax.grid(axis='x', alpha=.13); ledger_ax.set_axisbelow(True); ledger_ax.tick_params(axis='y', length=0)
    convergence = fig.add_axes([.61, .32, .35, .17])
    convergence.set_title(f'C. Mesh refinement at R = {mesh[0]["R"]:g} d, mode cutoff = {mesh[0]["nmax"]}',
                          fontsize=10.7, pad=12)
    convergence.semilogy(h, numerical, 'o-', color='#b38516', label='Numerical loss', ms=4)
    convergence.semilogy(h, modal, 'v-', color='#7951b0', label='Omitted mode norm', ms=4)
    convergence.set(xlabel='Mesh scale h / d (finer →)', ylabel='Input fraction', xlim=(max(h)*1.05, min(h)*.8), xticks=h)
    convergence.spines[['top', 'right']].set_visible(False); convergence.grid(alpha=.13)
    convergence.legend(loc='best', frameon=False, fontsize=8.8)
    sensitivity = ' → '.join(f'{100*row["lower_zero_flux"]:.2f}%' for row in radius)
    resolution = f'fixed h = {radius[0]["h"]:g} d' if fixed_radius_mesh else 'reported meshes differ'
    fig.text(.05, .188, f'Outer-radius dependence: R/d = 2 → 4 → 8 gives other-mouth model flux {sensitivity} ({resolution}).',
             fontsize=9.8, color='#8a650f')
    fig.text(.05, .139, 'Specified product mouth leads and free outer map. Field intensity is not a particle path; exterior convergence remains unresolved.',
             fontsize=9.6, color='#475569')
    fig.text(.05, .09, 'These conditional fluxes do not establish a physical return/escape probability or a stable medium.',
             fontsize=10, color='#475569')
    fig.text(.05, .039, '13 Sep 2027 → 12 Sep 2027 remains unproved. This is a model calculation, not a new physical law.',
             fontsize=10.2, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'exterior_boundary_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def magnetic_tail_progress():
    """Plot saved finite magnetic-tail sensitivity with explicit boundary scope."""
    import hashlib
    from matplotlib.patches import Circle
    root = Path(__file__).resolve().parent
    record = json.loads((root/'exterior_tail_transport_results.json').read_text(encoding='utf-8'))
    free_path = root/'exterior_boundary_map_results.json'
    free = json.loads(free_path.read_text(encoding='utf-8'))
    assert record['passed'] and all(record['tests'].values()) and free['passed']
    assert record['input_source_sha256'][free_path.name] == hashlib.sha256(free_path.read_bytes()).hexdigest()
    assert record['parameters'] == dict(q=1, mu=.5, b=.2, d=1., energy=.4)
    assert not record['convergence']['infinity_convergence_certified']
    tail_cases = sorted(record['radius_refinement'], key=lambda row: row['R'])
    free_cases = sorted(free['radius_refinement'], key=lambda row: row['R'])
    radii = np.array([row['R'] for row in tail_cases])
    assert np.array_equal(radii, [2., 4., 8.]) and [r['R'] for r in free_cases] == list(radii)
    for row in tail_cases+free_cases:
        assert row['h'] == .05 and row['nmax'] == 12 and row['energy'] == .4
        amplitude = complex(row['lower_zero_amplitude']['real'], row['lower_zero_amplitude']['imag'])
        assert np.isclose(abs(amplitude)**2, row['lower_zero_flux'], rtol=1e-13, atol=0)
    for row in tail_cases:
        assert row['L'] == 160 and abs(row['reduced_balance_residual']) < 1e-8
        ledger = row['lower_zero_flux']+row['retained_tail_net_flux']+row['mouth_nonzero_net_flux']
        ledger += sum(row['omitted_outgoing_norm'].values())+row['ledger']['numerical_loss']
        assert abs(ledger-1) < 1e-8
        assert np.isclose(row['retained_tail_net_flux'], row['spherical_tail_net_flux']+
                          row['signed_tail_geometry_flux_difference'], rtol=1e-12, atol=0)
    old_flux = np.array([row['lower_zero_flux'] for row in free_cases])
    tail_flux = np.array([row['lower_zero_flux'] for row in tail_cases])
    old_changes, tail_changes = [], []
    for cases, changes, key in ((free_cases, old_changes, 'previous_free_matching_radius_changes'),
                                (tail_cases, tail_changes, 'complex_matching_radius_changes')):
        amplitudes = np.array([complex(r['lower_zero_amplitude']['real'], r['lower_zero_amplitude']['imag']) for r in cases])
        changes.extend(abs(np.diff(amplitudes)))
        assert np.allclose(changes, record['convergence'][key], rtol=1e-10, atol=0)
    assert np.all(np.array(tail_changes) < np.array(old_changes))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13.5, 7.4))
    fig.suptitle('Improved exterior treatment reduces matching-radius sensitivity',
                 fontsize=17, fontweight='bold', y=.955)
    fig.text(.5, .885, 'Conditional q = 1 wave model  |  Prescribed product mouth leads remain unchanged',
             ha='center', fontsize=11, color='#475569')
    left = fig.add_axes([.045, .365, .43, .425]); left.set_axis_off()
    left.set_title('A. Where each exterior approximation is used', fontsize=11.5, pad=16)
    left.set(xlim=(-.84, 1.70), ylim=(-.83, .83)); left.set_aspect('equal', adjustable='box')
    left.add_patch(Circle((0, 0), .73, facecolor='#fff5df', edgecolor='none'))
    left.add_patch(Circle((0, 0), .58, facecolor='#e6f4f3', edgecolor='#007e87', lw=1.5))
    left.add_patch(Circle((0, 0), .31, facecolor='#e8eef8', edgecolor='#326bb3', lw=1.5))
    for z in (.115, -.115):
        left.add_patch(Circle((-.085, z), .055, facecolor='white', edgecolor='#475569', lw=1.1))
    left.text(.10, .225, 'R', fontsize=11, color='#326bb3')
    left.text(-.43, .38, 'L', fontsize=11, color='#007e87')
    labels = [((.66, .29), (1.0, .54), 'Free approximation\nbeyond L', '#8a650f'),
              ((.46, .07), (1.0, .05), 'Magnetic continuation\nfrom R to L', '#007e87'),
              ((.24, -.16), (1.0, -.45), 'Numerical domain\ninside R', '#326bb3')]
    for point, label_point, words, color in labels:
        left.annotate(words, xy=point, xytext=label_point, fontsize=10, color=color, va='center',
                      arrowprops={'arrowstyle': '-', 'color': color, 'lw': 1.1})
    left.annotate('Two mouth cuts', xy=(-.085, -.115), xytext=(-.66, -.64), fontsize=9.5,
                  color='#475569', arrowprops={'arrowstyle': '-', 'color': '#64748b', 'lw': 1})
    fig.text(.26, .31, 'Schematic, not to scale. The continuation ends at L = 160 d.',
             ha='center', fontsize=9.8, color='#475569')
    right = fig.add_axes([.61, .385, .35, .405])
    right.set_title('B. Conditional flux vs matching radius', fontsize=11.5, pad=16)
    right.plot(radii, 100*old_flux, 'o-', color='#64748b', lw=1.8, ms=6, label='Previous free-map treatment')
    right.plot(radii, 100*tail_flux, 'D-', color='#007e87', lw=2, ms=6, label='Tail + coordinate matching')
    right.set(xlabel='Matching radius R / d', ylabel='Other-mouth mode flux / unit input (%)',
              xticks=radii, xlim=(1.6, 8.4), ylim=(58.4, 65.8))
    right.spines[['top', 'right']].set_visible(False); right.grid(alpha=.15); right.set_axisbelow(True)
    right.legend(loc='upper right', frameon=False, fontsize=9.5)
    fig.text(.785, .31, 'Fixed mesh h = 0.05 d; angular mode cutoff = 12.', ha='center', fontsize=9.8, color='#475569')
    fig.text(.05, .24, 'Complex-amplitude change when R is enlarged (previous treatment → tail + coordinate matching):', fontsize=10, color='#475569')
    fig.text(.05, .187, f'R/d: 2 → 4     {old_changes[0]:.3g} → {tail_changes[0]:.3g}       |       '
             f'4 → 8     {old_changes[1]:.3g} → {tail_changes[1]:.3g}', fontsize=10.7, color='#007e87', fontweight='bold')
    fig.text(.05, .124, 'Finite differences measure sensitivity, not error bounds. Infinity convergence and physical mouth matching remain unresolved.',
             fontsize=9.7, color='#475569')
    fig.text(.05, .06, 'This conditional field calculation establishes no physical medium or 13 Sep 2027 → 12 Sep 2027 passage.',
             fontsize=10.2, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'magnetic_tail_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


def imaginary_response_progress():
    """Plot same-resolution imaginary response and a declared loop determinant."""
    from matplotlib.patches import FancyBboxPatch
    root = Path(__file__).resolve().parent
    record = json.loads((root/'exterior_loop_closure_results.json').read_text(encoding='utf-8'))
    assert record['passed'] and all(record['tests'].values())
    cases = sorted(record['cases'], key=lambda row: row['kappa'])
    kappa = np.array([row['kappa'] for row in cases]); assert np.array_equal(kappa, [.1, .4, 1.])
    to_complex = lambda value: complex(value['real'], value['imag'])
    response = np.array([to_complex(row['open_response']) for row in cases])
    assert np.all(np.isfinite(response)) and np.max(abs(response.imag)) < 1e-11
    logs = {0.: [], np.pi: []}
    for case, transmission in zip(cases, response):
        assert (case['h'], case['R'], case['nmax'], case['L']) == (.08, 2., 8, 40.)
        for row in case['closures']:
            phase, length = row['endpoint_phase'], row['optical_throat_length']
            factor = np.exp(1j*phase-case['kappa']*length)
            determinant = to_complex(row['determinant'])
            assert np.isclose(factor, to_complex(row['throat_factor']), rtol=1e-13, atol=0)
            assert np.isclose(determinant, 1-factor*transmission, rtol=1e-13, atol=0)
            assert np.isclose(np.log(abs(determinant)), row['log_absolute_determinant'], rtol=1e-12, atol=0)
            assert abs(to_complex(row['direct_boundary_determinant_ratio'])-determinant) < 1e-11
            assert row['full_field_identity_relative_error'] < 1e-11 and row['sewn_equation_relative_residual'] < 1e-11
            if length == 1.:
                assert phase in logs
                logs[phase].append(row['log_absolute_determinant'])
    assert all(len(values) == len(kappa) for values in logs.values())
    assert np.all(np.array(logs[0.]) < 0) and np.all(np.array(logs[np.pi]) > 0)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(13.5, 8.0))
    fig.suptitle('From an exterior field response to a conditional loop determinant',
                 fontsize=17, fontweight='bold', y=.96)
    fig.text(.5, .905, 'Imaginary energy is a calculation variable — these values are not probabilities or particle trajectories',
             ha='center', fontsize=10.5, color='#475569')
    diagram = fig.add_axes([.06, .69, .90, .16]); diagram.set_axis_off()
    boxes = [(0., .26, '#e6f4f3', '#007e87', 'Computed exterior', r'$T(i\kappa)$', 'Same numerical resolution below'),
             (.37, .27, '#fff5df', '#a77612', 'Assumed throat connection', r'$z=e^{i\theta-\kappa\ell}$', 'Length ℓ and phase θ are inputs'),
             (.75, .25, '#f0eaf8', '#7951b0', 'Conditional determinant', r'$D=1-zT$', 'Finite boundary-matrix ratio')]
    for x, width, fill, color, title, formula, detail in boxes:
        diagram.add_patch(FancyBboxPatch((x+.008, .09), width-.016, .82,
                          boxstyle='round,pad=.008,rounding_size=.025',
                          facecolor=fill, edgecolor=color, lw=1.1))
        diagram.text(x+width/2, .76, title, ha='center', va='center', fontsize=10.5, color=color, fontweight='bold')
        diagram.text(x+width/2, .48, formula, ha='center', va='center', fontsize=16, color=color)
        diagram.text(x+width/2, .22, detail, ha='center', va='center', fontsize=9.4, color='#475569')
    for start, finish in ((.275, .355), (.655, .735)):
        diagram.annotate('', xy=(finish, .50), xytext=(start, .50),
                         arrowprops={'arrowstyle': '->', 'color': '#64748b', 'lw': 1.6})
    left = fig.add_axes([.085, .32, .345, .28])
    left.set_title('A. Sampled exterior response', fontsize=11.5, pad=13)
    left.plot(kappa, response.real, 'o--', color='#007e87', lw=1.5, ms=6)
    for x, value in zip(kappa, response.real):
        left.annotate(f'{value:.3f}', xy=(x, value), xytext=(0, 10), textcoords='offset points',
                      ha='center', color='#007e87', fontsize=9.7)
    left.set(xlabel=r'Imaginary-energy parameter $\kappa d$', ylabel=r'Real response $T(i\kappa)$',
             xlim=(.02, 1.08), ylim=(0, 1), xticks=kappa)
    left.spines[['top', 'right']].set_visible(False); left.grid(alpha=.15); left.set_axisbelow(True)
    right = fig.add_axes([.615, .32, .345, .28])
    right.set_title('B. Same response, two assumed phases', fontsize=11.5, pad=13)
    right.axhline(0, color='#94a3b8', lw=1)
    for phase, color, marker, label in ((0., '#7951b0', 'o', r'$\theta=0$'), (np.pi, '#b38516', 'D', r'$\theta=\pi$')):
        right.plot(kappa, logs[phase], marker+'--', color=color, lw=1.5, ms=6, label=label)
    right.set(xlabel=r'Imaginary-energy parameter $\kappa d$', ylabel=r'$\log |D|$',
              xlim=(.02, 1.08), ylim=(-1.40, .73), xticks=kappa)
    right.spines[['top', 'right']].set_visible(False); right.grid(alpha=.15); right.set_axisbelow(True)
    right.legend(loc='lower right', frameon=False, fontsize=10)
    fig.text(.2575, .237, 'Three computed samples; dashed segments guide the eye.', ha='center', fontsize=9.7, color='#475569')
    fig.text(.7875, .237, 'Chosen optical throat length: ℓ / d = 1.', ha='center', fontsize=10, color='#475569')
    fig.text(.05, .17, 'Both panels: q = 1, h/d = 0.08, R/d = 2, L/d = 40, mode cutoff 8; fixed product mouth leads.',
             fontsize=10, color='#475569')
    fig.text(.05, .113, 'The endpoint phase is not physically derived. A full vacuum stress tensor and a stable medium remain unproved.',
             fontsize=10, color='#475569')
    fig.text(.05, .05, 'No physical time-travel route or 13 Sep 2027 → 12 Sep 2027 passage demonstrated.',
             fontsize=10.5, fontweight='bold', color='#007e87')
    outputs = tuple(root/f'imaginary_response_progress.{suffix}' for suffix in ('png', 'pdf'))
    for path in outputs: fig.savefig(path, dpi=180, facecolor='white')
    plt.close(fig); print(outputs[0]); return outputs


if __name__ == '__main__':
    main()
    source_progress()
    response_progress()
    probe_future_progress()
    proton_ring_progress()
    focusing_feedback_progress()
    return_loss_progress()
    offset_rf_progress()
    balance_kicker_progress()
    driven_outgoing_progress()
    coupled_coil_progress()
