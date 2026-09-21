# Directory guide

Current increment159: The [localized-collar comparison](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#localized-collar-comparison-and-a-physical-decision-criterion) reduces squared auxiliary-solver error by **6.10%**, from 1,187.04 to 1,114.61; all44 provider, calculation and independent-audit checks pass. Almost all improvement available in this finite pool is already captured. The physical supporting stress remains unknown. The next target is a **coarse absolute bound on the full renormalized source**, which could reject this one-field candidate without resolving a tiny signal to high relative accuracy. Further basis growth must be justified by that physical error budget. See [results](results/dihole_collar_selection_results.json), [audit](results/dihole_collar_audit_results.json), and [visual](figures/dihole_collar_selection_progress.png).

Current increment158: Selecting [fields by the error they remove](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#residual-localization-and-optimal-candidate-combinations) reduced the squared auxiliary-solver error by **69.39%**, from 3,878.29 to 1,187.05, using 64 added combinations. The tested radial additions produced almost all of the improvement; extra angular modes 6/7 contributed very little in this test. All 45 checks across the new components pass. A [direct spatial check](results/dihole_selected_profile_results.json) places 99.64% of the remaining bridge error in its collars, identifying localized radial profiles as the next target. The supporting source remains unproved. See [results](results/dihole_residual_selection_results.json) and [visual](figures/dihole_residual_selection_progress.png).

Current increment157: The [opening-field test](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#fields-crossing-the-openings-and-regional-residuals) adds twelve fields that extend continuously through the bridge openings. The squared auxiliary-solver residual falls by **10.336%**, with34 checks passing across the new components. A failed integration check was preserved and corrected by finer exterior integration. The source interval still includes both signs: no supporting medium or new physics law is established. See [results](results/dihole_mouth_enrichment_results.json) and [visual](figures/dihole_mouth_enrichment_progress.png).

Current increment156: [proper heat and common-filter tests](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#proper-heat-response-and-a-common-filter-residual-test), [heat calculation](scripts/dihole_proper_heat.py), [results](results/dihole_proper_heat_results.json), [common-filter result](results/dihole_common_pv_residual_results.json), and [visual](figures/dihole_proper_heat_progress.png). Eighteen checks pass. The combined344-field finite response changes by0.00847%, but short-scale modes and a large first-order residual still prevent a physical supporting-source value.

Current increment155: [retained-source error test](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#retained-source-error-test-on-the-saved-sewn-forms), [calculation](scripts/dihole_retained_resolvent.py), [results](results/dihole_retained_resolvent_results.json), and [visual](figures/dihole_retained_resolvent_progress.png). Seven checks pass, but the separated-shift error estimate is far too wide to determine the retained source sign. This rules out a numerical shortcut, not the physical candidate; the full supporting source remains unknown.

Current increment154: [frequency and regulator limits](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#frequency-and-regulator-limits-of-the-compact-response), [calculation](scripts/dirac_regulator_limits.py), [results](results/dirac_regulator_limits_results.json), and [visual](figures/dirac_regulator_limits_progress.png). The ordered compact first-variation limit exists in136's declared scheme; thirteen controls pass. A numerical source value with quantitative continuum and heat-correction errors remains the next target.

Current increment153: [actual zero-state confinement proof](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#uniform-zero-state-confinement-on-the-sewn-geometry), [calculation](scripts/dihole_zero_tail.py), [results](results/dihole_zero_tail_results.json), and [visual](figures/dihole_zero_tail_progress.png). The weighted-zero condition is proved for the actual sewn geometry, uniformly in angular sector. Eight controls pass, including interval-certified analytic constants. The full high-energy/regulator limit and physical source remain open.

Current increment152: [joint bounded-energy limit and zero-state test](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#uniform-low-energy-limit-and-zero-state-test), [controls](scripts/dirac_infrared_corner.py), [results](results/dirac_infrared_corner_results.json), and [visual](figures/dirac_infrared_corner_progress.png). Conditional local trace/derivative convergence is proved, and an abstract counterexample excludes an argument from symmetry alone. Eleven controls pass. Actual weighted-zero density and the full physical source remain open.

Current increment151: [localized Green identity and proper-mass limit](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md#localized-green-identity-and-proper-mass-limit), [algebra controls](scripts/dirac_green_mass_limit.py), [results](results/dirac_green_mass_limit_results.json), and [visual](figures/dirac_green_mass_limit_progress.png). Positive-mass and fixed-nonzero-frequency correspondence is proved with the actual unbounded-lapse assumptions. Ten finite algebra controls pass; the joint infrared/stress limit and physical supporting source remain unresolved.

Current increment150: [causal restriction and coupled principal test](docs/TILTED_CLOCK_ADMISSIBILITY.md), [calculator](scripts/tilted_clock_admissibility.py), [results](results/tilted_clock_admissibility_results.json), and [visual](figures/tilted_clock_admissibility_progress.png). Eleven saved-state checks pass. The same globally regular real timelike clock excludes closed metric-causal loops; local principal health does not establish global stability. Stop tuning this branch for backward travel. The existing quantum-source route's physical prescription remains unresolved.

Current increment149: [tilted-clock throat and horizon](docs/TILTED_CLOCK_THROAT.md), [local driver](scripts/tilted_clock_throat.py), [extension driver](scripts/tilted_clock_extension.py), [diagnostic driver](scripts/tilted_clock_diagnostics.py), and [visual](figures/tilted_clock_diagnostics_progress.png). Local and separate short horizon checks pass; the [long extension fails its accuracy target](results/tilted_clock_extension_results.json). The report links all saved states and source hashes. No stable global passage is established.

Current increment148: [throat source budget and static-branch rejection](docs/THROAT_SOURCE_BUDGET.md), [calculator](scripts/throat_source_budget.py), [results](results/throat_source_budget_results.json), and [visual](figures/throat_source_budget_progress.png). Eleven checks pass using saved states; a published no-go theorem and the specific density bound exclude the static aligned completion. A tilted-clock source is the next hypothesis, not an established medium.

Current increment147: [nonlinear evolution and null stress](docs/NONLINEAR_PACKET_EVOLUTION.md), [calculation](scripts/nonlinear_packet_evolution.py), [results](results/nonlinear_packet_evolution_results.json), [saved states](data/nonlinear_packet_evolution_state.npz), and [visual](figures/nonlinear_packet_evolution_progress.png). Eighteen checks pass through Ht=0.5; no stable bridge or backward passage is established. [Failed numerical pilots](results/nonlinear_packet_evolution_fd_trial.json) and their provenance sources are retained.

Current increment146: [nonlinear spherical initial constraints](docs/NONLINEAR_PACKET_CONSTRAINTS.md), [calculation](scripts/nonlinear_packet_constraints.py), [results](results/nonlinear_packet_constraints_results.json), [initial fields](data/nonlinear_packet_constraints_state.npz), and [visual](figures/nonlinear_packet_constraints_progress.png). Fourteen checks pass for one initial slice; nonlinear evolution, stability and time passage remain unproved.

Current increment145: [spherical packet and nonlinear-geometry limitation](docs/LOCALIZED_EXPANDING_PACKET.md), [packet calculation](scripts/localized_expanding_packet.py), [results](results/localized_expanding_packet_results.json), [saved fields](data/localized_expanding_packet_state.npz), [geometric check](scripts/localized_packet_geometry.py), [rejection result](results/localized_packet_geometry_results.json), and [visual](figures/localized_expanding_packet_progress.png). Linear negative stress is localized; a controlled nonlinear gravitational source is still missing.

Current increment144: [coupled expanding completion](docs/EXPANDING_FIELD_COMPLETION.md), [calculation](scripts/expanding_field_completion.py), [results](results/expanding_field_completion_results.json), and [visual](figures/expanding_field_completion_progress.png). The added energy source is included in the background; conditional mode amplification is bounded. A localized negative-null-energy source is still required.

Current increment143: [localized-field hypothesis and gravitational rejection](docs/LOCALIZED_FIELD_HYPOTHESIS.md), [calculation](scripts/localized_field_hypothesis.py), [results](results/localized_field_hypothesis_results.json), and [visual](figures/localized_field_hypothesis_progress.png). The completed test rejects the minimal static completion; a changed-background hypothesis is assessed separately.

Current increment142: [angular cross-check](scripts/dihole_angular_crosscheck.py), [results](results/dihole_angular_crosscheck_results.json), [saved matrices](data/dihole_angular_crosscheck_state.npz), and [visual](figures/dihole_angular_crosscheck_progress.png). Thirteen checks pass at radial resolution12.

Increment141: [exterior feedback calculation](scripts/dihole_exterior_feedback.py), [saved results](results/dihole_exterior_feedback_results.json), and [visual](figures/dihole_exterior_feedback_progress.png). The angular cross-check above is now complete.

Current increment140: [radial refinement report](docs/RADIAL_REFINEMENT.md), [driver](scripts/dihole_radial_refinement.py), [saved results](results/dihole_radial_refinement_results.json), and [visual progress](figures/dihole_radial_refinement_progress.png). Earlier evidence remains in place.

Files are grouped by purpose. Run commands from the project root and consult
the reports before starting potentially expensive calculations. Every chat
session must follow [AGENTS.md](AGENTS.md) when creating or relocating files.

| Folder | Contents |
| --- | --- |
| [scripts/](scripts/) | Python calculations, shared modules, and plot generators |
| [results/](results/) | Saved JSON results and numerical reports |
| [data/](data/) | NPZ states, initial coefficients, and resumable checkpoints |
| [figures/](figures/) | PNG previews and PDF exports |
| [docs/](docs/) | Research notes, derivations, and companion reports |
| [maintenance/](maintenance/) | Verification and cleanup/migration records |
| [provenance/](provenance/) | Preserved historical source artifacts |

Use `project_path` from `scripts/project_paths.py` to locate project assets in
Python. Direct commands use `python scripts/name.py`. For `python -c` imports,
set `$env:PYTHONPATH = (Join-Path (Get-Location) 'scripts')` in PowerShell.

The [layout migration](maintenance/layout-migration.json) records every old/new
path and hash. Saved JSON, NPZ, and figure bytes were preserved. Historical JSON
recipes still describe their original layout; translate their paths using this
mapping. Source-only path edits are recorded explicitly, and new runs record
actual new hashes. Existing Ori checkpoints retain strict source checks:
use their original source revision for exact historical resumption.

Run `python -B maintenance/verify_layout.py` for the lightweight layout check.

## Start here

| Material | Location |
| --- | --- |
| Current work and continuity | [TASK_STATE.md](TASK_STATE.md) |
| Research results and reproduction notes | [docs/RESEARCH.md](docs/RESEARCH.md) |
| Concise professor-assessment note | [docs/PROFESSOR_BRIEF.md](docs/PROFESSOR_BRIEF.md) |
| Independent bridge modes and finite-matrix frequency tails | [derivation](docs/INDEPENDENT_BRIDGE_FIELDS.md), [driver](scripts/dihole_bridge_enrichment.py), [bridge fields](scripts/dihole_enriched_bridge.py), [tail bounds](scripts/dihole_finite_frequency_tail.py), [results](results/dihole_bridge_enrichment_results.json), [matrix state](data/dihole_bridge_enrichment_state.npz), [figure](figures/dihole_bridge_enrichment_progress.png) |
| Globally sewn physical trial fields and basis sensitivity | [derivation](docs/GLOBAL_SEWN_FORM.md), [driver](scripts/dihole_sewn_form.py), [round quadrature](scripts/dihole_sewn_quadrature.py), [results](results/dihole_sewn_form_results.json), [figure](figures/dihole_sewn_form_progress.png) |
| Common physical Euclidean form and metric derivative | [operator](docs/PHYSICAL_EUCLIDEAN_OPERATOR.md), [implementation](scripts/dihole_euclidean_form.py), [results](results/dihole_euclidean_form_results.json), [figure](figures/dihole_euclidean_form_progress.png) |
| Complete finite angular multiplets and explicit action reference | [sectors](scripts/dihole_pv_sectors.py), [basis](scripts/dihole_azimuthal_bridge.py), [results](results/dihole_pv_sectors_results.json), [figure](figures/dihole_pv_sectors_progress.png), [scheme](docs/DIRAC_SCHEME_MATCHING.md), [counterterms](scripts/dihole_pv_counterterms.py) |
| Physical-mass Pauli--Villars connected pilot | [bridge](scripts/dihole_massive_bridge.py), [driver](scripts/dihole_pv_source.py), [results](results/dihole_pv_source_results.json), [figure](figures/dihole_pv_source_progress.png), [physical prescription and limits](docs/PHYSICAL_MASS_REGULATOR.md) |
| Direct bridge determinant and finite-regulator diagnosis | [driver](scripts/dihole_direct_source.py), [results](results/dihole_direct_source_results.json), [figure](figures/dihole_direct_source_progress.png), [normalization proof and limits](docs/RESEARCH.md) |
| Sewn ground state and local renormalization responses | [state/domain](docs/SEWN_SOURCE_STATE.md), [driver](scripts/dihole_local_source.py), [results](results/dihole_local_source_results.json), [figure](figures/dihole_local_source_progress.png), [limits](docs/RESEARCH.md) |
| Connected finite-window energy-and-pressure contribution | [driver](scripts/dihole_source_window.py), [results](results/dihole_source_window_results.json), [figure](figures/dihole_source_window_progress.png), [normalization and source comparison](docs/RESEARCH.md) |
| Actual curved exterior return sewn to the bridge | [curved helper](scripts/dihole_curved_exterior.py), [driver](scripts/dihole_return_pilot.py), [results](results/dihole_return_pilot_results.json), [figure](figures/dihole_return_pilot_progress.png), [derivation and limits](docs/RESEARCH.md) |
| Actual two-collar quantum response and compact source variation | [driver](scripts/dihole_bridge_scattering.py), [results](results/dihole_bridge_scattering_results.json), [figure](figures/dihole_bridge_scattering_progress.png), [derivation](docs/RESEARCH.md) |
| Actual central-cylinder zero profile and coupled propagation | [driver](scripts/dihole_cylinder_transfer.py), [results](results/dihole_cylinder_transfer_results.json), [figure](figures/dihole_cylinder_transfer_progress.png) |
| Smooth bridge with actual cut data and required null source | [driver](scripts/dihole_stretched_collar.py), [results](results/dihole_stretched_collar_results.json), [figure](figures/dihole_stretched_collar_progress.png), [derivation](docs/RESEARCH.md) |
| Phase-sensitive return-force bound and support sensitivity | [driver](scripts/quantum_return_support.py), [results](results/quantum_return_support_results.json), [figure](figures/quantum_return_support_progress.png), [assumptions](docs/RESEARCH.md) |
| Certified central Dirichlet spectral level | [driver](scripts/dihole_sine_certificate.py), [analytic remainder](scripts/dihole_sine_remainder.py), [results](results/dihole_sine_certificate_results.json), [saved matrices](data/dihole_sine_certificate_state.npz), [figure](figures/dihole_sine_certificate_progress.png), [proof and assumptions](docs/RESEARCH.md) |
| Positivity-compatible Schur spectrum and conditional matrix proof | [driver](scripts/dihole_schur_spectrum.py), [results](results/dihole_schur_spectrum_results.json), [saved trials](data/dihole_schur_spectrum_state.npz), [figure](figures/dihole_schur_spectrum_progress.png), [proof and premises](docs/RESEARCH.md) |
| Enclosed spatial-complement gap and coupling pilot | [driver](scripts/dihole_spatial_complement.py), [results](results/dihole_spatial_complement_results.json), [saved trials](data/dihole_spatial_complement_state.npz), [figure](figures/dihole_spatial_complement_progress.png), [proof](docs/RESEARCH.md) |
| Continuum-box angular tail and interval geometry | [driver](scripts/dihole_continuum_heat_tail.py), [results](results/dihole_continuum_heat_tail_results.json), [figure](figures/dihole_continuum_heat_tail_progress.png), [proof](docs/RESEARCH.md) |
| Angular completion of the fixed spatial heat model | [driver](scripts/dihole_heat_angular.py), [results](results/dihole_heat_angular_results.json), [saved matrices](data/dihole_heat_angular_state.npz), [figure](figures/dihole_heat_angular_progress.png), [proof](docs/RESEARCH.md) |
| Finite curved heat diagonal and angular sums | [driver](scripts/dihole_finite_heat.py), [results](results/dihole_finite_heat_results.json), [figure](figures/dihole_finite_heat_progress.png), [scope](docs/RESEARCH.md) |
| Optical heat subtraction and lapse counterterms | [driver](scripts/dihole_optical_heat.py), [results](results/dihole_optical_heat_results.json), [figure](figures/dihole_optical_heat_progress.png), [derivation](docs/OPTICAL_HEAT_SUBTRACTION.md) |
| Local lapse susceptibility and specified exterior state | [driver](scripts/dihole_lapse_susceptibility.py), [results](results/dihole_lapse_susceptibility_results.json), [figure](figures/dihole_lapse_susceptibility_progress.png), [state and stress bridge](docs/EXTERIOR_STRESS_BRIDGE.md) |
| Complete compact-probe certificate in the exact exterior | [driver](scripts/dihole_probe_certificate.py), [interval arithmetic](scripts/dihole_interval_arithmetic.py), [validated grid](scripts/dihole_interval_quadrature.py), [results](results/dihole_probe_certificate_results.json), [figure](figures/dihole_probe_certificate_progress.png) |
| Analytic compact-probe quadrature remainder only | [bound](scripts/dihole_quadrature_remainder.py), [interval results](results/dihole_quadrature_remainder_results.json); this historical remainder component is now combined with validated sums in the119 certificate |
| Natural quantum domain on the uncut exterior | [proof](docs/DIHOLE_DIRAC_DOMAIN.md), [consistency checks](scripts/dihole_domain_check.py), [results](results/dihole_domain_check_results.json), [figure](figures/dihole_domain_progress.png) |
| Compact response in the actual curved coefficients | [calculation](scripts/dihole_probe_residual.py), [results](results/dihole_probe_residual_results.json), [corrected fields](data/dihole_probe_residual_fields.npz), [figure](figures/dihole_probe_residual_progress.png) |
| Whole-plane saddle response with continuum bounds | [solver](scripts/null_saddle_resolvent.py), [independent interval checker](scripts/null_saddle_residual_interval.py), [results](results/null_saddle_resolvent_results.json), [fields](data/null_saddle_resolvent_fields.npz), [figure](figures/null_saddle_resolvent_progress.png) |
| Local frequency response and corner estimate | [scripts/null_saddle_spectrum.py](scripts/null_saddle_spectrum.py), [frequency arms](scripts/null_saddle_frequency_tail.py), [spectrum](results/null_saddle_spectrum_results.json), [corner proof](scripts/null_saddle_corner.py), [corner results](results/null_saddle_corner_results.json), [figure](figures/null_saddle_spectrum_progress.png) |
| Varying bilinear arms and matched local node | [scripts/null_saddle_tail.py](scripts/null_saddle_tail.py), [driver](scripts/null_saddle_matched.py), [results](results/null_saddle_matched_results.json), [figure](figures/null_saddle_matched_progress.png) |
| Local quantum scattering and computed field | [scripts/null_saddle_scattering.py](scripts/null_saddle_scattering.py), [analytic leads](scripts/null_saddle_leads.py), [results](results/null_saddle_scattering_results.json), [field](data/null_saddle_fields.npz), [figure](figures/null_saddle_scattering_progress.png) |
| Full meridional slab bound and local quantum mixing | [scripts/dihole_slab_bound.py](scripts/dihole_slab_bound.py), [results](results/dihole_slab_bound_results.json); [scripts/dihole_null_scaling.py](scripts/dihole_null_scaling.py), [results](results/dihole_null_scaling_results.json), [figure](figures/dihole_null_scaling_progress.png) |
| Return topology, central quantum barrier and conditional support threshold | [scripts/dihole_return_geometry.py](scripts/dihole_return_geometry.py), [results](results/dihole_return_geometry_results.json), [figure](figures/dihole_return_progress.png); exact assumptions in [research notes](docs/RESEARCH.md) |
| Quantum channels and curved axial response | [scripts/melvin_leads.py](scripts/melvin_leads.py), [results](results/melvin_leads_results.json), [figure](figures/melvin_channel_progress.png); [scripts/dihole_axial_tail.py](scripts/dihole_axial_tail.py) and [results](results/dihole_axial_tail_results.json) |
| Gravitating exterior and attachment data | [scripts/balanced_dihole.py](scripts/balanced_dihole.py), [results](results/balanced_dihole_results.json), [figure](figures/balanced_exterior_progress.png); [RN compatibility gate](scripts/rn_exterior_gate.py) and [results](results/rn_exterior_gate_results.json) |
| Smooth collar, magnetic response and same-geometry support target | [scripts/smooth_magnetic_collar.py](scripts/smooth_magnetic_collar.py), [results](results/smooth_magnetic_collar_results.json), [figure](figures/smooth_collar_progress.png); conditional RN exterior gate and coupled Dirac operator in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Finite magnetic matching on the actual mouth data | [scripts/magnetic_mouth_matching.py](scripts/magnetic_mouth_matching.py), [results](results/magnetic_mouth_matching_results.json), [figure](figures/magnetic_matching_progress.png); source definition and the remaining smooth-collar requirements in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Actual conforming field and failed residual target | [scripts/conforming_exterior_trial.py](scripts/conforming_exterior_trial.py), [saved field/results](results/conforming_exterior_trial_results.json), [figure](figures/conforming_trial_progress.png); [round-domain quadrature](scripts/physical_exterior_quadrature.py) and [controls](results/physical_exterior_quadrature_results.json) |
| Two-leg angular feedback and continuum-residual target | [scripts/feedback_residual_bound.py](scripts/feedback_residual_bound.py), [results](results/feedback_residual_bound_results.json), [diagram](figures/tail_completion_progress.png); finite mixing witness plus requirements for conforming physical trial fields |
| Complete retained low-frequency enclosure | [scripts/complete_threshold_bound.py](scripts/complete_threshold_bound.py), [results](results/complete_threshold_bound_results.json); all9 retained sectors, conditional finite N4/J9/DG result |
| Magnetic-tail angular quadrature error | [scripts/tail_quadrature_bound.py](scripts/tail_quadrature_bound.py), [results](results/tail_quadrature_bound_results.json); exact-GL truncation versus exact finite-J projection, separate from roundoff |
| All-angular high-frequency source bound | [scripts/ultraviolet_source_bound.py](scripts/ultraviolet_source_bound.py), [results](results/ultraviolet_source_bound_results.json), [frequency ledger](figures/source_remainder_progress.png) |
| One-sector zero-frequency enclosure | [scripts/threshold_return_bound.py](scripts/threshold_return_bound.py), [results](results/threshold_return_bound_results.json); finite N4/J9/DG model and non-interval floating evaluation |
| Buffered angular-feedback estimate | [scripts/buffered_angular_bound.py](scripts/buffered_angular_bound.py), [results](results/buffered_angular_bound_results.json); explicit decreasing bound, numerically too loose |
| Long-throat finite-window integral | [scripts/long_throat_integral.py](scripts/long_throat_integral.py), [results](results/long_throat_integral_results.json), [figure](figures/long_throat_integral_progress.png) |
| Infinite-angular direct-trace bound | [scripts/radial_angular_envelope.py](scripts/radial_angular_envelope.py), [results](results/radial_angular_envelope_results.json); conditional collar-feedback and threshold arguments in [docs/RESEARCH.md](docs/RESEARCH.md) |
| First massive mode on the saved smooth collar | [scripts/smooth_collar_mode.py](scripts/smooth_collar_mode.py), [results](results/smooth_collar_mode_results.json); reuses [scripts/dirac_return_closure.py](scripts/dirac_return_closure.py) through an optional profile argument |
| Longer-throat sampled source | [scripts/long_throat_source.py](scripts/long_throat_source.py), [results](results/long_throat_source_results.json), [figure](figures/long_throat_source_progress.png) |
| Defect-weighted angular remainder | [scripts/angular_defect_bound.py](scripts/angular_defect_bound.py), [results](results/angular_defect_bound_results.json); proof and assumptions in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Complete finite angular-source calculation | [scripts/complete_angular_source.py](scripts/complete_angular_source.py), [pilot and integration evidence](results/complete_angular_source_results.json), [figure](figures/complete_angular_source_progress.png); continuum-state and angular-cutoff arguments in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Three-sector reference sum and complete angular bases | [scripts/angular_reference_sum.py](scripts/angular_reference_sum.py), [saved results](results/angular_reference_sum_results.json), [scripts/complete_angular_modes.py](scripts/complete_angular_modes.py), [figure](figures/angular_reference_sum_progress.png) |
| Invariant AdS reference and interior source | [scripts/ads_reference_source.py](scripts/ads_reference_source.py), [saved results](results/ads_reference_source_results.json), [figure](figures/ads_reference_source_progress.png) |
| Geometric angular connection and missing sectors | [scripts/azimuthal_modes.py](scripts/azimuthal_modes.py), [connected source calculation](scripts/azimuthal_source.py), [saved results](results/azimuthal_source_results.json), [figure](figures/azimuthal_source_progress.png) |
| Opposite-chirality boundary matching and combined response | [scripts/dirac_sewing_audit.py](scripts/dirac_sewing_audit.py), [saved physical and adjoint controls](results/dirac_sewing_audit_results.json), [figure](figures/dirac_sewing_progress.png) |
| Throat pressure and energy-flow diagnostics | [scripts/throat_null_response.py](scripts/throat_null_response.py), [saved metric derivatives and checks](results/throat_null_response_results.json), [figure](figures/throat_null_progress.png) |
| Compact local source diagnostic | [scripts/local_lapse_response.py](scripts/local_lapse_response.py), [saved responses and direct perturbation checks](results/local_lapse_response_results.json), [figure](figures/local_lapse_progress.png) |
| Computed phase-dependent spectral diagnostic | [scripts/phase_vacuum_contrast.py](scripts/phase_vacuum_contrast.py), [saved determinants and integrals](results/phase_vacuum_contrast_results.json), [figure](figures/phase_vacuum_progress.png) |
| Curved modes and junction support | [scripts/curved_throat_modes.py](scripts/curved_throat_modes.py), [scripts/curved_exterior_closure.py](scripts/curved_exterior_closure.py), [closure results](results/curved_exterior_closure_results.json); [scripts/throat_junction_audit.py](scripts/throat_junction_audit.py), [required-stress results](results/throat_junction_audit_results.json) |
| Imaginary response and throat connection | [scripts/imaginary_weyl_leads.py](scripts/imaginary_weyl_leads.py), [scripts/imaginary_exterior.py](scripts/imaginary_exterior.py), [connected results](results/imaginary_exterior_results.json); [scripts/exterior_loop_closure.py](scripts/exterior_loop_closure.py), [closure results](results/exterior_loop_closure_results.json), [curved-mode derivation](docs/CURVED_THROAT_MATCHING.md) |
| Magnetic exterior continuation | [scripts/magnetic_tail_map.py](scripts/magnetic_tail_map.py), [tail checks](results/magnetic_tail_map_results.json); [scripts/exterior_tail_transport.py](scripts/exterior_tail_transport.py), [connected results](results/exterior_tail_transport_results.json), [figure](figures/magnetic_tail_progress.png) |
| Connected two-mouth exterior | [scripts/exterior_boundary_map.py](scripts/exterior_boundary_map.py), [saved results](results/exterior_boundary_map_results.json); [lead modes](scripts/spherical_weyl_leads.py), [mesh](scripts/two_mouth_mesh.py) |
| Incoming spinor modes | [docs/INCOMING_MODES.md](docs/INCOMING_MODES.md), [scripts/incoming_modes.py](scripts/incoming_modes.py), [saved results](results/incoming_modes_results.json) |
| Exact angular collar modes | [scripts/dipole_collar_modes.py](scripts/dipole_collar_modes.py), [saved results](results/dipole_collar_modes_results.json) |
| Local radial correction bound | [scripts/collar_dressing_bound.py](scripts/collar_dressing_bound.py), [saved results](results/collar_dressing_bound_results.json); derivation in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Open-boundary scattering | [docs/OPEN_BOUNDARY_SCATTERING.md](docs/OPEN_BOUNDARY_SCATTERING.md), [scripts/open_boundary_spinor.py](scripts/open_boundary_spinor.py), [saved results](results/open_boundary_spinor_results.json) |
| Quantum-source analysis | [docs/QUANTUM_COMPLETION.md](docs/QUANTUM_COMPLETION.md) |
| Source-energy analysis | [docs/SOURCE_FEASIBILITY.md](docs/SOURCE_FEASIBILITY.md) |
| Fiction-inspired diagnostics | [docs/FICTION_PHYSICS.md](docs/FICTION_PHYSICS.md) |

## Research families

| Family | Source and retained artifacts |
| --- | --- |
| Coupled models | `scripts/coupled_misner.py`, `scripts/coupled_ori.py`, corresponding `*_results.json` |
| Ori constraints and evolution | `scripts/disk_basis.py`, `scripts/ori_disk_constraints.py`, `scripts/ori_evolution.py`, `scripts/ori_large_analysis.py`; `ori_*`, `smooth_weak_*`, `polynomial_weak_*` data and comparisons |
| Casimir and quantum source | `scripts/casimir_cavity.py`, `scripts/ori_quantum_source.py`, corresponding `*_results.json` |
| Pulses and transport | `scripts/mmp_pulse.py`, `scripts/future_time.py`, corresponding `*_results.json` |
| Response and moving exterior | `scripts/polyakov_response.py`, `scripts/moving_exterior_maps.py`, `scripts/moving_exterior_nonlinear.py`; `polyakov_*`, `distributed_*`, `moving_exterior_*` results |
| Magnetic and Dirac calculations | `magnetic_*.py`, `dirac_*.py`, corresponding `*_results.json` |
| Spinful return | [scripts/spinful_return.py](scripts/spinful_return.py), [saved results](results/spinful_return_results.json), [figure](figures/spinful_return_progress.png) |
| Figures | Root `*_progress.png` / `*_progress.pdf` pairs and Ori evolution figures; generated by `scripts/quantum_horizon_plot.py`, `scripts/casimir_source_plot.py`, and `scripts/ori_large_analysis.py` |

The PNG and PDF versions serve different viewing/export purposes; both are
retained. JSON results, NPZ states, coefficient snapshots, failed-case records,
and checkpoint metadata are research evidence, not disposable temporary files.

## Existing directories

- [chronology_onset_audit/](chronology_onset_audit/): original audit sources,
  saved computations, report, run log, and SHA-256 manifest; kept together.
- [provenance/](provenance/): historical evolution source used for provenance,
  plus `ori_initial_data.cpython-312.pyc`, an unidentified historic
  **CPython 3.12-only** artifact without a matching source. It was archived
  from `__pycache__/ori_initial_data.cpython-312.pyc` with its SHA-256 unchanged;
  retain it for provenance, not as a supported executable entry point.
- [maintenance/](maintenance/): cleanup records. The
  [20 September 2026 manifest](maintenance/cleanup-2026-09-20.json) records exact
  deleted paths, byte counts, hashes, retained-file inventory, and verification.
- `.venv/`: local Python environment; retained untouched.
- `.cache/` and `__pycache__/`: runtime caches may reappear automatically.

## Earlier cleanup record

On 20 September 2026, cleanup removed 16 CPython bytecode files with retained
source counterparts and one generated Matplotlib font cache: **17 files,
2,284,220 bytes**. One orphan bytecode file was moved into `provenance/` as
described above. No scientific source, result, checkpoint, report, or figure
was moved or deleted during that earlier cache cleanup. The later layout
migration above relocates these files without deleting scientific artifacts.
No scientific calculations were rerun for cache cleanup.

Accepted coefficient snapshots and the audit run log were retained even where
content is redundant: their filenames participate in provenance or references.
Checkpoint files remain necessary to saved audits and resumption workflows.
The later migration updates and verifies those consumers together.

Initial verification matched all 369 retained files to their pre-cleanup
SHA-256 hashes. The subsequent archival check confirmed the bytecode's hash;
active `scripts/incoming_modes.py` changed concurrently and was left untouched by
cleanup. The final manifest records observed concurrent changes separately.
These seven files already differed from `chronology_onset_audit/sha256_manifest.json`:

- `chronology_onset_audit/symbolic_results.json`
- `chronology_onset_audit/mode_cases.json`
- `chronology_onset_audit/sample_trajectories.json`
- `chronology_onset_audit/kinetic_repair_cases.json`
- `chronology_onset_audit/initial_energy_results.json`
- `chronology_onset_audit/numerical_results.json`
- `chronology_onset_audit/kinetic_repair_results.json`

The cause was not investigated during cleanup. Both the files and historical
manifest were preserved without alteration; no integrity repair is claimed.
