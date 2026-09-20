# Project layout and collaboration

These instructions apply to every chat session and agent working in this project.
Read `TASK_STATE.md` for current research and `DIRECTORY_GUIDE.md` for navigation.

## Put new files in the right place

| Content | Directory |
| --- | --- |
| Active Python modules, calculations, and plot generators | `scripts/` |
| Saved JSON calculation results and reports | `results/` |
| NPZ states, initial coefficients, and resumable checkpoints | `data/` |
| Generated PNG and PDF figures | `figures/` |
| Research reports, derivations, and supporting Markdown | `docs/` |
| Layout verification and cleanup/migration records | `maintenance/` |
| Frozen historical source and provenance artifacts | `provenance/` |
| Original self-contained chronology audit | `chronology_onset_audit/` (preserve its internal layout) |

Keep the root limited to README, directory guide, task state, project instructions,
and setup/version-control files. Do not create new scientific `.py`, `.json`,
`.npz`, `.png`, or `.pdf` files at the root. A downloaded paper is a reference,
not a generated figure: use `docs/references/` and record its source/license.

## Paths and commands

- In `scripts/`, import `project_path` from `project_paths` and use
  `project_path('example_results.json')`, `project_path('example.npz')`, or
  `project_path('example_progress.png')`. It routes filenames by type.
- Do not use `Path(__file__).with_name(...)` for results/data/figures, or assume
  that the current working directory contains scientific assets.
- Run scripts from the root as `.\.venv\Scripts\python.exe -B scripts/example.py`.
  For `python -c` or archived-source imports, set
  `$env:PYTHONPATH = (Join-Path (Get-Location) 'scripts')` first.
- Use relative Markdown links: root README links to `figures/name.png`;
  reports in `docs/` link to `../figures/name.png` and `../scripts/name.py`.
- Use `$...$` for inline README mathematics and `$$...$$` for display mathematics;
  keep LaTeX compatible with KaTeX/GitHub and equations out of code spans.

## Preserve research and provenance

- Coordinate changes with active sessions; never revert their unrelated edits.
- Saved JSON/NPZ/figures and failed cases are evidence, not temporary files.
  Do not delete them just because a calculation finished.
- `maintenance/layout-migration.json` maps original paths and hashes to the new
  layout. Historical JSON source hashes and embedded recipes remain unchanged.
- Newly generated results must record actual current source hashes.
  `source_matches` accepts only exact bytes or the explicitly recorded layout
  migration; do not weaken it to accept later unrecorded edits.
- Ori checkpoint source validation remains strict. Exact historical resumption
  requires the original source revision recorded by the checkpoint; consult the
  migration's baseline commit. Do not bypass checkpoint validation.
- Preserve the original `chronology_onset_audit/` manifest and provenance sources.
  Existing historical checksum discrepancies are documented in the guide.

After path changes, run `.\.venv\Scripts\python.exe -B maintenance/verify_layout.py`.
It checks navigation, imports, migration integrity, and provenance without
running historical simulations or rewriting scientific results.
