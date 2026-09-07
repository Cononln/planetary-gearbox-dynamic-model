# Repository collaboration instructions

## Before working

Read `README.md`, `docs/CURRENT_STATE.md`, and the source files involved in the
requested change. This is an unpublished research repository. Existing results
and code belong to the user. Preserve unrelated edits and keep legacy runners
available unless their removal is explicitly requested.

## Model entry points

- Current three-sensor model: `src/pg_parameters_paper_v3.m`.
- Mechanism checks: `scripts/validate_paper_v3_model.m`.
- Response generator: `scripts/simulate_paper_v3_three_channel_response.m`.
- Matrices/paths: `src/pg_coupled_matrices.m`, `src/pg_path_frf.m`,
  `src/pg_fixed_path_kinematics.m`, `src/pg_sensor_observation.m`.
- TVMS/fault timing: `src/pg_tvms.m`, `src/pg_fault_event_truth.m`.
- Current sun signed-waveform pilot: `scripts/test_sun_signed_waveform.py`,
  `src_py/pearson_scope.py`, `src_py/tidal_waveform.py`,
  `src_py/signed_waveform.py`.

## Research integrity

- Distinguish source fault-contact time, sensor arrival time, local ringing
  phase, and common rotational phase. Improvement of one is not proof of all.
- State simulation assumptions, experimental observations, numerical checks,
  and untested hypotheses separately. Report weak or negative findings.
- The available ring asset contains frequencies only. Analytical mode shapes
  and masses are effective priors, not measured FE mode-shape reductions.
- Current signed-waveform matching, TSA and channel fusion retain signal sign;
  envelope-based routes are separate historical baselines. Do not silently
  switch the representation in a comparison.
- Fixed-rank SVD pilots are not the original PSI-weighted mode-selection method.
- Keep encoder information outside tacholess estimation and parameter selection;
  use it only for explicitly separated posterior evaluation.
- Freeze preprocessing, windows, lag limits and evaluation support when comparing
  methods. Use held-out checks; do not select parameters on the final test score.
- Confirm sampling rate, column mapping, sensor polarity and sensitivity for
  each recording family. A sensitivity difference in mV/g is not itself a ratio.
- Keep units explicit. The current 21/31/84 tooth counts need geometry review:
  the standard equal-module relation gives 21 + 2*31 = 83, not 84. Do not silently
  alter these counts, frequencies, recurrence definitions or fault labels.

## Change and test policy

- Make scoped, reviewable changes. Use `codex/<task-name>` branches by default.
- Prefer independent clone/worktree directories for concurrent writers.
- Do not reset, overwrite or force-push another contributor's changes.
- Run `python scripts/run_smoke_tests.py` for relevant Python changes. These are
  standalone test scripts, not conventional pytest tests.
- For MATLAB model changes, run `setup_paths; validate_paper_v3_model` where
  MATLAB is available. Report checks you could not execute.
- Full simulations and measured-data experiments require the active user's
  request; do not launch them automatically during setup, upload or review.
- Store new research runs in unique timestamped directories under `results/`.
  Log parameters, input provenance and revision with new experiments.
- Update `docs/CURRENT_STATE.md` when the main route or supported conclusions
  change. Do not turn a proposal into a reported finding.

## Data and publication

Commit source, small required model assets, templates and technical notes.
Keep raw measurements, CAD, papers, manuscripts, generated results, credentials
and local environments out of Git. Inspect staged paths and sizes before a
commit; do not use `git add -f` to bypass these boundaries. An exceptional new
data asset needs an explicit scope review. Never create credentials or grant
external AI services access as an inferred part of a code change.
