# Manuscript story and evidence map

## Working title

Tacholess multi-channel phase decoupling with three-planet spatial constraints for planetary gearbox fault diagnosis

## One-sentence story

Planetary-gearbox channels cannot be fused safely until the common speed-induced phase and the channel-dependent moving-path phase are separated; multi-harmonic common-phase estimation, three-planet spatial regularization, and held-out topology selection enable an amplitude-preserving correction that improves phase consistency without forcing unsupported sensor paths.

## Compact outline

1. Introduction: expose destructive multi-channel fusion as a phase-identification problem rather than a weighting problem.
2. Problem formulation: derive the common mesh phase and sensor-path phase observation model and state its gauge ambiguity.
3. Method: estimate common phase, fit pairwise three-planet path phasors, project them onto a sensor graph, select a safe topology by held-out carrier cycles, and retain common/differential outputs.
4. Validation: use a lumped-parameter/ring-mode model to demonstrate angle-frequency path phase; then use encoder-hidden measured records for phase scoring, severity tests, ablation, and fault-feature extraction.
5. Discussion: distinguish operational phase equalization from physical path recovery, compare with waveform-matching synchronization and amplitude-weighted fusion, and state the present generalization boundary.
6. Conclusions: summarize the physical insight, phase result, safeguard behavior, and diagnostic gain.

## Paragraph roles

- Introduction P1 (opening): planetary gearboxes generate spatially modulated vibration and require reliable weak-fault extraction.
- Introduction P2 (challenge): speed fluctuation and moving transfer paths create two different phase errors that become destructive during fusion.
- Introduction P3 (prior synchronization): tachometer-based and tacholess TSA address cycle timing but do not explicitly identify channel-dependent path phase.
- Introduction P4 (prior fusion): multi-channel fusion improves information coverage but amplitude weighting alone cannot prevent phase cancellation.
- Introduction P5 (gap): speed phase and path phase share the same observations and are only identifiable under cross-channel, cross-harmonic, and spatial constraints.
- Introduction P6 (solution): introduce the proposed joint phase-decoupling pipeline.
- Introduction P7 (evidence): summarize the most important model and measured-data results.
- Introduction P8 (contributions): list the three technical contributions and one validation contribution.

## Major claims and evidence

| Claim | Evidence | Status |
|---|---|---|
| Moving paths require an angle-frequency phase field rather than one constant delay. | The 18+7-DOF path model gives dual-channel phase spans of 181.37 deg at 168 Hz and 196.22 deg at 1848 Hz; path-Fourier sideband correlations are 0.9986 and 0.99995. | Supported as mechanistic simulation evidence. |
| Multi-harmonic common phase improves tacholess speed tracking. | On the 50% sun-fault record, mesh-frequency RMSE decreases from 1.518 Hz for single-channel H1 to 0.421 Hz; carrier-speed RMSE is 0.301 rpm. | Supported for one measured operating condition. |
| The differential path can be estimated accurately in the dual-channel sun-fault record. | Encoder-angle path-phase RMSE is 4.72 deg; coherence increases from 0.563 to 0.791. | Supported for one dual-channel record only. |
| A fixed three-sensor topology is unsafe. | Cross-validation rejects the full graph in healthy and PF50 records and selects reliable pairs; P025/P050/P075 also retain only the 0-120 deg edge. | Supported across the available 600-rpm records. |
| The full correction improves operational phase consistency while preserving amplitude. | PF50 HMPC increases from 0.252 to 0.502, PCC from 0.557 to 0.728, and normalized mean cycle SD decreases from 0.715 to 0.676; amplitude error is below 3.4e-16. | Supported as within-record ablation. |
| The spatial constraint prevents spurious correction. | Removing the spatial constraint raises healthy HMPC to 0.592, whereas the safe full method returns no path and healthy HMPC remains 0.136. | Supported as within-record ablation. |
| Phase correction can recover features lost by naive fusion. | On the sun-fault record, phase-aligned fusion gives 3.48 dB fault-line-family gain over naive fusion and raises TPSVD rank-3 energy fraction from 0.864 to 0.884. | Supported for one record; it does not exceed the best single sensor. |
| The method accurately recovers physical path phase on all three-channel records. | P025/P050/P075 encoder-path RMSE is 83.3-103.6 deg. | Not supported; must not be claimed. |
| The method generalizes across speeds and loads or improves classification accuracy. | No independent repeated-condition statistics or classifier experiment is available. | Needs evidence. |

## Reviewer-risk register

1. Only one 50-s record is available per condition in the controlled ablation, giving three complete 31-carrier-turn repeats.
2. The current evidence is concentrated at 600 rpm; cross-speed and cross-load generalization remains untested.
3. The third sensor is frequently rejected for differential-path correction. This is a safeguard result, not evidence that three-channel correction always succeeds.
4. The dual-channel sun-fault record supports accurate differential-path recovery, whereas the three-channel planet-fault records support only operational phase equalization.
5. Fault-feature gains are deterministic within-record results and must not be described as population-level diagnostic accuracy.
6. The ring-path model contains calibrated modal parameters and supports the mechanism, not absolute response prediction.

## Figure and table plan

- Fig. 1: dynamic moving-path mechanism and method workflow (new schematic still required).
- Fig. 2: modeled angle-dependent AM/PM and phase truth.
- Fig. 3: encoder-hidden phase/path validation on the sun-fault record.
- Fig. 4: RF-style before/after cycle comparison.
- Fig. 5: six-configuration ablation envelope spectra.
- Fig. 6: sun-fault envelope/TPSVD diagnosis comparison.
- Table 1: rig, gearbox, sensor, and acquisition parameters.
- Table 2: phase-estimation baseline comparison.
- Table 3: P025/P050/P075 topology-safeguard results.
- Table 4: six-configuration controlled ablation.

