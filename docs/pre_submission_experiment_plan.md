# Pre-submission experiment plan

## Objective

Upgrade the current mechanism-focused evidence into a reviewer-resistant Q1 evaluation without changing the core method. Every comparison must use the same raw records, channel calibration, analytic bands, record splits and encoder policy.

## P0: experiments required before submission

### E1. Cross-speed and cross-load robustness

- Conditions: retain 600 r/min and add 900 r/min; add any lower speed already available if acquisition quality is sufficient.
- Loads: retain the existing zero-load condition and add at least one nonzero load available on the test bench.
- Fault states: healthy, sun fault and planet fault; include at least one severity level for each fault family.
- Repetition: target at least five independently acquired 50-s records per condition. Do not split one record into pseudo-replicates for statistical inference.
- Frozen settings: select harmonics, bandwidths, spatial order and thresholds on training conditions only.
- Outputs: carrier-speed RMSE, HMPC, PCC, residual/coherence, common/differential fault-family SNR, TPSVD rank-3 energy and no-path return rate.

### E2. Fair strong-baseline comparison

Implement the following pipelines with identical calibration and frequency bands:

1. Naive equal-weight multi-channel fusion.
2. Cross-correlation/TSA alignment.
3. Conventional vibration-only tacholess TSA.
4. D-norm TSA.
5. Optimal phase synchronous averaging.
6. Recent multi-channel weighted fusion.
7. Proposed common phase only.
8. Proposed full phase decoupling.
9. Proposed full phase decoupling followed by the same weighted-fusion rule used in baseline 6.
10. Encoder-based TSA as an external upper reference, not as a vibration-only competitor.

Fairness rules:

- The encoder cannot be used for hyperparameter selection in vibration-only methods.
- All methods receive the same record duration and channel set.
- Weighted fusion must be evaluated both before and after phase equalization so that phase and amplitude contributions are not conflated.
- Report method failures and no-path returns; do not silently remove them.

### E3. Healthy-to-fault frozen-path transfer

- Fit common-phase and path-model settings on healthy records only.
- Freeze the spatial penalty, topology thresholds, harmonic set and filter bands.
- Apply the frozen model to sun- and planet-fault records at the same speed/load.
- Repeat with leave-one-speed-out and leave-one-load-out transfer.
- Primary question: does fault visibility improve without reducing the diagnostic differential residual or fitting fault morphology?
- Failure criterion: a gain in HMPC accompanied by a decrease in fault-family contrast or a large increase on healthy records is treated as over-compensation.

### E4. Independent-record statistics

- Unit of replication: independently acquired record, not carrier cycle.
- Report mean, standard deviation and 95% confidence interval for each main metric.
- Use paired comparisons because methods are applied to the same records.
- Select a paired t-test only after checking approximate difference normality; otherwise use a paired Wilcoxon test.
- Correct for multiple comparisons when ranking several baselines.
- Keep the current three-cycle RF displays as within-record visualizations, not statistical replicates.

## P1: robustness and reviewer questions

### E5. Parameter sensitivity

Vary one item at a time around the frozen default:

- retained harmonic set;
- analytic-band width;
- smoothing length for phase increments;
- maximum spatial order Q;
- ridge and off-grid penalty;
- topology admission gain and adverse-fold thresholds;
- sensor subset and one-channel dropout;
- additive noise/SNR.

Plot performance plateaus rather than selecting the best value independently for every record.

### E6. Runtime and scaling

- Report CPU model, RAM, software versions and whether processing is single-threaded.
- Measure wall time for band extraction, common-phase estimation, each path fit, graph selection and downstream diagnosis.
- Report total runtime per 50-s record and real-time factor.
- Scale record length, channel count and harmonic count separately.

### E7. Third-channel failure analysis

- Compare raw SNR, harmonic prominence, sensor direction and coherence for 0°, 120° and 240° channels.
- Inspect whether the rejected 240° edges coincide with local anti-resonance or housing-mode nodes.
- Repeat at a changed load or remounting if possible.
- Report whether the topology safeguard consistently rejects the same edge or adapts to the condition.

## Planned paper outputs

### Main-text additions

- Table: cross-speed/cross-load mean ± variability for the proposed method and strongest baselines.
- Figure: paired record-level baseline comparison for HMPC, path residual and fault-family SNR.
- Figure: healthy-to-fault transfer matrix across speed/load.
- Figure: compact parameter-sensitivity and runtime panel.

### Supplementary additions

- Exact preprocessing and hyperparameter table.
- All record identifiers and train/test assignments.
- Per-condition no-path decisions and failure cases.
- Source-data CSV files for every quantitative panel.

## Stop/go criteria

Proceed to submission only if:

1. The full method improves held-out phase metrics across speeds/loads without inflating healthy fault contrast.
2. At least one downstream fault metric improves consistently over both naive fusion and phase-only baselines.
3. The result remains competitive with Optimal-PSA and D-norm TSA under the same vibration-only inputs.
4. Healthy-trained path models transfer to fault records or safely return no path.
5. Gains are not driven by one record or one sensor, and runtime is practical for offline condition monitoring.
