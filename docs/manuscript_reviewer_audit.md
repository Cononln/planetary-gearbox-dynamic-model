# Reviewer-style audit for the phase-decoupling manuscript

## Overall verdict

The manuscript now has a coherent method paper structure and a defensible core contribution: it separates common speed-induced mesh phase from sensor-dependent moving-path phase, constrains the latter with three-planet spatial structure, and rejects unsupported sensor graphs before amplitude-preserving correction. The present evidence is sufficient for a mechanism-focused draft, but not yet sufficient for a strong MSSP/Q1 claim of general diagnostic superiority. The main remaining risk is experimental scope, not the internal storyline.

## 1. Novelty and contribution

**Status: conditional pass**

### Evidence already present

- The observation model explicitly distinguishes common phase and differential path phase rather than treating alignment as one cycle shift.
- Pairwise cross-phasors cancel the common phase algebraically and expose an angle- and harmonic-dependent path field.
- Three-planet spatial regularization is coupled to a sensor-graph projection and a no-path return, rather than forcing all measured channels into fusion.
- Unit-modulus correction separates phase equalization from amplitude weighting and preserves analytic amplitude to numerical precision.
- The healthy-record ablation provides a useful negative result: an unconstrained model can increase phase concentration without corresponding fault evidence.

### Reviewer risk

- Without fair implementation of Optimal-PSA, D-norm TSA and recent multi-channel weighted fusion, a reviewer may regard the work as a combination of tacholess TSA, Fourier path fitting and channel selection.
- The strict 3k prior is physically motivated, but the adaptive non-3k residual requires sensitivity and transfer evidence to establish that it is more than record-specific flexibility.

### Required action

- Add strong baselines under identical preprocessing, record splits and encoder policy.
- Add healthy-to-fault transfer to show that the path model does not learn fault morphology.

## 2. Writing clarity and reproducibility

**Status: pass with minor revision**

### Evidence already present

- Figure 1 maps the full method from synchronized vibration to common/differential diagnostic outputs.
- Equations define the observation model, common-phase estimator, spatial fit, graph projection, topology score and unit-modulus correction.
- Algorithm 1 provides an implementation sequence.
- The discussion clearly distinguishes operational equalization from absolute mechanical-zero recovery and physical path identification.
- Negative results and the single-record evidence boundary are stated explicitly.

### Remaining reproducibility gaps

- The final paper should tabulate the exact harmonic set, analytic-band width, smoothing length, maximum spatial order, ridge grid, off-grid penalty grid and topology thresholds.
- Runtime, hardware and software versions are not yet reported.
- The encoder decoding and sensitivity normalization rules should be summarized in a short reproducibility subsection or supplement.

## 3. Experimental strength

**Status: needs new experiments**

### Evidence already present

- Dynamic path truth supports the moving-path mechanism at 168 and 1848 Hz.
- Encoder-hidden dual-channel validation provides an external score for common and differential phase.
- Three-channel records exercise topology selection and expose an unreliable channel.
- The six-configuration ablation isolates common phase, spatial structure, topology protection and cross-fitting.
- Envelope and TPSVD analyses show downstream diagnostic consequences without claiming that aligned fusion beats the best sensor.

### Main weakness

- Most measured evidence is at 600 r/min and zero-load labels, with one 50-s record per condition and three complete mechanical repeats. This does not establish cross-speed, cross-load or population-level robustness.

### Required action

- Add repeated records at multiple speeds and at least one nonzero load.
- Freeze parameters before test records and report record-level variation.

## 4. Evaluation completeness

**Status: needs new experiments**

### Metrics already present

- Mesh-phase and carrier-speed RMSE against hidden encoder references.
- Differential-path RMSE where the oracle is repeatable.
- Inter-channel phase residual and coherence.
- HMPC and PCC as complementary phase/waveform scores.
- Fault-family signal-to-background ratio, cycle dispersion and TPSVD energy.
- Amplitude-preservation error.

### Missing comparisons and analyses

- Optimal phase synchronous averaging under the same vibration-only input.
- D-norm tacholess TSA and a conventional tacholess TSA baseline.
- Recent multi-channel weighted fusion, both before and after phase equalization.
- Encoder-based TSA as an upper-reference pipeline, clearly excluded from the vibration-only comparison ranking.
- Parameter sensitivity, runtime and memory scaling.
- Independent-record uncertainty and paired statistical comparison.

## 5. Method-design soundness

**Status: conditional pass**

### Strengths

- The identifiability boundary is stated: only differential path phase is recoverable from synchronized vibration without an additional reference.
- Common-phase estimation is frozen before path-topology selection, so a rejected path cannot corrupt the tacholess speed estimate.
- Held-out topology selection and the no-path return reduce the risk of automatic over-compensation.
- Common and differential outputs preserve localized evidence that coherent averaging could remove.

### Unresolved risks

- Within-record path fitting may absorb repeatable fault-related differential content.
- Adaptive non-3k content may be condition-specific when mounting or load changes.
- The third-channel rejection mechanism has not been separated among sensor direction, installation, housing mode shape and local SNR.
- The dynamic model validates normalized phase and sidebands, but not absolute acceleration.

### Required action

- Train the path model on healthy records and apply it unchanged to fault records.
- Test mounting/load/speed transfer and report when the safeguard returns no path.
- Include failure cases as part of the method result, not as discarded runs.

## Submission decision gate

### Ready now

- Method formulation and narrative.
- Mechanistic simulation.
- Encoder-hidden phase validation.
- Controlled within-record ablation.
- Publication-style core figures and bounded claims.

### Must be added before a strong Q1 submission

1. Multi-speed and multi-load repeated-record evaluation.
2. Fair comparison with Optimal-PSA, D-norm TSA, conventional tacholess TSA and multi-channel weighted fusion.
3. Healthy-to-fault frozen-model transfer.
4. Record-level uncertainty/statistics.
5. Runtime and parameter-sensitivity analysis.

Until these five items are complete, the paper should be presented as a validated phase-equalization mechanism rather than a generally superior planetary-gearbox diagnostic framework.
