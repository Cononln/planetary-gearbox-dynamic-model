# Core figure contracts

## Figure 1: method workflow

- Core conclusion: the method separates common speed phase from differential path phase, validates the path topology, and corrects phase without rescaling amplitude.
- Archetype: schematic-led composite.
- Output: double-column vector PDF/SVG plus 600-dpi TIFF.
- Evidence role: method comprehension, not quantitative validation.
- Reviewer risk: appearing as a generic pipeline; controlled by showing the physical sensor geometry, exact phase quantities, no-path return and common/differential outputs.

## Figure 2: dynamic moving-path field

- Core conclusion: a moving planet-to-housing path creates angle- and frequency-dependent magnitude and phase that one constant delay cannot represent.
- Archetype: quantitative grid.
- Hero evidence: dual-sensor differential phase spans at 168 and 1848 Hz.
- Validation evidence: path-Fourier/time-response sideband correlations.
- Reviewer risk: wrapped phase being mistaken for missing data; handled by explicit curve breaks and caption text.

## Figure 3: encoder-hidden validation

- Core conclusion: vibration-only estimation recovers a useful common coordinate and, on the dual-channel record, a repeatable differential path that agrees with an encoder-angle oracle.
- Archetype: quantitative grid.
- Hero evidence: 0.301 r/min carrier-speed RMSE and 4.72° differential-path RMSE.
- Controls: no-path, unconstrained, strict-3k and adaptive-3k variants.
- Reviewer risk: encoder leakage; handled by frozen choices and explicit scoring-only labels.

## Figure 6: sun-fault diagnosis and TPSVD

- Core conclusion: phase alignment recovers fault-family energy lost by naive fusion and modestly increases periodic-mode compactness, while remaining below the best single sensor.
- Archetype: asymmetric quantitative composite.
- Hero evidence: envelope spectrum and four fault-line amplitudes.
- Supporting evidence: shared-scale carrier-cycle matrices, cumulative TPSVD energy and leading-pattern order spectrum.
- Statistics: one 50-s record; no error bars or significance claims.
- Reviewer risk: display-scale bias; controlled by one common heatmap colour scale and identical spectral preprocessing.
