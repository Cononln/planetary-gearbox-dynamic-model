# Paper-facing three-channel model (v3.0)

The paper-facing simulation uses the physical three-sensor layout:

| Sensor | Circumferential installation angle | Sensitive direction |
|---|---:|---:|
| S1 | 0 degrees | outward radial, 0 degrees |
| S2 | 120 degrees | outward radial, 120 degrees |
| S3 | 240 degrees | outward radial, 240 degrees |

Its main route is

\[
\text{18-DOF loaded mesh dynamics}
\rightarrow
\text{seven-coordinate low-order elastic-ring path prior}
\rightarrow
\text{three fixed sensor responses}.
\]

The moving ring--planet mesh force couples directly to the elastic-ring
coordinates.  The complex path field is evaluated for each source planet and
each sensor using the same coupled matrices.  The supplied modal frequencies
are used as resonance-frequency priors.  Since modal shapes, modal masses,
and measured damping are unavailable, the model is an effective transfer-path
model: it supports normalized path amplitude, relative phase, and effective
group-delay studies, but does not claim an experimentally calibrated absolute
FRF.

`pg_motion_state` defines the carrier trajectory.  If a controlled speed
variation is enabled, the integrated carrier angle drives mesh position,
TVMS, transmission-error phase, and fault-event timing together.  This gives
a fair controlled benchmark for encoder-free phase recovery.

For a partial-width broken tooth, `pg_tvms` now records the stiffness loss
caused by the defective contact.  `pg_fault_event_truth` takes the maximum
loss instant in every complete contact interval as the source event time
\(T_n\).  It is saved with every simulated response and is independent of
sensor observations.  It therefore provides the ground truth required for
the fault-event arrival-time alignment experiment.

Run the light mechanism checks with:

```matlab
setup_paths
report = validate_paper_v3_model
```

Generate 3-second, 51.2-kHz three-channel responses with:

```matlab
output = simulate_paper_v3_three_channel_response
outputVariation = simulate_paper_v3_three_channel_response( ...
    'includeSpeedVariation',true)
```

The old v2 two-channel scripts remain as legacy reproductions.  Only the v3
workflow should be used for the paper's three-channel mechanism and
fault-event-alignment validation figures.
