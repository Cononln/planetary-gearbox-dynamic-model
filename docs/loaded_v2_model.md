# Loaded 18+7 DOF response model (v2.1)

The v2 workflow retains the 18 rigid planar coordinates of the sun, ring,
carrier, and three planets, and adds seven frequency-calibrated elastic-ring
coordinates for observation at the fixed 0 and 90 degree sensors.

Six coordinates retain the supplied mode pairs assigned to circumferential
orders 2, 3, and 4. A seventh coordinate provisionally assigns supplied
mode 15 (8711.641 Hz) to the axisymmetric radial-breathing order 0. This
coordinate preserves the mesh-frequency center line that otherwise cancels
exactly in an ideal three-planet spatial sum. Because only frequencies were
supplied, the order-0 identity must be verified from an FE mode-shape image.
Its provisional effective modal mass is 0.08 times the full-ring value; this
is an exposed path-calibration parameter and must ultimately be replaced by
an FE effective modal mass or fitted to the healthy measured 168-Hz line.

Unlike the original zero-load response generator, the v2 generalized force
contains a sun drive torque and a resisting carrier torque.  The current
illustrative carrier load is 100 N m.  The corresponding 20 N m sun torque
follows ideal power balance at the confirmed 5:1 speed ratio.  These values
are provisional and must be replaced by measured rig torque before absolute
amplitudes are compared with experiments.

The initial displacement is the loaded equilibrium at the first time step.
TVMS then modulates this nonzero mesh preload.  The healthy and 25-percent
cases remain fully compressive.  At 51.2 kHz the 50-percent case exhibits
incipient separation for about 0.010 percent of samples (roughly an 8-N
linear-contact overshoot).  This event is reported rather than hidden;
quantitative separation and re-impact require a unilateral-contact model.

The common deterministic transmission-error amplitude is provisionally
5 micrometres (10 micrometres peak-to-peak). The previous 0.08-micrometre
placeholder unrealistically suppressed the 168-Hz center line. This value
must be replaced by a measured no-load transmission-error trace or an
estimate tied to the manufactured gear-accuracy grade.

The tooth-pair load-sharing profile tends continuously to zero at contact
entry and exit. Because both contact ratios exceed one, the overlapping
pair retains a positive total mesh stiffness. A 25 or 50 percent
broken-width case removes the corresponding fraction of the 22-mm face
width only while tooth 1 is active. This avoids the nonphysical harmonic
comb created by a positive per-pair endpoint floor, but it is not yet a
full potential-energy tooth compliance calculation.

The former one-way SDOF path observers and empirical planet gains are not
used.  Ring-planet mesh forces couple directly into the seven elastic-ring
coordinates, and fixed sensor acceleration is projected from rigid-ring and
modal acceleration.  Because only modal frequencies, rather than measured
mode shapes and modal masses, are available, this observation model remains
qualitative.

Run:

    setup_paths
    output = simulate_loaded_v2_responses
    report = validate_loaded_v2_response

The principal results use a 3-s steady record at 51.2 kHz, providing about
0.333-Hz spectral resolution.  The first 0.2 s are discarded.
