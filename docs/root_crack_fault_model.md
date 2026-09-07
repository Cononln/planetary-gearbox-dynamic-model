# Energy-based single-sun-tooth root-crack model

The crack model is a separate fault option and leaves the partial-width
broken-tooth baseline reproducible. It follows the modeling assumptions in
Shen (2023): a single sun tooth has a root crack of extension `q`, the crack
angle is 65 degrees, and only the tooth bending and shear compliances are
changed. Hertz contact and axial compression remain healthy.

For each contact position, the sun and planet teeth are represented as
tapered involute cantilevers. The confirmed addendum/root diameters and the
nominal base-circle geometry determine the local tooth thickness. Numerical
integration gives bending, shear, and axial compliances. The cracked sun
section subtracts the transverse crack penetration `q*sin(v)` over the
root-side reach `q*cos(v)`. The resulting cracked/healthy pair-stiffness
ratio scales the calibrated tooth-pair TVMS only while sun tooth 1 is active.
Overlapping healthy pairs remain unchanged.

Create the reference crack levels with:

    c05 = pg_root_crack_case(0.5,p);
    c15 = pg_root_crack_case(1.5,p);
    c25 = pg_root_crack_case(2.5,p);

The optional first-order speed ripple follows

    n(t) = n0 + am*cos(2*pi*fs*t + phase).

Because Shen (2023) does not report `am`, the default ripple fraction is
zero. A controlled value can be supplied through `SpeedRippleFraction` for
a with/without-ripple ablation. It is never silently fitted.

The paper uses module 2 mm, whereas the present gearbox uses module 1.5 mm.
Consequently the literal 0.5/1.5/2.5-mm cases correspond to larger
dimensionless crack depths here (`q/m = 0.33/1.00/1.67`) than in the paper
(`q/m = 0.25/0.75/1.25`). The 2.5-mm case reaches the numerical minimum
ligament and should be interpreted as a near-severed limiting case, not as a
quantitatively matched experiment. Module-scaled reference equivalents are
0.375, 1.125, and 1.875 mm. Actual crack dimensions should replace these
study levels before experimental comparison.

Run `validate_root_crack_model` to export the contact-position stiffness
ratios, tooth-indexed TVMS, fault-order spectrum, and numerical metrics.

Run `simulate_root_crack_responses(20)` for one-second 0/90-degree
acceleration responses under a 20-N-m sun input and the ideal 100-N-m
carrier resisting torque. The solver initializes the static loaded
equilibrium, so the crack TVMS acts on a nonzero mesh preload. Passing zero
retains the former no-load comparison. Outputs include the 0-240-Hz
mesh/sideband spectra and resonance-band envelope spectra. The corrected
paper-comparison cases use a 5-um transmission-error carrier, 50-percent
crack face coverage, and module-scaled crack depths 0.375/1.125/1.875 mm.
The runner asserts that 168 Hz exceeds every 6-Hz-grid modulation line.
