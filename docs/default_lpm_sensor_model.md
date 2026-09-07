# Default 18-DOF lumped-parameter sensor model

## Gearbox dynamics

The sun, rigid ring, carrier, and three planets each retain planar
`[x,y,theta]` motion.  The 18 generalized coordinates satisfy

    M qdd + C(t) qd + K(t) q = F_TE(t).

The time dependence is caused by carrier-angle-dependent line-of-action
vectors and tooth-level TVMS.  For mesh `i`, the signed contact force
recovered from the solved displacement and velocity is

    F_i = k_i(t) [e_i(t)-b_i(t)'q]
        + c_i(t) [edot_i(t)-b_i(t)'qd].

Thus the sensor path is driven by the dynamic contact force, rather than by
an independently synthesized fault impulse.

## Broken-tooth TVMS

Each tooth pair remains active over its contact-ratio interval.  A smooth
`sin^2` load-sharing profile is integrated over 80 face-width slices.  The
25 and 50 percent cases remove 5.5 and 11 mm, respectively, from the 22-mm
sun-planet face width when the designated sun tooth is in contact.

## Fixed sensor observation

The three ring-planet mesh positions are

    psi_i(t) = phi_c(t) + 2*pi*(i-1)/3.

For sensor `s` at angle `theta_s`, the moving-source path weight is

    w_si(t) = w_floor + (1-w_floor)
              exp{kappa[cos(psi_i-theta_s)-1]}.

The nonzero floor preserves the 168-Hz mesh carrier.  The third spatial
harmonic of the three equally spaced moving paths produces the expected
`168 +/- 3*fc = 162/174 Hz` sidebands.

An ideal `[1,1,1]` planet-source gain makes the `fc` and `2*fc` spatial
orders cancel exactly.  The default `[1.00,0.96,1.04]` factors introduce a
zero-mean 4-percent load/path asymmetry, representing unequal planet load
sharing, pin-position tolerances, and non-axisymmetric casing transmission.
This restores the `fm +/- fc` and `fm +/- 2*fc` components.  The factors are
explicit calibration parameters and must be identified from measured load
sharing or sideband amplitudes before quantitative use.

Each sensor location is represented by a local SDOF observer:

    m_h ydd_s + c_h yd_s + k_h y_s
      = sum_i w_si(t) F_ring-planet,i(t).

The current observer resonance is 1839.17 Hz, calibrated from the supplied
first elastic ring-frequency pair, with damping ratio 0.02 and effective
mass equal to half the ring mass.  This is a lumped path approximation, not
an FE ring or a claim of measured modal mass.

The reported fixed-sensor acceleration is the local observer acceleration
plus the rigid-ring radial acceleration projected at 0 or 90 degrees.

## Interpretation boundary

This model is suitable for mechanism and phase-alignment validation.  The
absolute acceleration level, observer mass, path-window shape, support
stiffness, and damping remain calibration parameters until impact-test or
operating-response data are used for identification.
