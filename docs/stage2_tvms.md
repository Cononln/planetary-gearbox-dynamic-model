# Stage 2 — Partial-width broken-tooth TVMS

The fault is a face-width-direction partial tooth break cut to the root circle.
For the 22 mm active gear face width:

- 25% severity removes 5.5 mm;
- 50% severity removes 11.0 mm;
- the radial removal depth is 3.375 mm.

For a contact ratio epsilon between one and two, a tooth pair entering at mesh
counter j stays active for epsilon mesh periods. Therefore one or two tooth
pairs contribute at each instant. A sin-squared load-sharing profile is zero
at pair entry and exit, avoiding artificial stiffness jumps. Its scale is
calibrated so that the healthy time average equals the Stage-1 mean mesh
stiffness.

Every active tooth pair is identified explicitly. For a damaged tooth, the
local inverse compliance is integrated only across retained face-width slices.
Because the cut reaches the root circle, removed slices contribute zero tooth
stiffness. At double-pair contact, the overlapping healthy tooth remains
unchanged. Consequently the complete mesh stiffness is not a constant
(1-severity) multiple of the healthy TVMS.

Phase zero is defined as sun tooth 1 entering the planet-1 reference
sun-planet mesh. The other two sun meshes use the physically compatible
seven-tooth offsets. Their aggregate fault pattern repeats every seven mesh
cycles, giving 168/7 = 24 Hz. Planet-1 tooth 1 is tracked in both its
sun-planet and ring-planet mesh. The opposite-side ring contact is offset by
floor(31/2) = 15 tooth entries, producing two fault-contact events per
planet-relative revolution, approximately 10.839 Hz.

The current contact-path profile and mean mesh stiffness are calibration
parameters. The tooth identity, retained-width integral, contact ratios, fault
dimensions, and event kinematics are fixed physical parts of the model.
