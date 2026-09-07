function p = pg_parameters_loaded_v2
%PG_PARAMETERS_LOADED_V2 Provisional loaded model used for the v2 figures.
% The torque values are engineering assumptions until measured rig loads are
% supplied.  Geometry, CAD mass properties, and kinematics remain unchanged.

p = pg_parameters;
p.model.name = 'ZLS160 loaded 18+7 DOF v2.1';

% Fixed-ring speed-controlled operating point.  For the present 5:1 speed
% ratio, ideal power balance gives Tin = Tout*omega_c/omega_s.
p.operating.outputTorque = 100; % N m, provisional carrier resisting torque
p.operating.inputTorque = p.operating.outputTorque* ...
    abs(p.kin.omega_c/p.kin.omega_s);
p.operating.torqueRippleFraction = 0.01;
p.operating.torqueRippleFrequencyHz = p.kin.f_s;

% The prescribed-speed drive needs finite torsional impedance in the
% perturbation model; otherwise the ideal planetary rolling mode leaves the
% static stiffness singular.  This provisional value represents the input
% shaft/controller stiffness and must ultimately be identified.
p.support.sun.ktheta = 3.0e3;
p.support.sun.ctheta = 2.0;
p.support.carrier.ktheta = 0;
p.support.carrier.ctheta = 2.0;

% Each individual tooth pair must enter and leave contact continuously.
% A positive per-pair floor creates a stiffness jump when that pair is
% added or removed and injects a nonphysical broadband harmonic comb.  The
% contact ratios exceed one, so overlapping pairs keep the total mesh
% stiffness positive even though each pair tends smoothly to zero.
p.tvms.profileFloor = 0;
p.tvms.profileSin2 = 1;

% A common mesh-frequency component is required in addition to loaded TVMS
% modulation. The former 0.08-um value was below a realistic micrometre
% scale and made fault sidebands dominate the common mesh line. A 5-um
% amplitude (10-um peak-to-peak) is used as an explicit provisional static
% transmission-error value. Replace it with a measured no-load TE trace or
% gear-accuracy-based estimate before quantitative experimental comparison.
p.mesh.sunPlanet.teAmplitude = 5.0e-6;
p.mesh.ringPlanet.teAmplitude = 5.0e-6;
p.mesh.sunPlanet.zeta = 0.045;
p.mesh.ringPlanet.zeta = 0.045;
p.ringModal.zeta = 0.020;

% Preserve the mean (spatial order n=0) ring response that is cancelled by
% an ideal three-planet model containing only n=2,3,4 flexural pairs.  The
% supplied frequency list has no mode shapes, so mode 15 is provisionally
% assigned as the axisymmetric radial-breathing mode.  This assignment must
% be checked against an exported FE mode-shape image before quantitative
% comparison; it is not a fitted 168-Hz sinusoid.
p.ringModal.sourceModes = [7:12,15];
p.ringModal.waveNumbers = [2,2,3,3,4,4,0];
p.ringModal.family = {'cos','sin','cos','sin','cos','sin','cos'};
p.ringModal.axisymmetricSourceMode = 15;
% Frequencies alone do not determine modal mass. The first six analytical
% flexural modes retain their circumference-based masses. For the mean
% ring-to-housing path, 0.08 of the full ring mass is a provisional local
% participating-mass estimate. It is selected as an explicit calibration
% parameter, not as a hidden spectrum gain, and must be replaced by an FE
% effective modal mass or calibrated from the healthy 168-Hz response.
p.ringModal.modalMassScale = [ones(1,6),0.08];
p.ringModal.axisymmetricAssumption = [ ...
    'Mode 15 at 8711.641 Hz is provisionally treated as the n=0 radial ' ...
    'breathing mode until its FE mode shape is available.'];

p.v2.assumption = ['Provisional 100 N m carrier load; replace with the ' ...
    'measured operating torque before quantitative comparison.'];
p.v2.contactModel = ['Loaded linear contact with positive mean preload. ' ...
    'Incipient separation in the severe case is diagnosed but requires ' ...
    'a future unilateral-contact revision.'];
p.v2.spectralBalance = [ ...
    'A provisional n=0 breathing coordinate retains the physical mesh ' ...
    'center line; no harmonic is injected into the sensor signal.'];
end
