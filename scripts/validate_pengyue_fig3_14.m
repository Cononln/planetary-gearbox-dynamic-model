function report = validate_pengyue_fig3_14
%VALIDATE_PENGYUE_FIG3_14 Fast parameter/fault checks; no full simulation.

setup_paths;
p = pg_parameters_pengyue_fig314;
angles = p.reproduction.fig314.crackAnglesDeg;

expectedMesh = (13.3333333333333-2.09876543209877)*p.gear.zs;
expectedSunFault = expectedMesh/p.gear.zs*p.model.nPlanet;
assert(abs(p.kin.f_mesh-expectedMesh)<1e-9);
assert(abs(p.kin.f_sun_fault-expectedSunFault)<1e-9);
assert(abs(p.kin.f_mesh-191.111111111111)<1e-9);
assert(abs(p.kin.f_sun_fault-33.7254901960784)<1e-9);
assert(all([p.gear.zs p.gear.zp p.gear.zr]==[17 37 91]));
assert(p.reproduction.fig314.crackDepthMm==4);
assert(isequal(angles,[15 30 45 60 75]));

minRatio = zeros(size(angles));
for i = 1:numel(angles)
    f = pg_pengyue_root_crack_case(4,angles(i),p);
    minRatio(i) = min(f.energyLookup.stiffnessRatio);
end
assert(all(diff(minRatio)<0), ...
    'Crack stiffness reduction should increase with crack angle.');

report.meshFrequencyHz = p.kin.f_mesh;
report.sunFaultFrequencyHz = p.kin.f_sun_fault;
report.carrierFrequencyHz = p.kin.f_c;
report.planetRelativeFrequencyHz = p.kin.f_p_rel_c;
report.crackAnglesDeg = angles;
report.minimumPairStiffnessRatio = minRatio;
report.reported = p.reproduction.reported;
report.assumption = p.reproduction.assumption;

fprintf('Peng Yue Fig. 3-14 parameter validation PASSED\n');
fprintf('  fm       = %.6f Hz (paper: 191.1 Hz)\n',p.kin.f_mesh);
fprintf('  fs_fault = %.6f Hz (paper: 33.7 Hz)\n',p.kin.f_sun_fault);
fprintf('  fc       = %.6f Hz (paper: 2.1 Hz)\n',p.kin.f_c);
fprintf('  q        = %.3f mm\n',p.reproduction.fig314.crackDepthMm);
fprintf('  gamma    = 15 / 30 / 45 / 60 / 75 deg\n');
fprintf('  WARNING: face width, E/nu, absolute mean TVMS, Tin, numerical record\n');
fprintf('           length/tolerances and dimensional backlash scale are not\n');
fprintf('           fully reported by the thesis and remain explicit assumptions.\n');
end
