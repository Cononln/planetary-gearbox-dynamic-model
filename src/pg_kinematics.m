function kin = pg_kinematics(zs, zp, zr, sunRpm, encoderPpr)
%PG_KINEMATICS Nominal Willis kinematics for a fixed-ring planetary stage.

kin.omega_s = 2*pi*sunRpm/60;
kin.omega_r = 0;
kin.omega_c = (zs*kin.omega_s + zr*kin.omega_r)/(zs + zr);
kin.omega_p_rel_c = -(zs/zp)*(kin.omega_s - kin.omega_c);
kin.omega_p = kin.omega_c + kin.omega_p_rel_c;

kin.f_s = kin.omega_s/(2*pi);
kin.f_r = 0;
kin.f_c = kin.omega_c/(2*pi);
kin.f_p_rel_c = abs(kin.omega_p_rel_c)/(2*pi);
kin.f_p = kin.omega_p/(2*pi);
kin.f_mesh = zs*abs(kin.omega_s - kin.omega_c)/(2*pi);
kin.f_sun_fault = 3*kin.f_mesh/zs;
kin.f_planet_fault_double = 2*kin.f_mesh/zp;
kin.encoder_ppr = encoderPpr;
kin.f_encoder = encoderPpr*kin.f_c;
end
