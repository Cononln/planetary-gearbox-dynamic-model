function report = step1_validate_18dof
%STEP1_VALIDATE_18DOF Validate dimensions, mechanics, and nominal frequencies.

setup_paths;
p = pg_parameters;
angles = [0, pi/7, 2*pi/3, 2*pi];
symmetryK = zeros(size(angles));
symmetryC = zeros(size(angles));
minEigenK = zeros(size(angles));
naturalHz = cell(size(angles));

for ia = 1:numel(angles)
    [M,C,K,meta] = pg_assemble_baseline(angles(ia),p);
    symmetryK(ia) = norm(K-K.','fro')/max(norm(K,'fro'),eps);
    symmetryC(ia) = norm(C-C.','fro')/max(norm(C,'fro'),eps);
    minEigenK(ia) = min(eig(K));
    lambda = eig(K,M);
    lambda = sort(real(lambda(lambda > 0)));
    naturalHz{ia} = sqrt(lambda)/(2*pi);

    assert(isequal(size(M),[18,18]), 'Mass matrix is not 18-by-18.');
    assert(all(diag(M)>0), 'Mass matrix contains a non-positive diagonal.');
    assert(symmetryK(ia)<1e-12 && symmetryC(ia)<1e-12, ...
        'Assembled matrices lost symmetry.');
    assert(minEigenK(ia)>0, 'Stiffness matrix is not positive definite.');

    % Mesh projections must be invariant to a uniform x/y translation.
    rigidX = zeros(18,1);
    rigidY = zeros(18,1);
    allBodies = [{p.map.sun},{p.map.ring},{p.map.carrier},p.map.planet];
    for ib = 1:numel(allBodies)
        rigidX(allBodies{ib}(1)) = 1;
        rigidY(allBodies{ib}(2)) = 1;
    end
    for ip = 1:p.model.nPlanet
        assert(abs(meta.mesh(ip).bSunPlanet.'*rigidX)<1e-12);
        assert(abs(meta.mesh(ip).bSunPlanet.'*rigidY)<1e-12);
        assert(abs(meta.mesh(ip).bRingPlanet.'*rigidX)<1e-12);
        assert(abs(meta.mesh(ip).bRingPlanet.'*rigidY)<1e-12);
    end
end

[M0,C0,K0] = pg_assemble_baseline(0,p);
[M2,C2,K2] = pg_assemble_baseline(2*pi,p);
periodicResidual = max([norm(M2-M0,'fro')/norm(M0,'fro'), ...
    norm(C2-C0,'fro')/norm(C0,'fro'), norm(K2-K0,'fro')/norm(K0,'fro')]);
assert(periodicResidual < 1e-12, 'Matrices are not 2-pi carrier-periodic.');

assert(abs(p.kin.f_s-10)<1e-12);
assert(abs(p.kin.f_c-2)<1e-12);
assert(abs(p.kin.f_mesh-168)<1e-12);
assert(abs(p.kin.f_sun_fault-24)<1e-12);
assert(abs(p.kin.f_planet_fault_double-10.8387096774194)<1e-10);
assert(abs(p.kin.f_encoder-2048)<1e-12);
assert(max(abs(p.mesh.sunPlanet.tePhase))<1e-12);
assert(max(abs(p.mesh.ringPlanet.tePhase))<1e-12);

% Confirm that the first-order interface returns exactly 36 finite states.
dx0 = pg_rhs_baseline(0,zeros(36,1),p);
assert(isequal(size(dx0),[36,1]) && all(isfinite(dx0)));

report.nDof = p.map.n;
report.symmetryResidualK = max(symmetryK);
report.symmetryResidualC = max(symmetryC);
report.minStiffnessEigenvalue = min(minEigenK);
report.periodicResidual = periodicResidual;
report.firstTenNaturalHzAtPhi0 = naturalHz{1}(1:min(10,numel(naturalHz{1}))).';
report.kinematics = p.kin;

fprintf('Stage 1 validation PASSED\n');
fprintf('  DOF count             : %d\n',report.nDof);
fprintf('  max K symmetry error  : %.3e\n',report.symmetryResidualK);
fprintf('  max C symmetry error  : %.3e\n',report.symmetryResidualC);
fprintf('  2pi periodic residual : %.3e\n',report.periodicResidual);
fprintf('  fs, fc, fm            : %.3f, %.3f, %.3f Hz\n', ...
    p.kin.f_s,p.kin.f_c,p.kin.f_mesh);
fprintf('  sun/planet fault freq : %.3f, %.3f Hz\n', ...
    p.kin.f_sun_fault,p.kin.f_planet_fault_double);
fprintf('  first ten f_n (Hz)    :');
fprintf(' %.2f',report.firstTenNaturalHzAtPhi0);
fprintf('\n');
end
