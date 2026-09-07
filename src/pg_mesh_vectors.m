function mesh = pg_mesh_vectors(phi_c, p)
%PG_MESH_VECTORS Build line-of-action projection vectors for all meshes.
% The vectors b satisfy delta_mesh = b.'*q - e(t).

nq = p.map.n;
mesh = repmat(struct('planetAngle', 0, 'bSunPlanet', zeros(nq,1), ...
    'bRingPlanet', zeros(nq,1), 'nSunPlanet', zeros(2,1), ...
    'nRingPlanet', zeros(2,1)), 1, p.model.nPlanet);

for i = 1:p.model.nPlanet
    psi = phi_c + p.gear.planetPhase0(i);
    er = [cos(psi); sin(psi)];
    et = [-sin(psi); cos(psi)];

    % External sun-planet and internal ring-planet working lines of action.
    nsp = sin(p.gear.alphaSunPlanet)*er + ...
          cos(p.gear.alphaSunPlanet)*et;
    nrp = -sin(p.gear.alphaRingPlanet)*er + ...
           cos(p.gear.alphaRingPlanet)*et;

    bsp = zeros(nq,1);
    bsp(p.map.sun(1:2)) = nsp;
    bsp(p.map.planet{i}(1:2)) = -nsp;
    bsp(p.map.sun(3)) = p.gear.baseRadiusSun;
    bsp(p.map.planet{i}(3)) = p.gear.baseRadiusPlanet;

    brp = zeros(nq,1);
    brp(p.map.ring(1:2)) = nrp;
    brp(p.map.planet{i}(1:2)) = -nrp;
    brp(p.map.ring(3)) = p.gear.baseRadiusRing;
    brp(p.map.planet{i}(3)) = -p.gear.baseRadiusPlanet;

    mesh(i).planetAngle = psi;
    mesh(i).bSunPlanet = bsp;
    mesh(i).bRingPlanet = brp;
    mesh(i).nSunPlanet = nsp;
    mesh(i).nRingPlanet = nrp;
end
end
