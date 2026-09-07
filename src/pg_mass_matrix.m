function M = pg_mass_matrix(p)
%PG_MASS_MATRIX Assemble the diagonal 18-by-18 mass matrix.

M = zeros(p.map.n);
M = putBody(M, p.map.sun, p.body.sun.mass, p.body.sun.J);
M = putBody(M, p.map.ring, p.body.ring.mass, p.body.ring.J);
M = putBody(M, p.map.carrier, p.body.carrier.mass, p.body.carrier.J);
for i = 1:p.model.nPlanet
    M = putBody(M, p.map.planet{i}, p.body.planet.mass, p.body.planet.J);
end
end

function M = putBody(M, idx, mass, inertia)
M(idx, idx) = diag([mass, mass, inertia]);
end
