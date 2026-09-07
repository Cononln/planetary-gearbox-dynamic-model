function map = pg_dof_map(nPlanet)
%PG_DOF_MAP Return the index map for the 18-DOF model.
% Each rigid body has [x, y, theta].

if nargin < 1
    nPlanet = 3;
end
if nPlanet ~= 3
    error('pg_dof_map:PlanetCount', ...
        'This study is frozen to three planets; received %d.', nPlanet);
end

map.sun = 1:3;
map.ring = 4:6;
map.carrier = 7:9;
map.planet = cell(1, nPlanet);
for i = 1:nPlanet
    map.planet{i} = 9 + (3*(i-1) + (1:3));
end
map.n = 9 + 3*nPlanet;
map.names = {'sun-x','sun-y','sun-theta', ...
             'ring-x','ring-y','ring-theta', ...
             'carrier-x','carrier-y','carrier-theta', ...
             'planet1-x','planet1-y','planet1-theta', ...
             'planet2-x','planet2-y','planet2-theta', ...
             'planet3-x','planet3-y','planet3-theta'};
end
