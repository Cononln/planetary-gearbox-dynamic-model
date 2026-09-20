function response = pg_simulate_pengyue_dte(duration,p,fault,options)
%PG_SIMULATE_PENGYUE_DTE ODE45 DTE response for the Fig. 3-14 reproduction.
% Output is the sun-planet relative displacement x_sp^i of Eq. (2-1), which
% is the response quantity used in Sec. 3.4.1 rather than casing acceleration.

arguments
    duration (1,1) double {mustBePositive}
    p struct
    fault struct
    options.fs (1,1) double {mustBePositive} = p.model.fs
    options.RelTol (1,1) double {mustBePositive} = 1e-7
    options.AbsTol (1,1) double {mustBePositive} = 1e-9
    options.MaxStep (1,1) double {mustBePositive} = 1/(30*p.kin.f_mesh)
end

fs = options.fs;
time = (0:1/fs:duration).';
n = p.map.n;
[M0,C0,K0,meta0] = pg_assemble_pengyue_system(0,p,fault); %#ok<ASGLU>
F0 = pg_pengyue_excitation_force(0,p,meta0);
% Static loaded equilibrium is a reproducible initial condition; the runner
% later discards several carrier cycles before spectral analysis.
q0 = K0\F0;
x0 = [q0;zeros(n,1)];

odeOptions = odeset('RelTol',options.RelTol,'AbsTol',options.AbsTol, ...
    'MaxStep',options.MaxStep);
[time,state] = ode45(@rhs,time,x0,odeOptions);
q = state(:,1:n);
qd = state(:,n+1:end);

mesh = pg_mesh_vectors(0,p);
dte = zeros(numel(time),p.model.nPlanet);
transmissionError = zeros(size(dte));
for it = 1:numel(time)
    motion = pg_motion_state(time(it),p);
    for ip = 1:p.model.nPlanet
        phase = motion.meshPhaseRad+p.mesh.sunPlanet.tePhase(ip);
        e = p.mesh.sunPlanet.teAmplitude*sin(phase);
        transmissionError(it,ip) = e;
        dte(it,ip) = mesh(ip).bSunPlanet.'*q(it,:).'-e;
    end
end

response.time = time;
response.fs = fs;
response.q = q;
response.qd = qd;
response.dteSunPlanet = dte;
response.transmissionErrorSunPlanet = transmissionError;
response.fault = fault;
response.integration = 'MATLAB ode45 (Runge-Kutta 4-5)';
response.coordinateFrame = 'carrier-fixed reproduction assembly';
response.note = ['Friction and dimensionless backlash are logged but not ' ...
    'applied because the thesis does not report the displacement scale bc.'];

    function dx = rhs(t,x)
        qNow = x(1:n);
        qdNow = x(n+1:end);
        [M,C,K,meta] = pg_assemble_pengyue_system(t,p,fault);
        F = pg_pengyue_excitation_force(t,p,meta);
        qdd = M\(F-C*qdNow-K*qNow);
        dx = [qdNow;qdd];
    end
end
