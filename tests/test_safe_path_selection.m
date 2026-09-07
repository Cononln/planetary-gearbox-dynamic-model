function tests = test_safe_path_selection
tests = functiontests(localfunctions);
end

function testReliablePairIsSelected(testCase)
addProjectPaths;
rng(7);
fs = 2000;
t = (0:15999)'/fs;
carrierPhase = 2*pi*2*t;
meshPhase = 84*carrierPhase;
harmonicOrder = [1,2];
component = complex(zeros(numel(t),3,2));
pathPhase = [0.45*sin(3*carrierPhase), ...
    -0.45*sin(3*carrierPhase),zeros(size(carrierPhase))];
for ih = 1:2
    for is = 1:3
        phaseNoise = 0.025*randn(size(t));
        if is==3
            phaseNoise = 1.1*randn(size(t));
        end
        component(:,is,ih) = exp(1i*( ...
            harmonicOrder(ih)*meshPhase+pathPhase(:,is)+phaseNoise));
    end
end
options.cvStride = 5;
options.minCvGain = 0.01;
result = pg_select_safe_path_correction(component,carrierPhase, ...
    -6:6,ones(1,13),options);
verifyEqual(testCase,result.selectedCandidate,"pair_1_2");
verifyGreaterThan(testCase,result.fullCoherenceGain,0.10);
verifyLessThan(testCase,result.amplitudePreservationError,1e-12);
end

function testNoPathFallback(testCase)
addProjectPaths;
rng(11);
fs = 2000;
t = (0:11999)'/fs;
carrierPhase = 2*pi*2*t;
meshPhase = 84*carrierPhase;
component = complex(zeros(numel(t),3,2));
constantPhase = [0,0.7,-0.4];
for ih = 1:2
    for is = 1:3
        component(:,is,ih) = exp(1i*(ih*meshPhase+ ...
            constantPhase(is)+0.02*randn(size(t))));
    end
end
result = pg_select_safe_path_correction(component,carrierPhase, ...
    -6:6,ones(1,13),struct('cvStride',5));
verifyEqual(testCase,result.selectedCandidate,"no_path");
verifyEqual(testCase,result.correctedComponent,component);
verifyEqual(testCase,result.fullCoherenceGain,0,'AbsTol',eps);
end

function addProjectPaths
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'src'));
end
