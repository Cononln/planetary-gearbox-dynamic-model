function setup_paths
%SETUP_PATHS Add this project to the MATLAB path.

rootDir = fileparts(mfilename('fullpath'));
% Keep the bootstrap function itself reachable when MATLAB's unit-test
% runner temporarily changes the current directory.
addpath(rootDir);
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'scripts'));
addpath(fullfile(rootDir, 'tests'));
% Generated results are intentionally absent from a clean Git clone.
resultDir = fullfile(rootDir, 'results');
if ~isfolder(resultDir)
    mkdir(resultDir);
end
end
