function [f, a] = dte_spectrum(csvFile, colDyn)
%DTE_SPECTRUM single-sided amplitude spectrum from DTE csv
%   csv columns: 1 t_s, 2 xsp1dyn_um, 3 xsp1raw_um
    if nargin < 2, colDyn = 2; end
    d = readmatrix(csvFile);
    t = d(:, 1);
    x = d(:, colDyn);
    dtm = median(diff(t));
    x = x - mean(x);
    N = numel(x);
    X = fft(x)/N;
    nH = floor(N/2) + 1;
    a = 2*abs(X(1:nH));
    f = (0:nH-1)'/(N*dtm);
end
