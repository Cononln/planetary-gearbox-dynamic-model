function out = pg_extract_analytic_harmonics(signal,fs,centerFrequencyHz,halfBandwidthHz)
%PG_EXTRACT_ANALYTIC_HARMONICS Tacholess fixed-band analytic components.
% The input is an N-by-S real signal.  No shaft angle or pulse train is
% used.  Each requested positive-frequency band is converted to a complex
% analytic component by a tapered FFT mask.

if isvector(signal)
    signal = signal(:);
end
centerFrequencyHz = centerFrequencyHz(:).';
if isscalar(halfBandwidthHz)
    halfBandwidthHz = repmat(halfBandwidthHz,size(centerFrequencyHz));
else
    halfBandwidthHz = halfBandwidthHz(:).';
end
if numel(centerFrequencyHz) ~= numel(halfBandwidthHz)
    error('pg_extract_analytic_harmonics:BandSize', ...
        'One half bandwidth is required for each center frequency.');
end

[nSample,nSensor] = size(signal);
nBand = numel(centerFrequencyHz);
signal = signal-mean(signal,1);
spectrum = fft(signal,[],1);
frequency = (0:nSample-1)'*(fs/nSample);
component = complex(zeros(nSample,nSensor,nBand));
maskBank = zeros(nSample,nBand);

for ib = 1:nBand
    center = centerFrequencyHz(ib);
    outer = halfBandwidthHz(ib);
    inner = 0.80*outer;
    distance = abs(frequency-center);
    mask = zeros(nSample,1);
    mask(distance<=inner) = 1;
    transition = distance>inner & distance<outer;
    mask(transition) = 0.5*(1+cos(pi*(distance(transition)-inner)/ ...
        max(outer-inner,eps)));
    % A real sinusoid places half its amplitude in the positive-frequency
    % half-plane.  Doubling this band gives its analytic representation.
    mask = 2*mask;
    maskBank(:,ib) = mask;
    component(:,:,ib) = ifft(spectrum.*mask,[],1);
end

out.component = component;
out.fs = fs;
out.centerFrequencyHz = centerFrequencyHz;
out.halfBandwidthHz = halfBandwidthHz;
out.frequency = frequency;
out.mask = maskBank;
out.definition = ['Positive-frequency analytic bands extracted from ', ...
    'the vibration channels only; no encoder samples are used.'];
end
