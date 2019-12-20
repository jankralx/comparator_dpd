function [psd_dBm_Hz, freq_axis] = PSD(sig, hplot, fs, decim_fact)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R = 50;

global Fs;
if nargin < 3
    fs = Fs;
end
if nargin < 4
    decim_fact = 1;
end

R = 50;

spect = abs(fft(sig));
N = length(spect);
spect_1N = spect/N;
fft_W = spect_1N.^2 / R;
psd_W_Hz = fft_W / (fs / N);

% if decimation and filtration is needed it belongs to psd_W_Hz
if decim_fact > 1
    psd_W_Hz_decim = decimate(psd_W_Hz, decim_fact, 'fir');
else
    psd_W_Hz_decim = psd_W_Hz;
end

psd_dBm_Hz = 10*log10(psd_W_Hz_decim) + 30;     % +30 because of dBW to dBm


freq_axis = fftfreq(length(psd_dBm_Hz), fs, true);

if isa(hplot, 'matlab.graphics.chart.primitive.Line')
    set(hplot, 'XData', fftshift(freq_axis), 'YData',fftshift(psd_dBm_Hz));
end

