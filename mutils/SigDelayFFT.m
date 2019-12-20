function [ sig_d ] = SigDelayFFT( sig, delay )
%sig_delay_fft Delays signal through frequency domain using FFT.
%   sig - is a signal to be delayed
%   delay - is a required fractional delay in number of samples.
%           It can be negative and positive.

% check if the delay is non-zero, otherwise return the original signal
if delay == 0
    sig_d = sig;
    return
end

% sig has to be column vector for calculation so convert it to column
% vector
transposed = 0;
if size(sig, 1) ~= 1
    if size(sig,2) ~= 1
        % this is an error, this function cannot calculate with matrices
        error('Input parameter sig has to be 1-D vector.');
    end
    sig = sig.';
    transposed = 1;
end

% check if the input signal is complex
if isreal(sig)
    is_real = 1;
else
    is_real = 0;
end


Sig_fft = fft(sig);     % convert signal into frequency domain

% delay the signal in frequency domain
Sig_fft = ifftshift(fftshift(Sig_fft) ...
    .* exp(-1j*2*pi*delay/length(Sig_fft)*...
    (-length(Sig_fft)/2+1:length(Sig_fft)/2)));

sig_d = ifft(Sig_fft);  % inverse FFT to get the original signal


% if input was real, make output real as well
if is_real
    sig_d = real(sig_d);
end

% if input was transposed, transpose output as well
if transposed
    sig_d = sig_d.';
end

end

