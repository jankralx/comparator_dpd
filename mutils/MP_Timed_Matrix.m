function [U, xout] = MP_Timed_Matrix(x, P, M, ts, vicinity)
% MP_Timed_Matrix function returns matrix U for MP model calculation created
%   from given signal x. The matrix is created with rows for samples at given times
%   ts. Matrix is created for Memory Polynomial model
%   of nonlinearity order P and memory length M.
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the MP matrix U. Time of the x first sample 
%        is extected to be 1.
%
%   P  - DPD nonlinearity order
%
%   M  - DPD memory length
%
%   ts - vector of times which rows of the output matrix are created at.
%        The unit of ts is sampling period of given signal x.
%
%   [vicinity] - if time shifting signal in FFT domain it takes enormous
%                amount of time when the signal length is high. Vicinity is
%                provided to reduce calculation time of FFT/IFFT. With
%                given vicinity the signal shifting is done only using a
%                needed region plus minus given vicinity. If the vicinity
%                is not provided, the whole signal is going to be shifted.
%
%   Returns:
%   ========
%
%   U  - MP matrix created from signal x


% Authors: Jan Kral <kral.j@lit.cz>
% Date: 17.5.2018

coef_num = P*(M+1);
U = zeros(length(ts),coef_num);
xout = zeros(length(ts),1);

% if vicinity is not given, make it infinity so the whole signal is taken
% for shifting in FFT domain
if nargin < 5
    vicinity = inf;
end

for n = 1:length(ts)

    % shift input signal x in time such that first sample has time ts(n)
    
    ts_n = floor(ts(n));
    ts_fract = ts(n) - ts_n;
   
    % limit the signal only to the given vicinity
    start_i = ts_n - vicinity;
    stop_i = ts_n + vicinity;
    if start_i < 1
        start_i = 1;
    end
    if stop_i > length(x)
        stop_i = length(x);
    end
    x_sel = x(start_i:stop_i);      % make selection of vicinity
    ts_n = ts_n - start_i + 1;      % update ts_n to point into selected samples only    
    
    x_shifted = SigDelayFFT(x_sel, -ts_fract);
    xout(n) = x_shifted(ts_n);
    x_shifted = [zeros(M,1); x_shifted(1:ts_n)];
    
    for k = 1:P
        for q = 0:M
            coef_ind = k + (q * P);
            U(n,coef_ind) = x_shifted(end-q)*abs(x_shifted(end-q))^(k-1);
        end
    end

end

