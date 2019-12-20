function y = MP_Output(x, P, M, coef)

% MP_Output function returns signal y created by
%   transformation of signal x in memory polynomial model
%   with maximum order P and memory length M. Model is
%   characterised by coeficients coef.
%
%   Inputs:
%   =======
%   x  - signal to be transformed by MP model
%
%   P  - nonlinearity order
%
%   M  - memory length
%
%   coef - coefficients of the MP model
%
%   Returns:
%   ========
%
%   y - output of MP model - transformed signal x

% Authors: Jan Kral <kral.j@lit.cz>, Tomas Gotthans
% Date: 13.3.2018

%% Fast variant - precalculated values
y = zeros(size(x));
abs_x = abs(x);
for k = 1:P
    abs_x_k = abs_x.^(k-1);  % precalculate for faster operation
    temp_x = abs_x_k .* x;   % calculate non-shifted version of abs(x).^(k-1).*x
    for q = 0:M
        % coeficient index is calculated same way as by Tom Gotthans
        coef_ind = k + (q * P);
        shifted = circshift(temp_x, q);
        shifted(1:q) = 0;
        y = y + coef(coef_ind) * shifted;
    end
end

return;

%% Slow variant
y = zeros(size(x));
x_shifted = zeros(length(x),M);
for q = 0:M
    x_shifted(:,q+1) = circshift(x, q);
    x_shifted(1:q,q+1) = 0;
end

for k = 1:P
    for q = 0:M
        coef_ind = k + (q * P);
        y = y + coef(coef_ind) * ( ...
            x_shifted(:,q+1)*abs(x_shifted(:,q+1))^(k-1));
    end
end
