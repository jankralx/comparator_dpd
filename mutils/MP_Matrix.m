function U = MP_Matrix(x, P, M)

% MP_Matrix function returns matrix U created from the given signal x.
%   Matrix is created for Memory Polynomial model of nonlinearity order P
%   and memory length M.
%
%   Inputs:
%   =======
%   x  - signal to be ordered into the MP matrix U
%
%   P  - nonlinearity order
%
%   M  - memory length
%
%   Returns:
%   ========
%
%   U  - MP matrix created from signal x

% Authors: Jan Kral <kral.j@lit.cz>
% Date: 13.3.2018


% the MP is defined in the article
% 
% Equation 


coef_num = P*(M+1);
U = zeros(length(x),coef_num);

%% Fast variant - precalculated values

abs_x = abs(x);
for k = 1:P
    abs_x_k = abs_x.^(k-1);  % precalculate for faster operation
    temp_x = abs_x_k .* x;   % calculate non-shifted version of abs(x).^(k-1).*x
    for q = 0:M
        % coeficient index is calculated same way as by Tom Gotthans
        coef_ind = k + (q * P);
        U(:, coef_ind) = circshift(temp_x, q);
        U(1:q, coef_ind) = 0;
    end
end

return;

%% Slow variant
x_shifted = zeros(length(x),M);
for q = 0:M
    x_shifted(:,q+1) = circshift(x, q);
    x_shifted(1:q,q+1) = 0;
end

for k = 1:P
    for q = 0:M
        coef_ind = k + (q * P);
        U(:,coef_ind) = x_shifted(:,q+1)*abs(x_shifted(:,q+1))^(k-1);
    end
end


