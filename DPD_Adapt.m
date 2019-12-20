function [state] = DPD_Adapt(state, x, y, z, par)
% DPD_ADAPT(STATE, PARS) 
% x - PA input
% y - PA output
% z - desired PA output

%% DPD_Adapt initilization
if isempty(fieldnames(state))
    % register parameters with default values if no values are provided
    par = RegisterParam(par, 'dpd');
    par.dpd = RegisterParam(par.dpd, 'model', 'MP');
    par.dpd = RegisterParam(par.dpd, 'P', 5);
    par.dpd = RegisterParam(par.dpd, 'M', 1);
    par.dpd = RegisterParam(par.dpd, 'def_coefs', DPD_Default_Coefs(0.5, par));
    par.dpd = RegisterParam(par.dpd, 'mu0', 0.05);
    par.dpd = RegisterParam(par.dpd, 'mu_change_after', 30);
    par.dpd = RegisterParam(par.dpd, 'mu_zero_after', 80);
    par.dpd = RegisterParam(par.dpd, 'mu_min', par.dpd.mu0);
    par.dpd = RegisterParam(par.dpd, 'in_scaler', 0.5);
    par.dpd = RegisterParam(par.dpd, 'out_scaler', 2.3);
    
    state.iter = 0;
    
    switch par.dpd.algorithm
        case 'ILA'
            state.dpd_algorithm = 0;
        case 'DLA'
            state.dpd_algorithm = 1;
        case 'DLA-COMP'
            state.dpd_algorithm = 2;
        otherwise
            error('%s is unsupported DPD algorithm', par.dpd.algorithm);
    end        
    
    switch par.dpd.model
        case 'MP'
            state.dpd_model = 0;
        case 'DDR2'
            state.dpd_model = 1;
        case 'MP-Timed'
            state.dpd_model = 2;
        case 'DDR2-Timed'
            state.dpd_model = 3;
        otherwise
            error('%s is unsupported DPD model', par.dpd.model);
    end
    state.P = par.dpd.P;
    state.M = par.dpd.M;
    state.coefs = par.dpd.def_coefs;

    state.mu0 = par.dpd.mu0;
    state.mu_change_after = par.dpd.mu_change_after;
    state.mu_zero_after = par.dpd.mu_zero_after;
    state.mu_min = par.dpd.mu_min;
    state.mu_gamma = (state.mu0/state.mu_min).^...
        (1/(par.avg_res_after-state.mu_change_after));

    state.in_scaler = par.dpd.in_scaler;
    state.out_scaler = par.dpd.out_scaler;
end

%% DPD_Adapt execution

state.iter = state.iter + 1;        % update iteration
iter = state.iter;

% calculate new mu
state.mu = state.mu0 / state.mu_gamma.^(...
        (iter - state.mu_change_after)*...
        (iter > state.mu_change_after));
if state.mu < state.mu_min
    state.mu = state.mu_min;
end
if iter > state.mu_zero_after 
    state.mu = 0;
end

% scale signals for DPD calculation
norm_one = 0;
if norm_one
    z = z / max(abs(z));
    x = x / max(abs(x));
    y = y / max(abs(y));
else
    z = z / state.in_scaler;
    x = x / state.in_scaler;
    y = y / state.out_scaler;
end

state.z_max = max(abs(z));
state.y_max = max(abs(y));

% select signal for Matrix creation
switch state.dpd_algorithm
    case 0 % ILA
        u = y;
    case 1  % DLA
        u = z;
    case 2  % DLA with comparator
        u = z;
        [y_times, y_sel] = Comparator_Model(y);
end

% create DPD matrix
switch state.dpd_model
    case 0  % MP
        U = MP_Matrix(u, state.P, state.M);
    case 1  % DDR2
        U = DDR2_Matrix(u, state.P, state.M);
    case 2  % MP-Timed
        [U, z_sel] = MP_Timed_Matrix(u, state.P, state.M, y_times, 1024);
    case 3  % DDR2-Timed
        [U, z_sel] = DDR2_Timed_Matrix(u, state.P, state.M, y_times, 1024);
end

% update DPD coefficients
switch state.dpd_algorithm
    case 0  % ILA
        state.coefs = ((U'*U+0.001*eye(length(state.coefs)))\(U'*x));
        % without lambda = 0.001 ILA does not converge
%         state.coefs = ((U'*U)\(U'*x));
    case 1  % DLA
        % use Gauss-Newton method to update coefficients
        state.coefs = state.coefs - state.mu*((U'*U)\(U'*(y-z)));
    case 2  % DLA with comparator
        M = [real(U) -imag(U)];
        b = [real(state.coefs); imag(state.coefs)];
        b = b - state.mu*((M'*M)\(M'*(y_sel-real(z_sel))));
        state.coefs = b(1:(size(b,1)/2)) + 1j*b((size(b,1)/2)+1:end);
        
end






