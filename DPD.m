function [state, x] = DPD(state, z, par)
% DPD(STATE, Z, PAR)

if isempty(fieldnames(state))
    % register parameters with default values if no values are provided
    par = RegisterParam(par, 'dpd');
    par.dpd = RegisterParam(par.dpd, 'model', 'MP');
    par.dpd = RegisterParam(par.dpd, 'P', 5);
    par.dpd = RegisterParam(par.dpd, 'M', 1);
    par.dpd = RegisterParam(par.dpd, 'def_coefs', DPD_Default_Coefs(0.5, par));
    par.dpd = RegisterParam(par.dpd, 'in_scaler', 0.5);
    
    switch par.dpd.model
        case {'MP', 'MP-Timed'}
            state.dpd_model = 0;
        case {'DDR2', 'DDR2-Timed'}
            state.dpd_model = 1;
        otherwise
            error('%s is unsupported DPD model', par.dpd.model);
    end            
    state.P = par.dpd.P;
    state.M = par.dpd.M;
    state.coefs = par.dpd.def_coefs;
    state.in_scaler = 0.5;
    state.out_scaler = 1;
end
    
% adjust magnitude for DPD calculation
norm_one = 0;
if norm_one % normalize to one
    in_scaler = max(abs(z));
    z = z / in_scaler;
else
    z = z / state.in_scaler;
end

state.z_max = max(abs(z));

switch state.dpd_model
    case 0
        x = MP_Output(z, state.P, state.M, state.coefs);
    case 1
        x = DDR2_Output(z, state.P, state.M, state.coefs);
end


state.x_max = max(abs(x));

% recover magnitude
if norm_one
    x = x * in_scaler * state.out_scaler;
else
    x = x * state.in_scaler * state.out_scaler;
end