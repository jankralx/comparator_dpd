function [state, y] = PA_Model(state, x, par)
% PA_MODEL(STATE, X) returns output X of selected PA model

if isempty(fieldnames(state))
    load(par.coef_fname);
    state.model_name = model_name;
    state.pa_coefs = pa_coefs;
    state.in_scaler = in_scaler;
    state.out_scaler = out_scaler;
    state.mpars = mpars;
    
    % if measurement of a real Power amplifier is required, create a
    % measurement object
    if strcmp(state.model_name, 'SMU200A_FSVR')
        state.pa_meas = par.pa_meas;
        state.pa_meas.pa_exp_region = pa_exp_region;
        state.max_pa_peak = par.max_pa_peak;
    end

end

switch state.model_name
    case 'GMP'
        % scale the input
        x = x / state.in_scaler;
        % calculate model output
        y = GMP_Output(x, state.pa_coefs, state.mpars.Ka, ...
            state.mpars.La, state.mpars.Kb, state.mpars.Lb, state.mpars.Mb, ...
            state.mpars.Kc, state.mpars.Lc, state.mpars.Mc);
        % scale output
        y = y * state.out_scaler;
    case 'DDR2'
        % scale the input
        x = x / state.in_scaler;
        % calculate model output
        y = DDR2_Output(x, state.mpars.P, state.mpars.M, state.pa_coefs);
        % scale output
        y = y * state.out_scaler;
    case 'MP'
        % scale the input
        x = x / state.in_scaler;
        % calculate model output
        y = MP_Output(x, state.mpars.P, state.mpars.M, state.pa_coefs);
        % scale output
        y = y * state.out_scaler;
    case 'Saleh'
        y = state.alpha * x.^state.eta ./ (1 + state.beta*x.^2).^state.nu;
        
    case 'SMU200A_FSVR'
        if ~isvector(x) || nnz(isnan(x)) || nnz(abs(x) == Inf)
            y = [];
            return;
        end

        % limit input to specified dBm value
        R = 50;
        vlim = sqrt(10^(state.max_pa_peak/10)/1000*R);
        over = abs(x) > vlim;
        if nnz(over)
            % clip values greater than defined limit, preserve phase
            x(over) = vlim * exp(1j*angle(x(over)));
            warning('Clipping values greater than %.1f dBm', state.max_pa_peak);
        end

        y = state.pa_meas.GetOutput(x);
end
