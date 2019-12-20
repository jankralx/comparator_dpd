% PA coefficients selection
global Fs;
        
pa_st = struct();

switch par.pa.sel
    case 1
        % HonDri, 30mA quiscient current, fc = 1.6 GHz, GMP coefficients,
        % signal bandwidth 2.5 MHz
        pa_par.coef_fname = 'amplifiers/HonDri30_QAM_1G6Hz_B2M5Hz_hiPwr_GMP_coefs.mat';
        pa_st = PA_Model(pa_st, zeros(100,1), pa_par);
        par.sg.pep = -3;                   % signal generator - output power
        par.dpd.in_scaler = 0.39;
        par.dpd.out_scaler = 18;
        
    otherwise
        error('Unsupported PA model');
end
