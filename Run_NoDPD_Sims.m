function [par] = Run_NoDPD_Sims(par)
% RUN_NODPD_SIMS(PAR) runs simulations for selected PA amplifier and saves
% results into matlab file

% load results for PA without DPD if they exist
par.no_dpd_file = 'pa_no_dpd_cache.mat';
if exist(par.no_dpd_file, 'file')
    load(par.no_dpd_file);
%     if isequal(par, par_no_dpd)
        fprintf('Loading no DPD results from %s\n', par.no_dpd_file);
        par.sa_no_dpd = sa_no_dpd;
        return;
%     end
end

% we will use DPD to change output power, however, only one coefficient and
% therefore linear operation are permitted
par_sim = par;
par_sim.dpd.P = 1;
par_sim.dpd.M = 0;
par_sim.is_plot = 0;
par_sim.arch_ind = 1;

Prepare_PA_Model;           % prepare PA model

% initialise state variables
sg_st = struct();       
dpd_st = struct();
adapt_st = struct();
sa_st = struct();

fprintf('Running simulations for the selected PA without DPD\n');

% signal generator
[sg_st, z] = Signal_Generator(sg_st, par_sim);

% sweep over 40 dB in 200 points
pa_in_sweep = linspace(pa_st.in_scaler/100, pa_st.in_scaler, 100);

par_sim.iter_num = length(pa_in_sweep);
par_sim.sa.store_psd = 1;

% sweep the signal power
for iter = 1:par_sim.iter_num
    % adjust amplitude of the PA input signal
    x = z / max(abs(z)) * pa_in_sweep(iter);
    
    % PA model
    [pa_st, y] = PA_Model(pa_st, x);

    % signal analyser
    sa_st = Signal_Analyser(sa_st, x, y, z, sg_st, dpd_st, adapt_st, par_sim);
    drawnow();
end

fprintf('Saving no DPD results into %s\n', par.no_dpd_file);

sa_no_dpd.mchpwr = sa_st.mchpwr;
sa_no_dpd.nmse = sa_st.nmse;
sa_no_dpd.evm = sa_st.evm;
sa_no_dpd.acpr1 = sa_st.acpr1;
sa_no_dpd.acpr2 = sa_st.acpr2;
sa_no_dpd.psd = sa_st.psd;

par_no_dpd = par;
save(par.no_dpd_file, 'sa_no_dpd', 'par_no_dpd');

par.sa_no_dpd = sa_no_dpd;
