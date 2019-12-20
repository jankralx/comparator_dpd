addpath('mutils');

clear par;

global Fs;
Fs = 30e6;

par.is_plot = 1;

% select PA model
par.pa.sel = 1;

% set number of iterations
par.iter_num = 50;
par.avg_res_after = 60;

% signal generator
par.sg.type = 'QAM';

% signal analyser
par.sa.store_psd = 1;

% parameters of DPD
par.dpd.algorithm = 'DLA';
par.dpd.model = 'MP';
par.dpd.P = 7;
par.dpd.M = 3;
par.dpd.def_coefs = DPD_Default_Coefs(0.5, par);

Prepare_PA_Model;               % prepare PA model
par = Run_NoDPD_Sims(par);      % get simulation results for amplifier without DPD

% number of different architectures
par.num_arch = 2;
par = repmat(par, par.num_arch, 1);

% set different architectures
i = 1;
par(1).arch_ind = i;
par(1).arch_name = 'DLA';
i = i + 1;
par(i).arch_ind = i;
par(i).arch_name = 'DLA-COMP';
par(i).dpd.algorithm = 'DLA-COMP';
par(i).dpd.model = 'MP-Timed';


% initialise state variables
sg_st = struct();
sa_st = struct();
cp_st = struct();
dpd_st = cell(length(par),1);
adapt_st = cell(length(par),1);
for i = 1:length(par)
    dpd_st{i} = struct();
    adapt_st{i} = struct();
end

tic;

% signal generator
[sg_st, z] = Signal_Generator(sg_st, par(1));
    
for iter = 1:par(1).iter_num
    % signal generator
    [sg_st, z] = Signal_Generator(sg_st, par(1));

    % go through all architectures
    for i = 1:length(par)
        % DPD
        [dpd_st{i}, x] = DPD(dpd_st{i}, z, par(i));

        % PA model
        [pa_st, y] = PA_Model(pa_st, x);

        % DPD adaptation
        adapt_st{i} = DPD_Adapt(adapt_st{i}, x, y, z, par(i));
        dpd_st{i}.coefs = adapt_st{i}.coefs;              % copy coefficients to DPD

        % signal analyser
        sa_st = Signal_Analyser(sa_st, x, y, z, sg_st, dpd_st{i}, adapt_st{i}, par(i));
    end
    
    % control panel
    [cp_st, dpd_st, adapt_st] = Control_Panel(cp_st, par, dpd_st, adapt_st, sa_st);

    drawnow();
    
end

toc;