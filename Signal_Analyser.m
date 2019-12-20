function [state] = Signal_Analyser(state, x, y, z, sg_st, dpd_st, adapt_st, par)
% SIGNAL_ANALYSER(STATE, X, Y, Z, SG_STATE)

global Fs;

if isempty(fieldnames(state))
    par = RegisterParam(par, 'num_arch', 1);
    par = RegisterParam(par, 'arch_ind', 1);
    
    par = RegisterParam(par, 'sa', struct());
    par.sa = RegisterParam(par.sa, 'store_psd', 0);
    state.store_psd = par.sa.store_psd;
    
    Narch = par.num_arch;
    
    % initialization
    state.iter = 0;        % iteration
    state.iter_num = par.iter_num;
    
    % initialize variables for NMSE, EVM, ACPR
    state.mchpwr = zeros(Narch, state.iter_num, 1);
    state.nmse = zeros(Narch, state.iter_num, 1);
    state.evm = zeros(Narch, state.iter_num, 1);
    state.acpr1 = zeros(Narch, state.iter_num, 1);
    state.acpr2 = zeros(Narch, state.iter_num, 1);
    
    % reference values
    state.nmse_ref = zeros(state.iter_num, 1);
    state.evm_ref = zeros(state.iter_num, 1);
    state.acpr1_ref = zeros(state.iter_num, 1);
    state.acpr2_ref = zeros(state.iter_num, 1);
    
    % FFT spectra
    if state.store_psd
        state.psd_spect = zeros(Narch, state.iter_num, length(y));
    end
%     state.x = zeros(Narch, state.iter_num, length(x));
%     state.y = zeros(Narch, state.iter_num, length(y));
%     state.z = zeros(Narch, state.iter_num, length(z));
    
    % parameters for ACPR calculation
    state.acrp_channels = [1/2+0.1, 3/2+0.1, 3/2+0.2, 5/2+0.2];
    if isfield(sg_st, 'type')
        if strcmp(sg_st.type, 'OFDM')
            state.acpr_osf = sg_st.fft_len / (sg_st.num_rb * sg_st.rb_size);
        elseif strcmp(sg_st.type, 'QAM')
            state.acpr_osf = sg_st.oversample_fact;
        end
    else
        state.acpr_osf = 6;
    end

    % -------------------------------------------------------------------------
    % initialize figures
    % -------------------------------------------------------------------------
    
    if par.is_plot
        % Main Channel Power at PA output
        state.mchpwr_fn = 1;
        state.mchpwr_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.mchpwr_fn);
        if ~isempty(state.mchpwr_fh); clf(state.mchpwr_fh); end
        figure(state.mchpwr_fn);
        hold on;
        for i = 1:Narch
            state.mchpwr_ph{i} = plot(0,0);
        end
        hold off;
        title('Main channel power at PA output');
        xlabel('Iteration cycle');
        ylabel('Power (dBW)');
        axis([0 state.iter_num, -inf, inf]);

        % NMSE
        state.nmse_fn = 2;
        state.nmse_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.nmse_fn);
        if ~isempty(state.nmse_fh); clf(state.nmse_fh); end
        figure(state.nmse_fn);
        state.nmse_ref_ph = plot(0,0);
        hold on;
        for i = 1:Narch
            state.nmse_ph{i} = plot(0,0);
        end
        hold off;
        title('Normalised mean square exrror');
        xlabel('Iteration cycle');
        ylabel('NMSE (dB)');
        axis([0 state.iter_num, -inf, inf]);

        % EVM
        state.evm_fn = 3;
        state.evm_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.evm_fn);
        if ~isempty(state.evm_fh); clf(state.evm_fh); end
        figure(state.evm_fn);
        state.evm_ref_ph = plot(0,0);
        hold on;
        for i = 1:Narch
            state.evm_ph{i} = plot(0,0);        
        end
        hold off;
        title('Error Vector Magnitude');
        xlabel('Iteration cycle');
        ylabel('EVM (%)');
        axis([0 state.iter_num, -inf, inf]);

        % ACPR1
        state.acpr1_fn = 4;
        state.acpr1_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.acpr1_fn);
        if ~isempty(state.acpr1_fh); clf(state.acpr1_fh); end
        figure(state.acpr1_fn);
        state.acpr1_ref_ph = plot(0,0);
        hold on;
        for i = 1:Narch
            state.acpr1_ph{i} = plot(0,0);
        end        
        hold off;
        title('Adjacent channel power ratio - 1st adj. channel');
        xlabel('Iteration cycle');
        ylabel('ACPR1 (dB)');
        axis([0 state.iter_num, -inf, inf]);

        % ACPR2
        state.acpr2_fn = 5;
        state.acpr2_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.acpr2_fn);
        if ~isempty(state.acpr2_fh); clf(state.acpr2_fh); end
        figure(state.acpr2_fn);
        state.acpr2_ref_ph = plot(0,0);
        hold on;
        for i = 1:Narch
            state.acpr2_ph{i} = plot(0,0);
        end
        hold off;
        title('Adjacent channel power ratio - 2nd adj. channel');
        xlabel('Iteration cycle');
        ylabel('ACPR2 (dB)');
        axis([0 state.iter_num, -inf, inf]);
        
        % AM/AM characteristics
        state.amam_fn = 6;
        state.amam_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.amam_fn);
        if ~isempty(state.amam_fh); clf(state.amam_fh); end
        state.amam_fh = figure(state.amam_fn);
        switch Narch
            case 1
                sp_layout = [1 1];
            case 2
                sp_layout = [1 2];
            case {3, 4}
                sp_layout = [2 2];
            case {5,6}
                sp_layout = [3 2];
            otherwise
                error('Unsupported number of architectures');
        end
        for i = 1:Narch
            subplot(sp_layout(1), sp_layout(2), i);
            state.amam_pa_ph{i} = plot(0,0,'.');
            hold on;
            state.amam_lin_ph{i} = plot(0,0,'.');
            state.amam_dpd_ph{i} = plot(0,0,'.');
            hold off;
            title('AM/AM characteristics');
            xlabel('Input Magnitude (-)');
            ylabel('Output Magnitude (-)');
            axis([0 1, 0, 1]);
        end
        
        % constelation diagram
        state.iqdiag_fn = 7;
        state.iqdiag_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.iqdiag_fn);
        if ~isempty(state.iqdiag_fh); clf(state.iqdiag_fh); end
        state.iqdiag_fh = figure(state.iqdiag_fn);
        state.iqdiag_tx_ph = plot(0,0,'x');
        hold on;
        for i = 1:Narch
            state.iqdiag_rx_ph{i} = plot(0,0,'o');
        end
        hold off;
        title('IQ Constelation Diagram');
        xlabel('I (-)');
        ylabel('Q (-)');
        axis([-1.2 1.2 -1.2 1.2]);
        
        % frequency spectrum
        state.fspect_fn = 8;
        state.fspect_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.fspect_fn);
        if ~isempty(state.fspect_fh); clf(state.fspect_fh); end
        state.fspect_fh = figure(state.fspect_fn);
        state.fspect_ref_ph = plot(0,0);
        hold on;
        for i = 1:Narch
            state.fspect_ph{i} = plot(0,0);
        end
        hold off;
        title('Power Spectrum Density');
        xlabel('Frequency (Hz)');
        ylabel('PSD (dBm/Hz)');
    else
        for i = 1:Narch
            state.fspect_ph{i} = 0;
        end
    end
end

i = par.arch_ind;

% increment iteration only once per architectures
if i == 1
    state.iter = state.iter + 1;
end
iter = state.iter;


% state.x(i,iter,:) = x;
% state.y(i,iter,:) = y;
% state.z(i,iter,:) = z;

% calculate main channel power and ACPR
[acpr_mat, mchan_pwr] = ACPR(y, state.acpr_osf, state.acrp_channels, 0);
state.mchpwr(i,iter) = mchan_pwr;
state.acpr1(i,iter) = Avg_dB(acpr_mat(1,:), 2, 10);
state.acpr2(i,iter) = Avg_dB(acpr_mat(2,:), 2, 10);

% calculate NMSE
y_norm = lscov(y, z) * y;
state.nmse(i,iter) = NMSE(z, y_norm);

% calculate EVM
[~, ~, rx_symbs] = Signal_Demodulator(sg_st, y, z);
state.evm(i,iter) = EVM(rx_symbs,sg_st.symbs);

if par.is_plot
    % update figures
    set(state.mchpwr_ph{i}, 'XData', 1:iter, 'YData', state.mchpwr(i,1:iter));
    set(state.nmse_ph{i}, 'XData', 1:iter, 'YData', state.nmse(i,1:iter));
    set(state.evm_ph{i}, 'XData', 1:iter, 'YData', 100*state.evm(i,1:iter));
    set(state.acpr1_ph{i}, 'XData', 1:iter, 'YData', state.acpr1(i,1:iter));
    set(state.acpr2_ph{i}, 'XData', 1:iter, 'YData', state.acpr2(i,1:iter));

    % update AM/AM figures
    set(state.amam_pa_ph{i}, 'XData', abs(x)/dpd_st.in_scaler, 'YData', abs(y)/adapt_st.out_scaler);
    set(state.amam_lin_ph{i}, 'XData', abs(z)/dpd_st.in_scaler, 'YData', abs(y)/adapt_st.out_scaler);
    set(state.amam_dpd_ph{i}, 'XData', abs(z)/dpd_st.in_scaler, 'YData', abs(x)/dpd_st.in_scaler);
    
    % update constelation diagram
    set(state.iqdiag_rx_ph{i}, 'XData', real(rx_symbs), 'YData', imag(rx_symbs));
    
    if i == 1   % update reference values only once
        % find reference values for actual main channel power
        % for PA without DPD (consequently with power back-off)
        [state.nmse_ref(iter), state.evm_ref(iter), state.acpr1_ref(iter), state.acpr2_ref(iter), psd_ref] = ...
            Get_Reference_Values(par.sa_no_dpd, mchan_pwr);

        % update reference values
        set(state.nmse_ref_ph, 'XData', 1:iter, 'YData', state.nmse_ref(1:iter));
        set(state.evm_ref_ph, 'XData', 1:iter, 'YData', 100*state.evm_ref(1:iter));
        set(state.acpr1_ref_ph, 'XData', 1:iter, 'YData', state.acpr1_ref(1:iter));
        set(state.acpr2_ref_ph, 'XData', 1:iter, 'YData', state.acpr2_ref(1:iter));

        % reference spectrum
        set(state.fspect_ref_ph, 'XData', fftshift(fftfreq(length(psd_ref), Fs, true)), 'YData', fftshift(psd_ref));

        % reference constelation
        set(state.iqdiag_tx_ph, 'XData', real(reshape(sg_st.symbs,numel(sg_st.symbs),1)), 'YData', imag(reshape(sg_st.symbs,numel(sg_st.symbs),1)));
    end
end

% calculate and store psd if required
if par.is_plot || state.store_psd
    % update spectrum
    [psd_spect, ~] = PSD(y, state.fspect_ph{i});
    
    if state.store_psd
        state.psd(i,iter,:) = psd_spect;
    end
end