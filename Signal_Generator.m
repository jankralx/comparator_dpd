function [state, z] = Signal_Generator(state, par)
% SIGNAL_GENERATOR(SG_STATE, PARS) 

if isempty(fieldnames(state))
    
    par = RegisterParam(par, 'sg');
    par.sg = RegisterParam(par.sg, 'type', 'OFDM');
    par.sg = RegisterParam(par.sg, 'pep', 0);
    
    switch par.sg.type
        case 'OFDM'
            par.sg.typen = 0;
            par.sg = RegisterParam(par.sg, 'num_frames', 12);
            par.sg = RegisterParam(par.sg, 'imod');
            par.sg.imod = RegisterParam(par.sg.imod, 'M', 64);
            par.sg = RegisterParam(par.sg, 'num_rb', 57);
            par.sg = RegisterParam(par.sg, 'rb_size', 6);
            par.sg = RegisterParam(par.sg, 'fft_len', 4096);
            par.sg = RegisterParam(par.sg, 'filt');
            par.sg.filt = RegisterParam(par.sg.filt, 'L', par.sg.fft_len/2+1);
            par.sg.filt = RegisterParam(par.sg.filt, 'ena', true);
            
            state = par.sg;
            state.out_level = par.sg.pep - 11;  % [dBm] = pep - max(papr)

        case 'QAM'
            par.sg.typen = 1;
            
            % QAM Modulator
            par.sg = RegisterParam(par.sg, 'M', 64);
            par.sg = RegisterParam(par.sg, 'oversample_fact', 12);
            par.sg = RegisterParam(par.sg, 'filt', 'sqrt');
            par.sg = RegisterParam(par.sg, 'num_bits', 2*1024*10*log2(par.sg.M));
            
            state = par.sg;
            state.out_level = par.sg.pep - 3;  % [dBm] = pep - max(papr)

        otherwise
            error('%s is unsupported modulation', par.sg.type);
    end
    
    state.papr = [];
    state.par = par;
    
end

switch state.typen
    case 0
        bits_per_frame = log2(state.imod.M)*state.rb_size*state.num_rb;
        papr = 0;
        watchdog = 0;
        while (papr < 10 || papr > 11) && watchdog < 10
            tx_bin = randi([0, 1], state.num_frames*bits_per_frame, 1);           % generate random binary data
            [z, state.symbs] = OFDM_Modulator(tx_bin, state);
            papr = PAPR(z);
            watchdog = watchdog + 1;
        end
        if (papr < 10 || papr > 11)
            warning('Could not generate signal with set PAPR');
        end
    case 1
        papr = 0;
        watchdog = 0;
        while (papr < 6.5 || papr > 7.5) && watchdog < 10
            binin = rand(state.num_bits, 1) > 0.5;
            [z, state.symbs] = QAM_Modulator(binin, state);
            papr = PAPR(z);
            watchdog = watchdog + 1;
        end
        if (papr < 6.5 || papr > 7.5)
            warning('Could not generate signal with set PAPR');
        end

end

% adjust generator output power
% par.out_level = 6 - 11; % [dBm] = pep - max(papr) = 6 - 11 

R = 50;
z_pow = 10*log10(mean(abs(z).^2)/R)+30; % [dBm]

% adjust output power
z = z * 10^((state.out_level-z_pow)/20);
