function [state, rx_bin, rx_symbs] = Signal_Demodulator(sg_st, y, z)
% SIGNAL_DEMODULATOR(SG_STATE) 

switch sg_st.typen
    case 0
        y = y * lscov(y,z);
        [rx_bin, rx_symbs] = OFDM_Demodulator(y, sg_st, sg_st.symbs);   % demodulation
        sg_st.symbs = reshape(sg_st.symbs,numel(sg_st.symbs),1);
        rx_symbs = reshape(rx_symbs,numel(rx_symbs),1);
    case 1
        y = y * lscov(y,z);
        [rx_bin, rx_symbs, ~] = QAM_Demodulator(y, sg_st);
end

% normalize rx symbols to generated symbols by terms of MSE
rx_symbs = rx_symbs * lscov(rx_symbs, sg_st.symbs);

state = 0;