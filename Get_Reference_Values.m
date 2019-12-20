function [nmse, evm, acpr1, acpr2, psd] = Get_Reference_Values(ref, mchpwr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% find main channel power in reference results
sel_ind = find(mchpwr < ref.mchpwr, 1);
if isempty(sel_ind)
    sel_ind = length(ref.mchpwr);
    warning('Requesting reference values for mchpwr = %.2f dBW, reference available only up to %.2f dBW', ...
        mchpwr, ref.mchpwr(end));
end

nmse = ref.nmse(1,sel_ind);
evm = ref.evm(1,sel_ind);
acpr1 = ref.acpr1(1,sel_ind);
acpr2 = ref.acpr2(1,sel_ind);
psd = ref.psd(1,sel_ind,:);