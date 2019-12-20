function [papr] = PAPR(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

papr = 10*log10(max(abs(x)).^2/mean(abs(x).^2));

end

