function [Y] = phase_noise_Wiener(var_pn,L)
%PHASE_NOISE_WIENER generate Wiener (random walk, Lorentzian) phase noise
%using two input parameters:
% param.L: length of the sequence
% param.var_pn: variance of the phase noise
Y = cumsum(sqrt(var_pn)*randn(L,1));
end

