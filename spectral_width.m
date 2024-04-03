function [spect_width] = spectral_width(param,X)
%Estimate the spectral linewidth by integration of the phase noise power spectral density
N = length(X);
delta_f = param.Fs/N;
freq = [0 : N-1]*delta_f;
spect_width = sum(X(param.int_start : round(N/2)))*delta_f^2;
%spect_width = trapz(freq(param.int_start : round(N/2)),X(param.int_start : round(N/2)))*delta_f;
end

