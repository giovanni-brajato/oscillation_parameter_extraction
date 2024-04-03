function [spect_width,g_t_norm,tau] = coherence_prop(param,X)
%COHERENCE_PROP Calculate the spectral width of the phase noise
%autocorrelation function. Returns also the normalized autocorrelation

% X: the phase noise time sequence
% param.Ts: sampling time
% param.int_method: integration method to calculate the integral of the
% squared module of the autocorrelation function

% spect_width: spectral width of the phase noise ?
% g_t_norm: normalized autocorrelation of the phase noise
% tau: time delay of g_t_norm. Both should be plotted together, with tau on
% x-axis and g_t_norm on the y_axis

if  ~iscolumn(X)
    X = X';
end
E = cos(X) + 1i*sin(X); % compose the complex signal E with the provided phase noise X
g_t = xcorr(E,E,'coeff'); % calculate the complex autocorrelation
[~,idx]=max(g_t);
g_t_norm = g_t./g_t(idx); % normalize the autocorrelation (it is still complex)
N = length(g_t_norm); % calculate the dealy vector matching in size the autocorrelation vector
tau = [-N/2:N/2-1]*(param.Ts);
g_t_norm_abs = abs(g_t_norm).^2; % take the square of the autocorrelation magnitude
switch param.int_method % integrate the above to get an estimate of the coherence time
    case 'trapz'
        tau_c = trapz(tau,g_t_norm_abs);
    case 'simple'
        tau_c = sum(g_t_norm(1:end-1))*(1/param.Fs);
    otherwise
        disp('Integration method not defined!')
end
spect_width = 1/tau_c; % the spectral width is the inverse of the coherence time
end

