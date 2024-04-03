function [phase] = phase_estimation_VV(param,X)
% phase estimation from a beatnote signal	

% X: the detected single quadrature signal
% param.PE_type: phase estimation type: choose 'I' if we already know the
% beating frequency of the phase. It basically shifts the beating frequency
% to the origin, apply a lowpass filter (using the movmean) to the signal
% and then it extract the phase by taking the complex argument of it
% choose 'II' if we don't know the beating frequency and it needs to be
% estimated: it first extract the phase, then using a polynomial fitting
% tool it locates the beating frequency. The final step is the removal of
% the linear trend (beating frequency) from the phase noise.
% param.Fs: sampling frequency
% param.Ts: sampling time
% param.beat_freq: beating frequency
% param.L_av: smoothing parameter for the lowpass filter in param.PE_type =
% 'I'
% choose 'III' if the removal of the beating frequency is due to the
% detrend.m Matlab algorithm

% phase: the output phase noise


if iscolumn(X)
    X = X';
end
L = length(X);
switch param.PE_type
    case 'I'
        k = 1 : L;
        Y = hilbert(X(1:L)).*exp(-1i*(2*pi*param.beat_freq*param.Ts*k));
        Y = movmean(Y,param.L_av);
        phase = unwrap(angle(Y));
        
    case 'II'
        X = hilbert(X(1:L));
        t = (1 : length(X))*(1/param.Fs);
        phase = unwrap(angle(X));
        P = polyfit(t, phase, 1);
        phase = phase - P(1)*t;
        
    case 'III'
        X = hilbert(X(1:L));
        phase = detrend(unwrap(angle(X)));
end
end