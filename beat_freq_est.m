function [freq,val,pwr_X_fft,beat_freq] = beat_freq_est(param,X)
%BEAT_FREQ_EST: estimation of the beat (main) frequency of a signal
%   The way it works is very simple: it takes the FFT of the signal X and
%   then by inspection of the spectrum it locates the main peak of the
%   signal. 

% X : the signal (one quadrature)
% param: parameters of the signal, the only necessary one is param.Fs, that
% is the sampling frequency

% freq: vector of frequency (negative and positives) which correspond to
% the same frequencies one would obtain by calculating the double side
% power spectral density
% val: the power of the beat frequency peak
% pwr_X_fft: unnormalized double side power spectral density of the signal
% beat_freq: beat frequency in Hertz

if  ~iscolumn(X)
    X = X';  % the time direction is along the first dimension
end

X = X - mean(X); % removing DC component (if not done before). in this way we are sure the following max function is going to locate the beating frequency

L = length(X); 
df = param.Fs./L; % frequency bin resolution

X_fft=fft(X)./L;
pwr_X_fft = abs(X_fft).^2;
freq = (0 : L - 1)*df;
[val,idx] = max(pwr_X_fft); % the peak is located using the max function
beat_freq = freq(idx);

pwr_X_fft = fftshift(pwr_X_fft);
freq = [-L/2 : L/2-1]*df;
end

