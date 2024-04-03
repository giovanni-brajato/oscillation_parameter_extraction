function Y = comb_bp_filter(param,X)
%COMB_BP_FILTER 
% bandpass brickwall filter of frequency comb
f_ext = (0:param.L_trace-1).'./param.Ts/param.L_trace; % extended frequency vector

% prepare brickwall filter - lowpass
SQ_filter = zeros(param.L_trace,1); % in frequency domain
filter_taps = sum(f_ext < param.comb.bandpass_bw);
SQ_filter(1:floor(filter_taps)) = 1;

SQ_filter = [SQ_filter(1:param.L_trace/2);flipud(SQ_filter(1:param.L_trace/2))];
% force to even
SQ_filter = real(fft(real(ifft(SQ_filter))));

% brickwall filter - bandpass


[~,fc_index ] = min((abs(f_ext - param.comb.bandpass_cf)));
SQ_filter_bandpass = circshift(SQ_filter,fc_index);
SQ_filter_bandpass = [SQ_filter_bandpass(1:param.L_trace/2);flipud(SQ_filter_bandpass(1:param.L_trace/2))];
SQ_filter_bandpass = (fft(real(ifft(SQ_filter_bandpass)))); % force to real
Y = (ifft(fft(X).*SQ_filter_bandpass));

end

