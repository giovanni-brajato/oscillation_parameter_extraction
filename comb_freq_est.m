function freq_vector = comb_freq_est(param,X)
%COMB_FREQ_EST estimate the frequency of the comb, by locating the peaks of the power spectra of the comb
f_ext = (0:param.L_trace-1).'./param.Ts/param.L_trace; % extended frequency vector
f_res = f_ext(2);
% peak detection
pw_X = abs(fft(X));

[power_peaks,locations] = findpeaks(pw_X,'MinPeakDistance',param.comb.delta_f./f_res);
freq_locations = f_ext(locations);

frequency_candidates = freq_locations(freq_locations > param.comb.f_low);
power_peaks_candidates = power_peaks(freq_locations > param.comb.f_low);
freq_vector = frequency_candidates(frequency_candidates<param.comb.f_up);
power_vector = power_peaks_candidates(frequency_candidates<param.comb.f_up);

peak_level_dB = -20*log10(power_vector./max(power_vector));

frequency_candidates = frequency_candidates(peak_level_dB < param.comb.max_difference_dB);

end

