function [av,tau] = allan_variance(phi,Ts,m)
% Function ALLAN_VARIANCE
% phi: the phase noise signal over time
% Ts: sampling time of the phase noise signal
% m: sequence of integers that corresponds to the number of samples average
% on which the allan variance is calculated

% av: the allan variance sequence, correspondend to the observation time
% tau
L = length(phi);
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.

tau = m*Ts;

av = zeros(length(m), 1);
for i = 1:length(m)
    mi = m(i);
    av(i,:) = sum( ...
        (phi(1+2*mi:L) - 2*phi(1+mi:L-mi) + phi(1:L-2*mi)).^2, 1);
end
av = av ./ (2*tau.^2 .* (L - 2*m));
end