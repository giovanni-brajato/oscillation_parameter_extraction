function [amp] = amp_estimation_DD(X)
%amp_estimation_DD: estimate the power of the envelope of a signal using
%its Hilbert transform
% X: signal (single quadrature) over time
% amp: power of the signal over time
amp = abs(hilbert(X)).^2;
end

