function [dmf] = derivativeMf(mf, L)
%% Compute the derivative matched filter (dMF)
%
% [dmf] = derivativeMf(mf) returns the derivative matched filter (dMF)
% corresponding to a given matched filter.
%
% Args:
%    mf -> Matched filter taps.
%    L  -> Oversampling factor.

% First central difference based on Eq.(3.61) and using T=1/L
%
% Eq. (3.61) divides by "2*T", where T is the sampling interval. Here, we
% don't have "T", but we assume "T = 1/L", such that the denominator of
% (3.61) becomes "(2/L)" instead.
h = L * [0.5 0 -0.5];
central_diff_mf = conv(h, mf);

% Skip the tail and head so that the dMF length matches the MF length
dmf = central_diff_mf(2:end-1);

end