function [ polyFiltBank ] = polyDecomp(filt, L)
% Polyphase Decomposition
%
% Decompose a given FIR filter into L polyphase subfilters. For symbol
% timing recovery, the given filter is expected to be an interpolating
% filter, such as the one produced by function "intfilt()". Namely, "filt"
% is expected to be an L-band filter with zero-crossings after every L
% samples, except for the central value. This function decomposes the
% L-band filter into L subfilters, and each subfilter can be used
% independently to obtain a particular phase of the output sequence
% according to the estimated symbol timing offset.
%
%  Input Arguments:
% filt  -> Interpolating filter to be decomposed into polyphase branches.
% L     -> Interpolating filter's intrinsic upsampling factor.
%
%  Output
% polyFiltBank -> Polyphase interpolation filter bank

% First zero-pad the FIR filter to an integer multiple of L if necessary:
if (mod(length(filt), L) == 0)
    paddedFilt = filt;
else
    nZerosToPad = L - mod(length(filt), L);
    paddedFilt = [filt zeros(1, nZerosToPad)];
end

% Next, split the padded sequence into "L" branches/subfilters:
lenSubfilt = length(paddedFilt) / L;
polyFiltBank = reshape(paddedFilt, L, lenSubfilt);
end