function [ E ] = polyphaseFilterBank(L, I, P, Blim)
% Polyphase Interpolator Filter Bank
%
% polyphaseFilterBank(L, I, P, Blim) returns a matrix containing the
% subfilters (polyphase branches) of a polyphase interpolator.
%
%  Input Arguments:
% L    -> MF oversampling ratio (used to align the filter delays).
% I    -> Interpolator upsampling factor.
% P    -> Neighbor samples weighted by the interp filter.
% Blim -> Bandlimitedness of the interpolated sequence.
%
%  Output
% E    -> Polyphase interpolation filter bank
%
%% Design of L-band interpolation filter
% Interpolation is a process that combines upsampling with a subsequent
% low-pass filter to remove spectral images (anti-imaging filter).
% Typically, the adopted filter is the so-called L-band filter (or L-th
% Band Filter, or Nyquist Filter), namely a filter whose zero-crossings are
% located at integer multiples of L (the interpolation factor). Since the
% zero-interpolation applied to a given sequence (i.e., the upsampling
% operation) introduces L-1 zeros in between every sample of the original
% sequence, the subsequent filtering by an L-band filter (with
% zero-crossings spaced by L) preserves the input samples in the output and
% fills in the zeros based neighbor samples of the original sequence.
%
% Consider the following example sequence:
%
%   x = [1 2 3],
%
% and an L-band filter for L = 2:
%
%  h = [-0.1, 0, 0.6, 1, 0.6, 0, -0.1].
%
% The upsampled by 2 sequence is:
%
%   x_up = [ 1 0 2 0 3 0].
%
% Now consider the result when the peak of the flipped h (center value,
% equal to 1), is aligned with the third sample in x_up (equal to 2)
% through the convolution. The inner product becomes:
%
%    [ 0,  0.6,   1,  0.6,   0,  -0.1] --> (h)
%  .*[ 1,    0,   2,    0,   3,     0] --> (x_up)
%    ----------------------------------------
%   = 2
%
% Namely, the input sample is preserved in the output, because the
% zero-crossings in h coincide with the non-zero samples of x_up other than
% the one of interest (the one aligned with the peak of h).
%
% In contrast, consider the result when the peak of h is aligned with the
% fourth sample in the upsampled sequence (a zero that must be "filled"):
%
%    [ -0.1,  0,  0.6,   1,  0.6,   0] --> (h)
%  .*[    1,  0,    2,   0,    3,   0] --> (x_up)
%    ----------------------------------------
%   = -0.1 * 1 + 0.6 * 2 + 0.6 * 3
%   = 2.9
%
% In this case, the result consists of a weighted sum based primarily on
% the two adjacent neighbor non-zero samples of x_up (2 and 3).
%
% The number of non-zero samples to be weighted in the interpolation is
% determined by the parameter "P". From another perspective, P can be used
% to tune the filter length, which is "2*I*P - 1".
%
% Another parameter of interest is the so-called "bandlimitedness factor"
% (parameter "Blim"), which specifies the bandwidth occupied by the signal
% to be interpolated. More specifically, Blim indicates the signal is
% mostly contained within the bandwidth from "-Blim*pi" to "Blim*pi" (in
% two-sided representation). In case the signal occupies the full
% normalized spectrum from -pi to pi , then Blim=1. Since the sequence to
% be interpolated by the interpolator of the timing recovery scheme is
% actually the fractionally-spaced sequence at the receiver, we can safely
% assume its spectrum does not occupy the full bandwidth, so Blim can be
% less than 1 to allow for a smooth filter transition bandwidth.

% Finally, design the anti-imaging interpolation filter:
interpFilt = intfilt(I, P, Blim);
% By plotting interpFilt, it can be seen that the designed L-band filter
% has a peak at index P*I (considering MATLAB indexing), which implies a
% delay of "P*I - 1". However, since the interpolated sequence is
% immediately fully downsampled by the same factor of the interpolation,
% the delay becomes approximately P, aside from a fractional term. Another
% interpretation is that delay becomes the delay in each polyphase branch
% of the filter's polyphase decomposition.

%% Polyphase realization of the L-band interpolation filter
% The goal of a polyphase realization is to perform filtering before
% upsampling, at a lower rate, rather than after upsampling. The rationale
% is that upsampling interpolates the sequence with zeros, whose filtering
% is unnecessary and irrelevant. We only need to filter the original
% sequence and rearrange the outputs to produce the sample result that
% would be obtained through regular upsampling + filtering.
%
% A polyphase interpolator consists essentially of "I" parallel filters
% (also referred to as subfilters or polyphase branches), which operate in
% the same original lower-rate sequence of the input sequence to be
% interpolated (say Fs). Each subfilter is followed by an independent
% upsampler and then a delay-chain in the output. As explained in multirate
% signal processing books, this is equivalent to having the branches
% followed by a commutator that operates at the higher rate (I*Fs) and
% picks the output of the independent subfilters sequentially.
%
% The following lines obtain the polyphase branches:

% First zero-pad the impulse response of the anti-imaging FIR filter to an
% integer multiple of I:
nZerosToPad = I - mod(length(interpFilt), I);
paddedInterpFilt = [interpFilt zeros(nZerosToPad, 1)];

% Now, split the padded sequence into "I" branches. First, recognize that
% each subfilter has a length equivalent to:
lenSubfilt = length(paddedInterpFilt) / I;

% Then, save the individual subfilter responses in a matrix:
E = reshape(paddedInterpFilt, I, lenSubfilt);

end

