function [ polyFiltBank ] = polyInterpFilt(I, P, Blim)
% Polyphase Interpolator Filter Bank
%
% polyInterpFilt(I, P, Blim) returns a matrix containing the subfilters
% (polyphase branches) of a polyphase interpolator.
%
%  Input Arguments:
% I    -> Interpolation/upsampling factor.
% P    -> Neighbor samples weighted by the interpolation filter.
% Blim -> Bandlimitedness of the interpolated sequence.
%
%  Output"
% polyFiltBank -> Polyphase interpolation filter bank.
%
%  Example:
%
% Design a polyphase interpolator to upsample by a factor of two, while
% weighting 2*P=4 neighbor samples for each interpolant (P on each side):
%
%   polyFiltBank = polyInterpFilt(2, 2, 0.5)
%
% ---------
%
% Interpolation is a process that combines upsampling with a subsequent
% low-pass filter to remove spectral images (anti-imaging filter).
% Typically, the adopted filter is the so-called L-band filter (or L-th
% Band Filter, or Nyquist Filter), namely a filter whose zero-crossings are
% located at integer multiples of L (the interpolation factor). Since the
% zero-interpolation applied to a given sequence (i.e., the upsampling
% operation) introduces L-1 zeros in between every sample of the original
% sequence, the subsequent filtering by an Lth-band filter (with
% zero-crossings spaced by L) preserves the input samples in the output and
% fills in the zeros based neighbor samples of the original sequence.
%
% Consider the following example sequence:
%
%   x = [1 2 3],
%
% and an Lth-band filter for L=2:
%
%  h = [-0.1, 0, 0.6, 1, 0.6, 0, -0.1].
%
% The upsampled-by-2 sequence is:
%
%   x_up = [ 1 0 2 0 3 0 4 0 ].
%
% Now, consider the result when the peak of the flipped h (center value
% equal to 1) is aligned with the fifth sample in x_up (equal to 3) through
% the convolution. In this case, the inner product becomes:
%
%    [    -0.1,   0,   0.6,   1,   0.6,   0,  -0.1 ] --> (h)
%  .*[ 1,   0,    2,    0,    3,    0,    4,    0  ] --> (x_up)
%    ------------------------------------------------
%   = 3
%
% Namely, the input sample is preserved in the output since the
% zero-crossings in h coincide with the non-zero samples of x_up other than
% the one of interest (the one aligned with the peak of h).
%
% In contrast, consider the result when the peak of h is aligned with the
% fourth sample in the upsampled sequence (a zero that must be filled in):
%
%    [ -0.1,   0,   0.6,   1,   0.6,   0,  -0.1       ] --> (h)
%  .*[    1,   0,    2,    0,    3,    0,    4,    0  ] --> (x_up)
%    ------------------------------------------------
%   = (-0.1 * 1) + (0.6 * 2) + (0.6 * 3) + (-0.1 * 4)
%   = 2.5
%
% In this case, the result consists of a weighted sum based primarily on
% the two adjacent neighbor non-zero samples of x_up on each side. Namely,
% based on 2*P samples in total, P on each side of the target index.
%
% That is, the number of non-zero samples to be weighted in the
% interpolation is determined by the parameter "P". From another viewpoint,
% P can be used to tune the filter length, which is "2*I*P - 1". In the
% given example, I=2 and P=2, so the filter length is 7. Besides, "P"
% determines the filter delay. The interpolation filter has a peak at index
% P*I (considering MATLAB indexing), which implies a delay of "P*I - 1".
% However, when the interpolated sequence is immediately fully downsampled
% by the same factor I, the delay becomes approximately P, aside from a
% fractional term. Another interpretation is that delay becomes the delay
% in each polyphase branch of the filter's polyphase decomposition.
%
% Another parameter of interest is the so-called bandlimitedness factor
% ("Blim"), which specifies the bandwidth occupied by the signal to be
% interpolated. More specifically, Blim indicates the signal is mostly
% contained within the bandwidth from "-Blim*pi" to "Blim*pi" (in two-sided
% representation). In case the signal occupies the full normalized spectrum
% from -pi to pi, then Blim=1. Since the sequence to be interpolated by the
% interpolator of the timing recovery scheme is actually the
% fractionally-spaced sequence at the receiver, we can safely assume its
% spectrum does not occupy the full bandwidth, so Blim can be less than 1
% to allow for a smooth filter transition bandwidth.
%
% This function returns specifically the polyphase realization of the
% Lth-band filter. The goal of a polyphase realization is to perform
% filtering before upsampling, at a lower rate, rather than after
% upsampling. The rationale is that the upsampling step I-1 zeros between
% each sample of the input sequence, but the filtering of such zero values
% is unnecessary and irrelevant. It is only necessary to filter the
% original sequence and rearrange the outputs to produce the sample result
% that would be obtained through regular upsampling plus filtering.
%
% A polyphase interpolator consists essentially of "I" parallel filters
% (also called subfilters or polyphase branches), which operate in the same
% original low-rate input sequence to be interpolated (say, at rate Fs).
% The result of each subfilter can be serialized in sequence on the final
% output. This serialization is carried out by a "commutator", which
% operates at the higher rate (say, I*Fs) by picking the output of the
% independent subfilters sequentially. In many occasions, such as silicon
% implementations, it is significantly easier to implement a commutator at
% rate "I*Fs" while filtering at rate "Fs", than to implement a full filter
% at rate "I*Fs". Hence, the polyphase realization becomes handy.
%
% Nevertheless, for symbol timing recovery, the polyphase interpolation
% filter is used differently. It is not used to increase the rate of an
% input signal. Instead, it is used to compute the interpolated symbol
% between each selection of input samples. Namely, instead of a
% rate-increasing use case, the interpolator is used for rate-decreasing
% (downsampling) case. In this context, there is no commutator block
% following the polyphase interpolator. The polyphase interpolator takes an
% oversampled sequence on its input and filters I distinct phases of this
% sequence in parallel. Meanwhile, a symbol timing recovery loop selects
% the output of a single polyphase subfilter at a time whenever it decides
% it is time to obtain a new interpolated symbol. Each subfilter processes
% a particular phase of the oversampled signal, and the timing recovery
% loop can pick the appropriate phase according to its estimate of the
% symbol timing offset. Besides, note the input sequence can be oversampled
% by any arbitrary factor, not necessary by a factor of I. The polyphase
% interpolation factor I only determines how many phases of the input
% sequence are observed in parallel.

% Anti-imaging interpolation filter:
interpFilt = intfilt(I, P, Blim);
% Polyphase decomposition
polyFiltBank = polyDecomp(interpFilt, I);

end

