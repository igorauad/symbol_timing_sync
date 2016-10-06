function [ E ] = polyphaseFilterBank(L, I, P, Blim)
% Polyphase Interpolator Filter Bank
%   Returns a matrix containing the subfilters (polyphase branches) for a
%   polyphase interpolator.
%
%  Input Arguments:
% L    -> MF oversampling ratio (used to align the filter delays)
% I    -> Interpolator upsampling factor
% P    -> Neighbor samples weighted by the interp filter
% Blim -> Bandlimitedness of the interpolated sequence
%
%  Output
% E    -> Polyphase interpolation filter bank
%% Design of L-band interpolation filter
% Interpolation is a process that combines upsampling with a subsequent
% low-pass filter for removing spectral images (anti-imaging filter). A
% filter that is often used is the so-called L-band filter (or L-th Band
% Filter or Nyquist Filter), namely a filter whose zero-crossings are
% located at integer multiples of L (the interpolation factor). Since the
% zero-interpolation applied to a given sequence (i.e. the upsampling
% operation) introduces L-1 zeros in between every sample of the original
% sequence, the subsequent filtering by such an L-band filter (with zeros
% crossings spaced by L) preserves the input samples in the output and
% "fills" the zeros based on the values at a given number of neighbor
% samples of the original sequence. To fix ideas, consider an example
% sequence:
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
% Now consider the result when the peak of flipped h (center value, equal
% to 1), is aligned with the third sample in x_up (equal to 2) through the
% convolution. The inner product becomes:
%
%    [ 0,  0.6,   1,  0.6,   0,  -0.1] --> (h)
%  .*[ 1,    0,   2,    0,   3,     0] --> (x_up)
%    ----------------------------------------
%   = 2
% Namely, the input sample is preserved in the output, because
% zero-crossings in h coincide with the non-zero samples of x_up other than
% the one of interest (the one aligned with the peak of h).
%
% In constrast, consider the result when the peak of h is aligned with the
% fourth sample in the upsampled sequence (a zero that must be "filled"):
%
%    [ -0.1,  0,  0.6,   1,  0.6,   0] --> (h)
%  .*[    1,  0,    2,   0,    3,   0] --> (x_up)
%    ----------------------------------------
%   = -0.1 * 1 + 0.6 * 2 + 0.6 * 3
%   = 2.9
% In this case, the result consist of a weighted sum considering primarily
% the two adjacent neighbor non-zero samples of x_up.
%
% The number of non-zero samples to be weighted in the interpolation is
% mandated by the variable "Pintfilt" (in the parameters section). From
% another perspective, Pintfilt can be used to tune the filter length,
% which is "2*interpFactor*Pintfilt - 1".
%
% Another parameter of interest is the so-called "Bandlimitedness factor"
% (also in the parameters section), which informs the bandwidth that the
% signal to be interpolated occupies. More specifically, factor
% "bandlimFactor" informs that the signal is mostly contained within the
% bandwidth from -bandlimFactor*pi to bandlimFactor*pi (in its two-sided
% representation). In case it occupies the full normalized spectrum from
% -pi to pi , then the factor should be 1. Since sequence to be
% interpolated by this interpolator of the timing recovery scheme is
% actually the fractionally-spaced sequence at the receiver, we can safely
% assume its spectrum does not occupy the full bandwidth, bandlimFactor can
% be less than 1 to allow smoother filter transition bandwidth.

% Finally, design the anti-imaging filter:
B_antiImaging = intfilt(I, P, Blim);
% If B_antiImaging is plotted, it can be seen that the designed L-band
% filter has a peak at index Pintfilt*interpFactor (considering MATLAB
% indexing), which implies a delay of "Pintfilt*interpFactor - 1". However,
% since the interpolated sequence is immediately fully downsampled in the
% receiver (downsampled by the same factor of the interpolation), the delay
% becomes aproximately Pintfilt, aside from a fractional term. Another
% interpretation is that delay becomes the delay in each polyphase branch
% of the filter's polyphase decomposition. More details on that in the next
% code sections.

%% Polyphase realization of the L-band interpolation filter
% The goal of a polyphase realization of the interpolator is to perform
% filtering before upsampling, at a lower rate, rather than after
% upsampling. The rationale is that upsampling interpolates the sequence
% with zeros, whose filtering is unnecessary and irrelevant. We only need
% to filter the original sequence and rearrange the outputs to produce the
% sample result that would be obtained through regular upsampling +
% filtering.
%
% A polyphase interpolator consists essentially of "interpFactor" parallel
% filters or subfilters (also referred to as polyphase branches), which
% operate in the same original lower-rate sequence (at Fs). Each branch
% filter is followed by an independent upsampler and then a delay-chain in
% the output. As explained in most multirate signal processing books, this
% is equivalent as having the branches followed by a commutator that
% operates at the higher rate (interpFactor*Fs) and picks the output of the
% independent subfilters sequentially.
%
% The following lines obtain each of the subfilters or polyphase branches:

% First zero-pad the impulse response of the anti-imaging FIR filter to an
% integer multiple of interpFactor:
nZerosToPad = I - mod(length(B_antiImaging), I);
B_antiImaging_padded = [B_antiImaging zeros(nZerosToPad, 1)];

% Now, split the padded sequence into "interpFactor" branches. First,
% recognize that each subfilter has a length equivalent to:
Len_subfilt = length(B_antiImaging_padded)/I;
% Then, save the individual subfilter responses in a matrix:
E = reshape(B_antiImaging_padded, I, Len_subfilt);

end

