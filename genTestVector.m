function [x, y] = genTestVector(nSymbols, sps, rolloff, rcDelay, ...
    tedChoice, interpChoice, loopBw, dampingFactor)
% Generate input/output vectors for testing
%
% [x, y] = genTestVector(nSymbols, sps, rolloff, rcDelay, tedChoice, ...
%   interpChoice, loopBw, dampingFactor)
%
% Generates input and output test vectors to support the comparison of
% other symbol timing recovery implementations (e.g., in RTL or C/C++)
% against this project's implementation from function symbolTimingSync. The
% input test vector consists of the sample-spaced matched filter output due
% to a sequence of random unit-energy QPSK symbols. The output test vector
% is the corresponding symbol-spaced interpolated symbol sequence output by
% the symbol timing recovery loop.
%
% Args:
%   nSymbols -> Number of QPSK symbols to generate.
%   sps -> Oversampling ratio (samples per symbol).
%   rolloff -> Raised cosine pulse rolloff factor.
%   rcDelay -> Raised cosine pulse delay in symbols.
%   tedChoice -> Timing error detector to be used by symbolTimingSync.
%   interpChoice -> Interpolator method to be used by symbolTimingSync.
%   loopBw -> Normalized loop bandwidth.
%   dampingFactor -> Loop damping factor.
%
% Example with nSymbols=10, sps=2, rolloff=0.2, rcDelay=10,
% tedChoice='GTED' (Gardner TED), interpChoice=1 (linear interpolator),
% loopBw = 0.01, and dampingFactor=1:
%
%   genTestVector(10, 2, 0.2, 10, 'GTED', 1, 0.01, 1);
%
% NOTE: The input test vector comprises the "central part" of the
% convolution between the upsampled symbol sequence and the raised cosine
% pulse, after the pulse shaping filter's transitory. It is the result of
% calling "conv()" with "SHAPE=same".

% Test Constellation
M = 4; % Constellation size
Ex = 1; % Average symbol energy
const = qammod(0:(M-1), M);
Ksym = modnorm(const, 'avpow', Ex);
const = Ksym * const;

% Loop constants
Kp = calcTedKp(tedChoice, rolloff);
K0 = -1;
[ K1, K2 ] = piLoopConstants(Kp, K0, dampingFactor, loopBw, sps);

% Raised cosine pulse (Tx filter convolved with the matched filter)
h = rcosine(1, sps, 'normal', rolloff, rcDelay);

% Test symbols
data = randi([0 M-1], nSymbols, 1);
test_syms = Ksym * qammod(data, M);

% Matched filter output
test_syms_up = upsample(test_syms, sps);
y_mf = conv(test_syms_up, h, 'same');

% The matched filter output is the input to the symbol timing recovery loop
x = y_mf;

% Symbol timing recovery
y = symbolTimingSync(tedChoice, interpChoice, sps, x, [], K1, K2, ...
    const, Ksym);

printComplexVec(x, "x");
printComplexVec(y, "y");
end

function [] = printComplexVec(x, label)
    fprintf("%s = [", label);
    for i = 1:(length(x) - 1)
        fprintf('(%f%+fj), ', real(x(i)), imag(x(i)));
    end
    i = length(x);
    fprintf('(%f%+fj)]\n', real(x(i)), imag(x(i)));
end