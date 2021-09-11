function [ normTauE, g ] = simSCurve(TED, rollOff, rcDelay, nSymbols)
%% Compute the S-Curve for the data-aided TED operation
%
% Simulates the average error detected by a TED of choice for varying
% receiver timing offsets. The simulation sends a number of random 2-PAM
% symbols over the cascaded combination of the Tx and matched filters. Each
% transmitted symbol produces a timing error detection out of the TED, and
% the average of these detected errors becomes the S-curve value for that
% timing offset. This process repeats for all simulated timing offsets.
%
% [ normTauE, g ] = simSCurve(TED, rollOff, rcDelay, nSymbols) returns the
% S-curve g(normTauE) of the chosen TED obtained through simulation. It
% also returns the vector normTauE with the normalized time offset errors
% (within -0.5 to 0.5) where the S-curve is evaluated.
%    TED -> TED choice.
%    rollOff -> Rolloff factor.
%    rcDelay -> Raised cosine pulse delay (default: 10).
%    nSymbols -> Number of random 2-PAM symbols to simulate (default: 1e4).

if nargin < 3
    rcDelay = 10;
end

if nargin < 4
    nSymbols = 1e4;
end

% Oversampling (use a large value to observe the S-curve with enough
% resolution). Note the oversampling factor that is effectively used on a
% symbol timing recovery loop has nothing to do with the factor adopted
% here. Here, the only goal is to observe the S curve, and the S curve (a
% continuous time metric) should be independent of the oversampling factor.
L = 1e3;

% Root raised cosine Tx and Rx filters
htx = rcosdesign(rollOff, rcDelay, L);
mf = conj(fliplr(htx));  % Matched filter (MF)

% Corresponding raised cosine pulse shape
r_p = conv(htx, mf);

% Derivative matched filter (dMF)
h = L * [0.5 0 -0.5]; % first central difference (see Eq. 3.61)
dmf = conv(h, mf);
dmf = dmf(2:end-1); % Skip the tail and head samples

% Random data (note it must be random)
data = randi(2, 1, nSymbols) - 1;

% 2-PAM Tx symbols
txSym = real(pammod(data, 2));

% Tx upsampled sequence
tau = 0; % true time delay/offset
txUpSequence = upsample(txSym, L, tau);

% MF output sequence
mfOutput = conv(r_p, txUpSequence);

% dMF output sequence
dmfOutput = conv(conv(htx, txUpSequence), dmf);

% Vectors with the time offset estimates (tauEst) and the corresponding
% estimate errors (tauErr). Note "tauErr = tau - tauEst", as defined after
% Eq. 8.25. However, "tau = 0" here, so "tauErr = -tauEst". The intuition
% is as follows: suppose we are simulating the receiver signalling the
% strobe L/2 sample intervals earlier (before the ideal strobe location).
% In that case, the tauEst is "-L/2" (negative). However, the difference
% "tau - tauEst" is "+L/2" (positive).
tauEstVec = -(L/2):(L/2); % rx time delay/offset estimate
tauErrVec = tau - tauEstVec; % time offset estimate error
normTauE = tauErrVec / L; % normalized time offset estimate error
% NOTE: our tau metric (estimate or true) is given in terms of sample
% periods, not in seconds. Hence, the normalization is accordingly by "L",
% in units of sample periods, not by "Ts" as used in the book.

% S-Curve simulation
g = zeros(1, length(tauErrVec));
for i = 1:length(tauEstVec)
    tauEst = tauEstVec(i);
    tauErr = tau - tauEst;
    idealStrobeIdx = (rcDelay * L) + (1:L:(nSymbols * L));
    strobeIdx = idealStrobeIdx + tauEst;
    rxSym = mfOutput(strobeIdx);

    switch (TED)
    case 'MLTED' % Maximum Likelihood TED
        e = dmfOutput(strobeIdx) .* txSym;
    case 'ELTED' % Early-late - See Eq. (8.35)
        earlyIdx = strobeIdx + L/2;
        lateIdx = strobeIdx - L/2;
        e = txSym .* (mfOutput(earlyIdx) - mfOutput(lateIdx));
    case 'ZCTED' % Zero-crossing TED
        zcStrobeIdx = strobeIdx + L/2;
        zcSamples = mfOutput(zcStrobeIdx);
        e = zcSamples(1:end-1) .* (txSym(1:end-1) - txSym(2:end));
    case 'GTED' % Gardner TED
        zcStrobeIdx = strobeIdx(2:end) - L/2; % Zero-crossings
        prevStrobeIdx = strobeIdx(1:end-1); % Previous symbols
        currStrobeIdx = strobeIdx(2:end); % Current symbols
        e = mfOutput(zcStrobeIdx) .* ...
            (mfOutput(prevStrobeIdx) - mfOutput(currStrobeIdx));
        % NOTE: the GTED is purely data-aided. The other TEDs investigated
        % by this function are implemented in their decision-directed mode.
    case 'MMTED' % Mueller and Müller
        prevStrobeIdx = strobeIdx(1:end-1); % Previous symbols
        currStrobeIdx = strobeIdx(2:end); % Current symbols
        e = txSym(1:end-1) .* mfOutput(currStrobeIdx) - ...
            txSym(2:end) .* mfOutput(prevStrobeIdx);
    end

    % S-Curve value for this time offset estimation error
    idx = tauErrVec == tauErr;
    g(idx) = mean(e);
end

end