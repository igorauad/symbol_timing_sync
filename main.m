clearvars, clc, close all

%% Debug Configuration

debug_tl_static  = 0; % Show static debug plots after sync processing
debug_tl_runtime = 0; % Open scope for debugging of sync loop iterations

%% Parameters
L        = 32;         % Oversampling factor
M        = 4;          % Constellation order
N        = 2;          % Dimensions per symbol (1 for PAM, 2 for QAM)
nSymbols = 1e5;        % Number of transmit symbols
Bn_Ts    = 0.01;       % Loop noise bandwidth (Bn) times symbol period (Ts)
eta      = 1;          % Loop damping factor
rollOff  = 0.2;        % Pulse shaping roll-off factor
timeOffset = 25;       % Simulated channel delay in samples
fsOffsetPpm = 0;       % Sampling clock frequency offset in ppm
rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay
EsN0     = 20;         % Target Es/N0
Ex       = 1;          % Average symbol energy
TED      = 'ZCTED';    % TED (MLTED, ELTED, ZCTED, GTED, or MMTED)
intpl    = 2;          % 0) Polyphase; 1) Linear; 2) Quadratic; 3) Cubic
forceZc  = 0;          % Use to force zero-crossings and debug self-noise

%% System Objects
% Tx Filter
TXFILT = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', L, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);

% Rx Matched Filter (MF)
%
% NOTE: in most simulations, the decimation factor would be L below.
% However, here we want to process the filtered sequence as-is, without the
% downsampling stage. The symbol timing recovery loop processes the
% fractionally-spaced sequence and handles the downsampling process.
RXFILT = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...
    'DecimationFactor', 1, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);
mf = RXFILT.coeffs.Numerator;

% Derivative matched filter (dMF)
dmf = derivativeMf(mf, L);

% Digital Delay
DELAY = dsp.Delay(timeOffset);

% Reference constellation for MER measurement
if (N==2)
    const = qammod(0:M-1,M);
else
    const = pammod(0:M-1,M);
end
Ksym = modnorm(const, 'avpow', Ex);
const = Ksym * const;

% MER Meter
mer = comm.MER;
mer.ReferenceSignalSource = 'Estimated from reference constellation';
mer.ReferenceConstellation = const;

%% Loop Constants
% Time-error Detector Gain (TED Gain)
Kp = calcTedKp(TED, rollOff);

% Scale Kp based on the average symbol energy (at the receiver)
K  = 1; % Assume channel gain is unitary (or that an AGC is used)
Kp = K * Ex * Kp;
% NOTE: if using the GTED when K is not 1, note the scaling is K^2 not K.

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, but in the context of timing
% recovery loop.

% PI Controller Gains:
[ K1, K2 ] = piLoopConstants(Kp, K0, eta, Bn_Ts, L);

fprintf("Loop constants:\n");
fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);

% MATLAB's symbol synchronizer for comparison
tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
    {'Early-Late (non-data-aided)', ...
     'Zero-Crossing (decision-directed)', ...
     'Gardner (non-data-aided)', ...
     'Mueller-Muller (decision-directed)'
     });
if strcmp(TED, "MLTED")
    warning("MLTED not supported by MATLAB's synchronizer - using ZCTED");
    matlabTed = "ZCTED";
else
    matlabTed = TED;
end
SYMSYNC = comm.SymbolSynchronizer(...
    'TimingErrorDetector', tedMap(matlabTed), ...
    'SamplesPerSymbol', L, ...
    'NormalizedLoopBandwidth', Bn_Ts, ...
    'DampingFactor', eta, ...
    'DetectorGain', Kp);

%% Random Transmit Symbols
if (forceZc)
    % Force zero-crossings within the transmit symbols. Use to eliminate
    % the problem of self-noise and debug the operation of the loop
    data = zeros(nSymbols, 1);
    data(1:2:end) = M-1;
else
    data = randi([0 M-1], nSymbols, 1);
end

if (N==2)
    modSig = Ksym * qammod(data, M);
else
    modSig = real(Ksym * pammod(data, M));
end
% Ensure the average symbol energy is unitary, otherwise the loop constants
% must be altered (because Kp, the TED gain, must scale accordingly).

%% Simulation: Transmission -> Channel -> Reception (MF/dMF Filtering)
% Tx Filter
txSig = step(TXFILT, modSig);

% Sampling clock frequency offset
%
% The frequencies produced by the Tx and Rx sampling clocks are often
% significantly distinct, unless both sides adopt high-accuracy oscillators
% (e.g., atomic clocks) or clock disciplining mechanisms, such as with
% GPSDOs. Simulate this relative frequency offset by resampling the signal.
fsRatio = 1 + (fsOffsetPpm * 1e-6); % Rx/Tx clock frequency ratio
tol = 1e-9;
[P, Q] = rat(fsRatio, tol); % express the ratio as a fraction P/Q
txResamp = resample(txSig, P, Q);

% Channel
delaySig = step(DELAY, txResamp);
txSigPower = 1 / sqrt(L);
rxSig = awgn(delaySig, EsN0, txSigPower);

% Rx filter
rxSample = step(RXFILT, rxSig);

% dMF filtering
rxSampleDiff = filter(dmf, 1, rxSig);

%% Symbol Timing Recovery
% Downsampled symbols without symbol timing recovery
rxNoSync = downsample(rxSample, L);

% Downsampled symbols with perfect symbol timing recovery
rxPerfectSync = downsample(rxSample, L, timeOffset);

% Our symbol timing recovery implementation
[ rxSync1 ] = symbolTimingSync(TED, intpl, L, rxSample, rxSampleDiff, ...
    K1, K2, const, Ksym, debug_tl_static, debug_tl_runtime);
% MATLAB's implementation
rxSync2 = step(SYMSYNC, rxSample);

%% Plots and Measurement
skip = 0.2 * nSymbols; % skip the initial transitory when plotting

scatterplot(rxNoSync(skip:end))
title('No Timing Correction');

scatterplot(rxPerfectSync(skip:end))
title('Ideal Timing Correction');

scatterplot(rxSync1(skip:end))
title(sprintf('Our %s Timing Recovery', TED));

scatterplot(rxSync2(skip:end))
title(sprintf('MATLAB''s %s', matlabTed));

fprintf("\nMeasured MER:\n")
fprintf("No Timing Correction: %.2f dB\n", mer(rxNoSync(skip:end)))
fprintf("Ideal Timing Correction: %.2f dB\n", mer(rxPerfectSync(skip:end)))
fprintf("Our %s Timing Recovery: %.2f dB\n", TED, mer(rxSync1(skip:end)))
fprintf("MATLAB's %s Timing Recovery: %.2f dB\n", ...
    matlabTed, mer(rxSync2(skip:end)))