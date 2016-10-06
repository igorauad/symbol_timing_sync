clearvars, clc, %close all

%% Debug Configuration

debug_tl_static  = 0; % Show static debug plots after sync processing
debug_tl_runtime = 0; % Open scope for debugging of sync loop iterations

%% Parameters
L        = 32;         % Oversampling factor
M        = 4;          % Constellation order
nSymbols = 10000;      % Number of transmit symbols
Bn_Ts    = 0.01;       % PLL noise bandwidth (Bn) times symbol period (Ts)
eta      = 1/sqrt(2);  % PLL Damping Factor
rollOff  = 0.5;        % Pulse shaping roll-off factor
timeOffset = 5;        % Delay (in samples) added
rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay
SNR      = 25;         % Target SNR
Ex       = 1;          % Average symbol energy
TED      = 'MLTED';    % TED Type (currently only 'MLTED' is supported)
intpl    = 0;          % 0) Linear; 1) Polyphase Interpolator

%% System Objects

% Tx Filter
TXFILT  = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', L, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);

% Rx Filter (MF)
RXFILT  = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...
    'DecimationFactor',1, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);

% Digital Delay
DELAY   = dsp.Delay(timeOffset);

% Symbol Synchronizer
SYMSYNC = comm.SymbolSynchronizer('SamplesPerSymbol', L);

%% Matched Filter (MF)
mf  = RXFILT.coeffs.Numerator;

%% dMF
% IMPORTANT: use central-differences to match the results in the book
h = [0.5 0 -0.5]; % kernel function
central_diff_mf = conv(h, mf);
% Skip the filter delay
dmf = central_diff_mf(2:1+length(mf));

figure
plot(mf)
hold on, grid on
plot(dmf, 'r')
legend('MF', 'dMF')
title('MF vs. dMF')

%% PLL Design

% Time-error Detector Gain (TED Gain)
Kp = getTedKp(TED, L, rollOff, rcDelay);
% Scale Kp based on the average symbol energy (at the receiver)
K  = 1; % Assume channel gain is unitary
Kp = K*Ex*Kp;

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, but in the context of timing
% recovery loop.

% PI Controller Gains:
[ K1, K2 ] = timingLoopPIConstants(Kp, K0, eta, Bn_Ts, L)
% Note: MATLAB's implementations uses a default value of Kp = 2.7;

%% Random PSK Symbols
data    = randi([0 M-1], nSymbols, 1);
modSig  = real(modnorm(pammod(0:M-1,M), 'avpow', Ex) * pammod(data, M));
% Important, ensure to make the average symbol energy unitary, otherwise
% the PLL constants must be altered (because Kp, the TED gain, scales).

%%%%%%%%%%%%%%% Tx Filter  %%%%%%%%%%%%%%%
txSig    = step(TXFILT,modSig);

%%%%%%%%%%%%%%% Channel    %%%%%%%%%%%%%%%
delaySig = step(DELAY,txSig);
rxSig    = awgn(delaySig, SNR, 'measured');

%%%%%%%%%%%%%%% Rx filter  %%%%%%%%%%%%%%%
rxSample = step(RXFILT,rxSig);

%% dMF

rxSampleDiff = filter(dmf, 1, rxSig);

%% Decoder Inputs (MF Downsampled Output) without Timing Correction
scatterplot(downsample(rxSample, L), 2)
title('No Timing Correction');

%% Decoder Inputs after ML Timing Recovery

[ xx ] = symTimingLoop(intpl, L, rxSample, rxSampleDiff, K1, K2, ...
                       debug_tl_static, debug_tl_runtime);

% Scatter Plot
scatterplot(xx, 2)
title('Using MLTED Timing Recovery');

%% Decoder Inputs using MATLAB's Timing Error Correction

rxSync = step(SYMSYNC,rxSample);

% Scatter Plot
scatterplot(rxSync(1001:end),2)
title('Using MATLABs Zero-Crossing TED');