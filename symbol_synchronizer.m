clearvars, clc, %close all

%% Parameters
L        = 32;
nSymbols = 10000;
Bn_Ts    = 0.005;      % PLL noise bandwidth (Bn) times symbol period (Ts)
xi       = 1/sqrt(2);
rollOff  = 0.5;
rndDelay = 5;          % Random delay added to add a disturbance
rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay
SNR      = 25;
Ex       = 1;          % Average symbol energy
TED      = 'MLTED';    % TED Type

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
DELAY   = dsp.Delay(rndDelay);

% Symbol Synchronizer
SYMSYNC = comm.SymbolSynchronizer('SamplesPerSymbol', L);

%% MF
mf  = RXFILT.coeffs.Numerator;

%% dMF
% IMPORTANT: use central-differences to match the results in the book
h = (2)*[0.5 0 -0.5]; % kernel function
% It is the future symbol minus the past symbol, divided by the interval
% that corresponds to 2*T, where T is the sampling period. Here, we
% consider T = 1/L.
central_diff_mf = conv(h, mf);
% Skip the filter delay
dmf = central_diff_mf(2:1+length(mf));


figure
plot(mf)
hold on
plot(dmf, 'r')
legend('MF', 'dMF')

%% PLL Design

% Time-error Detector Gain (TED Gain)
Kp = getTedKp(TED, L, rollOff, rcDelay);
% Scale Kp based on the average symbol energy (at the receiver)
K  = 1; % Assume channel gain is unitary
Kp = K*Ex*Kp;

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, in the context of timing sync
% loop:

% PI Controller Gains:
[ K1, K2 ] = timingLoopPIConstants(Kp, K0, xi, Bn_Ts, L)

%% Random PSK Symbols
data    = randi([0 3], nSymbols, 1);
modSig  = real(modnorm(pammod(0:3,4), 'avpow', Ex) * pammod(data, 4));
% Important, ensure to make the average symbol energy unitary, otherwise
% the PLL constants must be altered (because Kp, the TED gain, scales).

%%%%%%%%%%%%%%% Tx Filter  %%%%%%%%%%%%%%%
txSig    = step(TXFILT,modSig);

%%%%%%%%%%%%%%% Channel    %%%%%%%%%%%%%%%
delaySig = step(DELAY,txSig);
rxSig    = awgn(delaySig, SNR, 'measured');

%%%%%%%%%%%%%%% Rx filter  %%%%%%%%%%%%%%%
rxSample = step(RXFILT,rxSig);

%% Decisions based on Timing Correction
scatterplot(downsample(rxSample, L), 2)

%% Decisions based on MLTED Timing Recovery
k         = 1;
underflow = 0;
mu_next   = 0;
CNT_next  = 1;
vi        = 0;

rIBuff = zeros(TXFILT.FilterSpanInSymbols*L, 1);
for n=2:length(rxSig)
    CNT = CNT_next;
    mu = mu_next;
    rI = mu * rxSig(n) + (1 - mu)*rxSig(n-1);
    % I believe "n-1" can be interpreted as the basepoint index
    % Note it has to be "rxSig", prior to the MF

    x = mf * [rI; rIBuff];
    xdot = dmf * [rI; rIBuff];

    if underflow == 1
        e = sign(x)*xdot;
        xx(k) = x;
        k = k+1;
    else
        % Upsample:
        e = 0;
    end

    vp = K1*e;
    vi = vi + K2*e;
    v(n) = vp + vi;
    W = 1/L + v(n);

    CNT_next = mod(CNT - W, 1);

    if (CNT - W < 0) % || (CNT - W > 1)
        underflow = 1;
        mu_next = CNT/W;
    else
        underflow = 0;
        mu_next = mu;
    end
    rIBuff = [rI; rIBuff(1:end-1)];
end

scatterplot(xx, 2)

figure
plot(v)

%% Decisions based on MATLAB Timing Error Correction

rxSync = step(SYMSYNC,rxSample);
scatterplot(rxSync(1001:end),2)