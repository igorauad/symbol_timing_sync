function [ Kp ] = getTedKp(TED, L, rollOff, delay)
% Returns the Kp corresponding to the timing error detector. This gain
% corresponds to the value of the TED S-curve evaluated at tau_e (the
% timing error) equal to 0. The value returned below assumes unitary
% average symbol energy at the receiver, that is, K*Ex = 0, where K is the
% gain from Tx to Rx and Ex is the original average symbol energy. The
% output of this function shoulds be scaled as needed.
%
%   Inputs
% TED      ->  TED Type
% L        ->  Oversampling factor
% rollOff  ->  Rolloff factor
% delay    ->  Delay for the cascaded composition between Tx and MF filter

%% Nyquist pulse

% Raised cosine (Tx filter convolved with MF), equivalent to the
% autocorrelation function of the pulse shaping filter
r_p = rcosine(1, L, 'normal', rollOff, delay);

%% Derivative of the Nyquist pulse
% IMPORTANT: use central-differences to match the results in the book
h = [0.5 0 -0.5]; % kernel function
% Central-differences:
r_p_diff = conv(h, r_p);
% Skip the filter delay
r_p_diff = r_p_diff(2:1+length(r_p));

%% TED Gain

% Assume constant gain and unitary average symbol energy
K  = 1;
Ex = 1;

% Zero-centered indices corresponding to one symbol interval
oneTsIndex = -L/2 : L/2 - 1;
ind = L*delay+1 + oneTsIndex; % Central indexes

%% S-Curve
switch (TED)
    case 'MLTED' % Maximum Likelihood
        % The S-Curve is "g = K*Ex*r_p_diff(-tau_e)", therefore:
        g = K*Ex*r_p_diff(fliplr(ind));
        % Note 1: The order is reversed (fliplr is used), because r_p_diff
        % in "g" is a function of "-tau_e"
        % Note 2: g(L/2), the central value, should be roughly 0

    case 'ELTED' % Early-late
        % The S-Curve is
        g = K*Ex*( r_p(fliplr(ind) + L/2) - r_p(fliplr(ind) - L/2) );

    case 'ZCTED'
        % The S-Curve is identical to the ELTED
        g = K*Ex*( r_p(fliplr(ind) + L/2) - r_p(fliplr(ind) - L/2) );
    case 'GTED'  % Gardner TED
    case 'MMTED' % Mueller-Muller
end

%% TED Gain
% The gain is the slope at tau_e = 0
% This slope is estimated as follows:
Kp = (g(L/2 + 1) - g(L/2 - 1))/(2/L);
% Note:
% Michael Rice's figure (e.g. 8.4.3) show the x-axis as tau_e/Ts,
% that is the timing error normalized by the symbol period. The two
% values surrounding the central value (at tau_e = 0) are at +T and
% -T. Therefore, at normalized time +T/Ts and -T/Ts. Since Ts =
% L*T, those are at normalized time +1/L and -1/L. The normalized
% interval between the two, therefore, corresponds to (2/L).

if (0)
    figure
    plot((1/L)*oneTsIndex, g)
end

