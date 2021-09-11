function [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay)
%% Compute the S-Curve of a given TED assuming data-aided operation
%
% [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay) returns the S-curve
% g(normTauE) of the chosen TED computed analytically based on closed-form
% expressions. It also returns the vector normTauE with the normalized time
% offset errors (within -0.5 to 0.5) where the S-curve is evaluated.
%    TED -> TED choice.
%    rollOff -> Rolloff factor.
%    rcDelay -> Raised cosine pulse delay (default: 10).
%
% Note: this function assumes constant channel gain (e.g., after automatic
% gain control) and unitary average symbol energy.

if nargin < 3
    rcDelay = 10;
end

%% Parameters and assumptions
% Oversampling factor
%
% Note the oversampling factor that is effectively used on a symbol timing
% recovery loop has nothing to do with the factor adopted here. Here, the
% only goal is to observe the S curve, and the S curve (a continuous time
% metric) should be independent of the oversampling factor. Hence, we use a
% large value to observe the S-curve with enough resolution.
L = 1e3;

% Assume constant gain and unitary average symbol energy
K  = 1;
Ex = 1;

%% Raised cosine pulse

% Raised cosine (Tx filter convolved with MF), equivalent to the
% autocorrelation function of the pulse shaping filter
r_p = rcosine(1, L, 'normal', rollOff, rcDelay);

%% Derivative of the raised cosine pulse
% NOTE: normally a receiver uses the dMF, which is the derivative matched
% filter (dMF) computed by differentiating the MF. However, the goal here
% is to compute the analytical S-curve, not to simulate an actual receiver.
% This computation requires the derivative of the entire raised cosine
% pulse (or pulse shape autocorrelation), as seen in Eq. (8.30). Hence, for
% convenience, we compute the latter directly, and not the dMF.
%
% IMPORTANT: use central-differences to match the results in the book. Note
% also that the first central difference from Eq. (3.61) divides by "2*T",
% where T is the sampling interval. Here, we don't have "T", but we assume
% "T = 1/L", such that the denominator of (3.61) becomes "(2/L)" instead.
h = L * [0.5 0 -0.5]; % kernel function
% Central-differences:
r_p_diff = conv(h, r_p);
% Skip the filter delay
r_p_diff = r_p_diff(2:1+length(r_p));

%% TED Gain

% Zero-centered indices corresponding to one symbol interval
tau_e = -L/2 : L/2; % timing error in units of sample periods
normTauE = tau_e / L; % normalized timing error
idx = L * rcDelay + 1 + fliplr(tau_e); % Central indexes
% Note: reverse the order using fliplr because r_p_diff in "g" is a
% function of "-tau_e".

%% S-Curve
switch (TED)
    case 'MLTED' % Maximum Likelihood - See Eq. (8.30)
        g = K * Ex * r_p_diff(idx);

    case 'ELTED' % Early-late - See Eq. (8.35)
        g = K * Ex * (r_p(idx + L/2) - r_p(idx - L/2));

    case 'ZCTED' % Zero-crossing - See Eq. (8.42)
        g = K * Ex * (r_p(idx + L/2) - r_p(idx - L/2));
        % NOTE: the the ZCTED and ELTED have identical S-curves

    case 'GTED'  % Gardner - See Eq. (8.45)
        C = sin(pi * rollOff / 2) / (4 * pi * (1 - (rollOff^2 / 4)));
        g = (4 * K^2 * Ex) * C * sin(2 * pi * tau_e / L);
        % NOTE:
        % - Eq. (8.45) normalizes tau_e by Ts inside the sine argument.
        %   Here, in contrast, we normalize it by L (oversampling factor).
        % - In Eq. (8.45), the scaling constant (4 * K^2 * Ex) is also
        %   divided by Ts. Here, we don't divide it, not even by L. This
        %   has been confirmed empirically using the "simSCurve" script.

    case 'MMTED' % Mueller and Müller - See Eq. (8.50)
        g = K * Ex * (r_p(idx + L) - r_p(idx - L));
end

end