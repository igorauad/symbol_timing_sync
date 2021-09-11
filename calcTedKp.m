function [ Kp ] = calcTedKp(TED, rollOff, method, plotSCurve, rcDelay)
% Returns the gain Kp corresponding to the chosen timing error detector
% (TED). The gain is computed as the slope of the TED's S-curve evaluated
% for tau_e (timing error) around 0. The returned value assumes constant
% channel gain (e.g., after using automatic gain control) and unitary
% average symbol. When this assumption does not hold, the Kp returned by
% this function must be scaled appropriately.
%
% [ Kp ] = calcTedKp(TED, rollOff, method, plotSCurve, rcDelay)
%     TED        -> TED choice (MLTED, ELTED, ZCTED, GTED, or MMTED).
%     rollOff    -> Rolloff factor.
%     method     -> 'analytic' or 'simulated' (default: 'analytic').
%     plotSCurve -> Whether to plot the S-Curve (default: false).
%     rcDelay    -> Raised cosine pulse delay (default: 10).

if nargin < 3
    method = 'analytic';
end

if nargin < 4
    plotSCurve = false;
end

if nargin < 5
    rcDelay = 10;
end

if (strcmp(method, 'simulated'))
    [ normTauE, g ] = simSCurve(TED, rollOff, rcDelay);
else
    [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay);
end

% Oversampling factor adopted for the evaluation. Assume the evaluation is
% over "-L/2:L/2", namely over L+1 timing offset values.
L = length(g) - 1;

% TED Gain: slope around the origin (tau_e=0), estimated as follows:
delta_y = g(L/2 + 1) - g(L/2 - 1);
delta_x = normTauE(L/2 + 1) - normTauE(L/2 - 1);
Kp = delta_y / delta_x;

if (plotSCurve)
    method(1) = upper(method(1));
    figure
    plot(normTauE, g)
    xlabel("Normalized timing error $(\tau_e / T_s)$", 'interpreter', 'latex')
    ylabel("$g(\tau_e)$", 'interpreter', 'latex')
    title(sprintf("%s S-Curve for the %s (rolloff=%.2f)", ...
        method, TED, rollOff), 'interpreter', 'latex')
    grid on
end