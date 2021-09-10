function [] = plotTedGain(tedChoice, rcDelay)
% Plot the TED gain (Kp) while varying the rolloff factor
%
% Args:
%    tedChoice -> TED to analyze
%    rcDelay -> RC pulse delay (defaults to 10)

if nargin < 2
    rcDelay = 10;
end

rollOff = 0:0.001:1;
Kp = zeros(length(rollOff),1);
for i = 1:length(rollOff)
    Kp(i) = getTedKp(tedChoice, rollOff(i), rcDelay);
end

figure
plot(rollOff, Kp)
grid on
xlabel('Excess Bandwidth')
ylabel('$K_p$', 'Interpreter', 'latex')

end