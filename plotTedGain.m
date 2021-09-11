function [] = plotTedGain(TED)
% [] = plotTedGain(TED) plots the TED gain (Kp) evaluated for varying
% roll-off factor values.
%    TED -> MLTED, ELTED, ZCTED, GTED, or MMTED.

rollOff = 0:0.001:1;
Kp = zeros(length(rollOff),1);
for i = 1:length(rollOff)
    Kp(i) = calcTedKp(TED, rollOff(i));
end

figure
plot(rollOff, Kp)
grid on
xlabel('Roll-off factor', 'Interpreter', 'latex')
ylabel('$K_p$', 'Interpreter', 'latex')
title(sprintf("%s gain vs. roll-off", TED), 'interpreter', 'latex')

end