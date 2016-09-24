rollOff = 0:0.001:1;
L       = 32;
delay   = 10;
tedType = 1;

%% TED Type
switch (tedType)
    case 0
        TED = 'MLTED'; % Maximum Likelihood
    case 1
        TED = 'ELTED' % Early-late
    case 2
        TED = 'ZCTED'
    case 3
        TED = 'GTED'
    case 4
        TED = 'MMTED' % Mueller-Muller
end

%% Iterate over excess bandwidths and compute Kp for each
Kp = zeros(length(rollOff),1);

for i = 1:length(rollOff)
    [ Kp(i) ] = getTedKp(TED, L, rollOff(i), delay);
end

%% Plot
figure
plot(rollOff, Kp)
grid on
xlabel('Excess Bandwidth')
ylabel('$K_p$', 'Interpreter', 'latex')