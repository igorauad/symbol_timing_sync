function [ xx ] = symTimingLoop(TED, intpl, L, mfOut, dMfOut, K1, K2, const, Ksym, debug_s, debug_r)
% Symbol Timing Loop
% ---------------------
%
% Implements symbol timing recovery with a configurable timing error detector
% (TED) and interpolator, using a proportional-plus-integrator (PI) controller
% and a Modulo-1 counter to control the interpolator. The TED can be configured
% as a maximum-likelihood (ML) TED (ML-TED) or a zero-crossing TED (ZC-TED),
% whereas the interpolator can be a linear interpolator or a polyphase
% interpolator.
%
% Input Arguments:
% TED     -> Defines which TED should be used
% intpl   -> Defines which interpolator should be used
% L       -> Oversampling Factor
% mfOut   -> MF output sequence sampled at L samples/symbol
% dMfOut  -> Derivative MF output sequence sampled at L samples/symbol
%            (used only with the ML-TED)
% K1      -> Proportional Gain
% K2      -> Integrator Gain
% const   -> Symbol constellation
% Ksym    -> Symbol scaling factor that must be undone prior to slicing
% debug_s -> Show static debug plots after loop processing
% debug_r -> Open scope objects for run-time debugging of loop iterations
%
%
% References:
%   [1] Michael Rice, Digital Communications - A Discrete-Time Approach. New
%   York: Prentice Hall, 2008.

if (nargin < 8)
    debug_s = 0;
    debug_r = 0;
end

% Interpolator Choice
interpChoice = intpl;

% Modulation order
M = numel(const);

%% Optional System Objects for Step-by-step Debugging of the Loop

% Constellation Diagram
if (debug_r)
    hScope = comm.ConstellationDiagram(...
        'SymbolsToDisplaySource', 'Property',...
        'SamplesPerSymbol', 1, ...
        'MeasurementInterval', 256, ...
        'ReferenceConstellation', ...
        Ksym * const);
    hScope.XLimits = [-1.5 1.5]*max(real(const));
    hScope.YLimits = [-1.5 1.5]*max(imag(const));
end

% Time scope used to debug the fractional error
if (debug_r)
    hTScopeCounter = dsp.TimeScope(...
        'Title', 'Fractional Inverval', ...
        'NumInputPorts', 1, ...
        'ShowGrid', 1, ...
        'ShowLegend', 1, ...
        'BufferLength', 1e5, ...
        'TimeSpanOverrunAction', 'Wrap', ...
        'TimeSpan', 1e4, ...
        'TimeUnits', 'None', ...
        'YLimits', [-1 1]);
end

%% If applicable, design a polyphase filter bank
if (interpChoice == 1)
    % Fixed Parameters
    interpFactor  = 32;    % Interpolation factor
    Pintfilt      = 4;     % Neighbor samples weighted by the interp filter
    bandlimFactor = 0.5;   % Bandlimitedness of the interpolated sequence

    % Design Polyphase Interpolator Filter Bank
    [ E ] = polyphaseFilterBank(L, interpFactor, Pintfilt, bandlimFactor);

    % Buffers for MF and dMF output samples
    mfOutBuf  = zeros(size(E, 2), 1);
    dMfOutBuf = zeros(size(E, 2), 1);
end

%% Timing Recovery Loop

% Constants
nSamples = length(mfOut);
nSymbols = ceil(nSamples / L);

% Preallocate
v    = zeros(nSamples, 1);
xx   = zeros(nSymbols, 1);
ee   = zeros(nSymbols, 1);
mu_k = zeros(nSymbols, 1);

% Initialize
k      = 0;     % interpolant/symbol index
strobe = 0;     % strobe signal
mu     = 0;     % fractional symbol timing offset estimate
cnt    = 1;     % modulo-1 counter
W      = 1 / L; % modulo-1 counter's nominal step
vi     = 0;     % PI filter integrator

for n = 1:length(mfOut)
    % When using the polyphase interpolator keep track of the interpolator
    % input buffers at the sample rate (not only within strobes)
    if (interpChoice == 1) % Polyphase interpolator
        % Update interpolator input buffers
        mfOutBuf  = [mfOut(n); mfOutBuf(1:end-1)];
        dMfOutBuf = [dMfOut(n); dMfOutBuf(1:end-1)];
    end

    if strobe == 1
        % Update the interpolant Index
        k = k+1;

        % Parallel interpolators
        %
        % NOTE: there are two parallel interpolators for the MF output and the
        % dMF output. However, the dMF is only required when using the ML-TED.
        switch (interpChoice)
            case 0 % Linear Interpolator (See Eq. 8.61)
                xI    = mu * mfOut(m_k + 1) + (1 - mu) * mfOut(m_k);
                xdotI = mu * dMfOut(m_k + 1) + (1 - mu) * dMfOut(m_k);
            case 1 % Polyphase interpolator
                % All subfilters filter the same samples (low-rate input
                % sequence)

                % Chose the output of one out of L polyphase subfilters.
                % Use mu(k) (the k-th fractional interval) to pick the
                % appropriate subfilter.
                chosenBranch = round((L-1)*mu) + 1;

                % Interpolants (from the chosen subfilter):
                xI    = E(chosenBranch, :) * mfOutBuf;
                xdotI = E(chosenBranch, :) * dMfOutBuf;
        end

        % Timing Error Detector:
        a_hat_k = Ksym * slice(xI / Ksym, M); % Data Symbol Estimate
        switch (TED)
            case 'MLTED' % Maximum Likelihood TED
                e = real(a_hat_k) * real(xdotI) + ...
                    imag(a_hat_k) * imag(xdotI);
                % Note: the error could be alternatively computed by
                % "sign(xI)*xdotI". The difference is that this would scale
                % the S-Curve, as commented for Eq. 8.29.
            case 'ZCTED' % Zero-crossing TED
                if (k > 1)
                    % Previous Data Symbol Estimate
                    a_hat_prev = Ksym * slice(xx(k-1) / Ksym, M);
                    % Timing Error
                    e = real(mfOut(m_k - L/2)) * (real(a_hat_prev) - real(a_hat_k)) + ...
                        imag(mfOut(m_k - L/2)) * (imag(a_hat_prev) - imag(a_hat_k));
                    % m_k - L/2 is the midpoint between the current and previous
                    % symbols (i.e., the current and previous basepoint indexes)
                else
                    e = 0; % the ZC-TED needs at least two symbols to start
                end
        end

        % Save some metrics to plot them later:
        xx(k)   = xI; % interpolant
        ee(k)   = e;  % timing error
        mu_k(k) = mu; % fractional symbol timing offset estimate

        % Real-time debugging scopes
        if (debug_r)
            step(hScope, xI)
            step(hTScopeCounter, mu);
        end
    else
        % Make the error null on the iterations without a strobe. This is
        % equivalent to upsampling the TED output.
        e = 0;
    end

    % Loop Filter
    vp   = K1 * e;        % Proportional
    vi   = vi + (K2 * e); % Integral
    v(n) = vp + vi;       % PI Output

    % Adjust the step used by the modulo-1 counter
    W = 1/L + v(n);

    % Check whether the counter will underflow on the next cycle, i.e., whenever
    % "cnt < W". When that happens, the strobe signal must indicate the
    % underflow occurrence and trigger updates on:
    %
    % - The basepoint index: set to the index right **before** the underflow.
    %   When strobe=1, it means an underflow will occur on the **next** cycle.
    %   Hence, the index before the underflow is exactly the current index.
    % - The estimate of the fractional symbol timing offset: the estimate is
    %   based on the counter value **before** the underflow (i.e., on the
    %   current cycle) and the current counter step, according to Eq. (8.89).
    strobe = cnt < W;
    if (strobe)
        m_k = n; % Basepoint index (the index **before** the underflow)
        mu = cnt / W; % Equation (8.89)
    end

    % Next modulo-1 counter value:
    cnt = mod(cnt - W, 1);
end

% Trim the output vector
xx = xx(1:k);

%% Static Debug Plots
if (debug_s)
    figure
    plot(ee)
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')

    figure
    plot(v)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(mu_k, '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')
end

end

%% Function to map Rx symbols into constellation points
function [z] = slice(y, M)
if (isreal(y))
    % Move the real part of input signal; scale appropriately and round the
    % values to get ideal constellation index
    z_index = round( ((real(y) + (M-1)) ./ 2) );
    % clip the values that are outside the valid range
    z_index(z_index <= -1) = 0;
    z_index(z_index > (M-1)) = M-1;
    % Regenerate Symbol (slice)
    z = z_index*2 - (M-1);
else
    M_bar = sqrt(M);
    % Move the real part of input signal; scale appropriately and round the
    % values to get ideal constellation index
    z_index_re = round( ((real(y) + (M_bar - 1)) ./ 2) );
    % Move the imaginary part of input signal; scale appropriately and
    % round the values to get ideal constellation index
    z_index_im = round( ((imag(y) + (M_bar - 1)) ./ 2) );

    % clip the values that are outside the valid range
    z_index_re(z_index_re <= -1)       = 0;
    z_index_re(z_index_re > (M_bar-1)) = M_bar-1;
    z_index_im(z_index_im <= -1)       = 0;
    z_index_im(z_index_im > (M_bar-1)) = M_bar-1;

    % Regenerate Symbol (slice)
    z = (z_index_re*2 - (M_bar-1)) + 1j*(z_index_im*2 - (M_bar-1));
end
end
