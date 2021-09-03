function [ xI ] = symTimingLoop(TED, intpl, L, mfOut, dMfOut, K1, K2, const, Ksym, debug_s, debug_r)
% Symbol Timing Loop
% ---------------------
%
% Implements symbol timing recovery with a configurable timing error detector
% (TED) and interpolator, using a proportional-plus-integrator (PI) controller
% and a Modulo-1 counter to control the interpolator. The TED can be configured
% as a maximum-likelihood (ML) TED (ML-TED) or a zero-crossing TED (ZC-TED),
% whereas the interpolator can be a linear, polyphase, quadratic, or cubic
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

%% Interpolator Design

% Polyphase filter bank
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

% Quadratic and cubic interpolators
%
% Define the matrix bl_i with the Farrow coefficients to multiply the samples
% surrounding the desired interpolant. Each column of bl_i_mtx holds b_l(i) for
% a fixed l (exponent of mu(k) in 8.76) and for i (neighbor sample index) from
% -2 to 1. After the fliplr operations, the first column becomes the one
% associated with l=0 and the last with l=2. Each of those columns are filters
% to process the samples from x(mk-1) to x(mk+2), i.e., the sample before the
% basepoint index (x(mk-1)), the sample at the basepoint index (x(mk)), and two
% samples ahead (x(mk+1) and x(mk+2)). Before the flipud, the first row of bl_i
% would have the taps for i=-2, which would multiply x(mk+2). For convenience,
% however, the flipping ensures the first row has the taps for i=+1, which
% multiply x(mk-1). This order facilitates the dot product used later.
if (interpChoice == 2)
    % Farrow coefficients for alpha=0.5 (see Table 8.4.1)
    alpha = 0.5;
    bl_i_mtx = flipud(fliplr(...
        [+alpha, -alpha, 0; ...
        -alpha, (1 + alpha), 0; ...
        -alpha, (alpha - 1), 1; ...
        +alpha, -alpha     , 0]));
elseif (interpChoice == 3)
    % Table 8.4.2
    bl_i_mtx = flipud(fliplr(...
        [+1/6,    0, -1/6, 0; ...
         -1/2, +1/2,   +1, 0; ...
         +1/2, -1  , -1/2, 1; ...
         -1/6, +1/2, -1/3, 0]));
end

%% Timing Recovery Loop

% Constants
nSamples = length(mfOut);
nSymbols = ceil(nSamples / L);

% Preallocate
xI = zeros(nSymbols, 1); % Output interpolants
mu = zeros(nSymbols, 1); % Fractional symbol timing offset estimate
v  = zeros(nSamples, 1); % PI output
e  = zeros(nSamples, 1); % Error detected by the TED

% Initialize
k      = 0;     % interpolant/symbol index
strobe = 0;     % strobe signal
cnt    = 1;     % modulo-1 counter
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
        % Parallel interpolators
        %
        % NOTE: there are two parallel interpolators for the MF output and the
        % dMF output. However, the dMF is only required when using the ML-TED.
        switch (interpChoice)
            case 0 % Linear Interpolator (See Eq. 8.61)
                xI(k) = mu(k) * mfOut(m_k + 1) + (1 - mu(k)) * mfOut(m_k);
                xdotI = mu(k) * dMfOut(m_k + 1) + (1 - mu(k)) * dMfOut(m_k);
            case 1 % Polyphase interpolator
                % All subfilters filter the same samples (low-rate input
                % sequence)

                % Chose the output of one out of L polyphase subfilters.
                % Use mu(k) (the k-th fractional interval) to pick the
                % appropriate subfilter.
                chosenBranch = round((L-1) * mu(k)) + 1;

                % Interpolants (from the chosen subfilter):
                xI(k) = E(chosenBranch, :) * mfOutBuf;
                xdotI = E(chosenBranch, :) * dMfOutBuf;
            case 2 % Quadratic Interpolator
                % Recursive computation based on Eq. 8.77
                v_l = mfOut(m_k - 1 : m_k + 2).' * bl_i_mtx;
                xI(k) = (v_l(3) * mu(k) + v_l(2)) * mu(k) + v_l(1);
            case 3 % Cubic Interpolator
                % Recursive computation based on Eq. 8.78
                v_l = mfOut(m_k - 1 : m_k + 2).' * bl_i_mtx;
                xI(k) = ((v_l(4) * mu(k) + v_l(3)) * mu(k) + v_l(2)) * mu(k) + v_l(1);
        end

        % Timing Error Detector:
        a_hat_k = Ksym * slice(xI(k) / Ksym, M); % Data Symbol Estimate
        switch (TED)
            case 'MLTED' % Maximum Likelihood TED
                e(n) = real(a_hat_k) * real(xdotI) + ...
                       imag(a_hat_k) * imag(xdotI);
                % Note: the error could be alternatively computed by
                % "sign(xI)*xdotI". The difference is that this would scale
                % the S-Curve, as commented for Eq. 8.29.
            case 'ZCTED' % Zero-crossing TED
                if (k > 1)
                    % Previous Data Symbol Estimate
                    a_hat_prev = Ksym * slice(xI(k-1) / Ksym, M);

                    % Timing Error
                    e(n) = real(mfOut(m_k - L/2)) * (real(a_hat_prev) - real(a_hat_k)) + ...
                           imag(mfOut(m_k - L/2)) * (imag(a_hat_prev) - imag(a_hat_k));
                    % m_k - L/2 is the midpoint between the current and previous
                    % symbols (i.e., the current and previous basepoint indexes)
                else
                    e(n) = 0; % the ZC-TED needs at least two symbols to start
                end
        end

        % Real-time debugging scopes
        if (debug_r)
            step(hScope, xI(k))
            step(hTScopeCounter, mu(k));
        end
    else
        % Make the error null on the iterations without a strobe. This is
        % equivalent to upsampling the TED output.
        e(n) = 0;
    end

    % Loop Filter
    vp   = K1 * e(n);        % Proportional
    vi   = vi + (K2 * e(n)); % Integral
    v(n) = vp + vi;          % PI Output

    % Adjust the step used by the modulo-1 counter (see below Eq. 8.86)
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
        k = k + 1; % Update the interpolant Index
        m_k = n; % Basepoint index (the index **before** the underflow)
        mu(k) = cnt / W; % Equation (8.89)
    end

    % Next modulo-1 counter value:
    cnt = mod(cnt - W, 1);
end

% Trim the output vector
xI = xI(1:k);

%% Static Debug Plots
if (debug_s)
    figure
    plot(e)
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')

    figure
    plot(v)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(mu, '.')
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
