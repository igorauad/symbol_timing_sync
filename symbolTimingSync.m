function [ xI ] = symbolTimingSync(TED, intpl, L, mfOut, dMfOut, K1, ...
    K2, const, Ksym, debug_s, debug_r)
% Symbol Timing Loop
% ---------------------
%
% Implements symbol timing recovery with a configurable timing error
% detector (TED) and interpolator, using a proportional-plus-integrator
% (PI) controller and a Modulo-1 counter to control the interpolator. The
% TED can be configured as a maximum-likelihood (ML) TED (ML-TED) or a
% zero-crossing TED (ZC-TED), whereas the interpolator can be a linear,
% polyphase, quadratic, or cubic interpolator.
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
%   [1] Michael Rice, Digital Communications - A Discrete-Time Approach.
%   New York: Prentice Hall, 2008.

if (nargin < 10)
    debug_s = 0;
    debug_r = 0;
end

% Interpolator Choice
interpChoice = intpl;

% Modulation order
M = numel(const);

% Make sure the MF and dMF output arguments are column vectors
if (size(mfOut, 1) == 1)
    mfOut = mfOut(:);
end

if (size(dMfOut, 1) == 1)
    dMfOut = dMfOut(:);
end

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
if (interpChoice == 0)
    % Fixed Parameters
    interpFactor  = 32;  % Interpolation factor
    Pintfilt      = 4;   % Neighbor samples weighted by the interp filter
    bandlimFactor = 0.5; % Bandlimitedness of the interpolated sequence

    % Design Polyphase Interpolator Filter Bank
    [ E ] = polyphaseFilterBank(L, interpFactor, Pintfilt, bandlimFactor);

    % To facilitate the convolution computation using inner products, flip
    % all subfilters (rows of E) from left to right.
    E = fliplr(E);
else
    polyBranch = []; % dummy variable
end

% Quadratic and cubic interpolators
%
% Define the matrix bl_i with the Farrow coefficients to multiply the
% samples surrounding the desired interpolant. Each column of b_mtx holds
% b_l(i) for a fixed l (exponent of mu(k) in 8.76) and for i (neighbor
% sample index) from -2 to 1. After the fliplr operations, the first column
% becomes the one associated with l=0 and the last with l=2. Each of those
% columns are filters to process the samples from x(mk-1) to x(mk+2), i.e.,
% the sample before the basepoint index (x(mk-1)), the sample at the
% basepoint index (x(mk)), and two samples ahead (x(mk+1) and x(mk+2)).
% Before the flipud, the first row of bl_i would have the taps for i=-2,
% which would multiply x(mk+2). For convenience, however, the flipping
% ensures the first row has the taps for i=+1, which multiply x(mk-1). This
% order facilitates the dot product used later.
if (interpChoice == 2)
    % Farrow coefficients for alpha=0.5 (see Table 8.4.1)
    alpha = 0.5;
    b_mtx = flipud(fliplr(...
        [+alpha,      -alpha, 0; ...
         -alpha, (1 + alpha), 0; ...
         -alpha, (alpha - 1), 1; ...
         +alpha,      -alpha, 0]));
elseif (interpChoice == 3)
    % Table 8.4.2
    b_mtx = flipud(fliplr(...
        [+1/6,    0, -1/6, 0; ...
         -1/2, +1/2,   +1, 0; ...
         +1/2,   -1, -1/2, 1; ...
         -1/6, +1/2, -1/3, 0]));
else
    b_mtx = []; % dummy variable
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
k      = 0; % interpolant/symbol index
strobe = 0; % strobe signal
cnt    = 1; % modulo-1 counter
vi     = 0; % PI filter integrator

% NOTE: by starting cnt with value 1, the first strobe is asserted when
% "n=L+1" and takes effect on iteration "n=L+2", while setting the
% basepoint index to "m_k=L+1". Furthermore, because the counter step is
% "W=1/L" before the first strobe and "cnt=0" when the first strobe is
% asserted, the first fractional interval estimate from Eq. (8.89) is
% "mu=0". Consequently, the first interpolant tends to be closer (or equal
% to) the sample at the basepoint index m_k. In other words, the first
% interpolant is "x(L+1)". This is not strictly necessary, but is important
% to understand, e.g., when evaluating the loop in unit tests. Also, this
% strategy is useful to ensure there is enough "memory" when the time comes
% to compute the first interpolant, as the interpolation equations use
% samples from the past. Lastly, note "cnt=1" only in the first iteration.
% In all other iterations, it is always within [0, 1).

% Due to the look-ahead scheme in the early-late TED, process
% "length(mfOut) - L" samples only. Also, when using the quadratic or cubic
% interpolators (which take "x(m_k + 1)" into account), leave one extra
% sample in the end. In all other cases, process all samples.
if (strcmp(TED, 'ELTED'))
    n_end = length(mfOut) - L;
elseif (intpl > 1)
    n_end = length(mfOut) - 1;
else
    n_end = length(mfOut);
end

for n = 1:n_end
    if strobe == 1
        % Interpolation
        if (interpChoice == 0)
            % When using a polyphase interpolator, update the polyphase
            % filter branch to be used next. Use mu (the k-th fractional
            % interval) to pick the appropriate subfilter.
            chosenPolyBranch = floor(L * mu(k)) + 1;
            polyBranch = E(chosenPolyBranch, :);
        end
        xI(k) = interpolate(interpChoice, ...
            mfOut, m_k, mu(k), b_mtx, polyBranch);
        xdotI = interpolate(interpChoice, ...
            dMfOut, m_k, mu(k), b_mtx, polyBranch);

        % Timing Error Detector:
        a_hat_k = Ksym * slice(xI(k) / Ksym, M); % Data Symbol Estimate
        switch (TED)
            case 'MLTED' % Maximum Likelihood TED
                % Decision-directed version of Eq. (8.98), i.e., Eq. (8.27)
                % adapted to complex symbols:
                e(n) = real(a_hat_k) * real(xdotI) + ...
                    imag(a_hat_k) * imag(xdotI);
            case 'ELTED' % Early-late TED
                % Early and late interpolants
                x_early = interpolate(interpChoice, ...
                    mfOut, m_k + L/2, mu(k), b_mtx, polyBranch);
                x_late = interpolate(interpChoice, ...
                    mfOut, m_k - L/2, mu(k), b_mtx, polyBranch);
                % Decision-directed version of (8.99), i.e., (8.34)
                % adapted to complex symbols:
                e(n) = real(a_hat_k) * (real(x_early) - real(x_late)) + ...
                    imag(a_hat_k) * (imag(x_early) - imag(x_late));
            case 'ZCTED' % Zero-crossing TED
                if (k > 1)
                    % Estimate of the previous data symbol
                    a_hat_prev = Ksym * slice(xI(k-1) / Ksym, M);

                    % Zero-crossing interpolant
                    %
                    % NOTE: "m(k) - L/2 + mu(k)" is the estimated instant
                    % of the midpoint between the current and previous
                    % symbols/interpolants, where the zero-crossing should
                    % be located when the loop converges.
                    x_zc = interpolate(interpChoice, ...
                        mfOut, m_k - L/2, mu(k), b_mtx, polyBranch);

                    % Decision-directed version of (8.100), i.e., (8.37)
                    % adapted to complex symbols:
                    e(n) = real(x_zc) * ...
                        (real(a_hat_prev) - real(a_hat_k)) + ...
                        imag(x_zc) * (imag(a_hat_prev) - imag(a_hat_k));
                else
                    e(n) = 0; % needs at least two symbols to start
                end
            case 'GTED' % Gardner TED
                % Zero-crossing interpolant, same as used by the ZCTED
                x_zc = interpolate(interpChoice, ...
                        mfOut, m_k - L/2, mu(k), b_mtx, polyBranch);

                % Equation (8.101):
                if (k > 1)
                    e(n) = real(x_zc) * (real(xI(k - 1)) - real(xI(k))) ...
                        + imag(x_zc) * (imag(xI(k - 1)) - imag(xI(k)));
                else
                    e(n) = 0; % needs at least two symbols to start
                end
            case 'MMTED' % Mueller and Müller TED
                if (k > 1)
                    % Estimate of the previous data symbol
                    a_hat_prev = Ksym * slice(xI(k-1) / Ksym, M);

                    % Decision-directed version of (8.102), i.e., (8.49)
                    % adapted to complex symbols:
                    e(n) = real(a_hat_prev) * real(xI(k)) - ...
                        real(a_hat_k) * real(xI(k - 1)) + ...
                        imag(a_hat_prev) * imag(xI(k)) - ...
                        imag(a_hat_k) * imag(xI(k - 1));
                else
                    e(n) = 0; % needs at least two symbols to start
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
    % NOTE: since e(n)=0 when strobe=0, the PI output can be simplified to
    % "v(n) = vi" on iterations without a strobe. It is only when strobe=1
    % that "vp != 0" and that vi (integrator output) changes. Importantly,
    % note the counter step W below changes briefly to "1/L + vp + vi" when
    % strobe=1 and, then, changes back to "1/L + vi" when strobe=0. In the
    % meantime, "vi" remains constant until the next strobe.

    % Adjust the step used by the modulo-1 counter (see below Eq. 8.86)
    W = 1/L + v(n);

    % Check whether the counter will underflow on the next cycle, i.e.,
    % whenever "cnt < W". When that happens, the strobe signal must
    % indicate the underflow occurrence and trigger updates on:
    %
    % - The basepoint index: set to the index right **before** the
    %   underflow. When strobe=1, it means an underflow will occur on the
    %   **next** cycle. Hence, the index before the underflow is exactly
    %   the current index.
    % - The estimate of the fractional symbol timing offset: the estimate
    %   is based on the counter value **before** the underflow (i.e., on
    %   the current cycle) and the current counter step, according to
    %   equation (8.89).
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
if (strobe) % ended on a strobe (before filling the k-th interpolant)
    xI = xI(1:k-1);
else
    xI = xI(1:k);
end

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

%% Interpolation
function [xI] = interpolate(method, x, m_k, mu, b_mtx, poly_h)
% [xI] = interpolate(interpChoice, x, mu) returns the interpolant xI
% obtained from the vector of samples x.
%
% Args:
%     method -> Interpolation method: polyphase (0), linear (1), quadratic
%               (2), or cubic (3).
%     x      -> Vector of samples based on which the interpolant shall be
%               computed, including the basepoint and surrounding samples.
%     m_k    -> Basepoint index, the index preceding the interpolant.
%     mu     -> Estimated fractional interval between the basepoint index
%               and the desired interpolant instant.
%     b_mtx  -> Matrix with the coefficients for the polynomial
%               interpolator used with method=2 or method=3.
%     poly_h -> Polyphase subfilter that should process the input samples
%               when using the polyphase interpolator (method=0).
    switch (method)
    case 0 % Polyphase interpolator
        N = length(poly_h);
        xI = poly_h * x((m_k - N + 1) : m_k);
    case 1 % Linear Interpolator (See Eq. 8.61)
        xI = mu * x(m_k + 1) + (1 - mu) * x(m_k);
    case 2 % Quadratic Interpolator
        % Recursive computation based on Eq. 8.77
        v_l = x(m_k - 1 : m_k + 2).' * b_mtx;
        xI = (v_l(3) * mu + v_l(2)) * mu + v_l(1);
    case 3 % Cubic Interpolator
        % Recursive computation based on Eq. 8.78
        v_l = x(m_k - 1 : m_k + 2).' * b_mtx;
        xI = ((v_l(4) * mu + v_l(3)) * ...
            mu + v_l(2)) * mu + v_l(1);
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
