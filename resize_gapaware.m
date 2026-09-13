function B = resize_gapaware(A, ratio, minValidFraction)
%RESIZE_GAPAWARE Gap-aware downsampling/resizing of a 2-D array.
%
%   B = resize_gapaware(A, ratio)
%   B = resize_gapaware(A, ratio, minValidFraction)
%
% NaN pixels are treated as missing data and DO NOT contribute to the
% resampling. The result is a normalized convolution:
%
%          resize(A .* valid)
%     B = --------------------
%            resize(valid)
%
% INPUTS
%   A                2-D array, may contain NaN stripes/gaps
%   ratio            resize factor:
%                       1   = same dimensions
%                       0.5 = half rows and columns
%                       0.25 = quarter rows and columns
%
%   minValidFraction optional minimum fraction of valid pixels contributing
%                    to an output pixel. Default = 0.
%
%                    Examples:
%                       0    : keep pixel if ANY valid data contributed
%                       0.25 : require >=25% valid coverage
%                       0.5  : require >=50% valid coverage
%                       1    : require completely valid support
%
% OUTPUT
%   B                resized array; unsupported pixels remain NaN
%
% Requires Image Processing Toolbox.

    if nargin < 3
        minValidFraction = 0;
    end

    validateattributes(A, {'numeric'}, {'2d'});
    validateattributes(ratio, {'numeric'}, ...
        {'scalar','positive','finite'});
    validateattributes(minValidFraction, {'numeric'}, ...
        {'scalar','>=',0,'<=',1});

    A = double(A);

    % Validity/weight mask
    W = isfinite(A);

    % Replace missing values by zero ONLY for the weighted numerator.
    % They receive zero weight and therefore do not affect the result.
    X = A;
    X(~W) = 0;

    if ratio == 1
        B = A;
        return
    end

    % -------------------------------------------------------------
    % Gap-aware normalized area resampling
    % -------------------------------------------------------------
    %
    % Numerator = integrated signal from valid pixels
    num = imresize(X, ratio, 'box', ...
                   'Antialiasing', true);

    % Denominator = fraction of valid support
    den = imresize(double(W), ratio, 'box', ...
                   'Antialiasing', true);

    % Normalize by the actual amount of valid detector area
    B = num ./ den;

    % Pixels with insufficient valid support remain undefined
    B(den <= max(minValidFraction, eps)) = NaN;

end
