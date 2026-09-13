function B = resize_nanaware(A, ratio)
%RESIZE_NANAWARE NaN-aware resizing without NaN propagation.
%
% Missing pixels do not contribute to the resampling. The output is
% normalized by the amount of valid input data contributing to each pixel.

    arguments
        A (:,:) {mustBeNumeric}
        ratio (1,1) double {mustBePositive}
    end

    A = double(A);

    valid = isfinite(A);

    X = A;
    X(~valid) = 0;

    % Resize data and corresponding statistical weight
    num = imresize(X,            ratio, 'bilinear', 'Antialiasing', true);
    den = imresize(double(valid), ratio, 'bilinear', 'Antialiasing', true);

    % Normalize only by contributing valid pixels
    B = num ./ den;

    % No valid information at all
    B(den < 1e-12) = NaN;
end