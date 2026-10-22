function [s, info] = scv(I, I_clr, d)
    % scv  Stochastic coefficient of variation of an irradiance series (Voyant et al., Renewable Energy 2026).
    %   s = scv(I, I_clr)      with the diffuse floor d = 0.2 of the article
    %   s = scv(I, I_clr, d)   with another diffuse floor, 0 <= d < 1
    %
    %   sCV = sqrt(2) / |d - 1| * sqrt( mean( (I - I_clr).^2 ) ) / sqrt( sigma_clr^2 + mu_clr^2 )   (Eq. 10)
    %
    % The numerator is the RMSE of the clear-sky model used as a forecast of I. The denominator is the
    % root mean square of the clear-sky series, since sigma_clr^2 + mu_clr^2 = mean(I_clr.^2) for the
    % population variance. sCV = 0 when I follows clear sky. sCV is normalised so that a sky switching at
    % random, half of the time, between clear sky and the diffuse floor d * I_clr has mean(sCV^2) = 1:
    % 1 is the reference of strong variability, not a pointwise maximum (I = d * I_clr at all times gives
    % sqrt(2)). The value is not clipped. Night samples (I = I_clr = 0) leave sCV unchanged.
    % Samples where I or I_clr is NaN are excluded; Inf values and an all-zero clear sky are errors.

    if nargin < 3 || isempty(d), d = 0.2; end
    if ~isscalar(d) || ~isfinite(d) || d < 0 || d >= 1
        error('scv:d', 'd must be a finite scalar in [0, 1).');
    end
    I = I(:); I_clr = I_clr(:);
    if numel(I) ~= numel(I_clr)
        error('scv:size', 'I and I_clr must have the same length.');
    end
    if any(isinf(I)) || any(isinf(I_clr))
        error('scv:inf', 'I and I_clr must not contain Inf values.');
    end
    ok = ~isnan(I) & ~isnan(I_clr);
    if ~any(ok)
        error('scv:empty', 'No sample where both I and I_clr are available.');
    end
    rms_clr = sqrt(mean(I_clr(ok) .^ 2));
    if rms_clr == 0
        error('scv:clearsky', 'The clear-sky series is zero on all available samples: sCV is undefined.');
    end
    rmse_cs = sqrt(mean((I(ok) - I_clr(ok)) .^ 2));
    s = sqrt(2) / abs(d - 1) * rmse_cs / rms_clr;
    info = struct('rmse_clear_sky', rmse_cs, 'rms_clear_sky', rms_clr, 'd', d, 'n', nnz(ok));
end
