function [F, info] = forecastability(I, I_clr, T, n, d)
    % forecastability  Forecastability F of an irradiance series at lag horizon n (Voyant et al., 2026).
    %   F = forecastability(I, I_clr, T, n)       T: time steps per day, n: horizon in time steps
    %   F = forecastability(I, I_clr, T, n, d)    diffuse floor d (default 0.2)
    %
    %   F = 1 + sCV * (rho_max - 1)                                                   (Eq. 13)
    %   rho_max = max over lags i in [n, T/2] of |rho(i)|, the autocorrelation of the residual
    %             eps = I - I_clr (Eq. 12), computed on the pairs where both values are available.
    %
    % F lies in [1 - sCV, 1]: it equals 1 - sCV when the residual has no memory (rho_max = 0) and 1 when
    % it is perfectly correlated (rho_max = 1). Higher F: easier to forecast at that horizon.

    if nargin < 5, d = []; end
    if n < 1 || n ~= round(n) || T < 2 || n > floor(T / 2)
        error('forecastability:lags', 'Need an integer 1 <= n <= T/2.');
    end
    [s, info] = scv(I, I_clr, d);
    e = I(:) - I_clr(:);
    if s == 0                                     % no residual at all: perfectly forecastable by clear sky
        F = 1; info.scv = 0; info.rho_max = NaN; info.lag_of_max = NaN; info.rho = []; info.lags = [];
        return
    end
    lags = n:floor(T / 2);
    rho = zeros(size(lags));
    for k = 1:numel(lags)
        a = e(1:end - lags(k)); b = e(1 + lags(k):end);
        ok = ~isnan(a) & ~isnan(b);
        rho(k) = pearson(a(ok), b(ok));
    end
    [rho_max, k] = max(abs(rho));
    F = 1 + s * (rho_max - 1);
    info.scv = s; info.rho_max = rho_max; info.lag_of_max = lags(k); info.rho = rho; info.lags = lags;
end

function r = pearson(a, b)
    a = a - mean(a); b = b - mean(b);
    r = sum(a .* b) / sqrt(sum(a .^ 2) * sum(b .^ 2));
end
