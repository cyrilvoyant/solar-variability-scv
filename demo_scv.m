% demo_scv  Teaching demonstration of sCV (variability) and F (forecastability).
%
% Part 1  Four typical skies on the same clear-sky cycle, from no variability to the reference case
%         sCV = 1 (half of the hours at clear sky, half at the diffuse floor d * I_clr).
% Part 2  Monte Carlo with the generator of the article: 100 series, random period T and noise memory;
%         sCV, F and the error of persistence and smart persistence at horizons 1, 4, 7 and 10 steps.
%
% Outputs (folder results/): skies_series.csv, skies_metrics.csv, montecarlo.csv.
% Figures: python make_figures.py. Base MATLAB only (R2021a or later).

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end
T = 24; d = 0.2;

%% Part 1 - four skies over 60 days
rng(3);
n_days = 60;
t = (1:T * n_days)';
h_day = mod(t - 1, T);                              % hour of the day
cs = 1000 * sin(pi * (h_day - 6) / 12) .* (h_day >= 6 & h_day <= 18);   % clear sky, 6 h to 18 h
sky.clear     = cs;
sky.light     = cs .* clip(movmean(1 - 0.08 * abs(randn(size(t))), 3), d, 1);
sky.broken    = cs .* clip(movmean(1 - 0.45 * abs(randn(size(t))), 2), d, 1);
sky.reference = cs .* ((rand(size(t)) < 0.5) * (1 - d) + d); % each hour: clear sky or diffuse floor, 50/50
names = fieldnames(sky);
M = table();
for i = 1:numel(names)
    x = sky.(names{i});
    s = scv(x, cs, d);
    F = forecastability(x, cs, T, 1, d);
    cv = std(x) / mean(x);
    M = [M; table(string(names{i}), cv, s, F, 'VariableNames', {'sky', 'CV', 'sCV', 'F_1h'})]; %#ok<AGROW>
end
disp(M);
writetable(M, fullfile(out, 'skies_metrics.csv'));
days = (10 * T + 1:13 * T)';
writetable(table(days - days(1), cs(days), sky.clear(days), sky.light(days), sky.broken(days), sky.reference(days), ...
    'VariableNames', {'hour', 'clear_sky', 'clear', 'light', 'broken', 'reference'}), fullfile(out, 'skies_series.csv'));

%% Part 2 - Monte Carlo, generator of the article (paper_scripts/boxplot.m)
% T uniform in [20, 100] (the article: [12, 100]; 20 keeps lags [n, T/2] non-empty up to n = 10),
% filter length uniform in [1, T], noise 1 + randn bounded to [0.2, 1.1], causal moving average,
% clear sky max(0, 1000 sin(2 pi t / T)), 500 cycles, horizons 1, 4, 7 and 10 steps.
rng(4);
n_runs = 100; horizons = [1 4 7 10];
R = [];
for run = 1:n_runs
    Tr = randi([20, 100]);
    L = randi([1, Tr]);
    tr = (1:500 * Tr)';
    w = clip(1 + randn(size(tr)), 0.2, 1.1);
    noise = filter(ones(L, 1) / L, 1, w);
    csr = max(0, 1000 * sin(2 * pi * tr / Tr));
    x = max(0, csr .* noise);
    s = scv(x, csr, d);
    for n = horizons
        F = forecastability(x, csr, Tr, n, d);
        k0 = max(Tr + 1, n + 1);                       % scored samples, as in the article (day and night)
        y = x(k0:end);
        p  = x(k0 - n:end - n);                        % persistence
        sp = max(0, x(k0 - n:end - n) + csr(k0:end) - csr(k0 - n:end - n)); % smart persistence (additive)
        R = [R; run, Tr, L, n, s, F, sqrt(mean((y - p) .^ 2)) / mean(y), sqrt(mean((y - sp) .^ 2)) / mean(y), ...
             std(y) / mean(y)]; %#ok<AGROW>
    end
end
MC = array2table(R, 'VariableNames', {'run', 'T', 'filter_length', 'horizon', 'sCV', 'F', 'nRMSE_P', 'nRMSE_SP', 'CV'});
writetable(MC, fullfile(out, 'montecarlo.csv'));
fprintf('\nMonte Carlo: %d series x %d horizons written to results/montecarlo.csv\n', n_runs, numel(horizons));

function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end
