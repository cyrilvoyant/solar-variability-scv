% MAIN SCRIPT FOR MONTE CARLO SIMULATION
clear;
clc;
close all;

% Parameters
n_MC = 50; % Reduced number of Monte Carlo simulations
n_delta_t_values = [1, 4, 7, 10]; % Different horizons
results = struct(); % Structure to store results

% Monte Carlo simulation
for mc = 1:n_MC
    % Random generation of T and filter_length
    T = randi([12, 100]); % T uniform between 12 and 100
    filter_length = randi([1, T]); % filter_length uniform between 1 and T
    p = 1; % Model order for PAR(p)
    d = 0.2; % Diffuse fraction

    % Generate signal
    n_cycles = 500;
    n_points = n_cycles * T;
    time = (1:n_points)';

    % Generate bounded noise
    white_noise = 1 + randn(n_points, 1);
    white_noise = max(min(white_noise, 1.1), 0.2); % Bound noise between 0.2 and 1.1
    correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
    signal = 1000 * sin(2 * pi * time / T); % Clear-sky signal
    signal = max(0, signal); % Ensure positivity
    noisy_signal = max(0, signal .* correlated_noise); % Noisy signal

    % Predictions
    for idx_n = 1:numel(n_delta_t_values)
        n_delta_t = n_delta_t_values(idx_n);
        predictions_P = predict_P(noisy_signal, n_delta_t);
        predictions_P_smart = predict_P_smart(noisy_signal, signal, n_delta_t);
        predictions_PAR = predict_PAR(noisy_signal, T, n_delta_t, p);

        % Calculate nRMSE
        start_index = max([T + 1, n_delta_t + 1]);
        mean_noisy_signal = mean(noisy_signal(start_index:end));
        nRMSE_P = sqrt(mean((noisy_signal(start_index:end) - predictions_P(start_index:end)).^2)) / mean_noisy_signal;
        nRMSE_P_smart = sqrt(mean((noisy_signal(start_index:end) - predictions_P_smart(start_index:end)).^2)) / mean_noisy_signal;
        nRMSE_PAR = sqrt(mean((noisy_signal(start_index:end) - predictions_PAR(start_index:end)).^2)) / mean_noisy_signal;

        % Variability metrics
        std_noisy_signal = std(noisy_signal(start_index:end));
        CV = std_noisy_signal / mean_noisy_signal;
        deviations = signal - noisy_signal;
        RMSE = sqrt(mean(deviations.^2));
        alpha = sqrt(2) / ((1 - d) * sqrt((mean(signal).^2 + std(signal).^2)));
        sCV = alpha * RMSE;

        % Forecastability metrics
        n_lags = floor(T / 2);
        [acf, ~] = autocorr(deviations, 'NumLags', n_lags);
        rho_max = max(acf(1 + n_delta_t:end)); % Ignore lag-0
        F = 1 - sCV + sCV * rho_max;

        % Delta Kt
        Kt = noisy_signal ./ signal;
        Kt(signal == 0) = 1; % Avoid division by zero
        delta_Kt = diff(Kt);
        sigma_delta_Kt = std(delta_Kt, 'omitnan');

        % MALR
        log_ratios = log(Kt(2:end) ./ Kt(1:end-1));
        log_ratios(Kt(1:end-1) == 0) = 0; % Avoid division by zero
        MALR = mean(abs(log_ratios), 'omitnan');

        % Store results in structure
        results(mc, idx_n).T = T;
        results(mc, idx_n).filter_length = filter_length;
        results(mc, idx_n).n_delta_t = n_delta_t;
        results(mc, idx_n).nRMSE = [nRMSE_P, nRMSE_P_smart, nRMSE_PAR];
        results(mc, idx_n).CV = CV;
        results(mc, idx_n).sCV = sCV;
        results(mc, idx_n).F = F;
        results(mc, idx_n).sigma_delta_Kt = sigma_delta_Kt;
        results(mc, idx_n).MALR = MALR;
    end
end

% Plotting: 4 groups (n_delta_t = 1, 4, 7, 10), each with 5 subplots
figure('Name', 'Monte Carlo Results: Scatter Plots', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 800]);

for idx_n = 1:numel(n_delta_t_values)
    n_delta_t = n_delta_t_values(idx_n);

    % Extract data
    CV_data = [];
    sCV_data = [];
    F_data = [];
    sigma_delta_Kt_data = [];
    MALR_data = [];
    nRMSE_data = [];

    for mc = 1:n_MC
        CV_data = [CV_data; results(mc, idx_n).CV];
        sCV_data = [sCV_data; results(mc, idx_n).sCV];
        F_data = [F_data; results(mc, idx_n).F];
        sigma_delta_Kt_data = [sigma_delta_Kt_data; results(mc, idx_n).sigma_delta_Kt];
        MALR_data = [MALR_data; results(mc, idx_n).MALR];
        nRMSE_data = [nRMSE_data; results(mc, idx_n).nRMSE];
    end

    % Ensure valid dimensions for plotting
    min_length = min([length(CV_data), length(sCV_data), length(F_data), ...
                      length(sigma_delta_Kt_data), length(MALR_data), size(nRMSE_data, 1)]);
    CV_data = CV_data(1:min_length);
    sCV_data = sCV_data(1:min_length);
    F_data = F_data(1:min_length);
    sigma_delta_Kt_data = sigma_delta_Kt_data(1:min_length);
    MALR_data = MALR_data(1:min_length);
    nRMSE_data = nRMSE_data(1:min_length, :);

    % Subplots for CV, sCV, delta_Kt, MALR, F
    subplot_idx_offset = (idx_n - 1) * 5;

    for sub_idx = 1:5
        subplot(4, 5, subplot_idx_offset + sub_idx);

        % Select X-data based on subplot index
        switch sub_idx
            case 1, x_data = CV_data; x_label = '$\mathrm{CV}$';
            case 2, x_data = sCV_data; x_label = '$\mathrm{sCV}$';
            case 3, x_data = sigma_delta_Kt_data; x_label = '$\Delta K_t$';
            case 4, x_data = MALR_data; x_label = '$\mathrm{MALR}$';
            case 5, x_data = F_data; x_label = '$\mathrm{F}$';
        end

        % Scatter plot
        scatter(x_data, nRMSE_data(:, 1), 10, 'b', 'filled'); hold on; % P
        scatter(x_data, nRMSE_data(:, 2), 10, 'r', 'filled'); % P_S
        scatter(x_data, nRMSE_data(:, 3), 10, 'g', 'filled'); % PAR

        % Polynomial trend (degree 2) for all predictors
        predictors = {'b', 'r', 'g'}; % Colors for P, P_S, PAR
        for i = 1:3
            coeffs = polyfit(x_data, nRMSE_data(:, i), 2); % Polynomial fit
            trend_line = polyval(coeffs, x_data);
            plot(sort(x_data), polyval(coeffs, sort(x_data)), predictors{i}, 'LineWidth', 1.5); % Trend line
            % R² Calculation
            SS_res = sum((nRMSE_data(:, i) - trend_line).^2);
            SS_tot = sum((nRMSE_data(:, i) - mean(nRMSE_data(:, i))).^2);
            R_squared = 1 - (SS_res / SS_tot);
            text(mean(x_data), max(nRMSE_data(:, i)) * (1 - 0.1 * i), ...
                 sprintf('$R^2 = %.2f$', R_squared), 'Interpreter', 'latex', ...
                 'FontSize', 12, 'Color', predictors{i});
        end

        % Titles and labels
        title(sprintf('$n = %d$: %s', n_delta_t, x_label), ...
            'Interpreter', 'latex', 'FontSize', 12);
        xlabel(x_label, 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\mathrm{nRMSE}$', 'Interpreter', 'latex', 'FontSize', 12);
        set(gca, 'FontSize', 12); % Adjust font size
    end
end

% Add overall title for the scatter plot figure
%sgtitle('Monte Carlo Results: Scatter Plots', 'Interpreter', 'latex', 'FontSize', 16);
