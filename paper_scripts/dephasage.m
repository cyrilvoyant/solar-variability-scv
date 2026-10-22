%% MAIN SCRIPT: Robustness Analysis with Corrected Manual Interpolation
clear;
clc;
close all;

% Paramètres
T = 480;         % Période
n_cycles = 500; % Nombre de cycles
n_points = n_cycles * T; % Nombre total de points
n_delta_t = 20;  % Décalage temporel pour la corrélation cyclique
filter_length = 20; % Longueur du filtre de lissage
p = 1; % Ordre du modèle PAR(p)
d = 0.2; % Fraction diffuse

% Générer le signal
time = (1:n_points)';
white_noise = 1 + randn(n_points, 1); % Bruit gaussien centré autour de 1
white_noise(white_noise < 0.2) = 0.2; % Limite inférieure
white_noise(white_noise > 1.1) = 1.1; % Limite supérieure
correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
signal = 1000 * sin(2 * pi * time / T); % Signal clair sinusoïdal
signal = max(0, signal);
noisy_signal_base = signal .* correlated_noise;
noisy_signal_base = max(0, noisy_signal_base);

% Remplacer par le signal clair uniquement pour tester la périodicité (optionnel)
% noisy_signal_base = signal;

% Déphasages fins : fractions entre deux points
fine_shifts = 0:0.001:0.05; % Déphasages fins (entre 0 et 1 période complète)
results = struct('shift', fine_shifts, 'CV', [], 'sCV', [], 'F', [], ...
                 'sigma_delta_Kt', [], 'MALR', []);

for shift_fraction = fine_shifts
    % Appliquer un déphasage strictement périodique
    shifted_signal = circshift(noisy_signal_base, round(shift_fraction * T));
    

    % Calcul des métriques : CV, sCV, F
    start_index = T + 1;
    std_shifted_signal = std(shifted_signal(start_index:end));
    mean_shifted_signal = mean(shifted_signal(start_index:end));
    mean_signal = mean(signal(start_index:end));
    CV = std_shifted_signal / mean_shifted_signal;

    deviations = signal - shifted_signal; % Déviation stochastique
    RMSE = sqrt(mean(deviations.^2));
    alpha = sqrt(2) / ((1 - d) * sqrt((mean(signal))^2+(std(signal))^2));
    sCV = alpha * RMSE; % Normaliser pour éviter sCV > 1

    % Autocorrélation maximale
    n_lags = floor(T / 2);
    [acf, ~, ~] = autocorr(deviations, 'NumLags', n_lags);
    rho_max = max(acf(1 + n_delta_t:end));
    F = 1 - sCV + sCV * rho_max;

    % Calcul de Kt
    Kt = shifted_signal ./ signal;
    Kt(signal <= 0.00001) = 1;

    % Sigma(delta Kt)
    delta_Kt = Kt(2:end) - Kt(1:end-1);
    sigma_delta_Kt = std(delta_Kt, 'omitnan');

    % MALR
    %valid_indices = (Kt(1:end-1) > 0.1) & (Kt(2:end) > 0.001);
    log_ratios = log(Kt(2:end)) - log(Kt(1:end-1));
    log_ratios(log_ratios>50)=50;
    log_ratios(log_ratios<-50)=-50;
    MALR = mean(abs(log_ratios), 'omitnan');

    % Stocker les résultats
    results.CV = [results.CV, CV];
    results.sCV = [results.sCV, sCV];
    results.F = [results.F, F];
    results.sigma_delta_Kt = [results.sigma_delta_Kt, sigma_delta_Kt];
    results.MALR = [results.MALR, MALR];
end

% Lisser les métriques pour éviter les paliers
results.CV = movmean(results.CV, 4);
results.sCV = movmean(results.sCV, 4);
results.F = movmean(results.F, 4);
results.sigma_delta_Kt = movmean(results.sigma_delta_Kt, 4);
results.MALR = movmean(results.MALR, 4);

% Afficher les résultats
fprintf('Analyse de robustesse terminée.\n');
save('robustness_results_manual.mat', 'results', 'fine_shifts');

%% Tracé des résultats
figure;

% Tracé des métriques CV, sCV, et F
subplot(3, 1, 1);
plot(fine_shifts, results.CV, '-o', 'DisplayName', '$CV$', 'LineWidth', 1, 'MarkerSize', 4);
hold on;
plot(fine_shifts, results.sCV, '-s', 'DisplayName', '$sCV$', 'LineWidth', 1, 'MarkerSize', 4);
plot(fine_shifts, results.F, '-^', 'DisplayName', '$F$', 'LineWidth', 1, 'MarkerSize', 4);
xlabel('Fine Shift (fraction of T)', 'FontSize', 12);
ylabel('Metrics', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10);
grid on;
title('Variability and Forecastability vs. Fine Shift', 'FontSize', 14);

% Tracé de Sigma(Delta Kt)
subplot(3, 1, 2);
plot(fine_shifts, results.sigma_delta_Kt, '-o', 'DisplayName', '$\sigma(\Delta K_t)$', 'LineWidth', 1, 'MarkerSize', 4);
xlabel('Fine Shift (fraction of T)', 'FontSize', 12);
ylabel('$\sigma(\Delta K_t)$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Sigma(Delta Kt)', 'FontSize', 14);

% Tracé de MALR
subplot(3, 1, 3);
plot(fine_shifts, results.MALR, '-o', 'DisplayName', 'MALR', 'LineWidth', 1, 'MarkerSize', 4);
xlabel('Fine Shift (fraction of T)', 'FontSize', 12);
ylabel('MALR', 'FontSize', 12);
grid on;
title('Mean Absolute Logarithmic Return (MALR)', 'FontSize', 14);

%% super plot

% Tracé combiné des métriques
figure;
subplot(1,2,1)

% Tracé de chaque métrique
plot(fine_shifts, results.CV, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '$CV$');
hold on;
plot(fine_shifts, results.sCV, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '$sCV$');
plot(fine_shifts, results.F, '-^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '$F$');
plot(fine_shifts, results.sigma_delta_Kt, '-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '$\sigma(\Delta K_t)$');
plot(fine_shifts, results.MALR, '-x', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'MALR');

% Réglages esthétiques
xlabel('Fine Shift (fraction of T)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Metrics', 'FontSize', 14, 'Interpreter', 'latex');
title({'Combined Variability Metrics vs. Fine Shift', ...
       '\textit{Analyzing the robustness of variability metrics}'}, ...
      'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
lgd =legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest', 'Box', 'on');
set(lgd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8([255; 255; 204; 200]));
grid on;

% Réglages des axes
ax = gca;
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
xlim([0, max(fine_shifts)]);
ylim([0, max([results.CV, results.sCV, results.F, results.sigma_delta_Kt, results.MALR]) + 0.1]);

% Ajouter une boîte pour les paramètres
dim = [0.65, 0.7, 0.25, 0.2]; % Position ajustée

% Texte de l'annotation
str = sprintf('$\\mathbf{Parameters:}$\n$T$: %.0f\n$n\\Delta t$: %.0f\n$\\mathrm{filter\\_length}$: %.0f\n$p$: %.2f', ...
              T, n_delta_t, filter_length, p);

% Ajout de l'annotation
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', [1, 1, 0.9], ...
           'EdgeColor', [1, 1, 0.9], 'LineWidth', 1.5);




%% Export du graphique
exportgraphics(gcf, 'combined_variability_metrics_with_params.png', 'Resolution', 300);