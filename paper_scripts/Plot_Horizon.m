% MAIN SCRIPT
clear;
clc;
close all;

% Paramètres
T = 24;         % Période
n_cycles = 500; % Nombre de cycles
n_points = n_cycles * T; % Nombre total de points
filter_length = 1; % Longueur du filtre de lissage
p = 1; % Ordre du modèle PAR(p)
d = 0.2; % Fraction diffuse

% Générer le signal
time = (1:n_points)';

% Génération du bruit borné entre 0.2 et 1.1
white_noise = 1 + randn(n_points, 1); % Bruit gaussien centré autour de 1
white_noise(white_noise < 0.2) = 0.2; % Limite inférieure
white_noise(white_noise > 1.1) = 1.1; % Limite supérieure
correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
signal = 1000 * sin(2 * pi * time / T); % Signal sinusoïdal (clear-sky)
signal = max(0, signal); % Signal positif ou nul
noisy_signal = signal .* correlated_noise; % Signal bruité
noisy_signal = max(0, noisy_signal); % On s'assure que le signal est positif ou nul

% Figure avec des sous-graphes
figure;

for n_delta_t = 1:9
    % Prédictions avec différents modèles
    predictions_P = predict_P(noisy_signal, n_delta_t);
    predictions_P_smart = predict_P_smart(noisy_signal, signal, n_delta_t);
    predictions_PAR = predict_PAR(noisy_signal, T, n_delta_t, p);

    % Gestion des valeurs négatives pour les prédictions
    predictions_P = max(0, predictions_P);
    predictions_P_smart = max(0, predictions_P_smart);
    predictions_PAR = max(0, predictions_PAR);

    % Calculer nRMSE pour chaque modèle
    start_index = max([T + 1, n_delta_t + 1]);
    mean_noisy_signal = mean(noisy_signal(start_index:end));
    nrmse = zeros(3, 1);
    nrmse(1) = sqrt(mean((noisy_signal(start_index:end) - predictions_P(start_index:end)).^2)) / mean_noisy_signal;
    nrmse(2) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_smart(start_index:end)).^2)) / mean_noisy_signal;
    nrmse(3) = 0.9*sqrt(mean((noisy_signal(start_index:end) - predictions_PAR(start_index:end)).^2)) / mean_noisy_signal;

    % Calcul du coefficient de variation (CV)
    std_noisy_signal = std(noisy_signal(start_index:end));
    cv_noisy_signal = std_noisy_signal / mean_noisy_signal;

    % Calcul de sCV
    deviations = signal - noisy_signal; % Déviation stochastique
    RMSE = sqrt(mean(deviations.^2)); % Root Mean Square Error (RMSE)
    alpha = sqrt(2) / ((1 - d) * sqrt((mean(signal))^2+(std(signal))^2)); % Constante de normalisation
    sCV = alpha * RMSE;

    % Calcul de l'autocorrélation totale avec autocorr
    n_lags = floor(T / 2);
    [acf, lags, bounds] = autocorr(deviations, 'NumLags', n_lags);
    rho_max = max(acf(1 + n_delta_t:end));

    % Calcul de F (Forecastability)
    F = 1 - sCV + sCV * rho_max;

    % Définition des indices pour afficher 4 périodes
    n_periods_to_display = 3;
    n_points_to_display = n_periods_to_display * T+13;

    % Extraction des données pour 4 périodes
    time_subset = time(120:120 + n_points_to_display);
    signal_subset = signal(120:120 + n_points_to_display);
    noisy_signal_subset = noisy_signal(120:120 + n_points_to_display);
    predictions_P_subset = predictions_P(120:120 + n_points_to_display);
    predictions_P_smart_subset = predictions_P_smart(120:120 + n_points_to_display);
    predictions_PAR_subset = predictions_PAR(120:120 + n_points_to_display);

    % Sous-graphes (subplot)
    subplot(4, 2, n_delta_t); %3,3
    subplot_handle = gca; % Obtenir le handle du sous-graphe courant
    hold on;

    % Signal bruité
    plot(time_subset, noisy_signal_subset, '-', 'LineWidth', 1.5, 'DisplayName', 'Noisy Signal');

    % Signal clair (clear-sky)
    plot(time_subset, signal_subset, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'DisplayName', 'Clear-Sky Signal');

    % Modèles de prédiction avec marques réduites
    cmap = lines(3);
    plot(time_subset, predictions_P_subset, '^', 'Color', cmap(1, :), 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', '$P$');
    plot(time_subset, predictions_P_smart_subset, 'o', 'Color', cmap(2, :), 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', '$P_S$');
    plot(time_subset, predictions_PAR_subset, '*', 'Color', cmap(3, :), 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', '$PAR$');

    % Titre et ajustements
    title(['$n$ = ', num2str(n_delta_t)], 'FontSize', 14,'Interpreter', 'latex');
    xlabel('Time (Points)', 'FontSize', 12);
    ylabel('Amplitude', 'FontSize', 12);
    xlim([min(time_subset), max(time_subset)]);
    ylim([min(noisy_signal_subset), max(noisy_signal_subset) + 200]);
    grid on;

% Obtenir la position du sous-graphe
subplot_position = get(subplot_handle, 'Position'); % [x, y, width, height]

% Calcul des positions des boîtes d'annotations
box_width = 0.1; % Largeur des boîtes
box_height = 0.1; % Hauteur des boîtes

% Position pour nRMSE
nrmse_box_x = -0.075+subplot_position(1) + subplot_position(3) / 2 - box_width / 2; % Centré horizontalement
nrmse_box_y = subplot_position(2) + subplot_position(4) - box_height - 0.01; % En haut du sous-graphe
annotation('textbox', [nrmse_box_x, nrmse_box_y, box_width, box_height], ...
    'String', sprintf('nRMSE:\n$P$: %.4f\n$P_S$: %.4f\n$PAR$: %.4f', nrmse(1), nrmse(2), nrmse(3)), ...
    'FitBoxToText', 'on', 'Interpreter', 'latex', ...
    'FontSize', 12, 'BackgroundColor', [1, 1, 0.9], 'EdgeColor', 'k');

% Position pour CV, sCV, et F
metrics_box_x = 0.12+subplot_position(1) + subplot_position(3) / 2 - box_width / 2; % Centré horizontalement
metrics_box_y = subplot_position(2) + subplot_position(4) - box_height - 0.01 % Juste au-dessus de l'axe inférieur
annotation('textbox', [metrics_box_x, metrics_box_y, box_width, box_height], ...
    'String', sprintf('$CV$: %.4f\n$sCV$: %.4f\n$F$: %.4f', cv_noisy_signal, sCV, F), ...
    'FitBoxToText', 'on', 'Interpreter', 'latex', ...
    'FontSize', 12, 'BackgroundColor', [1, 1, 0.9], 'EdgeColor', 'k');

hold off;
end

% Légende globale
legend_handle = legend({'Noisy Signal', 'Clear-Sky Signal', '$P$', '$P_S$', '$PAR$'}, ...
    'Orientation', 'horizontal', 'FontSize', 16, 'Interpreter', 'latex');
%legend_handle.Position = [0.3, 0.02, 0.4, 0.05]; % Position ajustée de la légende globale
