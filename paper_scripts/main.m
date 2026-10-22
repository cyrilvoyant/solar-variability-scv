% MAIN SCRIPT
clear;
clc;
close all;

% Paramètres
T = 24;         % Période
n_cycles = 500; % Nombre de cycles
n_points = n_cycles * T; % Nombre total de points
n_delta_t = 2;  % horizon
filter_length = 5; % Longueur du filtre de lissage
p=1; %ordre du modèle PAR(p) 1*T
d=0.2 % fraction diffuse

% Générer le signal
time = (1:n_points)';

% Génération du bruit borné entre 0.2 et 1.1
white_noise = 1 + randn(n_points, 1); % Bruit gaussien centré autour de 1
white_noise(white_noise < 0.2) = 0.2; % Limite inférieure
white_noise(white_noise > 1.1) = 1.1; % Limite supérieure
correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
signal = 1000 * sin(2 * pi * time / T); % signal est la composante sinusoïdale (clear-sky)
signal = max(0, signal); % On s'assure que signal est positif ou nul
noisy_signal = signal .* correlated_noise; % noisy_signal est le signal final bruité
noisy_signal = max(0, noisy_signal); % On s'assure que noisy_signal est positif ou nul

% Prédictions avec différents modèles
predictions_P = predict_P(noisy_signal, n_delta_t);
predictions_P_smart = predict_P_smart(noisy_signal, signal, n_delta_t);
predictions_PAR = predict_PAR(noisy_signal, T, n_delta_t,p);

% Gestion des valeurs négatives pour les prédictions
predictions_P = max(0, predictions_P);
predictions_P_smart = max(0, predictions_P_smart);
predictions_PAR = max(0, predictions_PAR);

% Calculer nRMSE pour chaque modèle
start_index = max([T+1, n_delta_t+1]);
nrmse = zeros(3, 1); % Pré-allocation pour stocker les nRMSE

% Calculer la VALEUR MOYENNE du signal bruité sur la plage considérée
mean_noisy_signal = mean(noisy_signal(start_index:end));
mean_signal = mean(signal(start_index:end));
std_signal = std(signal(start_index:end));

% Calculer les nRMSE en divisant par la VALEUR MOYENNE
nrmse(1) = sqrt(mean((noisy_signal(start_index:end) - predictions_P(start_index:end)).^2)) / mean_noisy_signal;
nrmse(2) = sqrt(mean((noisy_signal(start_index:end) - predictions_P_smart(start_index:end)).^2)) / mean_noisy_signal;
nrmse(3) = sqrt(mean((noisy_signal(start_index:end) - predictions_PAR(start_index:end)).^2)) / mean_noisy_signal;

% Afficher les résultats (AFFICHAGE UNIQUE)
fprintf('nRMSE pour Simple: %.4f\n', nrmse(1));
fprintf('nRMSE pour Smart: %.4f\n', nrmse(2));
fprintf('nRMSE pour PAR: %.4f\n', nrmse(3));

% Définition des indices pour afficher 4 périodes
n_periods_to_display = 4;
n_points_to_display = n_periods_to_display * T;

% Extraction des données pour 4 périodes
time_subset = time(115:115+n_points_to_display);
signal_subset = signal(115:115+n_points_to_display);
noisy_signal_subset = noisy_signal(115:115+n_points_to_display);
predictions_P_subset = predictions_P(115:115+n_points_to_display);
predictions_P_smart_subset = predictions_P_smart(115:115+n_points_to_display);
predictions_PAR_subset = predictions_PAR(115:115+n_points_to_display);

% Charger la colormap (parula)
cmap = lines(3); % Utilisation de 9 couleurs distinctes

% Marges ajustées pour plus d'espace vertical
vertical_margin = 0.08;

%% calcul des coefficient de variabilité et de forecastabilité

% Calcul du coefficient de variation (CV) du signal bruité
std_noisy_signal = std(noisy_signal(start_index:end));
cv_noisy_signal = std_noisy_signal / mean_noisy_signal;
% Afficher les résultats
fprintf('Coefficient de Variation (CV) : %.4f\n', cv_noisy_signal);

% Calcul de sCV
deviations = signal - noisy_signal; % Déviation stochastique
RMSE = sqrt(mean(deviations.^2)); % Root Mean Square Error (RMSE)
alpha = sqrt(2) / ((1 - d) * sqrt((mean_signal)^2+std_signal^2)); % Constante de normalisation
sCV = alpha * RMSE; % Stochastic Coefficient of Variation

% Calcul de l'autocorrélation totale avec autocorr
n_lags = floor(T / 2); % Nombre de décalages temporels
[acf, lags, bounds] = autocorr(deviations, 'NumLags', n_lags);
% Extraire l'autocorrélation maximale (ignorer lag 0)
rho_max = max(acf(1+n_delta_t:end)); % Ne pas prendre en compte l'autocorrélation au lag 0

% Calcul de F (Forecastability)
F = 1-sCV + sCV * rho_max;
fprintf('Stochastic Coefficient of Variation (sCV): %.4f\n', sCV);

% Calcul de Kt (indice clair) : noisy_signal / signal
Kt = noisy_signal ./ signal;
Kt(signal == 0) = 1; % Fixer Kt = 1 lorsque le dénominateur tend vers 0

% Calcul des différences Kt(t) - Kt(t-1) (delta Kt)
delta_Kt = Kt(2:end) - Kt(1:end-1);

% Calcul de sigma(delta Kt)
sigma_delta_Kt = std(delta_Kt, 'omitnan'); % Écart type des différences
fprintf('Sigma(delta Kt): %.4f\n', sigma_delta_Kt);

% Calcul de MALR (Mean Absolute Logarithmic Return)
log_ratios = log(Kt(2:end) ./ Kt(1:end-1));
log_ratios(Kt(1:end-1) == 0) = 0; % Fixer MALR=0 si le dénominateur est 0
MALR = mean(abs(log_ratios), 'omitnan'); % Moyenne des ratios logarithmiques
fprintf('MALR: %.4f\n', MALR);


%% Partie droite : Profil des prédictions
figure;
hold on;

% Signal bruité
plot(time_subset, noisy_signal_subset, '-', 'LineWidth', 2, 'DisplayName', 'Noisy Signal');

% Signal clair (clear-sky)
plot(time_subset, signal_subset,'--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'Clear-Sky Signal');

% Modèles de prédiction avec les noms corrigés
plot(time_subset, predictions_P_subset, '^', 'Color', cmap(1, :), 'LineWidth', 1, 'DisplayName', '$P$');
plot(time_subset, predictions_P_smart_subset, 'o', 'Color', cmap(2, :), 'LineWidth', 1, 'DisplayName', '$P_S$');
plot(time_subset, predictions_PAR_subset, '*', 'Color', cmap(3, :), 'LineWidth', 1, 'DisplayName', '$PAR$');


% Titres et légendes
title('Prediction Models Comparison Over 3 Periods', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (Points)', 'FontSize', 14);
ylabel('Signal Amplitude', 'FontSize', 14);

% Ajustement de la légende à l'intérieur
legend('Location', 'northeast', 'FontSize', 14, 'Box', 'on', 'Interpreter', 'latex'); % Utilisation de LaTeX pour les symboles

% Grille et limites des axes
grid on;
xlim([min(time_subset), max(time_subset)]);
ylim([min(noisy_signal_subset) , max(noisy_signal_subset) + 200]);

% Ajouter des annotations pour nRMSE et CV
dim1 = [0.4, 0.8, 0.2, 0.1]; % Position pour nRMSE (x, y, largeur, hauteur)
dim2 = [0.5, 0.8, 0.3, 0.1]; % Position pour CV (x, y, largeur, hauteur)
dim3 = [0.3, 0.8, 0.2, 0.1]; % Position pour nRMSE (x, y, largeur, hauteur)
str1 = sprintf('nRMSE:\n$P$: %.4f\n$P_S$: %.4f\n$PAR$: %.4f', nrmse(1), nrmse(2), nrmse(3));
str2 = sprintf('$CV$: %.4f\n$sCV$: %.4f\n$F$: %.4f', cv_noisy_signal, sCV, F);
str3 = sprintf('T: %.0f\n$n\\Delta t$: %.0f\n$filter\\_length$: %.0f\n$d$: %.1f', T, n_delta_t, filter_length, d);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k','Interpreter', 'latex','FontSize', 14);
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k','Interpreter', 'latex','FontSize', 14);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k','Interpreter', 'latex','FontSize', 14);

% Augmenter la taille des ticks (graduations) des axes
ax = gca; % Obtenir les axes courants
ax.FontSize = 14; % Taille de la police des ticks
ax.FontWeight = 'bold'; % Optionnel : Mettre en gras pour plus de lisibilité

hold off;





