function predictions = predict_PAR(noisy_signal, T, n_delta_t, p)
    % Periodic Autoregressive (PAR) Model Prediction
    % Parameters:
    % - noisy_signal: Input time series
    % - T: Period of the signal
    % - n_delta_t: Prediction horizon
    % - p: Order of the AR model
    %
    % Output:
    % - predictions: Predicted values

    n_points = length(noisy_signal);
    predictions = zeros(n_points, 1);
    ar_coefficients = zeros(T, p); % Coefficients for each phase

    % Estimate AR coefficients for each phase
    for phase = 1:T
        indices = phase:T:n_points; % Indices for the current phase
        if length(indices) > p
            % Construct lagged matrix
            X = zeros(length(indices) - p, p);
            y = noisy_signal(indices(p+1:end));
            for lag = 1:p
                X(:, lag) = noisy_signal(indices(p-lag+1:end-lag));
            end

            % Estimate AR coefficients for this phase
            if rank(X) == p
                ar_coefficients(phase, :) = (X \ y)';
            else
                warning('Rank deficient matrix for phase %d. Setting coefficients to zero.', phase);
                ar_coefficients(phase, :) = zeros(1, p);
            end
        else
            warning('Not enough data for phase %d. Setting coefficients to zero.', phase);
            ar_coefficients(phase, :) = zeros(1, p);
        end
    end

    % Generate predictions
    % for t = p + n_delta_t + 1:n_points
    %     phase = mod(t - 1, T) + 1; % Determine the phase
    %     lagged_values = noisy_signal(t-p:t-1); % Last p values
    %     predictions(t) = ar_coefficients(phase, :) * lagged_values;
    % end

    % Fill initial values with noisy_signal as fallback
    %predictions(1:p + n_delta_t) = noisy_signal(1:p + n_delta_t);
    for t = p + 1 : n_points - n_delta_t
    phase = mod(t - 1, T) + 1;
    temp_values = noisy_signal(t - p : t - 1); % Initial lagged values

    % Prédiction récursive pour l’horizon n_delta_t
    for step = 1 : n_delta_t
        % Calcul de la prédiction à ce pas
        next_value = ar_coefficients(phase, :) * temp_values;
        % Décalage des valeurs pour préparer le prochain pas
        temp_values = [temp_values(2:end); next_value];
    end

    % Stocker la prédiction finale à horizon n_delta_t
    predictions(t + n_delta_t) = next_value;
end

end
