function predictions = predict_P_cyclic(noisy_signal, n_delta_t, T)
    % Cyclique Persistence Prediction
    % Parameters:
    % - noisy_signal: Input signal with noise
    % - n_delta_t: Prediction horizon
    % - T: Period of the signal
    
    % Validate inputs
    n_points = length(noisy_signal);
    if T <= 0 || T >= n_points
        error('Invalid T: must be a positive integer less than the signal length.');
    end
    if n_delta_t < 0 || n_delta_t >= n_points
        error('Invalid n_delta_t: must be between 0 and length(noisy_signal).');
    end
    
    % Initialize predictions
    predictions = zeros(n_points, 1);
    
    % Compute predictions with fallback
    for i = n_delta_t + 1:n_points
        if i - T > 0
            predictions(i) = noisy_signal(i - T);
        else
            predictions(i) = noisy_signal(i - n_delta_t); % Fallback to simple persistence
        end
    end
end
