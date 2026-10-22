function predictions = predict_P(noisy_signal, n_delta_t)
    % Simple Persistence Prediction
    % Parameters:
    % - noisy_signal: Input signal with noise
    % - n_delta_t: Prediction horizon
    
    % Validate inputs
    if n_delta_t < 0 || n_delta_t > length(noisy_signal)
        error('Invalid n_delta_t: must be between 0 and length(noisy_signal).');
    end
    
    % Compute predictions
    predictions = [zeros(n_delta_t, 1); noisy_signal(1:end - n_delta_t)];
end
