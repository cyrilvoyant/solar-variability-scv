function predictions = predict_P_smart(noisy_signal, clear_sky_signal, n_delta_t)
    % Smart Persistence Prediction using an additive adjustment
    % Parameters:
    % - noisy_signal: Input signal (time series)
    % - clear_sky_signal: Clear-sky model signal (time series)
    % - n_delta_t: Prediction horizon (corresponds to lag h)

    % Validate inputs
    n_points = length(noisy_signal);
    if length(clear_sky_signal) ~= n_points
        error('The clear-sky signal must have the same length as the noisy signal.');
    end
    if n_delta_t <= 0 || n_delta_t >= n_points
        error('Invalid prediction horizon n_delta_t. It must be between 1 and the length of the signal.');
    end

    % Compute the clear-sky index
    nonzero_indices = clear_sky_signal > 0;
    k_CS = noisy_signal ./ clear_sky_signal;
    k_CS(~nonzero_indices)=1;

    % Initialize predictions
    predictions = zeros(n_points, 1);

    % Generate predictions using the additive smart persistence model
    for t = n_delta_t + 1:n_points
        predictions(t) = noisy_signal(t - n_delta_t) + (clear_sky_signal(t) - clear_sky_signal(t - n_delta_t));
    end

    % Fill initial values with noisy_signal as fallback
    predictions(1:n_delta_t) = noisy_signal(1:n_delta_t);

    predictions = max(0, predictions); % Placement CORRECT ici
end