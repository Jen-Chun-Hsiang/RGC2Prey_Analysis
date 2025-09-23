function [rms_error, error_len] = calculateFixedShiftRMSError(true_path, pred_path, shift_val, real_dim)
    % Calculate RMS error with a fixed shift
    % Inputs:
    %   true_path: true path coordinates [time_steps, 2] (x, y coordinates)
    %   pred_path: predicted path coordinates [time_steps, 2] (x, y coordinates)
    %   shift_val: fixed shift value (positive = pred leads true)
    %   real_dim: scaling factor [x_scale, y_scale]
    % Output:
    %   rms_error: mean RMS error across all valid time steps after shift
    
    % Scale paths to real dimensions first
    true_scaled = true_path .* reshape(real_dim, [1 2]);
    pred_scaled = pred_path .* reshape(real_dim, [1 2]);
    
    % Apply the fixed shift
    if shift_val > 0
        % Predicted signal leads - shift pred signal backwards (or true forward)
        true_shifted = true_scaled((shift_val+1):end, :);
        pred_shifted = pred_scaled(1:(end-shift_val), :);
    elseif shift_val < 0
        % Predicted signal lags - shift pred signal forwards (or true backward)
        shift_val = abs(shift_val);
        true_shifted = true_scaled(1:(end-shift_val), :);
        pred_shifted = pred_scaled((shift_val+1):end, :);
    else
        % No shift
        true_shifted = true_scaled;
        pred_shifted = pred_scaled;
    end
    
    % Calculate RMS error at each time step after shift
    if size(true_shifted, 1) > 0 && size(pred_shifted, 1) > 0
        error_per_timestep = sqrt(sum((true_shifted - pred_shifted).^2, 2));
        rms_error = error_per_timestep; % Return mean RMS error
    else
        rms_error = NaN; % Return NaN if no valid time steps
    end
    error_len = length(rms_error);
end