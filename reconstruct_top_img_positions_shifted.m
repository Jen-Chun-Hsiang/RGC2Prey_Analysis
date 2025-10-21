function [top_img_positions_shifted, top_img_disparity_positions_shifted, distances, disparity_deg] = reconstruct_top_img_positions_shifted(path, scaling_factors, interocular_dist, bounds)
    % Reconstruct top_img_positions_shifted and top_img_disparity_positions_shifted from path, scaling_factors, and interocular_dist.
    %
    % Args:
    %     path (matrix): Average path positions (time_steps, 2).
    %     scaling_factors (vector): Scaling factors per time step.
    %     interocular_dist (double): Interocular distance in cm.
    %     bounds (vector, optional): [x_min, x_max, y_min, y_max] for bounds checking. If empty, skip adjustment.
    %
    % Returns:
    %     top_img_positions_shifted (matrix): Reconstructed positions for left eye (time_steps, 2).
    %     top_img_disparity_positions_shifted (matrix): Reconstructed positions for right eye (time_steps, 2).
    
    % Step 1: Compute disparity (in degrees, then pixels)
    % Note: disparity_from_scaling_factor and convert_deg_to_pix must be defined or available
    [disparity_deg, distances] = disparity_from_scaling_factor(...
        scaling_factors, ...
        21, ...  % start_distance (hardcoded from code)
        4, ...   % end_distance (hardcoded from code)
        interocular_dist, ...
        [] ...   % fix_disparity (use computed disparity)
    );
    disparity_pix = convert_deg_to_pix(disparity_deg);  % Convert to pixels
    disparity_deg = disparity_deg(end-length(path)+1:end);  % Row vector for right eye adjustment
    distances = distances(end-length(path)+1:end);  % Row vector for right eye adjustment
    disparity_pix = disparity_pix(end-length(path)+1:end);  % Row vector for right eye adjustment
    disparity_pix = disparity_pix(:)/120;  % Ensure column vector
    % Step 2: Reconstruct positions (reverse the averaging)
    % top_img_positions_shifted ≈ path - [disparity/2, 0]
    top_img_positions_shifted = path;

    top_img_positions_shifted(:, 1) = top_img_positions_shifted(:, 1) - disparity_pix / 2;
    
    % top_img_disparity_positions_shifted ≈ path + [disparity/2, 0]
    top_img_disparity_positions_shifted = path;
    top_img_disparity_positions_shifted(:, 1) = top_img_disparity_positions_shifted(:, 1) + disparity_pix / 2;
    % Step 3: Apply bounds adjustment if bounds are provided (mimic adjust_trajectories)
    if ~isempty(bounds)
        x_min = bounds(1);
        x_max = bounds(2);
        y_min = bounds(3);
        y_max = bounds(4);
        % Clamp positions to bounds (simple approximation; adjust_trajectories may do more)
        top_img_positions_shifted(:, 1) = max(min(top_img_positions_shifted(:, 1), x_max), x_min);
        top_img_positions_shifted(:, 2) = max(min(top_img_positions_shifted(:, 2), y_max), y_min);
        top_img_disparity_positions_shifted(:, 1) = max(min(top_img_disparity_positions_shifted(:, 1), x_max), x_min);
        top_img_disparity_positions_shifted(:, 2) = max(min(top_img_disparity_positions_shifted(:, 2), y_max), y_min);
    end
end


function [disparities, distances] = disparity_from_scaling_factor(scaling_factors, start_distance, end_distance, iod_cm, fix_disparity)
    % disparity_from_scaling_factor: Maps scaling factors to distances and computes binocular disparities.
    %
    % Args:
    %     scaling_factors: Array of scaling factors (e.g., [1, 1.1, 1.1, 1.2, 1.4, ...])
    %     start_distance: Starting distance (e.g., 21 cm)
    %     end_distance: Ending distance (e.g., 4 cm)
    %     iod_cm: Interocular distance in cm
    %     fix_disparity: Optional fixed disparity value (scalar or empty). If provided, overrides computed disparities.
    %
    % Returns:
    %     disparities: Array of disparities in degrees
    %     distances: Array of corresponding distances in cm
    
    % Identify the range in scaling_factors
    min_sf = min(scaling_factors);
    max_sf = max(scaling_factors);
    
    % Edge case: if min_sf == max_sf, treat all distances as start_distance
    if min_sf == max_sf
        distances = start_distance * ones(size(scaling_factors));
    else
        % Map scaling factor to distance via linear interpolation
        sf = scaling_factors;
        distances = start_distance + (end_distance - start_distance) * (sf - min_sf) / (max_sf - min_sf);
    end
    
    if isempty(fix_disparity)
        % Compute disparities using binocular_disparity
        disparities = arrayfun(@(d) binocular_disparity(iod_cm, d), distances);
    else
        % Override with fixed disparity
        disparities = fix_disparity * ones(size(scaling_factors));
    end
end

function disparity_deg = binocular_disparity(iod_cm, distance_cm)
    % binocular_disparity: Calculate binocular disparity in degrees for an object at distance_cm
    % from the midpoint between two eyes separated by iod_cm.
    %
    % Args:
    %     iod_cm: Interocular distance in cm
    %     distance_cm: Distance to object in cm
    %
    % Returns:
    %     disparity_deg: Disparity in degrees
    
    % theta = 2 * atand( (IOD / 2) / distance )
    disparity_deg = 2 * atand((iod_cm / 2) / distance_cm);
end

function result = convert_deg_to_pix(x, deg2um, pix2um, scaling)
    % convert_deg_to_pix: Converts degrees to pixels using scaling factors.
    %
    % Args:
    %     x: Array of values in degrees
    %     deg2um: Conversion factor from degrees to micrometers (default: 32.5)
    %     pix2um: Conversion factor from pixels to micrometers (default: 4.375)
    %     scaling: Additional scaling factor (default: 0.54)
    %
    % Returns:
    %     result: Array of values in pixels
    
    if nargin < 2, deg2um = 32.5; end
    if nargin < 3, pix2um = 4.375; end
    if nargin < 4, scaling = 0.54; end
    
    result = x * deg2um / (pix2um / scaling);
end