%% PathPredictDistExtraction.m - Modular Cricket Path Prediction Analysis
%
% This script analyzes cricket path prediction data and provides modular
% functions for comparing different experiments and noise levels.
%
% FEATURES:
% 1. Load and analyze individual datasets
% 2. Compare different experiments side-by-side at fixed noise level
% 3. Compare different noise levels for the same experiment
% 4. Configurable shared variables for easy modification
%
% USAGE EXAMPLES:
% 
% 1. Compare different experiments at fixed noise level:
%    plotComparison({'2025091506', '2025091502'}, 'blend', '0.016');
%
% 2. Compare different noise levels for same experiment:
%    plotComparison('2025091506', 'blend', {'0.0', '0.016', '0.256'});
%
% 3. Change background type:
%    bg_type = 'grass';  % or 'blend'
%
% NOTE: Either experiments OR noise levels can be multiple, not both.
%
%% Configuration - Shared Variables
mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats\';
fig_save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Summary\Illustrator\'; 
bg_type = 'blend'; % or 'grass'

exp_name = '2025092108';
noise_level = '0.016';
is_correct_object_zone = 1;
[all_paths_r, all_paths_pred_r, seqLen, is_simple_contrast, all_id_numbers, all_scaling_factors, all_path_cm] = loadDataset(mat_folder, exp_name, bg_type, noise_level);

load_mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\';  % folder to save output MAT file
coverage_mat_file = fullfile(load_mat_folder, 'processed_cover_radius.mat');
cover_radius = load(coverage_mat_file, 'file_index_list', 'processed_cover_radius');
cover_radius = [cover_radius.file_index_list(:) cover_radius.processed_cover_radius(:)];

%% Single trial visualization (original functionality)
trial_id = 45; % Change this to visualize different trials
%real_dim = [240 180]*4.375/0.54;
real_dim = [120 90]*4.375/0.54;
cm_dim_scale = 4.375/0.54; % Convert cm to um
figure; 
subplot(1, 2, 1)
plot(all_paths_r(trial_id, :, 1)*real_dim(1), all_paths_r(trial_id, :, 2)*real_dim(2), '-b.'); hold on;
plot(all_paths_pred_r(trial_id, :, 1)*real_dim(1), all_paths_pred_r(trial_id, :, 2)*real_dim(2), '-r.');
plot(all_path_cm(trial_id, :, 1)*cm_dim_scale, all_path_cm(trial_id, :, 2)*cm_dim_scale, '-g.');
xlabel('X Position'); ylabel('Y Position'); title('Cricket Location Prediction');
legend('True Path', 'Predicted Path');

subplot(1, 2, 2); hold on;
error_distance = sqrt(sum((all_paths_r(trial_id, :, :).*reshape(real_dim, [1 1 2]) - all_paths_pred_r(trial_id, :, :).*reshape(real_dim, [1 1 2])).^2, 3));

% Calculate cross-correlation for this trial
true_path_trial = squeeze(all_paths_r(trial_id, :, :)) .* reshape(real_dim, [1 2]);
pred_path_trial = squeeze(all_paths_pred_r(trial_id, :, :)) .* reshape(real_dim, [1 2]);
[max_xcorr, shift_value] = calculateMaxCrossCorrelation(true_path_trial, pred_path_trial);

% Calculate fixed-shift correlation for this trial
fixed_shift = -9;
fixed_corr = calculateFixedShiftCorrelation(true_path_trial, pred_path_trial, fixed_shift);

% Calculate fixed-shift RMS error for this trial
fixed_rms = calculateFixedShiftRMSError(squeeze(all_paths_r(trial_id, :, :)), squeeze(all_paths_pred_r(trial_id, :, :)), fixed_shift, real_dim);

shift_val = 0;
[velocity_true, velocity_pred] = calculateVelocity(squeeze(all_paths_r(trial_id, :, :)), squeeze(all_paths_pred_r(trial_id, :, :)), fixed_shift, real_dim);
plot(1:seqLen, error_distance, 'm-');
plot(1:numel(fixed_rms), fixed_rms, 'k-');
plot(1:numel(velocity_true), velocity_true, 'b-');
xlabel('Time steps'); ylabel('Error distance (um)'); 
title(sprintf('Single Trial Error\nXCorr: %.3f, Shift: %d, FixedCorr(%.0f): %.3f\nFixedRMS(%.0f): %.1f', max_xcorr, shift_value, fixed_shift, fixed_corr, fixed_shift, fixed_rms));
%%
is_visual_degree = 1;
vis_deg_to_cm = 32.5;
if is_visual_degree
    vis_scale = 1/vis_deg_to_cm; % Convert cm to visual degrees
else
    vis_scale = 1; % No scaling
end

baseColors = [ 120, 120, 120;   % red-ish
               255, 0, 255;   % green-ish
               0, 255, 255]/255; % blue-ish
true_path_trial = squeeze(all_paths_r(trial_id, :, :));
pred_path_trial = squeeze(all_paths_pred_r(trial_id, :, :));
pred_cm_path_trial = squeeze(all_path_cm(trial_id, :, :));
[fixed_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
[fixed_cm_rms, ~] = calculateFixedShiftRMSError(true_path_trial, pred_cm_path_trial* cm_dim_scale, fixed_shift, ones(1, 2));
[velocity_true, velocity_pred] = calculateVelocity(true_path_trial, pred_path_trial, fixed_shift, real_dim);
figure;
subplot(2, 1, 1)
t = (0:rms_len-1)/100; % Time vector in seconds

plot(t, fixed_cm_rms * vis_scale, '-', 'Color', baseColors(3,:), 'LineWidth', 1);
hold on;
plot(t, fixed_rms * vis_scale, '-', 'Color', baseColors(2,:), 'LineWidth', 1);
ylabel('Dist to cricket (deg)'); 
ylim([0 15])
yticks(0:6:12);
yticklabels(arrayfun(@(v) sprintf('%d', v), 0:6:12, 'UniformOutput', false));
xlabel('Time (s)'); 
xlim([0 t(end)]);
xticks(0:0.8:t(end));
xticklabels(arrayfun(@(v) sprintf('%.1f', v), 0:0.8:t(end), 'UniformOutput', false));
box off
legend({'Fixed CM RMS', 'Fixed RMS'});

subplot(2, 1, 2)
plot(t(1:end-1), velocity_true * 100 * vis_scale, '-', 'Color', baseColors(1,:), 'LineWidth', 1);
xlabel('Time (s)'); 
xlim([0 t(end)]);
xticks(0:0.8:t(end));
xticklabels(arrayfun(@(v) sprintf('%.1f', v), 0:0.8:t(end), 'UniformOutput', false));
ylabel('Velocity (deg/s)'); 
ylim([0 150])
yticks(0:60:120);
yticklabels(arrayfun(@(v) sprintf('%d', v), 0:60:120, 'UniformOutput', false));
box off
legend({'Velocity True'});

%%
% keyboard;
% save_file_name = fullfile(fig_save_folder, sprintf('PredictionErrorTrace_%s_%s_%d', exp_name, noise_level, trial_id));
% print(gcf, [save_file_name '.eps'], '-depsc', '-painters'); % EPS format
% print(gcf, [save_file_name '.png'], '-dpng', '-r300'); % PNG, 600 dpi
%% Mean error across all trials (original functionality)
% Calculate all error metrics separately for each trial
all_err_dist = zeros(size(all_paths_r, 1), seqLen);
all_xcorr = zeros(size(all_paths_r, 1), 1);
all_shifts = zeros(size(all_paths_r, 1), 1);
all_fixed_corr = zeros(size(all_paths_r, 1), 1);

% Fixed shift value for correlation coefficient calculation
fixed_shift = -9;

for i = 1:size(all_paths_r, 1)
    % Extract trial data
    true_path_trial = squeeze(all_paths_r(i, :, :));
    pred_path_trial = squeeze(all_paths_pred_r(i, :, :));
    pred_cm_path_trial = squeeze(all_path_cm(i, :, :));
    cut_off = acceptance_zone_radius(double(all_id_numbers(i)), all_scaling_factors(i, 50:end), cover_radius, 0);
    
    % Calculate RMS error using separate function
    all_err_dist(i, :) = calculateRMSError(true_path_trial, pred_path_trial, real_dim);
    if is_correct_object_zone
        all_err_dist(i, :) = max(0, all_err_dist(i, :)  - cut_off);
    end
    
    cut_off = acceptance_zone_radius(double(all_id_numbers(i)), all_scaling_factors(i, 50:end), cover_radius, fixed_shift);
    

    % Scale paths for correlation calculations
    true_path_scaled = true_path_trial .* reshape(real_dim, [1 2]);
    pred_path_scaled = pred_path_trial .* reshape(real_dim, [1 2]);
    
    % Calculate max cross-correlation and optimal shift
    [all_xcorr(i), all_shifts(i)] = calculateMaxCrossCorrelation(true_path_scaled, pred_path_scaled);
    
    % Calculate correlation coefficient with fixed shift
    all_fixed_corr(i) = calculateFixedShiftCorrelation(true_path_scaled, pred_path_scaled, fixed_shift);
    
    % Calculate RMS error with fixed shift
    [fixed_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
    if i == 1
        seqLen = rms_len; % Update seqLen based on actual RMS error length
        all_fixed_rms = zeros(size(all_paths_r, 1), seqLen);
        all_fixed_cm_rms = zeros(size(all_path_cm, 1), seqLen);
        all_fixed_rms_cutoff = zeros(size(all_paths_r, 1), seqLen);
    end
    all_fixed_rms(i, :) = fixed_rms;
    if is_correct_object_zone
        all_fixed_rms(i, :) = max(0, all_fixed_rms(i, :)  - cut_off);
    end
    all_fixed_rms_cutoff(i, :) = max(0, fixed_rms' - cut_off);


    [fixed_cm_rms, ~] = calculateFixedShiftRMSError(true_path_scaled, pred_cm_path_trial* cm_dim_scale, fixed_shift, ones(1, 2));
    all_fixed_cm_rms(i, :) = fixed_cm_rms;
    if is_correct_object_zone
        all_fixed_cm_rms(i, :) = max(0, all_fixed_cm_rms(i, :)  - cut_off);
    end

end

% Calculate statistics for each metric
mean_err = mean(all_err_dist, 1);
n_trials = size(all_err_dist, 1);
sem_err = std(all_err_dist, 0, 1) ./ sqrt(max(1, n_trials));

mean_xcorr = mean(all_xcorr);
sem_xcorr = std(all_xcorr) / sqrt(max(1, n_trials));

mean_shift = mean(all_shifts);
sem_shift = std(all_shifts) / sqrt(max(1, n_trials));

mean_fixed_corr = mean(all_fixed_corr);
sem_fixed_corr = std(all_fixed_corr) / sqrt(max(1, n_trials));

mean_fixed_rms = mean(all_fixed_rms, 1);
sem_fixed_rms = std(all_fixed_rms, [], 1) / sqrt(max(1, n_trials));

mean_fixed_rms_cutoff = mean(all_fixed_rms_cutoff, 1);
sem_fixed_rms_cutoff = std(all_fixed_rms_cutoff, [], 1) / sqrt(max(1, n_trials));

mean_fixed_cm_rms = mean(all_fixed_cm_rms, 1);
sem_fixed_cm_rms = std(all_fixed_cm_rms, [], 1) / sqrt(max(1, n_trials));

%%

colors = lines(4);
figure;
subplot(1, 3, 1); hold on;
x = (0:numel(mean_err)-1)/100;
shadePlot(x, mean_err, sem_err, colors(1, :));
x = (0:numel(mean_fixed_rms)-1)/100;
shadePlot(x, mean_fixed_rms, sem_fixed_rms, colors(2, :));
x = (0:numel(mean_fixed_rms_cutoff)-1)/100;
shadePlot(x, mean_fixed_rms_cutoff, sem_fixed_rms_cutoff, colors(3, :));
x = (0:numel(mean_fixed_cm_rms)-1)/100;
shadePlot(x, mean_fixed_cm_rms, sem_fixed_cm_rms, colors(4, :));

% ylim([0 1000])
xlabel('Time (s)');
ylabel('Error distance (um)');
title('Mean prediction error with SEM shaded');

subplot(1, 3, 2); hold on;
yr = mean(all_err_dist(is_simple_contrast, :), 1);
x = (0:numel(yr)-1)/100;
sr = std(all_err_dist(is_simple_contrast, :), 0, 1) ./ sqrt(max(1, numel(yr)));
shadePlot(x, yr, sr, colors(1, :));
yr = mean(all_err_dist(~is_simple_contrast, :), 1);
x = (0:numel(yr)-1)/100;
sr = std(all_err_dist(~is_simple_contrast, :), 0, 1) ./ sqrt(max(1, numel(yr)));
shadePlot(x, yr, sr, colors(2, :));
ylim([0 1000])
xlabel('Time (s)');
ylabel('Error distance (um)');
title('Background difference in prediction error with SEM shaded');

subplot(1, 3, 3); hold on;
yr = mean(all_fixed_rms(is_simple_contrast, :), 1);
x = (0:numel(yr)-1)/100;
sr = std(all_fixed_rms(is_simple_contrast, :), 0, 1) ./ sqrt(max(1, numel(yr)));
shadePlot(x, yr, sr, colors(1, :));
yr = mean(all_fixed_rms(~is_simple_contrast, :), 1);
x = (0:numel(yr)-1)/100;
sr = std(all_fixed_rms(~is_simple_contrast, :), 0, 1) ./ sqrt(max(1, numel(yr)));
shadePlot(x, yr, sr, colors(2, :));
ylim([0 1000])
xlabel('Time (s)');
ylabel('Error distance (um)');
title('Background difference in fixed-shift RMS error with SEM shaded');
%%
return


subplot(1, 2, 2)
plotWithSEM(x, mean_xcorr, sem_xcorr, 'r');
ylim([0 1000])
xlabel('Time (s)');
ylabel('Cross-Correlation');
title(sprintf('Mean Cross-Correlation with SEM shaded\nXCorr: %.3f±%.3f, Shift: %.1f±%.1f\nFixed Corr(-5): %.3f±%.3f', mean_xcorr, sem_xcorr, mean_shift, sem_shift, mean_fixed_corr, sem_fixed_corr));
legend({'Mean XCorr'});
grid on;

%% COMPARISON EXAMPLES - Uncomment to use
% Compare different experiments at fixed noise level
plotComparison({'2025091502', '2025091506'}, 'blend', '0.256');
% plotComparison({'2025091501', '2025091503'}, 'blend', '0.256');

% Compare different noise levels for the same experiment
% plotComparison('2025091503', 'blend', {'0.016', '0.064', '0.256'});
title('Mean prediction error with SEM shaded');
legend([hMean], {'Mean error'});
grid on;
hold off;

%% ========== FUNCTIONS ==========

function [all_paths_r, all_paths_pred_r, seqLen, is_simple_contrast, all_id_numbers, all_scaling_factors, all_path_cm] = loadDataset(mat_folder, exp_name, bg_type, noise_level)
    % Load and process a single dataset
    % Inputs:
    %   mat_folder: path to the folder containing .mat files
    %   exp_name: experiment name (e.g., '2025091506')
    %   bg_type: background type ('blend' or 'grass')
    %   noise_level: noise level as string (e.g., '0.016', '0.0', '0.256')
    % Outputs:
    %   all_paths_r: reshaped true paths
    %   all_paths_pred_r: reshaped predicted paths
    %   seqLen: sequence length
    
    filename = sprintf('%s_cricket_%s_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                      exp_name, bg_type, noise_level);
    filepath = fullfile(mat_folder, filename);
    
    if ~exist(filepath, 'file')
        error('File not found: %s', filepath);
    end
    
    fprintf('Loading: %s\n', filename);
    data = load(filepath);
    [all_paths_r, seqLen] = reshapeAllPaths(data.all_paths);
    all_path_cm = reshapeAllPaths(double(data.all_path_cm));
    all_paths_pred_r = squeeze(data.all_paths_pred);
    is_simple_contrast = cellfun(@(x) contains(x, 'gray_image'), data.all_bg_file);
    all_id_numbers = data.all_id_numbers;
    all_scaling_factors = data.all_scaling_factors;
    
    assert(isequal(size(all_paths_r), size(all_paths_pred_r)), ...
           'Reshaped paths and predicted paths must have the same dimensions');
end

function plotComparison(exp_names, bg_type, noise_levels)
    % Plot comparison of multiple datasets
    % Usage patterns:
    %   1. Compare multiple experiments at fixed noise level:
    %      plotComparison({'2025091506', '2025091502'}, 'blend', '0.016')
    %   2. Compare multiple noise levels for fixed experiment:
    %      plotComparison('2025091506', 'blend', {'0.0', '0.016', '0.256'})
    % 
    % Inputs:
    %   exp_names: single exp name (string) OR cell array of exp names
    %   bg_type: background type ('blend' or 'grass')
    %   noise_levels: single noise level (string) OR cell array of noise levels
    %   
    % Note: Either exp_names OR noise_levels should be multiple, not both
    
    mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats\';
    real_dim = [220 140]*4.375/0.54;
    
    % Define colors for different comparisons
    
    
    % Convert single strings to cell arrays for uniform handling
    if ischar(exp_names)
        exp_names = {exp_names};
    end
    if ischar(noise_levels)
        noise_levels = {noise_levels};
    end
    
    % Validate input: only one dimension should be multiple
    if length(exp_names) > 1 && length(noise_levels) > 1
        error('Either exp_names OR noise_levels should be multiple, not both');
    end
    
    % Determine comparison type
    if length(exp_names) > 1
        comparison_type = 'experiment';
        n_comparisons = length(exp_names);
    else
        comparison_type = 'noise';
        n_comparisons = length(noise_levels);
    end
    %colors = {'r', 'b', 'g', 'm', 'c', 'k', [0.8 0.5 0.2], [0.5 0.8 0.2]};

    colors = lines(n_comparisons); % Use distinct colors from colormap
    figure;
    hold on;
    legend_labels = {};
    xcorr_info = {}; % Store cross-correlation info for display
    
    for i = 1:n_comparisons
        if strcmp(comparison_type, 'experiment')
            % Multiple experiments, single noise level
            exp_name = exp_names{i};
            noise_level = noise_levels{1};
            label = sprintf('%s (noise=%s)', exp_name, noise_level);
        else
            % Single experiment, multiple noise levels
            exp_name = exp_names{1};
            noise_level = noise_levels{i};
            label = sprintf('%s (noise=%s)', exp_name, noise_level);
        end
        
        [all_paths_r, all_paths_pred_r, seqLen] = loadDataset(mat_folder, exp_name, bg_type, noise_level);
        
        % Calculate all error metrics separately for each trial
        all_err_dist = zeros(size(all_paths_r, 1), seqLen);
        all_xcorr = zeros(size(all_paths_r, 1), 1);
        all_shifts = zeros(size(all_paths_r, 1), 1);
        all_fixed_corr = zeros(size(all_paths_r, 1), 1);
        all_fixed_rms = zeros(size(all_paths_r, 1), 1);
        
        % Fixed shift value for correlation coefficient calculation
        fixed_shift = -5;
        
        for j = 1:size(all_paths_r, 1)
            % Extract trial data
            true_path_trial = squeeze(all_paths_r(j, :, :));
            pred_path_trial = squeeze(all_paths_pred_r(j, :, :));
            
            % Calculate RMS error using separate function
            all_err_dist(j, :) = calculateRMSError(true_path_trial, pred_path_trial, real_dim);
            
            % Scale paths for correlation calculations
            true_path_scaled = true_path_trial .* reshape(real_dim, [1 2]);
            pred_path_scaled = pred_path_trial .* reshape(real_dim, [1 2]);
            
            % Calculate max cross-correlation and optimal shift
            [all_xcorr(j), all_shifts(j)] = calculateMaxCrossCorrelation(true_path_scaled, pred_path_scaled);
            
            % Calculate correlation coefficient with fixed shift
            all_fixed_corr(j) = calculateFixedShiftCorrelation(true_path_scaled, pred_path_scaled, fixed_shift);
            
            % Calculate RMS error with fixed shift
            all_fixed_rms(j) = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
        end
        
        % Calculate statistics for each metric
        mean_err = mean(all_err_dist, 1);
        n_trials = size(all_err_dist, 1);
        sem_err = std(all_err_dist, 0, 1) ./ sqrt(max(1, n_trials));
        
        mean_xcorr = mean(all_xcorr);
        sem_xcorr = std(all_xcorr) / sqrt(max(1, n_trials));
        
        mean_shift = mean(all_shifts);
        sem_shift = std(all_shifts) / sqrt(max(1, n_trials));
        
        mean_fixed_corr = mean(all_fixed_corr);
        sem_fixed_corr = std(all_fixed_corr) / sqrt(max(1, n_trials));
        
        mean_fixed_rms = mean(all_fixed_rms);
        sem_fixed_rms = std(all_fixed_rms) / sqrt(max(1, n_trials));
        
        x = (0:seqLen-1)/100;
        
        % plotWithSEM(x, mean_err, sem_err, colors{mod(i-1, length(colors))+1});
        plotWithSEM(x, mean_err, sem_err, colors(i, :));
        legend_labels{end+1} = sprintf('%s (XCorr: %.3f±%.3f, FixedRMS: %.1f±%.1f)', label, mean_xcorr, sem_xcorr, mean_fixed_rms, sem_fixed_rms);
        xcorr_info{end+1} = sprintf('%s: XCorr = %.3f ± %.3f, Shift = %.1f ± %.1f, FixedCorr(-5) = %.3f ± %.3f, FixedRMS(-5) = %.1f ± %.1f', label, mean_xcorr, sem_xcorr, mean_shift, sem_shift, mean_fixed_corr, sem_fixed_corr, mean_fixed_rms, sem_fixed_rms);
    end
    
    ylim([0 1000]);
    xlabel('Time (s)');
    ylabel('Error distance (um)');
    if strcmp(comparison_type, 'experiment')
        title(sprintf('Experiment Comparison (bg=%s, noise=%s)', bg_type, noise_levels{1}));
    else
        title(sprintf('Noise Level Comparison (bg=%s, exp=%s)', bg_type, exp_names{1}));
    end
    legend(legend_labels, 'Location', 'best');
    
    % Display cross-correlation information in command window
    fprintf('\nCross-Correlation Results:\n');
    for i = 1:length(xcorr_info)
        fprintf('%s\n', xcorr_info{i});
    end
    fprintf('\n');
    
    % grid on;
    hold off;
end

function plotWithSEM(x, mean_err, sem_err, color)
    % Plot mean with SEM shading
    ci_upper = mean_err + sem_err;
    ci_lower = mean_err - sem_err;
    
    % Create shaded region for SEM with transparency
    if ischar(color)
        % Convert character color to RGB for shading
        rgb_color = getColorRGB(color);
    else
        rgb_color = color;
    end
    
    % Create fill object (stored for potential future use)
    fill([x fliplr(x)], [ci_upper fliplr(ci_lower)], rgb_color, ...
         'LineStyle', 'none', 'FaceAlpha', 0.3);
    
    % Plot mean line on top
    plot(x, mean_err, 'Color', color, 'LineWidth', 1.5);
end

function rgb = getColorRGB(color_char)
    % Convert character color codes to RGB values
    switch color_char
        case 'r', rgb = [1, 0, 0];
        case 'g', rgb = [0, 1, 0];
        case 'b', rgb = [0, 0, 1];
        case 'c', rgb = [0, 1, 1];
        case 'm', rgb = [1, 0, 1];
        case 'y', rgb = [1, 1, 0];
        case 'k', rgb = [0, 0, 0];
        otherwise, rgb = [0.5, 0.5, 0.5];
    end
end

function [max_xcorr, shift_value] = calculateMaxCrossCorrelation(true_path, pred_path)
    % Calculate maximum cross-correlation between true and predicted paths
    % Inputs:
    %   true_path: true path coordinates [time_steps, 2] (x, y coordinates)
    %   pred_path: predicted path coordinates [time_steps, 2] (x, y coordinates)
    % Outputs:
    %   max_xcorr: maximum cross-correlation coefficient (normalized)
    %   shift_value: lag/shift at which maximum correlation occurs (positive = pred leads true)
    
    % Concatenate x and y coordinates to create 1D signals
    true_signal = [true_path(:, 1); true_path(:, 2)];
    pred_signal = [pred_path(:, 1); pred_path(:, 2)];
    
    % Remove mean to center the signals (important for correlation coefficient)
    true_signal = true_signal - mean(true_signal);
    pred_signal = pred_signal - mean(pred_signal);
    
    % Calculate cross-correlation with limited lag range (±20 points)
    max_lag = 20;
    [xcorr_vals, lags] = xcorr(true_signal, pred_signal, max_lag, 'normalized');
    
    % Find maximum cross-correlation value and its position
    [max_xcorr, max_idx] = max(xcorr_vals);
    shift_value = lags(max_idx);
end

function rms_error = calculateRMSError(true_path, pred_path, real_dim)
    % Calculate RMS error between true and predicted paths
    % Inputs:
    %   true_path: true path coordinates [time_steps, 2] (x, y coordinates)
    %   pred_path: predicted path coordinates [time_steps, 2] (x, y coordinates)
    %   real_dim: scaling factor [x_scale, y_scale]
    % Output:
    %   rms_error: RMS error array [time_steps]
    
    % Scale paths to real dimensions
    true_scaled = true_path .* reshape(real_dim, [1 2]);
    pred_scaled = pred_path .* reshape(real_dim, [1 2]);
    
    % Calculate RMS error at each time step
    rms_error = sqrt(sum((true_scaled - pred_scaled).^2, 2));
end

function corr_coeff = calculateFixedShiftCorrelation(true_path, pred_path, shift_val)
    % Calculate correlation coefficient with a fixed shift
    % Inputs:
    %   true_path: true path coordinates [time_steps, 2] (x, y coordinates)
    %   pred_path: predicted path coordinates [time_steps, 2] (x, y coordinates)
    %   shift_val: fixed shift value (positive = pred leads true)
    % Output:
    %   corr_coeff: correlation coefficient at the specified shift
    
    % Concatenate x and y coordinates to create 1D signals
    true_signal = [true_path(:, 1); true_path(:, 2)];
    pred_signal = [pred_path(:, 1); pred_path(:, 2)];
    
    % Apply the fixed shift
    if shift_val > 0
        % Predicted signal leads - shift pred signal backwards (or true forward)
        true_shifted = true_signal((shift_val+1):end);
        pred_shifted = pred_signal(1:(end-shift_val));
    elseif shift_val < 0
        % Predicted signal lags - shift pred signal forwards (or true backward)
        shift_val = abs(shift_val);
        true_shifted = true_signal(1:(end-shift_val));
        pred_shifted = pred_signal((shift_val+1):end);
    else
        % No shift
        true_shifted = true_signal;
        pred_shifted = pred_signal;
    end
    
    % Calculate correlation coefficient
    if length(true_shifted) > 1 && length(pred_shifted) > 1
        corr_matrix = corrcoef(true_shifted, pred_shifted);
        corr_coeff = corr_matrix(1, 2);
        % Handle NaN case (when signals are constant)
        if isnan(corr_coeff)
            corr_coeff = 0;
        end
    else
        corr_coeff = 0;
    end
end



% ...existing code...
function [velocity_true, velocity_pred, error_len] = calculateVelocity(true_path, pred_path, shift_val, real_dim)
    % calculateVelocity - Compute velocity (step) vectors after optional fixed shift
    % Inputs:
    %   true_path: [T x 2] true coordinates
    %   pred_path: [T x 2] predicted coordinates
    %   shift_val: integer shift (positive = pred leads true)
    %   real_dim: [1x2] scaling factors for x and y
    % Outputs:
    %   rms_error: [(T_shifted-1) x 1] RMS error per timestep comparing step vectors
    %   error_len: length of rms_error
    %
    % This computes movement vectors between consecutive points (current - previous)
    % for true and predicted paths, then returns the euclidean norm of their difference
    % at each time point after applying the specified fixed shift.
    
    % Scale paths to real dimensions
    true_scaled = true_path .* reshape(real_dim, [1 2]);
    pred_scaled = pred_path .* reshape(real_dim, [1 2]);
    
    % Apply fixed shift (same convention as other functions)
    if shift_val > 0
        true_shifted = true_scaled((shift_val+1):end, :);
        pred_shifted = pred_scaled(1:(end-shift_val), :);
    elseif shift_val < 0
        s = abs(shift_val);
        true_shifted = true_scaled(1:(end-s), :);
        pred_shifted = pred_scaled((s+1):end, :);
    else
        true_shifted = true_scaled;
        pred_shifted = pred_scaled;
    end
    
    % Need at least two points to form step vectors
    if size(true_shifted,1) < 2 || size(pred_shifted,1) < 2
        velocity = NaN;
        error_len = 0;
        return;
    end
    
    % Compute step (delta) vectors: current - previous
    true_steps = diff(true_shifted, 1, 1);   % (T_shifted-1) x 2
    pred_steps = diff(pred_shifted, 1, 1);   % (T_shifted-1) x 2

    % If lengths mismatch (shouldn't after shifting) trim to common length
    n = min(size(true_steps,1), size(pred_steps,1));
    true_steps = true_steps(1:n, :);
    pred_steps = pred_steps(1:n, :);
    
    % RMS error per timestep comparing the step vectors
    velocity_true = sqrt(sum((true_steps).^2, 2));
    velocity_pred = sqrt(sum((pred_steps).^2, 2));
    error_len = numel(velocity_true);
end
% ...existing code...
