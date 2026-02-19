%% MultiExperimentComparison.m - Multi-Experiment Cricket Path Prediction Analysis
%
% This script compares cricket path prediction data across multiple experiments
% providing visual comparisons and verification of consistency.
%
% FEATURES:
% 1. Load and compare multiple experiments
% 2. Single trial visualization across experiments with different colors
% 3. Fixed RMS error comparison plots across time
% 4. Velocity comparison and path consistency verification
%
% USAGE:
% Set the experiment names in exp_names cell array and run the script
% Plot are learned from plot_three_traces_gradient.m
%% Configuration - Shared Variables
clear; clc; 
mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats\';
fig_save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Summary\Illustrator\'; 
bg_type = 'blend'; % or 'grass'
noise_level = '0.016'; % Fixed noise level for comparison
disparity_sets = {}; % e.g., {} to compare different experiments
%disparity_sets = {'0.0','3.0','6.0','12.0'}; % e.g., {'0.0','3.0','6.0','12.0'} (id 20) to compare fixed disparity runs

% Configure multiple experiments to compare
exp_name_tag = 'Demo_only_center';
exp_names = {'2026021102'}; %94 Add/modify experiment names here
% gain control (ON) ''2026021102', '2026021104', '2026021105' (OFF) '2026021106',  '2026021107',   '2026021108' (26)
% temporal shift (ON) '2026021401', '2026021404', '2026021405' (OFF) , '2026021406', '2026021408', '2026021410' (94)
% Fixed disparity training '2025100507', '2025100501',  '2025100502',  '2025100503', '2025100508','2025100504','2025100505','2025100506'
% Interocular distance '2025101402',  '2025101408', '2025101404', '2025101403',  '2025101409', '2025101406' (29)
% Interoccular distance '2025100604', '2025100605', '2025100602', '2025100606', '2025100607', '2025100603'
% Interoccular distance '2025092904', '2025092905', '2025092902', '2025092906', '2025092907', '2025092903'
% surround inhibition   '2025091805', '2025091802', '2025091807', '2025091808', '2025091803', '2025091810'
% varied coverage       '2025092109', '2025091802', '2025092107', '2025092110', '2025091803', '2025092108'
% varied density        '2025092105', '2025091802', '2025092102', '2025092106', '2025091803', '2025092104'
color_ids = [3]; %3, 2, 1, 8, 7, 6
trial_id = 26; % Trial ID to visualize across experiments 16 34 5
disp_trajectory_id = [1, 2];
is_y_axis_flip = true; % set true only for movie drawing
is_plot_ground_truth = true;
is_plot_pred_trace = true;

% Disparity shift trace control: 'none', 'left', or 'right'
% 'none': use all_paths_pred_r as usual
% 'left': replace predicted trace with top_img_positions_shifted (left eye)
% 'right': replace predicted trace with top_img_disparity_positions_shifted (right eye)
disparity_shift_trace = 'none'; % Options: 'none', 'left', 'right'
interocular_dist = 1.0; % Interocular distance in cm for reconstruction

% Visual settings
is_visual_degree = 1;
vis_deg_to_um = 32.5;
vis_scale = is_visual_degree / vis_deg_to_um + (1 - is_visual_degree); % Convert cm to visual degrees if enabled

% Physical dimensions and scaling
real_dim = [120 90]*4.375/0.54;
cm_dim_scale = 4.375/0.54; % Convert cm to um
fixed_shift = -9; % Fixed shift for correlation analysis


if isempty(disparity_sets)
    use_disparity_mode = false;
    comparison_labels = exp_names;
    base_exp_name = '';
else
    use_disparity_mode = true;
    comparison_labels = disparity_sets;
    base_exp_name = exp_names{1};
end

n_experiments = numel(comparison_labels);
if n_experiments == 0
    error('Either exp_names or disparity_sets must contain at least one entry.');
end

if isempty(color_ids) || numel(color_ids) < n_experiments
    color_ids = 1:n_experiments;
else
    color_ids = color_ids(1:n_experiments);
end
disp_trajectory_id = min(max(disp_trajectory_id, 1), n_experiments);

if use_disparity_mode
    comparison_mode_label = 'Disparity';
    comparison_title_suffix = 'Disparities';
    output_tag = sprintf('%s_%s_disp', exp_name_tag, base_exp_name);
else
    comparison_mode_label = 'Experiment';
    comparison_title_suffix = 'Experiments';
    output_tag = exp_name_tag;
end


% Load coverage data
load_mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\';
coverage_mat_file = fullfile(load_mat_folder, 'processed_cover_radius.mat');
cover_radius = load(coverage_mat_file, 'file_index_list', 'processed_cover_radius');
cover_radius = [cover_radius.file_index_list(:) cover_radius.processed_cover_radius(:)];

%% Load data for all comparisons
experiments = cell(n_experiments, 1);
comparison_names = cell(n_experiments, 1);
comparison_disp_vals = cell(n_experiments, 1);
bounds = [-100, 100, -60, 60];
interocular_dist = 1.0;

fprintf('Loading data for %d %s comparisons...\n', n_experiments, lower(comparison_mode_label));
for i = 1:n_experiments
    if use_disparity_mode
        current_exp_name = base_exp_name;
        current_disp_val = comparison_labels{i};
        current_label = sprintf('disp %s', current_disp_val);
    else
        current_exp_name = comparison_labels{i};
        current_disp_val = '';
        current_label = current_exp_name;
    end

    fprintf('Loading comparison %d/%d: %s\n', i, n_experiments, current_label);
    [all_paths_r, all_paths_pred_r, seqLen, is_simple_contrast, all_id_numbers, ...
     all_scaling_factors, all_path_cm, all_paths_bg] = loadDataset(mat_folder, current_exp_name, bg_type, noise_level, current_disp_val);
    
    experiments{i}.name = current_label;
    experiments{i}.source_exp = current_exp_name;
    experiments{i}.disp_value = current_disp_val;
    experiments{i}.all_paths_r = all_paths_r;
    experiments{i}.all_paths_pred_r = all_paths_pred_r;
    experiments{i}.all_path_cm = all_path_cm;
    experiments{i}.all_paths_bg = all_paths_bg;
    experiments{i}.seqLen = seqLen;
    experiments{i}.is_simple_contrast = is_simple_contrast;
    experiments{i}.all_id_numbers = all_id_numbers;
    experiments{i}.all_scaling_factors = all_scaling_factors;
    experiments{i}.n_trials = size(all_paths_r, 1);

    comparison_names{i} = current_label;
    comparison_disp_vals{i} = current_disp_val;
    fprintf('  - Loaded %d trials with sequence length %d\n', experiments{i}.n_trials, seqLen);
end

%% 1. Single trial visualization across experiments with different colors
baseColors = [
    60, 0, 60;     % 1 ON-Temporal
    120, 0, 120;   % 2 OFF-Temporal
    180, 0, 180;   % 3 OFF-Nasal
    240, 0, 240;   % 4
    0, 60, 0;      % 5
    0, 120, 0;     % 6
    0, 180, 0;     % 7
    0, 240, 0;     % 8
    120, 120, 120; % 9
    180, 120, 0;   % 10
    0, 120, 180;   %11
    120, 120, 180; %12
]/255;
%true_color = baseColors(10, :); % Brown for ground truth (for movie frame)
true_color = [0.4, 0.4, 0.4]; 

% Extend colors if needed
if n_experiments > size(baseColors, 1)
    additional_colors = lines(n_experiments - size(baseColors, 1));
    baseColors = [baseColors; additional_colors];
end

fprintf('\nCreating single trial visualization for trial %d...\n', trial_id);

% Figure 1: Single trial paths comparison
figure('Name', sprintf('Single Trial Comparison (%s) - Trial %d', comparison_mode_label, trial_id), 'Position', [100, 100, 1200, 800]);

% Subplot 1: Path trajectories with gradient colors
subplot(2, 2, 1); hold on;

% Parameters for gradient plotting
alpha_start = 0.25;             % how "light" the start is (0..1)
alpha_end   = 1.0;              % how "dark" the end is (0..1)



legend_handles = [];
legend_entries = {};

% Get ground truth data from first experiment for reference
if length(disp_trajectory_id) == 1
    ref_exp = experiments{disp_trajectory_id};
else
    ref_exp = experiments{1};
end
if is_plot_ground_truth && trial_id <= ref_exp.n_trials
    true_path_x = ref_exp.all_paths_r(trial_id, :, 1) * real_dim(1) * vis_scale;
    true_path_y = ref_exp.all_paths_r(trial_id, :, 2) * real_dim(2) * vis_scale;
    
    N = length(true_path_x);
    segAlpha = linspace(alpha_start, alpha_end, N-1);
    
    % Plot ground truth with gray gradient
    for j = 1:(N-1)
        alpha = segAlpha(j);
        color = (1 - alpha) * [1 1 1] + alpha * true_color;
        plot(true_path_x(j:j+1), true_path_y(j:j+1), '-', 'Color', color, 'LineWidth', 2.5);
    end
    % Add legend handle for ground truth
    legend_color = (1 - alpha_end) * [1 1 1] + alpha_end * true_color;
    legend_handles(end+1) = plot(nan, nan, '-', 'Color', legend_color, 'LineWidth', 2.5);
    legend_entries{end+1} = 'Ground Truth';
    
    % Mark start and end for ground truth
    markerSize = 60;
    startColor = (1 - segAlpha(1)) * [1 1 1] + segAlpha(1) * true_color;
    endColor = (1 - segAlpha(end)) * [1 1 1] + segAlpha(end) * true_color;
    scatter(true_path_x(1), true_path_y(1), markerSize/2, startColor, 'filled', 'MarkerEdgeColor','k');
    scatter(true_path_x(end), true_path_y(end), markerSize/2, endColor, 'filled', 'MarkerEdgeColor','k');
end

% Plot prediction traces only for experiments specified in disp_trajectory_id
if is_plot_pred_trace
    for idx = 1:length(disp_trajectory_id)
        i = disp_trajectory_id(idx);
        if i <= n_experiments
            exp = experiments{i};
            if trial_id <= exp.n_trials
                % Get trajectory data - apply disparity shift if requested
                if strcmp(disparity_shift_trace, 'none')
                    % Use standard predicted path
                    pred_path_normalized = squeeze(exp.all_paths_pred_r(trial_id, :, :));
                else
                    % Reconstruct shifted positions using disparity
                    path = squeeze(exp.all_paths_r(trial_id, :, :));
                    scaling_factors = exp.all_scaling_factors(trial_id, :);
                    bounds = []; % No bounds adjustment by default
                    
                    [top_img_positions_shifted, top_img_disparity_positions_shifted, disparity_deg, distances] = ...
                        reconstruct_top_img_positions_shifted(path, scaling_factors, interocular_dist, bounds);
                    
                    if strcmp(disparity_shift_trace, 'left')
                        pred_path_normalized = top_img_positions_shifted;
                    elseif strcmp(disparity_shift_trace, 'right')
                        pred_path_normalized = top_img_disparity_positions_shifted;
                    else
                        error('disparity_shift_trace must be ''none'', ''left'', or ''right''');
                    end
                end
                
                % Convert to visual degrees
                pred_path_x = pred_path_normalized(:, 1) * real_dim(1) * vis_scale;
                pred_path_y = pred_path_normalized(:, 2) * real_dim(2) * vis_scale;
                
                N = length(pred_path_x);
                segAlpha = linspace(alpha_start, alpha_end, N-1);
                
                % Plot predicted path with experiment-specific color gradient
                pred_color = baseColors(color_ids(i),:);
                for j = 1:(N-1)
                    alpha = segAlpha(j);
                    color = (1 - alpha) * [1 1 1] + alpha * pred_color;
                    plot(pred_path_x(j:j+1), pred_path_y(j:j+1), '-', 'Color', color, 'LineWidth', 2.5);
                end
                % Add legend handle for this experiment's prediction
                legend_color = (1 - alpha_end) * [1 1 1] + alpha_end * pred_color;
                legend_handles(end+1) = plot(nan, nan, '-', 'Color', legend_color, 'LineWidth', 2.5);
                legend_entries{end+1} = sprintf('%s Prediction', exp.name);
                
                % Mark start and end for predictions
                markerSize = 60;
                startColor = (1 - segAlpha(1)) * [1 1 1] + segAlpha(1) * pred_color;
                endColor = (1 - segAlpha(end)) * [1 1 1] + segAlpha(end) * pred_color;
                scatter(pred_path_x(1), pred_path_y(1), markerSize/2, startColor, 'filled', 'MarkerEdgeColor','k');
                scatter(pred_path_x(end), pred_path_y(end), markerSize/2, endColor, 'filled', 'MarkerEdgeColor','k');
            end
        end
    end
end

% Set axis limits based on twice the real_dim converted to visual degrees
x_lim_max = real_dim(1) * vis_scale;  % Half-width in visual degrees
y_lim_max = real_dim(2) * vis_scale;  % Half-height in visual degrees
% xlim([-x_lim_max, x_lim_max]);
% ylim([-y_lim_max, y_lim_max]);

% Add boundary rectangle
x_min = -real_dim(1) * vis_scale;
x_max =  real_dim(1) * vis_scale;
y_min = -real_dim(2) * vis_scale;
y_max =  real_dim(2) * vis_scale;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1);

% Add scale bars for 10 degrees
if is_visual_degree
    % Vertical scale bar (10 deg)
    plot(-20*ones(1, 2), [-20 -10], '-k', 'LineWidth', 2);
    text(-22, -15, '10 deg', 'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    % Horizontal scale bar (10 deg)
    plot([-20 -10], -20*ones(1, 2), '-k', 'LineWidth', 2);
    text(-15, -22, '10 deg', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

xlabel('X Position (deg)'); ylabel('Y Position (deg)'); 
title(sprintf('Cricket Location Prediction - Path Trajectories (%s)', comparison_mode_label));
if ~isempty(legend_handles)
    legend(legend_handles, legend_entries, 'Location', 'best');
else
    legend off;
end
if is_y_axis_flip
    set(gca, 'YDir', 'reverse');
end
axis off; box off;

% Subplot 2: Fixed RMS Error comparison
subplot(2, 2, 2); hold on;
for i = 1:n_experiments
    exp = experiments{i};
    if trial_id <= exp.n_trials
        true_path_trial = squeeze(exp.all_paths_r(trial_id, :, :));
        pred_path_trial = squeeze(exp.all_paths_pred_r(trial_id, :, :));
        [fixed_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
        
        t = (0:rms_len-1)/100; % Time vector in seconds
        plot(t, fixed_rms * vis_scale, '-', 'Color', baseColors(color_ids(i),:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('%s Fixed RMS', exp.name));
    end
end
xlabel('Time (s)'); ylabel('Distance to cricket (deg)'); 
ylim([0 35]); box off;
yticks(0:15:30);
yticklabels(arrayfun(@(y) sprintf('%d', y), 0:15:30, 'UniformOutput', false));
xlim([0 t(end)]);
xticks(0:0.9:1.8);
xticklabels(arrayfun(@(x) sprintf('%.1f', x), 0:0.9:1.8, 'UniformOutput', false));
title('Fixed RMS Error Comparison');
legend('Location', 'best');
% ylim([0 15]); box off;

% Subplot 3: Fixed CM RMS Error comparison  
% subplot(2, 2, 3); hold on;
% for i = 1:n_experiments
%     exp = experiments{i};
%     if trial_id <= exp.n_trials
%         true_path_trial = squeeze(exp.all_paths_r(trial_id, :, :));
%         pred_cm_path_trial = squeeze(exp.all_path_cm(trial_id, :, :));
%         [fixed_cm_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial .* reshape(real_dim, [1 2]), ...
%                                                               pred_cm_path_trial * cm_dim_scale, fixed_shift, ones(1, 2));
        
%         t = (0:rms_len-1)/100; % Time vector in seconds
%         plot(t, fixed_cm_rms * vis_scale, '-', 'Color', baseColors(i,:), 'LineWidth', 1.5, ...
%              'DisplayName', sprintf('%s CM RMS', exp.name));
%     end
% end
% xlabel('Time (s)'); ylabel('Distance to cricket (deg)'); 
% title('Fixed CM RMS Error Comparison');
% legend('Location', 'best');
% ylim([0 15]); box off;

% Subplot 4: Velocity comparison with validation
subplot(2, 2, 4); hold on;

% Validate velocity consistency across experiments and calculate single trial velocities
cricket_velocities = cell(n_experiments, 1);
bg_velocities = cell(n_experiments, 1);

for i = 1:n_experiments
    exp = experiments{i};
    if trial_id <= exp.n_trials
        true_path_trial = squeeze(exp.all_paths_r(trial_id, :, :));
        bg_path_trial = squeeze(exp.all_paths_bg(trial_id, :, :));
        [velocity_cricket, velocity_bg] = calculateVelocity(true_path_trial, bg_path_trial, fixed_shift, real_dim);
        
        cricket_velocities{i} = velocity_cricket;
        bg_velocities{i} = velocity_bg;
    end
end

% Validation: Check if velocities are identical across experiments
if n_experiments > 1
    fprintf('\n=== VELOCITY CONSISTENCY VALIDATION ===\n');
    for i = 2:n_experiments
        % Check cricket velocity consistency
        if ~isempty(cricket_velocities{1}) && ~isempty(cricket_velocities{i})
            cricket_diff = mean(abs(cricket_velocities{1} - cricket_velocities{i}));
            try
                if cricket_diff* 100 * vis_scale  < 1e-10
                    assert(cricket_diff* 100 * vis_scale  < 1e-10, sprintf('Cricket velocities differ between %s and %s (max diff: %.2e)', ...
                        comparison_names{1}, comparison_names{i}, cricket_diff));
                    fprintf('✓ Cricket velocities are IDENTICAL between %s and %s\n', comparison_names{1}, comparison_names{i});
                elseif cricket_diff* 100 * vis_scale  < 5
                    fprintf('⚠ WARNING: Cricket velocities are SIMILAR between %s and %s (max diff: %.2e)\n', ...
                        comparison_names{1}, comparison_names{i}, cricket_diff);
                else
                    error('Cricket velocities differ significantly');
                end
            catch ME
                fprintf('⚠ ERROR: %s\n', ME.message);
                error('Velocity validation failed for cricket velocities');
            end
        end
        
        % Check background velocity consistency
        if ~isempty(bg_velocities{1}) && ~isempty(bg_velocities{i})
            bg_diff = max(abs(bg_velocities{1} - bg_velocities{i}));
            try
          assert(bg_diff* 100 * vis_scale < 1e-10, sprintf('Background velocities differ between %s and %s (max diff: %.2e)', ...
              comparison_names{1}, comparison_names{i}, bg_diff));
          fprintf('✓ Background velocities are IDENTICAL between %s and %s\n', comparison_names{1}, comparison_names{i});
            catch ME
                fprintf('⚠ ERROR: %s\n', ME.message);
                error('Velocity validation failed for background velocities');
            end
        end
    end
    fprintf('✓ All velocity validations passed successfully\n\n');
end

% Plot single cricket and background velocity traces (since they should be identical across experiments)
if ~isempty(cricket_velocities{1})
    t = (0:length(cricket_velocities{1})-1)/100; % Time vector in seconds
    
    % Plot cricket velocity
    plot(t, cricket_velocities{1} * 100 * vis_scale, '-', 'Color', [0.2, 0.6, 0.2], 'LineWidth', 2, ...
         'DisplayName', 'Cricket Velocity');
    
    % Plot background velocity  
    plot(t, bg_velocities{1} * 100 * vis_scale, '-', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 2, ...
         'DisplayName', 'Background Velocity');
end

xlabel('Time (s)'); ylabel('Velocity (deg/s)'); 
title(sprintf('Velocity Traces (Trial %d)', trial_id));
xlim([0 t(end)]);
xticks(0:0.9:1.8);
xticklabels(arrayfun(@(x) sprintf('%.1f', x), 0:0.9:1.8, 'UniformOutput', false));
legend('Location', 'best');
ylim([0 150]); box off;
yticks(0:50:150);
yticklabels(arrayfun(@(y) sprintf('%d', y), 0:50:150, 'UniformOutput', false));

%%
save_file_name = fullfile(fig_save_folder, sprintf('PredictionTrace_%s_%s_%d', output_tag, noise_level, trial_id));
print(gcf, [save_file_name '.eps'], '-depsc', '-vector'); % EPS format
print(gcf, [save_file_name '.png'], '-dpng', '-r300'); % PNG, 600 dpi

keyboard
%%
%% 2. Fixed RMS and CM RMS as function of time frames across comparisons
fprintf('Creating fixed RMS error comparison plots across %s...\n', lower(comparison_title_suffix));

figure('Name', sprintf('Fixed RMS Error Comparison Across %s', comparison_title_suffix), 'Position', [200, 150, 1400, 600]);

% Calculate statistics for each experiment
exp_stats = cell(n_experiments, 1);

% Initialize matrix to store average fixed RMS across trials and experiments
% This will be size (max_trials, n_experiments) where each column is one experiment
max_trials = max(arrayfun(@(i) experiments{i}.n_trials, 1:n_experiments));
DataTab = NaN(max_trials, n_experiments);

for i = 1:n_experiments
    exp = experiments{i};
    fprintf('Processing experiment %s (%d trials)...\n', exp.name, exp.n_trials);
    
    % Initialize arrays for all trials
    all_fixed_rms = [];
    all_fixed_cm_rms = [];
    
    for j = 1:exp.n_trials
        true_path_trial = squeeze(exp.all_paths_r(j, :, :));
        pred_path_trial = squeeze(exp.all_paths_pred_r(j, :, :));
        pred_cm_path_trial = squeeze(exp.all_path_cm(j, :, :));
        
        % Calculate fixed RMS errors
        [fixed_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
        try 
            [fixed_cm_rms, ~] = calculateFixedShiftRMSError(true_path_trial .* reshape(real_dim, [1 2]), ...
                                                            pred_cm_path_trial * cm_dim_scale, fixed_shift, ones(1, 2));
        catch 
            fixed_cm_rms = NaN(size(fixed_rms));
        end

        if j == 1
            all_fixed_rms = zeros(exp.n_trials, rms_len);
            all_fixed_cm_rms = zeros(exp.n_trials, rms_len);
        end
        
        all_fixed_rms(j, :) = fixed_rms;
        all_fixed_cm_rms(j, :) = fixed_cm_rms;
        
        % Store average fixed RMS across time for this trial in the cross-experiment matrix
        DataTab(j, i) = mean(fixed_rms);
    end
    
    % Calculate mean and SEM
    exp_stats{i}.name = exp.name;
    exp_stats{i}.mean_fixed_rms = mean(all_fixed_rms, 1);
    exp_stats{i}.sem_fixed_rms = std(all_fixed_rms, [], 1) / sqrt(exp.n_trials);
    exp_stats{i}.mean_fixed_cm_rms = mean(all_fixed_cm_rms, 1);
    exp_stats{i}.sem_fixed_cm_rms = std(all_fixed_cm_rms, [], 1) / sqrt(exp.n_trials);
    exp_stats{i}.rms_len = rms_len;
    
    % Also store the full trial-by-time matrix for this experiment
    exp_stats{i}.all_fixed_rms = all_fixed_rms;
    exp_stats{i}.all_fixed_cm_rms = all_fixed_cm_rms;
end

DataTab = [(1:max_trials)' is_simple_contrast(:) DataTab*vis_scale];
DataTab = sortrows(DataTab, [2 3]) % Sort by simple contrast column
% Display information about the cross-experiment matrix
fprintf('\nCreated avg_fixed_rms_across_trials_and_exp matrix:\n');
fprintf('  Size: %d trials x %d experiments\n', size(DataTab, 1), size(DataTab, 2));
fprintf('  Each value is the time-averaged fixed RMS for that trial and experiment\n');
fprintf('  Experiment order: %s\n', strjoin(exp_names, ', '));
fprintf('  Use this matrix to compare trial-by-trial performance across experiments\n\n');

% Plot Fixed RMS Error
subplot(1, 2, 1); hold on;
for i = 1:n_experiments
    stats = exp_stats{i};
    t = (0:stats.rms_len-1)/100;
    
    % Plot with SEM shading
    plotWithSEM(t, stats.mean_fixed_rms * vis_scale, stats.sem_fixed_rms * vis_scale, baseColors(i,:));
end
xlabel('Time (s)'); ylabel('Distance to cricket (deg)');
title('Fixed RMS Error Across Experiments');
ylim([0 15]); box off;

% Create legend
legend_entries = comparison_names;
legend(legend_entries, 'Location', 'best');

% Plot Fixed CM RMS Error  
subplot(1, 2, 2); hold on;
for i = 1:n_experiments
    stats = exp_stats{i};
    t = (0:stats.rms_len-1)/100;
    
    % Plot with SEM shading
    plotWithSEM(t, stats.mean_fixed_cm_rms * vis_scale, stats.sem_fixed_cm_rms * vis_scale, baseColors(i,:));
end
xlabel('Time (s)'); ylabel('Distance to cricket (deg)');
title('Fixed CM RMS Error Across Experiments');
ylim([0 15]); box off;
legend(legend_entries, 'Location', 'best');

%%  Compare the results of both experiments
% need to save the first DataTab as DataTab2, and run the script with another experiment
if exist('DataTab2', 'var')
    DataTabAll = [sortrows(DataTab2, 1), sortrows(DataTab, 1)];
    DataTabAll= sortrows(DataTabAll, [2 6])
end


%% 3. Velocity comparison and path consistency verification
fprintf('Verifying path consistency across experiments...\n');

figure('Name', 'Velocity and Path Consistency Verification', 'Position', [300, 200, 1400, 800]);

% Velocity comparison
subplot(2, 2, 1); hold on;
velocity_stats = cell(n_experiments, 1);

for i = 1:n_experiments
    exp = experiments{i};
    all_velocity_true = [];
    all_velocity_bg = [];
    
    for j = 1:min(exp.n_trials, 50) % Limit to first 50 trials for efficiency
        true_path_trial = squeeze(exp.all_paths_r(j, :, :));
        bg_path_trial = squeeze(exp.all_paths_bg(j, :, :));
        
        [velocity_true, velocity_bg] = calculateVelocity(true_path_trial, bg_path_trial, fixed_shift, real_dim);
        
        if j == 1
            all_velocity_true = zeros(min(exp.n_trials, 50), length(velocity_true));
            all_velocity_bg = zeros(min(exp.n_trials, 50), length(velocity_bg));
        end
        
        all_velocity_true(j, :) = velocity_true;
        all_velocity_bg(j, :) = velocity_bg;
    end
    
    velocity_stats{i}.name = exp.name;
    velocity_stats{i}.mean_velocity_true = mean(all_velocity_true, 1);
    velocity_stats{i}.sem_velocity_true = std(all_velocity_true, [], 1) / sqrt(size(all_velocity_true, 1));
    velocity_stats{i}.mean_velocity_bg = mean(all_velocity_bg, 1);
    velocity_stats{i}.sem_velocity_bg = std(all_velocity_bg, [], 1) / sqrt(size(all_velocity_bg, 1));
end

% Plot True object velocity
for i = 1:n_experiments
    stats = velocity_stats{i};
    t = (0:length(stats.mean_velocity_true)-1)/100;
    plotWithSEM(t, stats.mean_velocity_true * 100 * vis_scale, stats.sem_velocity_true * 100 * vis_scale, baseColors(i,:));
end
xlabel('Time (s)'); ylabel('Velocity (deg/s)');
title('True Cricket Velocity Across Experiments');
ylim([0 150]); box off;
legend(legend_entries, 'Location', 'best');

% Plot Background velocity
subplot(2, 2, 2); hold on;
for i = 1:n_experiments
    stats = velocity_stats{i};
    t = (0:length(stats.mean_velocity_bg)-1)/100;
    plotWithSEM(t, stats.mean_velocity_bg * 100 * vis_scale, stats.sem_velocity_bg * 100 * vis_scale, baseColors(i,:));
end
xlabel('Time (s)'); ylabel('Velocity (deg/s)');
title('Background Velocity Across Experiments');
ylim([0 150]); box off;
legend(legend_entries, 'Location', 'best');

% Path consistency verification - compare first few trials
subplot(2, 2, 3); hold on;
n_verify_trials = min(5, min(arrayfun(@(i) experiments{i}.n_trials, 1:n_experiments)));

for trial = 1:n_verify_trials
    for i = 1:n_experiments
        exp = experiments{i};
        true_path = squeeze(exp.all_paths_r(trial, :, :));
        plot(true_path(:, 1), true_path(:, 2), '-', 'Color', baseColors(i,:), 'LineWidth', 1);
    end
end
xlabel('X Position (normalized)'); ylabel('Y Position (normalized)');
title(sprintf('True Paths Consistency (Trials 1-%d)', n_verify_trials));
legend(legend_entries, 'Location', 'best');
box off;

% Background paths consistency verification
subplot(2, 2, 4); hold on;
for trial = 1:n_verify_trials
    for i = 1:n_experiments
        exp = experiments{i};
        bg_path = squeeze(exp.all_paths_bg(trial, :, :));
        plot(bg_path(:, 1), bg_path(:, 2), '-', 'Color', baseColors(i,:), 'LineWidth', 1);
    end
end
xlabel('X Position (normalized)'); ylabel('Y Position (normalized)');
title(sprintf('Background Paths Consistency (Trials 1-%d)', n_verify_trials));
legend(legend_entries, 'Location', 'best');
box off;

%% Quantitative consistency verification
fprintf('\n=== PATH CONSISTENCY VERIFICATION ===\n');
reference_exp = experiments{1}; % Use first experiment as reference

for i = 2:n_experiments
    comp_exp = experiments{i};
    
    % Compare dimensions
    fprintf('Comparing %s vs %s:\n', reference_exp.name, comp_exp.name);
    fprintf('  Sequence lengths: %d vs %d\n', reference_exp.seqLen, comp_exp.seqLen);
    
    % Compare a subset of trials for true paths
    n_compare = min([reference_exp.n_trials, comp_exp.n_trials, 10]);
    path_diff_sum = 0;
    bg_diff_sum = 0;
    
    for j = 1:n_compare
        ref_path = squeeze(reference_exp.all_paths_r(j, :, :));
        comp_path = squeeze(comp_exp.all_paths_r(j, :, :));
        path_diff = mean(sqrt(sum((ref_path - comp_path).^2, 2)));
        path_diff_sum = path_diff_sum + path_diff;
        
        ref_bg = squeeze(reference_exp.all_paths_bg(j, :, :));
        comp_bg = squeeze(comp_exp.all_paths_bg(j, :, :));
        bg_diff = mean(sqrt(sum((ref_bg - comp_bg).^2, 2)));
        bg_diff_sum = bg_diff_sum + bg_diff;
    end
    
    avg_path_diff = path_diff_sum / n_compare;
    avg_bg_diff = bg_diff_sum / n_compare;
    
    fprintf('  Avg true path difference (first %d trials): %.6f\n', n_compare, avg_path_diff);
    fprintf('  Avg background path difference (first %d trials): %.6f\n', n_compare, avg_bg_diff);
    
    if avg_path_diff < 1e-10
        fprintf('  ✓ True paths are IDENTICAL\n');
    else
        fprintf('  ⚠ True paths DIFFER\n');
    end
    
    if avg_bg_diff < 1e-10
        fprintf('  ✓ Background paths are IDENTICAL\n');
    else
        fprintf('  ⚠ Background paths DIFFER\n');
    end
    fprintf('\n');
end

fprintf('Multi-experiment comparison complete!\n');

%% ========== FUNCTIONS ==========
% Most functions are copied from PathPredictDistExtraction.m

function [all_paths_r, all_paths_pred_r, seqLen, is_simple_contrast, all_id_numbers, all_scaling_factors, all_path_cm, all_paths_bg] = loadDataset(mat_folder, exp_name, bg_type, noise_level, disp_val)
    % Load and process a single dataset
    % Always expect disp_val as 5th argument; if empty or '', use standard filename
    if ~isempty(disp_val) && ~(ischar(disp_val) && isempty(strtrim(disp_val)))
        filename = sprintf('%s_cricket_%s_disp%s_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                          exp_name, bg_type, disp_val, noise_level);
    else
        filename = sprintf('%s_cricket_%s_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                          exp_name, bg_type, noise_level);
    end
    filepath = fullfile(mat_folder, filename);

    if ~exist(filepath, 'file')
        error('File not found: %s', filepath);
    end

    fprintf('  Loading: %s\n', filename);
    data = load(filepath);
    [all_paths_r, seqLen] = reshapeAllPaths(data.all_paths);
    if size(data.all_path_cm, 2) ==1
        all_path_cm = data.all_path_cm; % Already in correct shape
    else
        all_path_cm = reshapeAllPaths(double(data.all_path_cm));
    end
    all_paths_bg = reshapeAllPaths(double(data.all_paths_bg));
    all_paths_pred_r = squeeze(data.all_paths_pred);
    is_simple_contrast = cellfun(@(x) contains(x, 'gray_image'), data.all_bg_file);
    all_id_numbers = data.all_id_numbers;
    all_scaling_factors = data.all_scaling_factors;

    assert(isequal(size(all_paths_r), size(all_paths_pred_r)), ...
           'Reshaped paths and predicted paths must have the same dimensions');
end

function plotWithSEM(x, mean_err, sem_err, color)
    % Plot mean with SEM shading
    ci_upper = mean_err + sem_err;
    ci_lower = mean_err - sem_err;
    
    % Create shaded region for SEM with transparency
    if ischar(color)
        rgb_color = getColorRGB(color);
    else
        rgb_color = color;
    end
    
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

function [velocity_true, velocity_pred, error_len] = calculateVelocity(true_path, pred_path, shift_val, real_dim)
    % calculateVelocity - Compute velocity (step) vectors after optional fixed shift
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
        velocity_true = NaN;
        velocity_pred = NaN;
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
    
    % Velocity magnitude per timestep
    velocity_true = sqrt(sum((true_steps).^2, 2));
    velocity_pred = sqrt(sum((pred_steps).^2, 2));
    error_len = numel(velocity_true);
end


%%
N = size(path,1);
test_interocular_dist_set = [0.5, 1.0];
figure; hold on
for i = 1:length(test_interocular_dist_set)
    test_interocular_dist = test_interocular_dist_set(i);
    distances = linspace(21, 4, N);                       % linearly from 21 cm to 4 cm
    if exist('binocular_disparity','file')==2 || exist('binocular_disparity','builtin')==5
        disparity_deg = arrayfun(@(d) binocular_disparity(test_interocular_dist, d), distances);
    else
        disparity_deg = 2 * atand((test_interocular_dist/2) ./ distances);   % fallback formula
    end

    disparity_deg = disparity_deg(:);  
    plot(distances, disparity_deg');
end
xlim([4 21]);
xticks(4:8:20);
xticklabels(arrayfun(@(x) sprintf('%d', x), 4:8:20, 'UniformOutput', false));
yticks(0:7:14);
yticklabels(arrayfun(@(y) sprintf('%d', y), 0:7:14, 'UniformOutput', false));
ylim([0 15]);
xlabel('Distance to Object (cm)'); ylabel('Disparity (deg)');
title('Disparity vs Distance for Given Interocular Distance');
box off

save_file_name = fullfile(fig_save_folder, 'disparity_vs_distance');
print(gcf, [save_file_name '.eps'], '-depsc', '-vector'); % EPS format
print(gcf, [save_file_name '.png'], '-dpng', '-r300'); % PNG, 600 dpi
