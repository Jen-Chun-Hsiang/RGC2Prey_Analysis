clear; clc;
%% ——— User picks which experiment to plot ———
exp_id = 3;  
% 1: Noise contribution to ON/OFF grid
% 2: Temporal filter biphasic
% 3: Surround inhibition

%% ——— Common settings ———
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
fig_save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Summary\Illustrator\'; 
Tag         = '_cricket_location_prediction_200_prediction_error_with_path';
n_epoch     = 200;

switch exp_id
   
    case 1
        title_name    = '(Model) Fixed disparity';
        Dates         = {'2025071501'};
        Second_name   = 'disp';
        Second_level  = {'0.0', '1.0', '2.0', '4.0', '8.0'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = 'blend';
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '4.0', '16.0', '8.0',...
                         '1.0', '0.0', '1.5'};
        plot_line_ids = [1:5];
        fname_pattern = '%s_cricket_%s_disp%s_noise%s%s';
    case 2
        title_name    = '(Model) Fixed disparity';
        Dates         = {'2025092902'};
        Second_name   = 'disp';
        Second_level  = {'0.5', '1.0', '2.0', '4.0', '8.0'};
        Noise_level   = {'0.016','0.032','0.064'};
        BG_folder     = 'blend';
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '4.0', '16.0', '8.0',...
                         '1.0', '0.0', '1.5'};
        plot_line_ids = [1:5];
        fname_pattern = '%s_cricket_%s_disp%s_noise%s%s';
    case 3
        title_name    = '(Model) Fixed disparity test (trained dynamic)';
        Dates         = {'2025101403'};  % 2025100602 2025101403 2025101402
        Second_name   = 'disp';
        Second_level  = {'0.0', '3.0', '6.0', '12.0'};
        Noise_level   = {'0.016','0.032','0.064', '0.128','0.256'};
        BG_folder     = 'blend';
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '4.0', '16.0', '8.0',...
                         '1.0', '0.0', '1.5'};
        plot_line_ids = [1:4];
        fname_pattern = '%s_cricket_%s_disp%s_noise%s%s';
        exp_name_tag = 'fixed-disparity-test-OFF';
    
    case 4
        title_name    = '(Model) Fixed disparity test (trained fixed)';
        Dates         = {'2025100507'};
        Second_name   = 'disp';
        Second_level  = {'0.0', '3.0', '6.0', '12.0'};
        Noise_level   = {'0.016','0.032','0.064', '0.128','0.256'};
        BG_folder     = 'blend';
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '4.0', '16.0', '8.0',...
                         '1.0', '0.0', '1.5'};
        plot_line_ids = [1:4];
        fname_pattern = '%s_cricket_%s_disp%s_noise%s%s';
        exp_name_tag = sprintf('fixed-disparity-test_fix_training_%s_ON', Dates{1});
    otherwise
        error('exp_id must be 1, 2 or 3');
end

if ~exist('exp_name_tag','var') || isempty(exp_name_tag)
    exp_name_tag = sprintf('exp%d_%s', exp_id, Dates{1});
end

%% ——— Shared data‐loading and plotting code ———

N_days   = numel(Second_level);
N_levels = numel(Noise_level);

load_mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\';
coverage_mat_file = fullfile(load_mat_folder, 'processed_cover_radius.mat');
cover_radius_struct = load(coverage_mat_file, 'file_index_list', 'processed_cover_radius');
cover_radius = [cover_radius_struct.file_index_list(:) cover_radius_struct.processed_cover_radius(:)];

fixed_shift = -9;
sf_scale = 0.54;
pix_to_um = 4.375; % each pixel is 4.375 um
real_dim = [120 90] * pix_to_um / sf_scale;
is_correct_object_zone = true;
bg_type = 'grass';  % 'blend', 'grass', 'simple'
if ~exist('is_degree','var')
    is_degree = true;
end
num_sample = [];

if is_degree
    unit_factor = 1/32.5; % convert microns to visual degrees
else %#ok<UNRCH>
    unit_factor = 1;
end

Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
Data_n = nan(N_days, N_levels);
DataP_m = nan(N_days, N_levels);
DataP_s = nan(N_days, N_levels);
DataCM_m = nan(N_days, N_levels);
DataCM_s = nan(N_days, N_levels);
Data_v = [];
DataP_v = [];
DataCM_v = [];
D_validation = [];
KeyP = cell(N_days, N_levels);

for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf(fname_pattern, Dates{1}, BG_folder, Second_level{i}, Noise_level{j}, Tag);
        clear test_losses training_losses all_paths all_paths_pred all_id_numbers all_scaling_factors all_bg_file all_path_cm
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses',...
            'all_paths', 'all_paths_pred', 'all_id_numbers', 'all_scaling_factors', 'all_bg_file', 'all_path_cm');
        n_sample      = numel(test_losses);
        if isempty(D_validation)
            num_sample = n_sample;
            Data_v = nan(N_days, N_levels, num_sample);
            DataP_v = nan(N_days, N_levels, num_sample);
            DataCM_v = nan(N_days, N_levels, num_sample);
            D_validation = nan(N_days, N_levels, num_sample);
        else
            assert(n_sample == num_sample, 'Mismatch in number of test samples');
        end
        Data_m(i, j)  = mean(test_losses);
        Data_s(i, j)  = std(test_losses)/sqrt(n_sample);
        Data_v(i, j, 1:n_sample) = test_losses;
        Data_n(i, j)  = n_sample;

        all_paths_r = reshapeAllPaths(all_paths);
        D_validation(i, j, 1:n_sample) = all_paths_r(:, end-1, 1) * unit_factor;
        try
            all_path_cm = reshapeAllPaths(all_path_cm);
        catch
            disp('Center of mass path not calculated, using regular path instead');
            all_path_cm = all_paths_r;
        end
        all_path_cm = reshapeAllPaths(all_path_cm);
        all_paths_pred_r = squeeze(all_paths_pred);
        is_simple_contrast = cellfun(@(x) contains(x, 'gray_image'), all_bg_file);
        assert(n_sample == size(all_paths_r, 1), 'Mismatch in number of test samples and loaded paths');

        % Store per-sample subject keys for paired stats across disparity levels
        % (robust even if ordering differs across files).
        KeyP{i, j} = makePredationSubjectKeys(all_id_numbers, all_bg_file);

        for ii = 1:n_sample
            true_path_trial = squeeze(all_paths_r(ii, :, :));
            pred_path_trial = squeeze(all_paths_pred_r(ii, :, :));
            pred_cm_path_trial = squeeze(all_path_cm(ii, :, :));
            true_path_scaled = true_path_trial .* reshape(real_dim, [1 2]);

            cut_off = acceptance_zone_radius(double(all_id_numbers(ii)), all_scaling_factors(ii, 50:end), cover_radius, fixed_shift);
            cut_off_cm = acceptance_zone_radius(double(all_id_numbers(ii)), all_scaling_factors(ii, 50:end), cover_radius, fixed_shift);
            [fixed_rms, rms_len] = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim);
            [fixed_cm_rms, rms_len_cm] = calculateFixedShiftRMSError(true_path_scaled, pred_cm_path_trial * pix_to_um / sf_scale, fixed_shift, ones(1, 2));

            if ii == 1
                all_fixed_rms = zeros(n_sample, rms_len);
                all_fixed_cm_rms = zeros(n_sample, rms_len_cm);
            end

            all_fixed_rms(ii, :) = double(fixed_rms(:))';
            all_fixed_cm_rms(ii, :) = double(fixed_cm_rms(:))';

            if is_correct_object_zone
                all_fixed_rms(ii, :) = max(0, all_fixed_rms(ii, :) - cut_off(:)');
                all_fixed_cm_rms(ii, :) = max(0, all_fixed_cm_rms(ii, :) - cut_off_cm(:)');
            end
        end

        all_fixed_rms = all_fixed_rms * unit_factor;
        all_fixed_cm_rms = all_fixed_cm_rms * unit_factor;

        switch bg_type
            case 'blend'
                DataP_m(i, j)  = mean(all_fixed_rms, 'all');
                DataP_s(i, j)  = std(mean(all_fixed_rms, 1))/sqrt(n_sample);
                DataP_v(i, j, 1:size(all_fixed_rms,1)) = mean(all_fixed_rms, 2);
                DataCM_m(i, j) = mean(all_fixed_cm_rms, 'all');
                DataCM_s(i, j) = std(mean(all_fixed_cm_rms, 1))/sqrt(n_sample);
                DataCM_v(i, j, 1:size(all_fixed_cm_rms,1)) = mean(all_fixed_cm_rms, 2);
            case 'grass'
                n_grass = sum(~is_simple_contrast);
                if n_grass == 0
                    DataP_m(i, j) = NaN; DataP_s(i, j) = NaN;
                    DataCM_m(i, j) = NaN; DataCM_s(i, j) = NaN;
                else
                    DataP_m(i, j)  = mean(all_fixed_rms(~is_simple_contrast, :), 'all');
                    DataP_s(i, j)  = std(mean(all_fixed_rms(~is_simple_contrast, :), 1))/sqrt(n_grass);
                    DataCM_m(i, j) = mean(all_fixed_cm_rms(~is_simple_contrast, :), 'all');
                    DataCM_s(i, j) = std(mean(all_fixed_cm_rms(~is_simple_contrast, :), 1))/sqrt(n_grass);
                end

                % Preserve original sample indices (important for paired stats)
                DataP_v(i, j, ~is_simple_contrast) = mean(all_fixed_rms(~is_simple_contrast, :), 2);
                DataCM_v(i, j, ~is_simple_contrast) = mean(all_fixed_cm_rms(~is_simple_contrast, :), 2);

                % Excluded samples stay as NaN
                DataP_v(i, j, is_simple_contrast) = NaN;
                DataCM_v(i, j, is_simple_contrast) = NaN;
            case 'simple'
                n_simple = sum(is_simple_contrast);
                if n_simple == 0
                    DataP_m(i, j) = NaN; DataP_s(i, j) = NaN;
                    DataCM_m(i, j) = NaN; DataCM_s(i, j) = NaN;
                else
                    DataP_m(i, j)  = mean(all_fixed_rms(is_simple_contrast, :), 'all');
                    DataP_s(i, j)  = std(mean(all_fixed_rms(is_simple_contrast, :), 1))/sqrt(n_simple);
                    DataCM_m(i, j) = mean(all_fixed_cm_rms(is_simple_contrast, :), 'all');
                    DataCM_s(i, j) = std(mean(all_fixed_cm_rms(is_simple_contrast, :), 1))/sqrt(n_simple);
                end

                % Preserve original sample indices (important for paired stats)
                DataP_v(i, j, is_simple_contrast) = mean(all_fixed_rms(is_simple_contrast, :), 2);
                DataCM_v(i, j, is_simple_contrast) = mean(all_fixed_cm_rms(is_simple_contrast, :), 2);

                % Excluded samples stay as NaN
                DataP_v(i, j, ~is_simple_contrast) = NaN;
                DataCM_v(i, j, ~is_simple_contrast) = NaN;
            otherwise
                error('bg_type must be blend or grass or simple');
        end
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

% Note: do not trim the 3rd dimension based on a single slice.
% For bg_type='grass' or 'simple', excluded samples are represented as NaN.
% Keeping the original indexing avoids accidental misalignment across
% noise levels or disparity conditions.

x      = 1:N_levels;
colors = lines(N_days);
%%
validation_level_id = min(2, N_levels);
validation_ids = plot_line_ids;
figure; hold on
for idx = 1:numel(validation_ids)
    id = validation_ids(idx);
    plot(squeeze(D_validation(id, validation_level_id, :)), 'Color', colors(id, :), 'LineWidth', 1.0);
end
xlabel('Sample index');
if is_degree
    ylabel('Validation distance (degrees)');
else %#ok<UNRCH>
    ylabel('Validation distance (cm)');
end
legend(Second_level(validation_ids), 'Location', 'best');
title(sprintf('Final validation distances | noise %s', Noise_level{validation_level_id}));

hFig = figure;

% ——— Right: RMS Error vs Noise ———
clear legs
subplot(1, 3, 3); hold on
e = errorbar(x, mean(DataCM_m, 1), mean(DataCM_s, 1), 'CapSize', 0, 'LineWidth',1.5);
e.Color   = 0.4*ones(1,3);
legs{1}   = 'Center of Mass';
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, DataP_m(i, :), DataP_s(i, :), 'CapSize', 0, 'LineWidth',1.5);
    e.Color   = colors(i, :);
    legs{end+1}   = sprintf('%s %s', Second_name, Second_level{i});
end

legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(1:numel(Noise_level)); xticklabels(Noise_level);
xlim([0.5 numel(Noise_level)+0.5]);
ylim([0 12]);
yticks(0:6:12);
if is_degree
    ylabel('Error dist. (degrees)');
else %#ok<UNRCH>
    ylabel('Error dist. (cm)');
end

% ——— Middle: Test‐loss vs Noise ———
clear legs
subplot(1, 3, 2); hold on
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0, 'LineWidth',1.5);
    e.Color   = colors(i, :);
    legs{k}   = Second_level{i};
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x); xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left: Training‐loss vs Epoch ———
subplot(1, 3, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(i, :), 'LineWidth', 1.2);
    legs{k} = Second_level{i};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))

%%
if is_correct_object_zone
    object_tag = 'correctZone';
else
    object_tag = 'noCorrection';
end

save_file_name = fullfile(fig_save_folder, sprintf('SummaryPredResultsDisparity_%s_%s_%s', exp_name_tag, bg_type, object_tag));
print(gcf, [save_file_name '.eps'], '-depsc', '-vector'); % EPS format
print(gcf, [save_file_name '.png'], '-dpng', '-r300'); % PNG, 600 dpi

%% --- Stats: per-noise paired nonparam across disparity levels (key-matched) ---
% Mirrors PlotSimulationResults_case.m usage, but matches samples via keys to
% ensure the paired assumption holds even if file ordering differs.
try
    res = runFDR_ANOVA_nonparam_matched(DataP_v(plot_line_ids, :, :), KeyP(plot_line_ids, :));
    sig = find(res.signif_levels);
    if ~isempty(sig)
        fprintf('\nFDR-significant noise levels (alpha=%.3f):\n', res.alpha);
        for kk = 1:numel(sig)
            j = sig(kk);
            fprintf('  noise=%s | p=%.3g | q=%.3g | nMatched=%d\n', Noise_level{j}, res.pvals(j), res.qvals(j), res.n_matched_per_level(j));
        end
    else
        fprintf('\nNo FDR-significant noise levels (alpha=%.3f).\n', res.alpha);
    end
catch ME
    warning('Stats analysis failed: %s', ME.message);
    res = [];
end

