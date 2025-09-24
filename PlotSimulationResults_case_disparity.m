clear; clc;
%% ——— User picks which experiment to plot ———
exp_id = 1;  
% 1: Noise contribution to ON/OFF grid
% 2: Temporal filter biphasic
% 3: Surround inhibition

%% ——— Common settings ———
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag         = '_cricket_location_prediction_200_prediction_error';
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
    otherwise
        error('exp_id must be 1, 2 or 3');
end

%% ——— Shared data‐loading and plotting code ———

N_days   = numel(Second_level);
N_levels = numel(Noise_level);

Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf(fname_pattern, Dates{1}, BG_folder, Second_level{i}, Noise_level{j}, Tag);
        clear test_losses training_losses
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses');
        n_sample      = numel(test_losses);
        Data_m(i, j)  = mean(test_losses);
        Data_s(i, j)  = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(N_days);
%%
hFig = figure;

% ——— Right: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
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
subplot(1, 2, 1); hold on
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

