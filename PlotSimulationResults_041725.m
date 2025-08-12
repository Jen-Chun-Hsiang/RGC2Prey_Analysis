%% Center-suround (seed: 251X)
title_name = '(Model) End2End training results';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error.mat';
Dates = {'2025041001', '2025041002','2025041003', '2025041004', '2025041005', '2025041006', '2025041007'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'blend_', 'blend_', 'blend_', 'blend_','blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'0.2 (1)','0.4 (2)', '0.8 (3)', '1.6(4)', '2.0 (5)', '1.2 (6)','2.4 (7)'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_' BG_folder{i} 'noise' Noise_level{j} Tag]))        
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end
x = 1:N_levels;
colors = lines(N_days);
figure; 
subplot(1, 2, 2); hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
legend(legs);
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.04])
yticks(0:0.02:0.04);
yticklabels({'0', '0.02', '0.04'})

subplot(1, 2, 1); hold on
for i = 1:N_days
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :));
    legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
sgtitle(sprintf('%s in cricket prediction', title_name))
xlabel('Training epochs')
xticks(0:50:200);
xticklabels({'0', '50', '100', '150', '200'})
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);
yticklabels({'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'})
%%
title_name = '(Model) Density and regional difference';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025022301', '2025022302', '2025022303', '2025022304', '2025022305'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'500 0.378','500 0.288', '142 0.54', '142 0.709', '142 0.54 2'};
n_epoch = 200;

% ——— New: only plot these day-indices ———
plot_line_ids = [3, 4];   % choose which Dates / LSTM_layer_n entries to plot

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = length(Dates);
N_levels = length(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load all data as before
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample       = length(test_losses);
        Data_m(i, j)   = mean(test_losses);
        Data_s(i, j)   = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right panel: test‐loss vs noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = sprintf('%s', unique_string{string_ids(i)});
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left panel: training‐loss vs epoch ———
subplot(1, 2, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = sprintf('%s', unique_string{string_ids(i)});
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200);
xticklabels({'0', '50', '100', '150', '200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))
