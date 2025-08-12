%% Boundary size & drift grating (seed: 2513)
close all; clear
title_name = '(Model)';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021001', '2025021002', '2025021003'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'boundary size', 'same boundary size', 'drift-grating'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

subplot(1, 2, 1); hold on
for i = 1:N_days
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :));
    legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
title(sprintf('%s in cricket prediction', title_name))
xlabel('Training epochs')
xticks(0:50:200);
xticklabels({'0', '50', '100', '150', '200'})
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);
yticklabels({'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'})


%% Last layer projection to encode relative drift grating (seed: 2514)
close all; clear
title_name = '(Model) Auxiliary cost ratio';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021101', '2025021102', '2025021103'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'0', '0.25', '0.5'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

subplot(1, 2, 1); hold on
for i = 1:N_days
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :));
    legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
title(sprintf('%s in cricket prediction', title_name))
xlabel('Training epochs')
xticks(0:50:200);
xticklabels({'0', '50', '100', '150', '200'})
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);
yticklabels({'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'})

%% Auxiliary lost with combined LSTM (feedback) to encode relative drift grating (seed: 2514)

title_name = '(Model) Location, Model Feedback, Background image';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021201', '2025021202', '2025021203', '2025021204'}; %, 
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'rloc | lstm-proj | drift', 'loc | lstm-proj | drift',...
    'rloc | one-proj | drift', 'rloc | lstm-proj | texture'};% , 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Temporal filter (seed: 251X)

title_name = '(Model) Temporal filter biphasic';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021301', '2025021302', '2025021303', '2025021204'}; %, 
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'Biphasic | drift', 'Monphasic | drift','Semibiphasic| drift','Biphasic | blend ' };% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Density filter (seed: 251X)
title_name = '(Model) Retina regional difference';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021402', '2025021401', '2025021403'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'Temporal', 'Middle','Nasal' };% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Coverage (seed: 251X)
title_name = '(Model) Coverage';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021405', '2025021401', '2025021404'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'0.378','0.459','0.54' };% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Coverage (seed: 251X)
title_name = '(Model) Density';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021407', '2025021401', '2025021406'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'142','442','1020' };% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Mask SF radius (seed: 251X)
title_name = '(Model) Masked SF radius';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021404', '2025021408'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'37','50'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses;
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Surround strength (seed: 251X)
title_name = '(Model) Surround strength';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021601', '2025021602', '2025021603'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'0.0','-0.4', '-0.8'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Alignment (seed: 251X)
title_name = '(Model) ON OFF grid alignment';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021801', '2025021802', '2025021803', '2025021804', '2025021805'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'0.0','1.0', '0.5', '1.0 norm', '0.0 norm'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Alignment (seed: 251X)
title_name = '(Model) Density and regional difference';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025022301', '2025022302', '2025022303', '2025022304', '2025022305'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'500 0.378','500 0.288', '142 0.54', '142 0.709', '142 0.54 2'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Binocular (seed: 251X)
title_name = '(Model) Binocular';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025022402', '2025022404','2025022401', '2025022403', '2025022405'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'No 500','No 1000', 'Bi 500', 'Bi 500 norm', 'Bi 500 0'};% 
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]))
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Center-suround (seed: 251X)
title_name = '(Model) Binocular';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
% Dates = {'2025022701', '2025022702','2025022703', '2025022704'};
Dates = {'2025030101', '2025030102','2025030103'}; 
% BG_folder = {'', '', '', ''};
% BG_folder = {'', '', ''};
BG_folder = {'blend_', 'blend_', 'blend_'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
% LSTM_layer_n = {'-0.0 75','-0.4 75', '-0.8 75', '-0.0 50'};
% LSTM_layer_n = {'-0.0 75 (t)','-0.4 75 (t)', '-0.8 75 (t)'};
LSTM_layer_n = {'-0.0 75 (2)','-0.4 75 (2)', '-0.8 75 (2)'};% 
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
colors = jet(N_days);
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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

%% Center-suround (seed: 251X)
title_name = '(Model) Noise contribution to different ON OFF grid';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025030201', '2025030202','2025030203', '2025030204', '2025030401', '2025030402', '2025030403'}; %, 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'', '', '', '','blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'n0.1','n0.3', 'n0.5', '0.1 s-0.8', 'n0.1 norm', 'n0.3 norm','n0.5 norm'};% 
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Center-suround (seed: 251X)
title_name = '(Model) Surround inhibition';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025030501', '2025030502','2025030503', '2025030504'}; 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
% BG_folder = {'blend_', 'blend_', 'blend_', 'blend_'};
BG_folder = {'texture-blend_', 'texture-blend_', 'texture-blend_', 'texture-blend_'};
% LSTM_layer_n = {'-1.0','-1.4', '-1.8', '0.0'};
LSTM_layer_n = {'-1.0 (t)','-1.4 (t)', '-1.8 (t)', '0.0 (t)'};
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Center-suround (seed: 251X)
title_name = '(Model) Surround inhibition';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025030701', '2025030702','2025030703', '2025030704',...
        '2025030801', '2025030802','2025030803', '2025030804'}; 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'blend_', 'blend_', 'blend_', 'blend_',...
              'blend_', 'blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'0.0','-0.4', '-0.8', '-1.2',...
                 '0.0','-0.4', '-0.8', '-1.2'};
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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


%% Center-suround (seed: 251X)
title_name = '(Model) Binocular grid configuration';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025031001', '2025031002','2025031003', '2025031004','2025031005', '2025031006',...
    '2025042201', '2025042202', '2025042203', '2025042205'}; 
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'blend_', 'blend_', 'blend_', 'blend_', 'blend_', 'blend_',...
    'blend_', 'blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'Bi (500)','Two grid Bi (500)', 'Bi (1000)', 'Bi Dist 0.0 (500)', 'Mo (1000)', 'Mo ON-OFF (500)',...
    'Mo ON-OFF 2 (500)', 'Mo ON-OFF 3 (500)', 'Mo ON-ON (500)', 'Mo ON-ON (500) aligned'};
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
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

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
%% Center-suround (seed: 251X)
title_name = '(Model) Binocular grid configuration';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025031001', '2025031002','2025031003', '2025031004','2025031005', '2025031006', ...
         '2025042201', '2025042202', '2025042203', '2025042205'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'blend_', 'blend_', 'blend_', 'blend_', 'blend_', 'blend_', ...
             'blend_', 'blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'Bi (500)','Two grid Bi (500)', 'Bi (1000)', 'Bi Dist 0.0 (500)', ...
                'Mo (1000)', 'Mo ON-OFF (500)', 'Mo ON-OFF 2 (500)', ...
                'Mo ON-OFF 3 (500)', 'Mo ON-ON (500)', 'Mo ON-ON (500) aligned'};
n_epoch = 200;

% ——— New: only plot these day‐indices ———
plot_line_ids = [1 3];   % example: will plot Dates 1, 3, 7, and 10

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = length(Dates);
N_levels = length(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load all data
for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf('%s_cricket_%snoise%s%s', Dates{i}, BG_folder{i}, Noise_level{j}, Tag);
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses')
        n_sample       = numel(test_losses);
        Data_m(i, j)   = mean(test_losses);
        Data_s(i, j)   = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = unique_string{string_ids(i)};
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
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))


%% Alignment (seed: 251X)
title_name = '(Model) ON OFF grid alignment';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021801', '2025021802', '2025021803', '2025021804', '2025021805'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
LSTM_layer_n = {'0.0','1.0', '0.5', '1.0 norm', '0.0 norm'};
n_epoch = 200;

% ——— New: only plot these day‐indices ———
plot_line_ids = [4, 5];   % e.g. only plot Dates{1}, Dates{2}, and Dates{5}

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = numel(Dates);
N_levels = numel(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load data
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_cricket_noise' Noise_level{j} Tag]), 'test_losses', 'training_losses')
        n_sample     = numel(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right panel: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x); xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left panel: Training‐loss vs Epoch ———
subplot(1, 2, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))

%% Center-suround (seed: 251X)
title_name = '(Model) Noise contribution to different ON OFF grid';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025030201', '2025030202','2025030203', '2025030204', '2025030401', '2025030402', '2025030403'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'', '', '', '','blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'n0.1','n0.3', 'n0.5', '0.1 s-0.8', 'n0.1 norm', 'n0.3 norm','n0.5 norm'};
n_epoch = 200;

% ——— New: only plot these day‐indices ———
plot_line_ids = [5, 6, 7];   % e.g. will plot Dates{1}, Dates{3}, and Dates{6}

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = numel(Dates);
N_levels = numel(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load all data
for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf('%s_cricket_%snoise%s%s', Dates{i}, BG_folder{i}, Noise_level{j}, Tag);
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses');
        n_sample     = numel(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right panel: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x); xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left panel: Training‐loss vs Epoch ———
subplot(1, 2, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))

%%
title_name = '(Model) Temporal filter biphasic';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025021301', '2025021302', '2025021303', '2025021204'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'Biphasic | drift', 'Monphasic | drift', 'Semibiphasic| drift', 'Biphasic | blend '};
n_epoch = 200;

% ——— New: only plot these day‐indices ———
plot_line_ids = [1, 2, 3];   % e.g. will only plot the 1st and 3rd entries from Dates/LSTM_layer_n

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = numel(Dates);
N_levels = numel(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load all data
for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf('%s_cricket_noise%s%s', Dates{i}, Noise_level{j}, Tag);
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses');
        n_sample     = numel(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right panel: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x); xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left panel: Training‐loss vs Epoch ———
subplot(1, 2, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))

%% Center-suround (seed: 251X)
title_name = '(Model) Surround inhibition';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025030701', '2025030702','2025030703', '2025030704', ...
         '2025030801', '2025030802','2025030803', '2025030804'};
Noise_level = {'0.0','0.002', '0.004', '0.008', '0.016', '0.032'};
BG_folder = {'blend_', 'blend_', 'blend_', 'blend_', ...
             'blend_', 'blend_', 'blend_', 'blend_'};
LSTM_layer_n = {'0.0','-0.4', '-0.8', '-1.2', ...
                '0.0','-0.4', '-0.8', '-1.2'};
n_epoch = 200;

% ——— New: only plot these day‐indices ———
plot_line_ids = [1, 2, 3, 4];   % e.g. will plot entries 1, 3, 5 & 7

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = numel(Dates);
N_levels = numel(Noise_level);

% Preallocate
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);

% Load all data
for i = 1:N_days
    for j = 1:N_levels
        fname = sprintf('%s_cricket_%snoise%s%s', Dates{i}, BG_folder{i}, Noise_level{j}, Tag);
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses');
        n_sample     = numel(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

figure;

% ——— Right panel: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(string_ids(i), :);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Noise levels')
xticks(x); xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);

% ——— Left panel: Training‐loss vs Epoch ———
subplot(1, 2, 1); hold on
clear legs
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    plot(Data_t(i, :), 'Color', colors(string_ids(i), :), 'LineWidth', 1.2);
    legs{k} = unique_string{string_ids(i)};
end
legend(legs, 'Location', 'best');
xlabel('Training epochs')
xticks(0:50:200); xticklabels({'0','50','100','150','200'});
ylabel('Training losses');
ylim([0 0.5])
yticks(0.0:0.1:0.5);

sgtitle(sprintf('%s in cricket prediction', title_name))



