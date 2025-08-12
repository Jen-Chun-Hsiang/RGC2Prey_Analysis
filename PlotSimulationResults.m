Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
Dates = {'1217202401', '1217202402', '1217202403'};
N_days = length(Dates);
A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A test_losses(:)];
end


[p, h] = ranksum(A(:, 2), A(:, 3))

%% Surround strength
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
Dates = {'1218202401', '1218202402', '1218202403'};
Surround_strength = [0, 0.4, 0.8];
N_days = length(Dates);

A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A training_losses(:)];
    legs{i} = sprintf('s:%0.2G (%s)', Surround_strength(i), Dates{i});
end

colors = parula(N_days);
figure; hold on
for i = 1:N_days
    plot(A(:, i), 'Color', colors(i, :));
end
legend(legs);
title('Impact of surround strength on trajectory prediction')
xlabel('Epoch (#)')
ylabel('Training losses')
%% Coordinate correction
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
Dates = {'1220202401', '1220202402', '1220202403', '1220202404', '1220202405', '1220202406'};
Coord_adj_dir = [-1 1 -1 1 0 -1];
Coord_adj_type = {'body', 'body', 'head', 'head', 'body', 'body'};
N_days = length(Dates);

A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A training_losses(:)];
    legs{i} = sprintf('%s (d=%1G)', Coord_adj_type{i}, Coord_adj_dir(i));
end

colors = parula(N_days);
figure; hold on
for i = 1:N_days
    plot(A(:, i), 'Color', colors(i, :));
end
legend(legs);
title('Coordinate correction on trajectory prediction')
xlabel('Epoch (#)')
ylabel('Training losses')

%% Quantization
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
Dates = {'1220202401', '1223202401', '1223202402', '1223202403'};
Quantization_fac = [100, 64, 32, 16];
N_days = length(Dates);

A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A training_losses(:)];
    legs{i} = sprintf('%0.3G (%s)', Quantization_fac(i), Dates{i});
end

colors = parula(N_days);
figure; hold on
for i = 1:N_days
    plot(A(:, i), 'Color', colors(i, :));
end
legend(legs);
title('Quantization on trajectory prediction')
xlabel('Epoch (#)')
ylabel('Training losses')


%% Denisty (close to real)
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
Dates = {'1223202403', '1224202401', '1224202402', '1224202403', '1225202401', '1225202402', '1225202403',...
    '1226202401', '1226202402', '1226202403'};
% Density = [176, 50, 113, 176, 50, 113, 176, 50, 113];
Num_target = [500 500, 142, 320, 500, 142, 320, 500, 142, 320];
Quantization_fac = [16 16 16 16 128 128 128 32 32 32];
Mask_radius = [15.4 8.33 8.33 8.33 8.33 8.33 8.33 17.7 17.7 17.7];
SF_scale = [0.2 0.54 0.54 0.54 0.54 0.54 0.54 0.54 0.54 0.54]; % 0.2 0.2 0.2

N_days = length(Dates);

A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A training_losses(:)];
    legs{i} = sprintf('N(%d) Q(%d) M(%0.3G) S(%0.3G)', Num_target(i), Quantization_fac(i),...
        Mask_radius(i), SF_scale(i));
end

colors = parula(N_days);
figure; hold on
for i = 1:N_days
    plot(A(:, i), 'Color', colors(i, :));
end
legend(legs);
title('Density on trajectory prediction')
xlabel('Epoch (#)')
ylabel('Training losses')

%% Surround strength with background grass "blend" no quantization but rectification
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_200_prediction_error';
is_different_seed = 0;
switch is_different_seed
    case 1
        Dates = {'1229202401', '1229202407', '2024010201', '1229202404', '1229202402', '1229202408', '2024010202',...
            '1229202405', '1229202403', '1229202409', '2024010203', '1229202406'};
        surround_fac = [0 0 0 0.2 0.4 0.4 0.4,...
            0.6 0.8 0.8 0.8 1.0];
    case 0
        Dates = {'1229202401', '1229202407', '1229202404', '1229202402', '1229202408',...
            '1229202405', '1229202403', '1229202409', '1229202406'};
        surround_fac = [0 0  0.2 0.4 0.4 ,...
            0.6 0.8 0.8 1.0];
end
N_days = length(Dates);

A = [];
for i = 1:N_days
    load(fullfile(Folder_Name, [Dates{i} Tag]))
    A = [A training_losses(:)];
    legs{i} = sprintf('%0.3G', surround_fac(i));
end

colors = parula(N_days);
figure; hold on
for i = 1:N_days
    plot(A(:, i), 'Color', colors(i, :));
end
legend(legs);
title('Quantization on trajectory prediction')
xlabel('Epoch (#)')
xticks(0:50:200);
xticklabels({'0', '50', '100', '150', '200'})
ylim([0.05 0.3])
ylabel('Training losses')
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})
%% Surround strength with noise
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2024010201', '1229202404', '1229202402', '2024010203'};
Noise_level = {'0', '0.005', '0.01', '0.02', '0.04'};
surround_fac = [0 0.2 0.4 0.8];
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
end
x = 1:N_levels;
colors = lines(N_days);
figure; hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(i, :);
    legs{i} = sprintf('%0.3G', surround_fac(i));
end
legend(legs);
title('Surround inhibition in cricket prediction')
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.05 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

%% Density with noise
clear
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025010301', '2025010302', '2025010303'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
Density = [500 320 142];
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
        n_sample = length(test_losses);
        Data_m(i, j) = mean(test_losses);
        Data_s(i, j) = std(test_losses)/sqrt(n_sample);
    end
end
x = 1:N_levels;
colors = lines(N_days);
figure; hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(i, :);
    legs{i} = sprintf('%0.3G', Density(i));
end
legend(legs);
title('Surround inhibition in cricket prediction')
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.05 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

%% Complex model
clear
title_name = '(Model) Different architectures of prediction models';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025010502', '2025010501', '2025010503', '2025010601', '2025010602', '2025010603'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
Exps = {'1 layer', '2 layers', '2 layers +noise', '2 layers +2 flats +noise', '3 layers +noise', '3 layers +2 flats +noise'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(Exps);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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

%% Complex model
close all; clear
title_name = '(Training) Training noise levels';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025010701', '2025010702', '2025010703','2025010801', '2025010802', '2025010803'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
Noise_trained = {'0.0', '0.005', '0.01', '0.0', '0.005', '0.01'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(Noise_trained);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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
%% Complex model
close all; clear
title_name = '(Model) Number of LSTM hidden layer';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025011001', '2025011002', '2025011003', '2025011101', '2025011102', '2025011103'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'3', '4', '5', '3', '4', '5'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
title('(Model) Number of LSTM hidden layers in cricket prediction')
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.05 0.3])
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


%% Complex model
close all; clear
title_name = '(Model) LSTM hidden layer size';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025011701', '2025011702', '2025011703', '2025011801', '2025011802', '2025011803'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'64', '128', '256', '64', '128', '256'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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

%% Complex model
close all; clear
title_name = '(Model) CNN feature size';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025011901', '2025011902', '2025011903', '2025012001', '2025012002', '2025012003'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'128', '256', '512', '128', '256', '512'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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

%% Complex model
close all; clear
title_name = '(Model) CNN kernel size';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012101', '2025012102', '2025012103'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'4-2 4-1', '3-1 4-1', '3-1 3-1'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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

%% Complex model
close all; clear
title_name = '(Model) Depth of CNN layer size';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012201', '2025012202', '2025012203'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'2 (1-1)', '3 (0.5-1-1)', '5 (0.1-0.25-0.5-1-2)'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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


%% Complex model
close all; clear
title_name = '(Model) CNN depth & layer channel size';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012201', '2025012202', '2025012301', '2025012203','2025012302', '2025012303'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'2 (1-1)', '3 (0.5-1-1)', '3 (0.5-1-1)', '5 (0.1-0.25-0.5-1-2)', '5 (0.25-0.25-0.5-1-2)', '5 (0.5-0.5-0.5-1-2)'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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


%% Complex model
close all; clear
title_name = '(Model) Minimum learning rate (CAWR)';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012401', '2025012402', '2025012403'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'1e-6', '5e-7', '1e-7'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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
ylim([0.05 0.3])
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


%% Complex model
close all; clear
title_name = '(Model) Bypassing RGCs';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012701', '2025012702', '2025012703'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'Direct image', 'RGC', 'half sf_mask'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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


%% Complex model
close all; clear
title_name = '(Model) Density & Spatial Receptive Field';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025012801', '2025012802', '2025012803', '2025012901', '2025012902', '2025012903'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'Temporal', 'Middle', 'Nasal', 'Temporal', 'Middle', 'Nasal'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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

%% Complex model
close all; clear
title_name = '(Model) Spatial Receptive Size (ratio) fixed cell number';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025013001', '2025012902', '2025013002'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'0.378', '0.459', '0.54'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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


%% Complex model
close all; clear
title_name = '(Model) Density (fixed RF size)';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025013003', '2025012902', '2025013004'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'1020', '442', '142'};
n_epoch = 200;
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days = length(Dates);
N_levels = length(Noise_level);
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
for i = 1:N_days
    for j = 1:N_levels
        load(fullfile(Folder_Name, [Dates{i} '_noise' Noise_level{j} Tag]))
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


%% Complex model
close all; clear
title_name = '(Model) Temporal filter';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025020301', '2025020302', '2025020303'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'median', 'more biphasic', 'less biphasic'};
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

%% ON OFF contribution
close all; clear
title_name = '(Model) ON and OFF alpha RGCs';
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag = '_cricket_location_prediction_200_prediction_error';
Dates = {'2025020501', '2025020502', '2025020503'};
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};
LSTM_layer_n = {'ON-OFF (442)', 'ON-OFF (221)', 'ON (442)'};
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