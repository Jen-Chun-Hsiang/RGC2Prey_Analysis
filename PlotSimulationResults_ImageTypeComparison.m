clear; clc;
% Define fixed folder path
folder_path = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats'; % Update with your actual directory

% Define timestamp and noise level variations
% timestamps = {'2025021401', '2025021402', '2025021403'};  % Modify as needed
timestamps = {'2025021401', '2025021402', '2025021403', '2025021404', '2025021405', '2025021406',...
    '2025021407', '2025021408'};  
Noise_level = {'0.0', '0.002', '0.004', '0.008', '0.016', '0.032'};  % Modify as needed
LSTM_layer_n = {'Temporal', 'Middle','Nasal' };% 
x_scale = 120;
y_scale = 90;

[unique_string, ~, string_ids] = unique(LSTM_layer_n);
Data = [];
N_days = length(timestamps);
N_levels = length(Noise_level);
for i = 1:N_days
    
    for j = 1:N_levels
        % Construct file name dynamically
        filename = sprintf('%s_cricket_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                            timestamps{i}, Noise_level{j});

        % Full file path
        full_path = fullfile(folder_path, filename);

        
        % Check if file exists before loading
        if isfile(full_path)
            fprintf('Loading file: %s\n', full_path);
            load(full_path);  % Load .mat file
            all_paths_pred = squeeze(all_paths_pred);
            n_tp = size(all_paths_pred, 2);
            n_sp = size(all_paths_pred, 1);

            pred_error = mean((reshape(permute(double(all_paths_pred), [1, 3, 2]), n_sp, [])-all_paths).^2, 2);
           
            pred_error = pred_error(:);
            bg_type = cellfun(@(x) contains(x, 'gray_image'), all_bg_file);
            bg_ids = cellfun(@(x) extractX(x) , all_bg_file);
            ob_ids = all_id_numbers;
        else
            fprintf('File not found: %s\n', full_path);
        end
        Data = [Data; [repmat([i, j], numel(pred_error), 1) double(bg_type) double(ob_ids) pred_error bg_ids]];
    end
end

% Helper function to extract x
function x = extractX(s)
    tokens = regexp(s, '(\d+)_0', 'tokens');
    if ~isempty(tokens)
        x = str2double(tokens{1}{1});
    else
        x = -1; % Handle invalid format
    end
end
%%
Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
for i = 1:N_days
    for j = 1:N_levels
        cids = Data(:, 1) == i & Data(:, 2) == j & Data(:, 3) == 1;
        Data_m(i, j) = mean(Data(cids, 5));
        Data_s(i, j) = std(Data(cids, 5))/sqrt(sum(cids));
    end
end
x = 1:N_levels;
colors = lines(N_days);
figure; 
subplot(1, 2, 1); hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    % e.Color = colors(string_ids(i), :);
    % legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
% legend(legs);
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})
title('Plain gary background')

Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
for i = 1:N_days
    for j = 1:N_levels
        cids = Data(:, 1) == i & Data(:, 2) == j & Data(:, 3) == 0;
        Data_m(i, j) = mean(Data(cids, 5));
        Data_s(i, j) = std(Data(cids, 5))/sqrt(sum(cids));
    end
end

subplot(1, 2, 2); hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    % e.Color = colors(string_ids(i), :);
    % legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
% legend(legs);
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'});
title('Naturalistic grass background')
%% Different cricket
ob_ids = unique(Data(:, 4));
N_ob = length(ob_ids);
Data_m = nan(N_days, N_ob);
Data_s = nan(N_days, N_ob);
for i = 1:N_days
    for j = 1:N_ob
        cids = Data(:, 1) == i & Data(:, 4) == ob_ids(j) & Data(:, 3) == 1;
        if sum(cids)>0
            Data_m(i, j) = mean(Data(cids, 5));
            Data_s(i, j) = std(Data(cids, 5))/sqrt(sum(cids));
        
        end
    end
end
x = 1:N_ob;
colors = lines(N_days);
figure; hold on
clear legs
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    % e.Color = colors(string_ids(i), :);
    % legs{i} = sprintf('%s', unique_string{string_ids(i)});
end
% legend(legs);
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Images')
% xticks(x);
% xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})

%% Different background
bg_ids = unique(Data(:, 6));
bg_ids(bg_ids == -1) = [];
N_bg = length(bg_ids);
Data_m = nan(1, N_bg);
Data_s = nan(1, N_bg);

for j = 1:N_bg
    cids = Data(:, 6) == bg_ids(j);
    Data_m(j) = mean(Data(cids, 5));
    Data_s(j) = std(Data(cids, 5))/sqrt(sum(cids));
end
x = 1:N_bg;
figure; 
subplot(1, 2, 1); 
e = errorbar(x, Data_m, Data_s, 'CapSize', 0);
    
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Images')
% xticks(x);
% xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.4])
yticks(0.1:0.1:0.4);
yticklabels({'0.1', '0.2', '0.3', '0.4'})
box off
xlim([1 N_bg])

Data_m = nan(N_levels, N_bg);
Data_s = nan(N_levels, N_bg);
colors = lines(N_levels);
[n_unique_string, ~, n_string_ids] = unique(Noise_level);
for i = 1:N_levels
    for j = 1:N_bg
        cids = Data(:, 1) == i & Data(:, 6) == bg_ids(j);
        Data_m(i, j) = mean(Data(cids, 5));
        Data_s(i, j) = std(Data(cids, 5))/sqrt(sum(cids));
    end
end
x = 1:N_bg;
subplot(1, 2, 2); hold on
for i = 1:N_levels
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(n_string_ids(i), :);
    legs{i} = sprintf('%s', n_unique_string{n_string_ids(i)});
end
legend(legs);
% title(sprintf('%s in cricket prediction', title_name))
xlabel('Images')
xlim([1 N_bg])
% xticks(x);
% xticklabels(Noise_level);
ylabel('Test losses');
ylim([0.0 0.4])
yticks(0.1:0.1:0.4);
yticklabels({'0.1', '0.2', '0.3', '0.4'})

