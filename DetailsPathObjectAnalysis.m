% Define fixed folder path
folder_path = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats'; % Update with your actual directory

% Define timestamp and noise level variations
timestamps = {'2025013001', '2025013002', '2025013003'};  % Modify as needed
noise_levels = {'0.0', '0.005'};  % Modify as needed

% Iterate over timestamps and noise levels
is_gray = [];
test_loss = [];
scaling_fac = [];
time_loss = [];
for i = 1:length(timestamps)
    for j = 1:length(noise_levels)
        % Construct file name dynamically
        filename = sprintf('%s_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                            timestamps{i}, noise_levels{j});

        % Full file path
        full_path = fullfile(folder_path, filename);
        
        % Check if file exists before loading
        if isfile(full_path)
            fprintf('Loading file: %s\n', full_path);
            load(full_path);  % Load .mat file
            all_paths_pred = squeeze(all_paths_pred);
            n_tp = size(all_paths_pred, 2);
            n_sp = size(all_paths_pred, 1);
            
            c_trial_loss = nan(n_sp, n_tp);
            for k = 1:n_sp
                c_path = all_paths(k, :);
                c_path = reshape(c_path, [2, n_tp]);
                c_path_bg = all_paths_bg(k, :);
                c_path_bg = reshape(c_path_bg, [2, n_tp]);
                c_path_pred = double(squeeze(all_paths_pred(k, :, :, :)))';
                c_trial_loss(k, :) = mean((c_path_pred-c_path).^2, 1);
            end

            gray_mask = contains(all_bg_file, 'gray_image');
            is_gray = [is_gray; gray_mask(:)'];
            all_scaling_factors = all_scaling_factors(:, (end-n_tp+1):end);
            test_loss = [test_loss; test_losses(:)'];
            scaling_fac = [scaling_fac; mean(all_scaling_factors, 1)'];
            time_loss = [time_loss; mean(c_trial_loss(is_gray(1, :)'==1, :), 1)];
        else
            fprintf('File not found: %s\n', full_path);
        end
    end
end
%%
Data_g = [];
Data_i = [];
for i = 1:size(test_loss, 1)
    Data_g = [Data_g; mean(test_loss(i, is_gray(i, :)==1))];
    Data_i = [Data_i; mean(test_loss(i, is_gray(i, :)==0))];
end

figure; hold on
bar(1:2, [mean(Data_g) mean(Data_i)], 'EdgeColor', 'w');
e=errorbar(1:2, [mean(Data_g) mean(Data_i)], ...
    [std(Data_g)/sqrt(numel(Data_g)) std(Data_i)/sqrt(numel(Data_i))], '_', 'CapSize', 0, ...
    'LineWidth',0.5);
e.Color = 'k';
xlabel('Background types')
xticks(1:2);
xlim([0.5 2.5])
xticklabels({'Gray', 'Image'});
ylabel('Test losses');
ylim([0.0 0.17])
yticks(0:0.05:0.15);
yticklabels({'0', '0.05', '0.1', '0.15'});

%% 

%%
test_id = 1;% 1 6
x_scale = 120;
y_scale = 90;

a = all_paths(test_id, :);
a = reshape(a, [2, 202]).*[x_scale y_scale]';
b = all_paths_bg(test_id, :);
b = reshape(b, [2, 202]).*[x_scale y_scale]';
c = double(squeeze(all_paths_pred(test_id, :, :, :)))'.*[x_scale y_scale]';
d = all_path_cm(test_id, :);
d = reshape(d, [2, 202]);
figure; 
subplot(1, 2, 1); hold on
plot(a(1, :), 'k');
plot(b(1, :), 'g');
plot(c(1, :), 'b');
plot(d(1, :), 'r');
legend({'Path','BG', 'Pred', 'CenterM'})
title('x coords');
ylabel('Coords');
xlabel('Time (#frame)');
subplot(1, 2, 2); hold on
plot(a(2, :), 'k');
plot(b(2, :), 'g');
plot(c(2, :), 'b');
plot(d(2, :), 'r');
legend({'Path','BG', 'Pred', 'CenterM'})
title('y coords');
ylabel('Coords');
xlabel('Time (#frame)');


%% Examine difference between model prediction and Center of mass of firing rates
% Define fixed folder path
folder_path = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats'; % Update with your actual directory

title_name = 'Simplest white spot prediction (Model vs center RFs)';
% Define timestamp and noise level variations
timestamps = {'2025020301'};  % Modify as needed
Noise_level = {'0.0', '0.005', '0.01', '0.02', '0.04'};  % Modify as needed
test_object_folder = 'white-spot';
unique_string = {'Model pred', 'Center FRs'};
% Iterate over timestamps and noise levels
x_scale = 120;
y_scale = 90;
pred_error = [];
centerRF_error = [];
for i = 1:length(timestamps)
    for j = 1:length(Noise_level)
        % Construct file name dynamically
        filename = sprintf('%s_%s_noise%s_cricket_location_prediction_200_prediction_error_with_path.mat', ...
                            timestamps{i}, test_object_folder, Noise_level{j});

        % Full file path
        full_path = fullfile(folder_path, filename);
        
        % Check if file exists before loading
        if isfile(full_path)
            fprintf('Loading file: %s\n', full_path);
            load(full_path);  % Load .mat file
            all_paths_pred = squeeze(all_paths_pred);
            n_tp = size(all_paths_pred, 2);
            n_sp = size(all_paths_pred, 1);
            
            
            a = mean((reshape(permute(all_paths_pred, [1, 3, 2]), n_sp, [])-all_paths).^2, 2);
            b = mean((reshape(all_path_cm, n_sp, 2, n_tp)./[x_scale y_scale]-...
                reshape(all_paths, n_sp, 2, n_tp)).^2, [2 3]);
            pred_error = [pred_error; a(:)'];
            centerRF_error = [centerRF_error; b(:)'];
        else
            fprintf('File not found: %s\n', full_path);
        end
    end
end

N_days = 2;
N_levels = length(Noise_level);
x = 1:N_levels;
colors = lines(N_days);
Data_m = [mean(pred_error, 2)'; mean(centerRF_error, 2)'];
Data_s = [std(pred_error, [], 2)'; std(centerRF_error, [], 2)'];
figure; hold on
for i = 1:N_days
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0);
    e.Color = colors(i, :);
    legs{i} = sprintf('%s', unique_string{i});
end
legend(legs);
title(sprintf('%s in cricket prediction', title_name))
xlabel('Noise levels')
xticks(x);
xticklabels(Noise_level);
ylabel('Prediction losses (+/-SD)');
ylim([0.0 0.3])
yticks(0.1:0.1:0.3);
yticklabels({'0.1', '0.2', '0.3'})