clear; clc;
%% ——— User picks which experiment to plot ———
exp_id = 16;  
% 1: Noise contribution to ON/OFF grid
% 2: Temporal filter biphasic
% 3: Surround inhibition

%% ——— Common settings ———
Folder_Name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\Results\Mats';
Tag         = '_cricket_location_prediction_200_prediction_error';
n_epoch     = 200;

switch exp_id
    case 1
        title_name    = '(Model) Noise contribution to different ON OFF grid';
        Dates         = {'2025030201','2025030202','2025030203','2025030204', ...
                         '2025030401','2025030402','2025030403'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = {'','','','', 'blend_','blend_','blend_'};
        LSTM_layer_n  = {'n0.1','n0.3','n0.5','0.1 s-0.8','n0.1 norm','n0.3 norm','n0.5 norm'};
        plot_line_ids = [5,6,7];
        fname_pattern = '%s_cricket_%snoise%s%s';
        
    case 2
        title_name    = '(Model) Temporal filter biphasic';
        Dates         = {'2025021301','2025021302','2025021303','2025021204'};
        Noise_level   = {'0.0','0.005','0.01','0.02','0.04'};
        BG_folder     = {};                          % not used in pattern
        LSTM_layer_n  = {'Biphasic | drift','Monphasic | drift','Semibiphasic| drift','Biphasic | blend '};
        plot_line_ids = [1,2,3];
        fname_pattern = '%s_cricket_noise%s%s';
        
    case 3
        title_name    = '(Model) Surround inhibition';
        Dates         = {'2025030701','2025030702','2025030703','2025030704', ...
                         '2025030801','2025030802','2025030803','2025030804'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,8);
        LSTM_layer_n  = {'0.0','-0.4','-0.8','-1.2','0.0','-0.4','-0.8','-1.2'};
        plot_line_ids = [1,2,3,4];
        fname_pattern = '%s_cricket_%snoise%s%s';

    case 4
        title_name    = '(Model) Repeat binocular test';
        Dates         = {'2025031001','2025042206','2025031004','2025042207','2025042208','2025042209',...
                         '2025042210','2025042211','2025042212','2025042213','2025042214','2025042215',...
                         '2025042216','2025042217'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'Bi (500) 1','Bi (500) 2','Bi dist0 (500) 1','Bi dist0 (500) 2', 'ON (500)', 'OFF (500)',...
                         'ON-OFF', 'ON (-0.1)-OFF', 'ON (-0.2)-OFF', 'ON (-0.3)-OFF', 'Bi ON-OFF', 'Bi OFF',...
                         'Bi ON(-0.2)', 'Bi ON(-0.1)'};
        plot_line_ids = [11:14];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 5
        title_name    = '(Model) ON surround inhibition';
        Dates         = {'2025042208', '2025050201','2025050202','2025050203','2025050204','2025050205',...
                         '2025050701','2025050702','2025050703'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'ON sur0.0', 'ON sur-0.4','ON sur-0.8','ON sur-1.2','ON sur-0.2', 'ON sur-0.6',...
                         'ON sur-0.4 (2)','ON sur-0.8 (2)','ON sur-1.2 (2)'};
        plot_line_ids = [7:9];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 6
        title_name    = '(Model) Contrast Offset & Range';
        Dates         = {'2025052101', '2025052102','2025052103','2025052104','2025052105','2025052106',...
                         '2025052301', '2025052302','2025052303','2025052304','2025052305','2025052306'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'Mean-Nan', 'Mean-1.0','Mean-2.0','Mean-0.0','Mean+1.0', 'Mean+2.0',...
                         'bContr1.0','bContr1.5','bContr2.0','bContr0.75','bContr0.5','bContr0.25'};
        plot_line_ids = [7:12];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 7
        title_name    = '(Model) Density';
        Dates         = {'2025021407', '2025021401', '2025021406'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({''},1,length(Dates));
        LSTM_layer_n  = {'142','442','1020'};
        plot_line_ids = [1:3];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 8
        title_name    = '(Model) Noise and training';
        Dates         = {'2025060401', '2025060402','2025060403','2025060404','2025060405'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'0.004', '0.000', '0.008', '0.016', '0.032'};
        plot_line_ids = [1:5];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 9
        title_name    = '(Model) Fixed disparity';
        Dates         = {'2025071501', '2025071502','2025071503', '2025071504', '2025071505','2025071506',...
                         '2025071507', '2025071508','2025071509'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '4.0', '16.0', '8.0',...
                         '1.0', '0.0', '1.5'};
        plot_line_ids = [1:4 6:9];
        fname_pattern = '%s_cricket_%snoise%s%s';
    case 10
        title_name    = '(Model) Fixed disparity';
        Dates         = {'2025072001', '2025072002','2025072003', '2025072004', '2025072005','2025072006',...
                         };
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'No fixed', '0.5', '2.0', '1.0', '4.0', '8.0'};
        plot_line_ids = [1:6];
        fname_pattern = '%s_cricket_%snoise%s%s';

    case 11
        title_name    = '(Model) obersered predation results';
        Dates         = {'2025080901', '2025080902', '2025080903', '2025080904', '2025081001','2025081002'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'ON-N', 'ON-T', 'OFF-N', 'OFF-T', 'ON-N (syn)', 'OFF-T (half)'};
        plot_line_ids = [1:6];
        fname_pattern = '%s_cricket_%snoise%s%s';

    case 12
        title_name    = '(Model) obersered predation results';
        Dates         = {'2025081201', '2025081202', '2025081203', '2025081204'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'ON-N', 'ON-T', 'OFF-N', 'OFF-T'};
        plot_line_ids = [1:4];
        fname_pattern = '%s_cricket_%snoise%s%s';
    
    case 13
        title_name    = '(Model) obersered predation results (incorrect number)';
        Dates         = {'2025082901', '2025082902', '2025082903', '2025082904'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'OFF-N', 'ON-T', 'OFF-T', 'ON-N'};
        plot_line_ids = [1:4];
        fname_pattern = '%s_cricket_%snoise%s%s';

    case 14
        title_name    = '(Model) obersered predation results';
        Dates         = {'2025083101', '2025083104'};
        Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'OFF-N', 'ON-T', 'OFF-T', 'ON-N'};
        plot_line_ids = [1];
        fname_pattern = '%s_cricket_%snoise%s%s';
    
    case 15
        title_name    = '(Model) obersered predation results';
        Dates         = {'2025091201', '2025091202', '2025091203', '2025091204', '2025091501', '2025091502', '2025091503',};  %, 
        % Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        Noise_level   = {'0.0','0.016','0.032','0.064','0.128','0.256'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'OFF-N', 'ON-T', 'OFF-T', 'ON-N', 'OFF-N (2)', 'ON-T (2)', 'OFF-T (2)'};
        plot_line_ids = [1:7];
        fname_pattern = '%s_cricket_%snoise%s%s';

    case 16
        title_name    = '(Model) obersered predation results';
        Dates         = {'2025091501', '2025091502', '2025091503', '2025091504', '2025091505', '2025091506', '2025091507'};  %, 
        % Noise_level   = {'0.0','0.002','0.004','0.008','0.016','0.032'};
        Noise_level   = {'0.016','0.032','0.064','0.128','0.256'};
        BG_folder     = repmat({'blend_'},1,length(Dates));
        LSTM_layer_n  = {'OFF-N', 'ON-T', 'OFF-T', 'ON-N', 'ON-T (-0.006)', 'ON-T (-0.18)', 'ON-T (-0.36)'};
        plot_line_ids = [2 5:7];
        fname_pattern = '%s_cricket_%snoise%s%s';
    
    otherwise
        error('exp_id must be 1, 2 or 3');
end

%% ——— Shared data‐loading and plotting code ———
[unique_string, ~, string_ids] = unique(LSTM_layer_n);
N_days   = numel(Dates);
N_levels = numel(Noise_level);

Data_m = nan(N_days, N_levels);
Data_s = nan(N_days, N_levels);
Data_t = nan(N_days, n_epoch);
Data_n = nan(N_days, N_levels);

for i = 1:N_days
    for j = 1:N_levels
        % choose bg prefix if provided
        if numel(BG_folder)==N_days
            bg = BG_folder{i};
        else
            bg = '';
        end
        fname = sprintf(fname_pattern, Dates{i}, bg, Noise_level{j}, Tag);
        load(fullfile(Folder_Name, fname), 'test_losses', 'training_losses');
        n_sample      = numel(test_losses);
        Data_n(i, j)  = n_sample;
        Data_m(i, j)  = mean(test_losses);
        Data_s(i, j)  = std(test_losses)/sqrt(n_sample);
    end
    Data_t(i, :) = training_losses(1:n_epoch);
end

x      = 1:N_levels;
colors = lines(numel(unique_string));

hFig = figure;

% ——— Right: Test‐loss vs Noise ———
subplot(1, 2, 2); hold on
for k = 1:numel(plot_line_ids)
    i = plot_line_ids(k);
    e = errorbar(x, Data_m(i, :), Data_s(i, :), 'CapSize', 0, 'LineWidth',1.5);
    e.Color   = colors(string_ids(i), :);
    legs{k}   = unique_string{string_ids(i)};
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
%%
print(hFig, "myHighResPlot.png", "-dpng", "-r300");
%%
print("myHighResPlot.png", "-dpng", "-r300");
