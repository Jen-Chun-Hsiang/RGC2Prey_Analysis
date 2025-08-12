close all; clc; clear
%%
FilNam = ['\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey',...
    '\RunningTimeEvaluation.xlsx'];
DocTable = readtable(FilNam, 'Sheet', 'CricketHunting', 'Range', 'B:AC');


%%
rmids = [18 19];
DocTable(rmids, :) = [];
data = DocTable(:, {'mov_avg', 'mov_avg_1', 'mov_avg_2', 'mov_avg_3','batch_size'});
data = table2array(data);
data = data(:, 1:4)./data(:, 5);
worker = table2array(DocTable(:, 'n_worker'));

n_data = size(data, 1);

% Define the colors for each column in the stack
columnColors = [
    159 2 81;
    236 121 154; 
    21 52 80;
    68 114 148]/255;

% Create the stacked bar chart
figure;
hold on;
h = bar(data, 'stacked');

% Apply the colors to each stack
for i = 1:length(h)
    h(i).FaceColor = 'flat'; % Set to allow individual bar color control
    h(i).CData = repmat(columnColors(i, :), size(data, 1), 1);
    h(i).EdgeColor = 'w';
end

xlabel('Exp (#)');
ylabel('Single data Time (s)');
xlim([0.5 n_data+0.5])

yyaxis right
plot(worker, 'k');
ylabel('number of workers')

title('Cricket Hunting');
legend({'Data loading', 'Data transfer', 'Forward pass', 'Backpropagation', 'Worker'}, 'Location', 'best');
hold off;

%%
DocTable = readtable(FilNam, 'Sheet', 'RetinalPerceiver', 'Range', 'B:V');
% Extract numeric and string parts if needed:
%DocNum = DocTable{:, :}; % Access numeric data

%ColName2Ind = @(String) find(strcmpi(DocStr(1, :), String));


%%
data = DocTable(:, {'mov_avg', 'mov_avg_1', 'mov_avg_2', 'mov_avg_3','batch_size'});
data = table2array(data);
data = data(:, 1:4)./data(:, 5);

worker = table2array(DocTable(:, 'n_worker'));
n_data = size(data, 1);

% Define the colors for each column in the stack
columnColors = [
    159 2 81;
    236 121 154; 
    21 52 80;
    68 114 148]/255;

% Create the stacked bar chart
figure;
hold on;
h = bar(data, 'stacked');

% Apply the colors to each stack
for i = 1:length(h)
    h(i).FaceColor = 'flat'; % Set to allow individual bar color control
    h(i).CData = repmat(columnColors(i, :), size(data, 1), 1);
    h(i).EdgeColor = 'w';
end
xlabel('Exp (#)');
ylabel('Single data Time (s)');
xlim([0.5 n_data+0.5])
yticks(0:0.1:0.5)
yticklabels({'0', '', '0.2', '', '0.4', ''})

yyaxis right
plot(worker, 'k');
ylabel('number of workers')
yticks(0:4:16)
yticklabels({'0', '', '8', '', '16'})
% Add labels and title

legend({'Data loading', 'Data transfer', 'Forward pass', 'Backpropagation', 'Worker'}, 'Location', 'best');

title('RetinalPerceiver (MEA recording)');

hold off;


%% relative improvement
SaveFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\VideoSpikeDataset\Summary\cosyne2025\';
FleNam = sprintf('%sFigure_speed_up_process', SaveFolder);
print('-depsc','-painters','-loose', '-r300',FleNam)
saveas(gcf,[FleNam '.png']);
