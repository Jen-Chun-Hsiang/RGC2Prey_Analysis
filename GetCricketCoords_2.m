% Script to process PNG files, collect coordinates, dimensions, and filenames
clc;
clear;

% Prompt user to select a folder
folderPath = uigetdir(pwd, 'Select Folder Containing PNG Files');
if folderPath == 0
    disp('No folder selected. Exiting script.');
    return;
end

% Get all PNG files in the folder
pngFiles = dir(fullfile(folderPath, '*.png'));
if isempty(pngFiles)
    disp('No PNG files found in the selected folder. Exiting script.');
    return;
end

% Initialize summary variables
fileNames = {};
coordinates = [];
imageDimensions = []; % To store [Height, Width] of each image

% Process each PNG file
for i = 1:length(pngFiles)
    fileName = pngFiles(i).name;
    filePath = fullfile(folderPath, fileName);
    img = imread(filePath);
    [height, width, ~] = size(img); % Get image dimensions
    
    selected = false;
    while ~selected
        % Display the image with the filename as the title
        figure;
        imshow(img);
        title(fileName, 'Interpreter', 'none');
        
        % Prompt user to select a point
        disp(['Processing file: ', fileName]);
        disp('Click on a point in the image to select its coordinates.');
        [x, y] = ginput(1);
        
        % Show the selected point on the image
        hold on;
        plot(x, y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        
        % Ask the user to confirm or reselect
        choice = questdlg('Do you want to proceed with this point or reselect?', ...
                          'Confirm Selection', ...
                          'Proceed', 'Reselect', 'Proceed');
        
        if strcmp(choice, 'Proceed')
            selected = true;
            % Store the file name, coordinates, and dimensions
            fileNames{end+1} = fileName; %#ok<SAGROW>
            coordinates(end+1, :) = [x, y]; %#ok<SAGROW>
            imageDimensions(end+1, :) = [height, width]; %#ok<SAGROW>
        else
            % Close the figure to reselect
            close(gcf);
        end
    end
    
    % Close the figure
    close(gcf);
end

%% Extract numbers from filenames
extractedNumbers = cellfun(@(name) str2double(regexp(name, '\d+', 'match', 'once')), fileNames);

% Summarize results
summary = table(fileNames', ...
                coordinates(:, 1), coordinates(:, 2), ...
                imageDimensions(:, 1), imageDimensions(:, 2), ...
                extractedNumbers', ...
                'VariableNames', {'FileName', 'X', 'Y', 'Height', 'Width', 'image_id'});
summary.coord_x = summary.X-summary.Width*0.5;
summary.coord_y = summary.Y-summary.Height*0.5;

summary_table = summary;
summary_struct = table2struct(summary);
summary_array = [summary.image_id summary.coord_x summary.coord_y];
summary_array_name = {'image_id', 'coord_x', 'coord_y'};
% Display the summary
disp('Summary of selected points:');
disp(summary);

%% 
% Save the results to a file in the same location as this script
save_file_name = 'selected_points_summary_body.mat';
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey';
save(fullfile(save_folder, save_file_name), 'summary', 'summary_table', 'summary_struct', ...
    'summary_array', 'summary_array_name');


