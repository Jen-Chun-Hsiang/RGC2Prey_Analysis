% Script to process PNG files, collect coordinates and filenames
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

% Process each PNG file
for i = 1:length(pngFiles)
    fileName = pngFiles(i).name;
    filePath = fullfile(folderPath, fileName);
    img = imread(filePath);
    
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
            % Store the file name and coordinates
            fileNames{end+1} = fileName; %#ok<SAGROW>
            coordinates(end+1, :) = [x, y]; %#ok<SAGROW>
        else
            % Close the figure to reselect
            close(gcf);
        end
    end
    
    % Close the figure
    close(gcf);
end

% Summarize results
summary = table(fileNames', coordinates(:, 1), coordinates(:, 2), ...
                'VariableNames', {'FileName', 'X', 'Y'});

% Display the summary
disp('Summary of selected points:');
disp(summary);

% Optionally save the results
savePath = fullfile(folderPath, 'selected_points_summary.mat');
save(savePath, 'summary');
disp(['Summary saved to: ', savePath]);
