
close all; clear; clc;
FolderName = uigetdir;
files = dir(FolderName);
expression = '_.mat'; % Grk4_GCamP6s
replace = '.mat'; % Grk4_GCaMP6
for i = 1:length(files)
    if strfind(files(i).name, expression)
        NewFileNames = regexprep(files(i).name, expression, replace);
        movefile(fullfile(FolderName, files(i).name), fullfile(FolderName, NewFileNames));
    end
end

disp('Finished');

%% (1) After copy
close all; clear all; clc;
FolderName = uigetdir;
files = dir(FolderName);
target = 'Ai148_AAV-Grm6Cre_NatMovPred_';
expression = '.mat.mat';
replace = '.mat'; 
for i = 1:length(files)
    if strfind(files(i).name, target)
        NewFileNames = regexprep(files(i).name, expression, replace);
        movefile(fullfile(FolderName, files(i).name), fullfile(FolderName, NewFileNames));
    end
end

disp('Finished');
%% (2) Remove copy
close all; clear all; clc;
FolderName = uigetdir;
files = dir(FolderName);
target = '2025091804';
expression = '2025091804';
replace = 'c2_2025091804'; 
for i = 1:length(files)
    if strfind(files(i).name, target)
        NewFileNames = regexprep(files(i).name, expression, replace);
        movefile(fullfile(FolderName, files(i).name), fullfile(FolderName, NewFileNames));
    end
end

disp('Finished');



