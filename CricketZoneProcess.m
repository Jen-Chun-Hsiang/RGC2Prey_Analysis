% filepath: \\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\load_and_plot_pngs.m
% Script: load PNG alpha channels, save into a MAT struct with names starting
%         with "processed_x" (x from filename processed_x.png), load a .mat
%         summary file and plot each alpha with a marker at the center
%         given by summary_struct(id).X and summary_struct(id).Y.
%
% Configure folders / file names below before running.

load_png_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\CricketDataset\Images\cropped\cricket';  % folder with processed_x.png files
load_mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\';  % folder containing the MAT file with summary_struct
load_mat_file_name = 'selected_points_summary_body.mat';   % MAT file that contains variable summary_struct
save_mat_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\RISserver\RGC2Prey\';  % folder to save output MAT file
target_cover_perc = 0.85;
% Output MAT that will contain variables/fields starting with 'processed_'
out_mat_file = fullfile(load_png_folder, 'processed_cover_radius.mat');

% Find png files matching processed_x.png
png_files = dir(fullfile(load_png_folder, 'processed_*.png'));
if isempty(png_files)
    warning('No files matching processed_*.png found in %s', load_png_folder);
    return;
end

% Load summary MAT
summary_mat_path = fullfile(load_mat_folder, load_mat_file_name);
if ~exist(summary_mat_path, 'file')
    error('Summary MAT file not found: %s', summary_mat_path);
end
S = load(summary_mat_path);
if ~isfield(S, 'summary_struct')
    error('Loaded MAT does not contain variable ''summary_struct''');
end
summary_struct = S.summary_struct;
nSumm = numel(summary_struct);
image_id_list = [summary_struct.image_id];

% Container to save
file_index_list = [];
processed_cover_radius = [];

for fi = 1:numel(png_files)
    fname = png_files(fi).name;
    fpath = fullfile(load_png_folder, fname);
    % extract integer x from filename processed_x.png using regexp
    tok = regexp(fname, '^processed_(\d+)\.png$', 'tokens', 'once');
    if isempty(tok)
        warning('Skipping file with unexpected name format: %s', fname);
        continue;
    end
    x = str2double(tok{1});
    file_index_list(end+1) = x; %#ok<SAGROW>
    
    % read image and alpha channel
    try
        [~, ~, alpha] = imread(fpath);  % PNG with alpha returns 3rd output
    catch
        % imread may error if file missing or not PNG; skip
        warning('Could not read alpha from %s; skipping.', fpath);
        continue;
    end
    % If no alpha channel returned, make an all-ones alpha
    if isempty(alpha)
        info = imfinfo(fpath);
        alpha = ones(info.Height, info.Width, 'uint8') * 255;
    end
    
    id = find(x == image_id_list, 1);
    % Validate id range
    if id >= 1 && id <= numel(summary_struct)
        % Extract coordinates directly
        Xc = summary_struct(id).X;
        Yc = summary_struct(id).Y;
        % Ensure numeric values
        if ischar(Xc) || isstring(Xc), Xc = str2double(char(Xc)); end
        if ischar(Yc) || isstring(Yc), Yc = str2double(char(Yc)); end
        have_xy = ~(isempty(Xc) || isnan(Xc) || isempty(Yc) || isnan(Yc));
    else
        have_xy = false;
        error('Index x=%d is out of range for summary_struct (numel=%d).', x, numel(summary_struct));
    end
    
    % Plot alpha and overlay center marker if available
    hfig = figure('Name', sprintf('Alpha: %s', fname), 'NumberTitle', 'off');
    imagesc(double(alpha)); colormap(gray); colorbar;
    axis image; title(sprintf('%s  (x=%d)', fname, x), 'Interpreter', 'none');
    hold on;
    if have_xy
        % summary_struct coordinates are used directly
        plot(Xc, Yc, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
        legend(sprintf('summary_struct(%d).X,Y = (%.1f, %.1f)', id, Xc, Yc), 'TextColor', 'r');
    end

    [reach_zone_radius, covered_frac] = compute_reach_zone_radius(double(alpha), Xc, Yc, target_cover_perc);

    processed_cover_radius(end+1) = reach_zone_radius; %#ok<SAGROW>
    % overlay circle at (Xc,Yc) showing reach_zone_radius
    if ~isempty(reach_zone_radius) && ~isnan(reach_zone_radius) && ~isnan(Xc) && ~isnan(Yc)
        theta = linspace(0, 2*pi, 360);
        xr = Xc + reach_zone_radius * cos(theta);
        yr = Yc + reach_zone_radius * sin(theta);
        % prefer viscircles if available, otherwise fall back to plot
        try
            viscircles([Xc, Yc], reach_zone_radius, 'Color', 'r', 'LineWidth', 1);
        catch
            plot(xr, yr, 'r-', 'LineWidth', 1.5);
        end
        % annotate achieved coverage percentage near the circle
        txt = sprintf('%.1f%% covered', covered_frac * 100);
        text(Xc + reach_zone_radius, Yc, [' ' txt], 'Color', 'r', ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    hold off;
    
    drawnow;

    keyboard
end


save(out_mat_file, 'file_index_list', 'processed_cover_radius');

fprintf('Saved alpha variables into %s (variables named processed_<x>)\n', out_mat_file);