% plot_three_traces_gradient.m
% Plot 3 same-length traces with a color gradient (light -> dark) along each trace.

% Parameters
N = 200;                        % number of points per trace (same length)
t = linspace(0, 4*pi, N);       % common parameter for generating traces
alpha_start = 0.25;             % how "light" the start is (0..1)
alpha_end   = 1.0;              % how "dark" the end is (0..1)
is_visual_degree = 1;
vis_deg_to_cm = 32.5;
legendLabels = {'Ground truth', 'Model prediction', 'Center of Mass'};
% Example traces (same length). Replace with your own data if needed.
x1 = all_paths_r(trial_id, :, 1)*real_dim(1);
y1 = all_paths_r(trial_id, :, 2)*real_dim(2);

x2 = all_paths_pred_r(trial_id, :, 1)*real_dim(1);
y2 = all_paths_pred_r(trial_id, :, 2)*real_dim(2);

x3 = all_path_cm(trial_id, :, 1)*cm_dim_scale;
y3 = all_path_cm(trial_id, :, 2)*cm_dim_scale;

if is_visual_degree
    vis_scale = 1/vis_deg_to_cm; % Convert cm to visual degrees
else
    vis_scale = 1; % No scaling
end
traces = { [x1; y1]*vis_scale, [x2; y2]*vis_scale, [x3; y3]*vis_scale };

% Base (target) colors for each trace (RGB)
baseColors = [ 120, 120, 120;   % red-ish
               255, 0, 255;   % green-ish
               0, 255, 255]/255; % blue-ish

figure('Color','w');
hold on;
axis tight;
box on;

% For legend proxies
legendHandles = gobjects(3,1);


% Precompute alpha ramp for segments (N-1 segments for N points)
segAlpha = linspace(alpha_start, alpha_end, N-1);

for k = 1:3
    XY = traces{k};
    x = XY(1,:);
    y = XY(2,:);
    base = baseColors(k,:);
    % Plot each segment with interpolated color = mix(base, white, alpha)
    % color = (1-alpha)*white + alpha*base  -> light->dark as alpha grows
    for i = 1:(N-1)
        alpha = segAlpha(i);
        color = (1 - alpha) * [1 1 1] + alpha * base;
        plot(x(i:i+1), y(i:i+1), '-', 'Color', color, 'LineWidth', 2.5);
    end
    % Add a proxy line (darkest color) for legend
    legendHandles(k) = plot(nan, nan, '-', 'Color', (1 - alpha_end)*[1 1 1] + alpha_end*base, 'LineWidth', 2.5);
end

% Optional: mark start and end points
markerSize = 60;
for k = 1:3
    XY = traces{k};
    startColor = (1 - segAlpha(1)) * [1 1 1] + segAlpha(1) * baseColors(k,:);
    endColor   = (1 - segAlpha(end)) * [1 1 1] + segAlpha(end) * baseColors(k,:);
    scatter(XY(1,1), XY(2,1), markerSize/2, startColor, 'filled', 'MarkerEdgeColor','k');
    scatter(XY(1,end), XY(2,end), markerSize/2, endColor,   'filled', 'MarkerEdgeColor','k');
end

x_min = -120 * cm_dim_scale * vis_scale;
x_max =  120 * cm_dim_scale * vis_scale;
y_min = -90  * cm_dim_scale * vis_scale;
y_max =  90  * cm_dim_scale * vis_scale;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1);
if is_visual_degree
    plot(-20*ones(1, 2), [-20 -10], '-k', 'LineWidth', 2);
    text(-22, -15, '10 deg', 'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    plot([-20 -10], -20*ones(1, 2), '-k', 'LineWidth', 2);
    text(-15, -15, '10 deg', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end


legend(legendHandles, legendLabels, 'Location', 'best');
xlabel('X');
ylabel('Y');
title('Three traces with color gradient (light -> dark)');
set(gca, 'FontSize', 12);
axis equal;
hold off;
box off
axis off
drawnow;

%%
% save_file_name = fullfile(fig_save_folder, sprintf('PredictionTrace_%s_%s_%d', exp_name, noise_level, trial_id));
% print(gcf, [save_file_name '.eps'], '-depsc', '-painters'); % EPS format
% print(gcf, [save_file_name '.png'], '-dpng', '-r300'); % PNG, 600 dpi