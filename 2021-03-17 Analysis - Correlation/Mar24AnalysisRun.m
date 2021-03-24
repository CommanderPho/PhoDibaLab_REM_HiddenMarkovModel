% 'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run'
data_config.output.March24OutputFilepath = '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/ResultsTemp/March 24/smoothed.mat';

temp.smoothedVariablesList = {'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run'};

if ~exist('gau_sdf','var')
    %% Load from the file:
    load(data_config.output.March24OutputFilepath, 'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run');
end

% plot(tbin_centers, gau_sdf,'g');
plotting_options.use_panel = true;

%% Customize the plotting command to use (stem, plot, area, etc):
% active_plot_cmd = @(ax,x,y) stem(ax, x, y);
% active_plot_cmd = @(ax,x,y) plot(ax, x, y);
% plotting_options.active_plot_cmd = @(x,y) plot(x, y);
plotting_options.active_plot_cmd = @(x,y) reduce_plot(x, y);


plotting_options.activeFigure = figure(19);
clf;

num_cells = size(gau_sdf,1);

subset_cells = 1:32;
% Max num is 32 for panel
num_subset_cells = length(subset_cells);

plotting_options.absolute_y_max = max(gau_sdf(subset_cells,:),[],'all');
plotting_options.absolute_y_min = min(gau_sdf(subset_cells,:),[],'all');
plotting_options.y_lims = [plotting_options.absolute_y_min, plotting_options.absolute_y_max];

if plotting_options.use_panel
    p = panel(plotting_options.activeFigure);
    p.margin = [13 10 2 2];
    p.pack(num_subset_cells, 1);
    p.de.margin = 2;
    
    plotting_results.axes_objects = gobjects(num_subset_cells,1);
    
    for i = 1:num_subset_cells
       p(i,1).select(); % Activate the panel (subplot)
       plotting_results.axes_objects(i) = gca;
       active_cell_index = subset_cells(i); % get the appropriate index from the list of subsets.
       plotting_options.active_plot_cmd(tbin_centers,  gau_sdf(active_cell_index,:));
       axis([tbin_centers(1) tbin_centers(end) plotting_options.y_lims(1) plotting_options.y_lims(2)]);
       set(plotting_results.axes_objects(i), 'xtick', [], 'ytick', []);
    end
    linkaxes(plotting_results.axes_objects, 'x');
    
%     p(1).xlabel('time (unitless)');
    sgtitle('Spike Density Function (SDF)');
else
    fh = splot(tbin_centers, 1:size(gau_sdf,1), gau_sdf');
end


% plotting_options.window_duration = 10;
% %% Build the scrollable interaction bar that sits below the main raster plot:
% scrollHandles = scrollplot(fh, 'WindowSizeX', plotting_options.window_duration);