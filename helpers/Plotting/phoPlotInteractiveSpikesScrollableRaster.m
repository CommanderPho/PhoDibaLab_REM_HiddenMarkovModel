% phoPlotInteractiveSpikesScrollableRaster.m
% Pho Hale, 02/05/2021

% Produces a raster plot that displays the spike trains for each unit in a window of customizable length.
%   Spawns an interactive slider that allows you to specify the current window to look at, acting as a paginated manner.

%% Filtering Options:
filter_config.filter_included_cell_types = {};
% filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

plotting_options.window_duration = 10; % 10 seconds


temp.curr_timesteps_array = across_experiment_results{1, 1}.timesteps_array;  
temp.curr_active_processing = across_experiment_results{1, 1}.active_processing;


% define cell type ({'pyramidal', 'contaminated', 'interneurons'}) colors: 
if ~isfield(temp.curr_active_processing.definitions.speculated_unit_info, 'classColors')
    temp.curr_active_processing.definitions.speculated_unit_info.classColors = [0.8, 0.5, 0.1
       0.5, 0.1, 0.1
       0.0, 0.7, 0.7
     ];
end

plotting_options.active_timesteps = temp.curr_timesteps_array{1};
plotting_options.total_duration = plotting_options.active_timesteps(end) - plotting_options.active_timesteps(1);
% plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;
plotting_options.showOnlyAlwaysStableCells = true;

plotting_options.num_windows = seconds(ceil(plotting_options.total_duration / plotting_options.window_duration));

% Build a new slider controller
iscInfo.slider_identifier = 'PhoScrollableSpikeRasterPlot';
iscInfo.curr_i = 1;
iscInfo.NumberOfSeries = plotting_options.num_windows;


extantFigH = figure(12);

%% Scrollplot mode:
curr_i = 1;
[extantFigH, currPlotHandle, currStateMapHandle, plot_outputs] = pho_plot_spikeRaster(temp.curr_active_processing, filter_config, plotting_options, extantFigH, curr_i);
scrollHandles = scrollplot(currPlotHandle, 'WindowSizeX', plotting_options.window_duration);

% linkaxes([currPlotHandle currStateMapHandle],'x'); 

%% INFO:
% Can get scroll handles (the blue adjustment handle positions) using
% [scrollHandles.ScrollMin, scrollHandles.ScrollMax];


% scrollHandles.ScrollPatchHandle: Patch (scrollPatch) object
% scrollHandles.ScrollPatchHandle.Vertices: [4x2 double] like:
%     [57     0
%     57    80
%     67    80
%     67     0]

% scrollHandles.ScrollSideBarHandles: [1x2 Line array]
% scrollHandles.ScrollSideBarHandles(1).Color = [1 0 0];


%% Set the current window to the specified range:
xlim(scrollHandles.ParentAxesHandle, [1800 1860])

% xlim(currStateMapHandle, xlim(scrollHandles.ParentAxesHandle))

%% Get the current window:
xlim(scrollHandles.ParentAxesHandle)



% %% Pho Scrollable Mode:
% slider_controller = fnBuildCallbackInteractiveSliderController(iscInfo, @(curr_i) (pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i)) );
% 
%% Plot function called as a callback on update
function [plotted_figH, rasterPlotHandle, stateMapHandle, plot_outputs] = pho_plot_spikeRaster(active_processing, filter_config, plotting_options, extantFigH, curr_windowIndex)
    

    %% Get filter info for active units
    [plot_outputs.filter_active_units, plot_outputs.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, plotting_options.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
        filter_config.filter_maximum_included_contamination_level);
    
    temp.num_active_units = sum(plot_outputs.filter_active_units, 'all');
    fprintf('Filter: Including %d of %d total units\n', temp.num_active_units, length(plot_outputs.filter_active_units));

    
    
    %% Build colors
%     active_processing.definitions.speculated_unit_info.classColors
    
    temp.active_units_speculated_type = double(active_processing.spikes.speculated_unit_type(plot_outputs.filter_active_units));
        
%     plotting_options.unitBackgroundColors = [1.0, 1.0, 1.0;  0.9, 0.9, 0.9]';
    plotting_options.unitBackgroundColors = active_processing.definitions.speculated_unit_info.classColors(temp.active_units_speculated_type,:)';
    
    if exist('extantFigH','var')
        plotted_figH = figure(extantFigH); 
    else
        plotted_figH = figure(12);
    end
    clf
    curr_rasterWindowOffset = double(curr_windowIndex - 1) * plotting_options.window_duration;
    curr_window = [0 plotting_options.window_duration];
    curr_window = curr_window + curr_rasterWindowOffset;
    % RelSpikeStartTime: x-axis offset from the left edge of the plot
    % rasterWindowOffset: x-axis window start time
%     hold off
    

    [plot_outputs.x_points, plot_outputs.y_points, rasterPlotHandle] = phoPlotSpikeRaster(active_processing.spikes.time(plot_outputs.filter_active_units),'PlotType','vertline','rasterWindowOffset', curr_rasterWindowOffset,'XLimForCell', curr_window, ...
        'TrialBackgroundColors', plotting_options.unitBackgroundColors);

    if plotting_options.showOnlyAlwaysStableCells
        numAlwaysStableCells = sum(plot_outputs.filter_active_units, 'all');
        ylabel('Stable Unit Index')
        title(sprintf('Spike Train for Always Stable Units (%d of %d total)', numAlwaysStableCells, length(active_processing.spikes.time)));
    else
        ylabel('Unit Index')
        title(sprintf('Spike Train for %d second window from [%d, %d]', plotting_options.window_duration, curr_window(1), curr_window(end)));
    end
    xlabel('Time [seconds]')
    
    ax = gca;
    ax.YGrid = 'on';
%     ax.YMinorGrid = 'on';
    yticks(ax, 1:temp.num_active_units);
    
    %% Add the state_map:
    state_statemapPlottingOptions.orientation = 'horizontal';
    state_statemapPlottingOptions.plot_variable = 'behavioral_state';
    state_statemapPlottingOptions.vertical_state_mode = 'combined';
    state_statemapPlottingOptions.x_axis = 'timestamp'; % Timestamp-appropriate relative bins
    
    % Get the position of the main raster plot axes:
    temp.updated_main_axes_pos = ax.Position;    
    
    %% Puts Above
    temp.statemap_pos = temp.updated_main_axes_pos;
    temp.statemap_pos(2) = temp.updated_main_axes_pos(2) + temp.updated_main_axes_pos(4);
    temp.statemap_pos(4) = 0.05;
    
    subplot('Position', temp.statemap_pos);
    [stateMapHandle] = fnPlotStateDiagram(active_processing, state_statemapPlottingOptions);
    
    % Link the state map to the main raster plot. Important for when the raster plot is scrolled after the scrollHandles are added.
    linkaxes([ax stateMapHandle],'x'); 
    
end
