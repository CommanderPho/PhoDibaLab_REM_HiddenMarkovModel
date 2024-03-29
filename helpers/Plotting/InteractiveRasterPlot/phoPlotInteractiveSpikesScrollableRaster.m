% phoPlotInteractiveSpikesScrollableRaster.m
% Pho Hale, 02/05/2021

% Produces a raster plot that displays the spike trains for each unit in a window of customizable length.
%   Spawns an interactive slider that allows you to specify the current window to look at, acting as a paginated manner.

%% Filtering Options:
% filter_config.filter_included_cell_types = {};
filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

% plotting_options.window_duration = 2; % 2 seconds
plotting_options.window_duration = 20; % 10 seconds


%% Sorting/Ordering Units:
plotting_options.sorting_config.unit_sort_indicies = sortedTuningCurveIndicies;


if exist('across_experiment_results','var')
    temp.curr_timesteps_array = across_experiment_results{1, 1}.timesteps_array;  
    temp.curr_active_processing = across_experiment_results{1, 1}.active_processing;
% elseif exist('spikeStruct','var')
%     temp.curr_timesteps_array = spikeStruct.t;  
%     temp.curr_active_processing = active_processing;

else
    temp.curr_timesteps_array = timesteps_array;  
    temp.curr_active_processing = active_processing;
end

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
[extantFigH, currRasterPlotHandles, currStateMapHandle, plot_outputs] = pho_plot_spikeRaster(temp.curr_active_processing, filter_config, plotting_options, extantFigH, curr_i);
currPlotHandles = currRasterPlotHandles;

%% Add Position Curves if possible:
plotting_options.subplots.position_curves = true;
if plotting_options.subplots.position_curves
    if exist('fileinfo','var')
        [t, t_rel, x, y, linearPos] = phoPlotInteractiveRasterExtras.processFileInfoPositionExtendedExtras(fileinfo);
        [~, ~] = phoPlotInteractiveRasterExtras.addPositionSubplot(seconds(t_rel), ...
            {linearPos}, ...
            currPlotHandles.axesHandle, currPlotHandles, ...
            plot_outputs.other_subplots_rect);
    else
        fprintf('variable "fileinfo" does not exist, skipping plots of position data\n');
    end
end

%% Add Ripple Periods if possible:
plotting_options.subplots.ripple_periods = false;
if plotting_options.subplots.ripple_periods
    if exist('source_data','var')
        % t_rel is aligned with the timestamps in the active_processing.position_table's timestamp column
        ripples_time_mat = source_data.ripple.RoyMaze1.time;
        warning('RoyMaze1 is hardcoded for ripple periods!')
    %     ripples_time_mat = seconds((ripples_time_mat - source_data.behavior.RoyMaze1.time(1,1)) ./ 1e6); % Convert to relative timestamps since start
        ripples_time_mat = (ripples_time_mat - source_data.behavior.RoyMaze1.time(1,1)) ./ 1e6; % Convert to relative timestamps since start
    
        [currPlotHandles.ripplesHandle] = phoPlotInteractiveRasterExtras.addRipplePeriodsSubplot(ripples_time_mat, ...
            currPlotHandles.axesHandle, ...
            plot_outputs.other_subplots_rect); % Overlay the subplots
    %         plot_outputs.mainplot_rect); % Overlay the existing main plot
    
        % [spike_linearPos, spike_unitSpikeLaps, spike_unitIsRippleSpike] = phoPlotInteractiveRasterExtras.processSpikeStructExtendedExtras(spikeStruct);
    
    else
        fprintf('variable "source_data" does not exist, skipping plots of position data\n');
    end
end


%% Build the scrollable interaction bar that sits below the main raster plot:
scrollHandles = scrollplot(currPlotHandles.linesHandle, 'WindowSizeX', plotting_options.window_duration);
% scrollHandles = scrollplot(currStateMapHandle, 'WindowSizeX', plotting_options.window_duration);
% scrollHandles = scrollplot({currStateMapHandle currPlotHandles.ripplesHandle, currPlotHandles.linesHandle}, 'WindowSizeX', plotting_options.window_duration);



%% Add Blurred Spike Overlays if possible:
plotting_options.subplots.blurred_spike_overlays = false;
if plotting_options.subplots.blurred_spike_overlays
    if exist('unitStatistics','var')
        phoPlotInteractiveRasterExtras.addBlurredSpikeOverlays(unitSpikeCells, scrollHandles, currPlotHandles, plotting_options);
    else
        warning('unitStatistics variable does not exist, so cannot overlay the blurred spike outputs. Continuing.')
    end % end if exist('unitStatistics','var')
end

%% INFO:
% Can get scroll handles (the blue adjustment handle positions) using
% [scrollHandles.ScrollMin, scrollHandles.ScrollMax];

% curr_window_length = scrollHandles.ScrollMax - scrollHandles.ScrollMin;


% scrollHandles.ScrollPatchHandle: Patch (scrollPatch) object
% scrollHandles.ScrollPatchHandle.Vertices: [4x2 double] like:

%% Set the current window to the specified range:
% xlim(scrollHandles.ParentAxesHandle, [1800 1860])
xlim(scrollHandles.ParentAxesHandle, [1800 (1800 + plotting_options.window_duration)]);



%% Rectangle Selection:
% plot_outputs.compOutputs.TrialBackgroundRects.handles % Handles to the Rectangle objects that represent each trial row and span all of the data

plotting_options.lineVisibleColor = [0.95 0.95 0.95];

plotting_options.trialSelection.RectangleProperties.normal.EdgeColor = [0.00,0.00,0.00];
plotting_options.trialSelection.RectangleProperties.normal.LineStyle = 'none';
plotting_options.trialSelection.RectangleProperties.normal.LineWidth = 0.5;

plotting_options.trialSelection.RectangleProperties.selected.EdgeColor = [0.00,1.00,0.00];
plotting_options.trialSelection.RectangleProperties.selected.LineStyle = '-.';
plotting_options.trialSelection.RectangleProperties.selected.LineWidth = 2;

% Set the actual background rectangle properties for access in the phoSelectionAnnotations(...) functions:
plotting_options.trialSelection.TrialBackgroundRects = plot_outputs.compOutputs.TrialBackgroundRects;

%% Add Annotations Testing if possible:
plotting_options.subplots.annotations_testing = false;
if plotting_options.subplots.annotations_testing
    [test_annotations] = testAnnotations_phoTrialRectangles();
    plotting_options.annotations = test_annotations;
end

% % Select:
% paramCell = struct2argsList(plotting_options.trialSelection.RectangleProperties.selected);
% 
% % Deselect:
% paramCell = struct2argsList(plotting_options.trialSelection.RectangleProperties.normal);
% 
% set(plot_outputs.compOutputs.TrialBackgroundRects.handles(1), paramCell{:});



%% Add enable interactive section selection:
plotting_options.subplots.interactive_section_selection = false;
if plotting_options.subplots.interactive_section_selection
    %enable interactive section selection
    phoSelectionAnnotations(currPlotHandles.linesHandle, plotting_options);
end


function [annotations] = testAnnotations_phoTrialRectangles()
    %% testAnnotations_phoTrialRectangles: build test annotations
    % annotations: cell array of RasterplotAnnotation objects.
    annotationLists.unitIDs = [{1, 5, 9}, {2, 3, 4, 7}]; % The absolute unit ID (original/unfiltered ID) of each unit included
    annotationLists.referenceTimes = [1080, nan; 1440, 1467]; % The index timestamps of which to make a line-annotation
    % the scalar is an example of a point annotation (to highlight an event) and the second a window annotation (to highlight a range of timestamps)
    annotationLists.notes = {'', ''};
    
    for i = 1:length(annotationLists.notes)
        startTimestamp = annotationLists.referenceTimes(i, 1);
        endTimestamp = annotationLists.referenceTimes(i, 2);
        comment = annotationLists.notes{i};
        typeName = 'Temp';
        annotations{i} = RasterplotAnnotation(typeName, startTimestamp, endTimestamp, comment, annotationLists.unitIDs{i});
    end

end

% Instructions for use:
% - Double click the left mouse button to place the first mark line to a 
%   new possition;
% - Double click the right mouse button to place the second mark line to a 
%   new possition;
% - Press "Enter" key in order to save the selected 2D plot section as a
%   mat. data file;
% - Press "End" key in order to close the interactive section selection
%   regime.


% %% Testing for figure sizing:
% % Get figure width:
% temp.fig_width = extantFigH.Position(3);
% 
% xlim(scrollHandles.ParentAxesHandle)
% 
% scrollHandles.Position
% % Postition: [0.13 0.11 0.775 0.0815]
% 
% scrollHandles.XAxis

% %% Pho Scrollable Mode:
% slider_controller = fnBuildCallbackInteractiveSliderController(iscInfo, @(curr_i) (pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i)) );
% 
%% Plot function called as a callback on update
function [plotted_figH, rasterPlotHandles, stateMapHandle, plot_outputs] = pho_plot_spikeRaster(active_processing, filter_config, plotting_options, extantFigH, curr_windowIndex)
    %% pho_plot_spikeRaster

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
    % Note that we set the unitBackgroundColorOpacity here by adding a row for the alpha value
    plotting_options.unitBackgroundColorAlpha = 0.0;
    plotting_options.unitBackgroundColors = [plotting_options.unitBackgroundColors; plotting_options.unitBackgroundColorAlpha .* ones([1 size(plotting_options.unitBackgroundColors, 2)])];

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
    
    % Can specify line colors here, default is black ('k'):
%     plotting_options.spikeLinesFormat.Color = 'k';

    plotting_options.spikeLinesFormat.Color = [1.0, 0.0, 0.0]';

    %% Test filtering spikes by ripple criteria
%     active_processing.spikes.time(plot_outputs.filter_active_units)
%     active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units)


%     fnFilterSpikesWithCriteria(active_processing.spikes, {}, {});
%     active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes


    active_spike_times = active_processing.spikes.time(plot_outputs.filter_active_units); %% Original
    
    %% Re-ordering/sorting
    % Re-order the units in the spike table
    active_spike_times = active_spike_times(plotting_options.sorting_config.unit_sort_indicies);
    plotting_options.unitBackgroundColors(:, plotting_options.sorting_config.unit_sort_indicies);
    
    [plot_outputs.x_points, plot_outputs.y_points, rasterPlotHandles.linesHandle, plot_outputs.compOutputs] = phoPlotSpikeRaster(active_spike_times, ...
        'PlotType','vertline', ...
        'rasterWindowOffset', curr_rasterWindowOffset, ...
        'XLimForCell', curr_window, ...
        'TrialBackgroundColors', plotting_options.unitBackgroundColors, ...
        'LineFormat', plotting_options.spikeLinesFormat);
%         'SpikeColors', [1.0, 0.0, 0.0]');

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
    rasterPlotHandles.axesHandle = ax;
    ax.YGrid = 'on';
%     ax.YMinorGrid = 'on';
    yticks(ax, 1:temp.num_active_units);
    yticklabels(ax, num2cellstr(plotting_options.sorting_config.unit_sort_indicies));
    

    % Resize the main plot to prepare for adding various subplots in the margins (such as the sleep_state plot, and the position plot
    [plot_outputs.mainplot_rect, plot_outputs.subplots_rect] = phoPlotInteractiveRasterExtras.reallocateForAddingSubplots(ax, 0.10);

    plot_outputs.state_subplots_rect = plot_outputs.subplots_rect;
    plot_outputs.other_subplots_rect = plot_outputs.subplots_rect;
    
    plot_outputs.state_subplots_rect(4) = plot_outputs.subplots_rect(4) * 0.5;
    plot_outputs.other_subplots_rect(4) = plot_outputs.subplots_rect(4) * 0.5;
    
    plot_outputs.state_subplots_rect(2) = plot_outputs.state_subplots_rect(2) + plot_outputs.other_subplots_rect(4); % Shift up the state subplot rect above the other one

    %% Add the state_map above the raster plot:
    plotting_options.subplots.behavioral_state_map = true;
    if plotting_options.subplots.behavioral_state_map
        [stateMapHandle] = phoPlotInteractiveRasterExtras.addStateMapSubplot(active_processing, ax, plot_outputs.state_subplots_rect);
    else
         stateMapHandle = [];
    end

    % Add extra interface elements
%     [extraControlsPanel] = phoPlotInteractiveRasterExtras.addExtraControls(plotted_figH, rasterPlotHandles.axesHandle);

%     [linearPos, unitSpikeLaps, unitIsRippleSpike] = phoPlotInteractiveRasterExtras.processSpikeStructExtendedExtras(spikeStruct)

%     [t, ~, ~, linearPos] = phoPlotInteractiveRasterExtras.processFileInfoPositionExtendedExtras(fileinfo);
%     [~, ~] = phoPlotInteractiveRasterExtras.addPositionSubplot(seconds(t), ...
%         {linearPos}, ...
%         ax, rasterPlotHandles, ...
%         plot_outputs.other_subplots_rect);

%     [~, ~] = phoPlotInteractiveRasterExtras.addPositionSubplot(seconds(active_processing.position_table.timestamp), ...
% %         {active_processing.position_table.x, active_processing.position_table.y}, ...
%         {active_processing.position_table., active_processing.position_table.y}, ...
%         ax, rasterPlotHandles, ...
%         other_subplots_rect);

end
