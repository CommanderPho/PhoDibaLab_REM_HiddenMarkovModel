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
[extantFigH, currPlotHandles, currStateMapHandle, plot_outputs] = pho_plot_spikeRaster(temp.curr_active_processing, filter_config, plotting_options, extantFigH, curr_i);


%% Build the scrollable interaction bar that sits below the main raster plot:
scrollHandles = scrollplot(currPlotHandles.linesHandle, 'WindowSizeX', plotting_options.window_duration);

%% Add Blurred Spike Overlays:
numBlurredSpikeOutputs = length(unitStatistics.blurredSpikeOutputs);
unitStatistics.blurredStats.max = cellfun(@max, unitStatistics.blurredSpikeOutputs);
unitStatistics.blurredStats.min = cellfun(@min, unitStatistics.blurredSpikeOutputs);
unitStatistics.blurredStats.range = unitStatistics.blurredStats.max - unitStatistics.blurredStats.min;

%% Axes:

% temp.currRasterAxisPosition = currPlotHandles.axesHandle.Position;
temp.currRasterAxisPosition = currPlotHandles.axesHandle.Position;
% Subdivide its height into equal rectangles
temp.currNumOfHeightSubdivisions = length(plotting_options.trialSelection.TrialBackgroundRects.pos);

temp.subplotHeight = temp.currRasterAxisPosition(3) ./ temp.currNumOfHeightSubdivisions;
% Get the subplot's y-offset for each subplot:



hold on;
for i = 1:numBlurredSpikeOutputs
%     ax(i) = subplot(numBlurredSpikeOutputs,1,i);
    % Need to convert to parent-space:
    
    currSubplotPositionRect = temp.currRasterAxisPosition;
    currSubplotPositionRect(3) = temp.subplotHeight; % Set to the common height
    currSubplotPositionRect(2) = temp.currRasterAxisPosition(2) + ((i-1) * temp.subplotHeight);
    
%     [figureXPoints, figureYPoints] = axescoord2figurecoord(figureXPoints, figureYPoints);
    ax(i) = axes('Position', currSubplotPositionRect,'Color','none');
    % Normalize the blurredSpikeOutputs down to unit height for plotting:
    h(i) = plot(ax(i), seconds(temp.curr_timesteps_array{2}), (unitStatistics.blurredSpikeOutputs{i} ./ unitStatistics.blurredStats.range(i)));
    xlabel(ax(i), [])
    xticks(ax(i), [])
    box off
    yticks(ax(i), [])
    ylabel(ax(i),'')
end

%% Set the current window to the specified range:
% xlim(ax, xlim(scrollHandles.ParentAxesHandle))
linkaxes([currPlotHandles.axesHandle ax],'x'); % Link all blurred axes to the main rasterplot axes


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
% xlim(scrollHandles.ParentAxesHandle)


%% Rectangle Selection:
% plot_outputs.compOutputs.TrialBackgroundRects.pos % [nTrials x 4]
% plot_outputs.compOutputs.TrialBackgroundRects.pos(:, 2) % y-offsets [0.5, 1.5, 2.5, ...]

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


%% Annotations Testing:
[test_annotations] = testAnnotations_phoTrialRectangles();
plotting_options.annotations = test_annotations;

% % Select:
% paramCell = struct2argsList(plotting_options.trialSelection.RectangleProperties.selected);
% 
% % Deselect:
% paramCell = struct2argsList(plotting_options.trialSelection.RectangleProperties.normal);
% 
% set(plot_outputs.compOutputs.TrialBackgroundRects.handles(1), paramCell{:});

% enable interactive section selection
phoSelectionAnnotations(currPlotHandles.linesHandle, plotting_options);

% i = 1;
% annotations{i}.unitIDs = {1, 5, 9}; % The absolute unit ID (original/unfiltered ID) of each unit included
% annotations{i}.referenceTimes = {1080, [1440 1467]}; % The index timestamps of which to make a line-annotation
% % the scalar is an example of a point annotation (to highlight an event) and the second a window annotation (to highlight a range of timestamps)
% annotations{i}.notes = {'', ''};
% 
% obj = RasterplotAnnotation(typeName, startTimestamp, endTimestamp, comment, unitIDs);
% 
% endTimestamp



% figure;
% hold on;





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



function performToggleRectangleSelection_phoTrialRectangles(rectangleHandle, plottingOptions)
    

end

function mouseDownCallback_phoTrialRectangles(src, eventdata)
    
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
    plotting_options.spikeLinesFormat.Color = 'k';


    [plot_outputs.x_points, plot_outputs.y_points, rasterPlotHandles.linesHandle, plot_outputs.compOutputs] = phoPlotSpikeRaster(active_processing.spikes.time(plot_outputs.filter_active_units),'PlotType','vertline','rasterWindowOffset', curr_rasterWindowOffset,'XLimForCell', curr_window, ...
        'TrialBackgroundColors', plotting_options.unitBackgroundColors, ...
        'LineFormat', plotting_options.spikeLinesFormat);

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
    
    %% Add the state_map:
    state_statemapPlottingOptions.orientation = 'horizontal';
    state_statemapPlottingOptions.plot_variable = 'behavioral_state';
    state_statemapPlottingOptions.vertical_state_mode = 'combined';
    state_statemapPlottingOptions.x_axis = 'timestamp'; % Timestamp-appropriate relative bins
    
    % Get the position of the main raster plot axes:
    temp.updated_main_axes_pos = ax.Position;    
    
    %% Puts Above the main raster plot box:
    temp.statemap_pos = temp.updated_main_axes_pos;
    temp.statemap_pos(2) = temp.updated_main_axes_pos(2) + temp.updated_main_axes_pos(4);
    temp.statemap_pos(4) = 0.05;
    
    subplot('Position', temp.statemap_pos);
    [stateMapHandle] = fnPlotStateDiagram(active_processing, state_statemapPlottingOptions);
    
    % Link the state map to the main raster plot. Important for when the raster plot is scrolled after the scrollHandles are added. 
    linkaxes([ax stateMapHandle],'x'); 
    
end
