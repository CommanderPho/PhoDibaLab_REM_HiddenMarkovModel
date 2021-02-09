% phoPlotInteractiveSpikesScrollableRaster.m
% Pho Hale, 02/05/2021

% Produces a raster plot that displays the spike trains for each unit in a window of customizable length.
%   Spawns an interactive slider that allows you to specify the current window to look at, acting as a paginated manner.

plotting_options.window_duration = 3; % 10 seconds

plotting_options.active_timesteps = timesteps_array{1};
plotting_options.total_duration = plotting_options.active_timesteps(end) - plotting_options.active_timesteps(1);
plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;


% temp.num_whole_windows = floor(plotting_options.total_duration / plotting_options.window_duration);
% temp.num_partial_

plotting_options.num_windows = seconds(ceil(plotting_options.total_duration / plotting_options.window_duration));


% Build a new slider controller
iscInfo.slider_identifier = 'PhoScrollableSpikeRasterPlot';
iscInfo.curr_i = 1;
iscInfo.NumberOfSeries = plotting_options.num_windows;

% timesteps_array = cellfun((@(dt) seconds(active_processing.behavioral_epochs.start_seconds(1):dt:active_processing.behavioral_epochs.end_seconds(end))), ...
%  processing_config.step_sizes, 'UniformOutput', false);






% plotFigureStates = {};

% plotFigureStates{end+1} = PlotFigureState('isc2DPlot', should_show_2d_plot, ...
%     @(extantFigH, curr_i) (pho_plot_2d(final_data_explorer_obj.dateStrings, final_data_explorer_obj.uniqueAmps, final_data_explorer_obj.uniqueFreqs, final_data_explorer_obj.finalOutPeaksGrid, final_data_explorer_obj.multiSessionCellRoi_CompListIndicies, extantFigH, curr_i)));
% 
% plotFigureStates{end+1} = PlotFigureState('iscStimulusTracesPlot', should_show_stimulus_traces_plot, ...
%     @(extantFigH, curr_i) (pho_plot_stimulus_traces(final_data_explorer_obj, extantFigH, curr_i)));

% plotFigureStates{end+1} = PlotFigureState('isc2DPlot', should_show_2d_plot, ...
%     @(extantFigH, curr_i) (plot_manager_cellRoiPlot.pho_plot_2d(curr_i)));

% slider_controller = PhoInteractiveScrollableSpikeRasterPlotManager.getInstance(iscInfo, plot_manager_cellRoiPlot, valid_only_quality');
% slider_controller = PhoInteractiveCallbackSliderDefault.getInstance(iscInfo, plot_manager_cellRoiPlot, valid_only_quality');


extantFigH = figure(3);


%% Scrollplot mode:
curr_i = 1;
[extantFigH, currPlot] = pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i);
scrollHandles = scrollplot(currPlot, 'WindowSizeX', plotting_options.window_duration);

% %% Pho Scrollable Mode:
% slider_controller = fnBuildCallbackInteractiveSliderController(iscInfo, @(curr_i) (pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i)) );
% 
%% Plot function called as a callback on update
function [plotted_figH, plotHandle] = pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_windowIndex)
    if exist('extantFigH','var')
        plotted_figH = figure(extantFigH); 
    else
        plotted_figH = figure(3);
    end
    clf
    curr_rasterWindowOffset = double(curr_windowIndex - 1) * plotting_options.window_duration;
    curr_window = [0 plotting_options.window_duration];
    curr_window = curr_window + curr_rasterWindowOffset;
    % RelSpikeStartTime: x-axis offset from the left edge of the plot
    % rasterWindowOffset: x-axis window start time
%     hold off
    
    if plotting_options.showOnlyAlwaysStableCells
        isAlwaysStable = (active_processing.spikes.stability_count == 3);
        numAlwaysStableCells = sum(isAlwaysStable, 'all');
        [~, ~, plotHandle] = plotSpikeRaster(active_processing.spikes.time(isAlwaysStable),'PlotType','vertline','rasterWindowOffset', curr_rasterWindowOffset,'XLimForCell', curr_window);
        ylabel('Stable Unit Index')
        title(sprintf('Spike Train for Always Stable Units (%d of %d total)', numAlwaysStableCells, length(active_processing.spikes.time)));
    else
        [~, ~, plotHandle] = plotSpikeRaster(active_processing.spikes.time,'PlotType','vertline','rasterWindowOffset', curr_rasterWindowOffset,'XLimForCell', curr_window);
        ylabel('Unit Index')
        title(sprintf('Spike Train for %d second window from [%d, %d]', plotting_options.window_duration, curr_window(1), curr_window(end)));
    end
    xlabel('Time [seconds]')
    
    
    set(gca,'XTick',[]);
    drawnow;
end
