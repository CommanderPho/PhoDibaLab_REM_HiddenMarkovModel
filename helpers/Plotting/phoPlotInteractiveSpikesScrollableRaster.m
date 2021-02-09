% phoPlotInteractiveSpikesScrollableRaster.m
% Pho Hale, 02/05/2021

% Produces a raster plot that displays the spike trains for each unit in a window of customizable length.
%   Spawns an interactive slider that allows you to specify the current window to look at, acting as a paginated manner.

plotting_options.window_duration = 10; % 10 seconds

plotting_options.active_timesteps = timesteps_array{1};
plotting_options.total_duration = plotting_options.active_timesteps(end) - plotting_options.active_timesteps(1);
plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;

plotting_options.num_windows = seconds(ceil(plotting_options.total_duration / plotting_options.window_duration));

% Build a new slider controller
iscInfo.slider_identifier = 'PhoScrollableSpikeRasterPlot';
iscInfo.curr_i = 1;
iscInfo.NumberOfSeries = plotting_options.num_windows;


extantFigH = figure(12);

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
        plotted_figH = figure(12);
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
    
    drawnow;
end
