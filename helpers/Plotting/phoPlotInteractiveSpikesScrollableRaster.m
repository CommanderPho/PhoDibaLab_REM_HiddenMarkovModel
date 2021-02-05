% phoPlotInteractiveSpikesScrollableRaster.m
% Pho Hale, 02/05/2021

plotting_options.window_duration = 10; % 10 seconds

plotting_options.active_timesteps = timesteps_array{1};
plotting_options.total_duration = plotting_options.active_timesteps(end) - plotting_options.active_timesteps(1);


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



% 
% figure(3)
% 
% % RelSpikeStartTime: x-axis offset from the left edge of the plot
% % rasterWindowOffset: x-axis window start time
% 
% 
% plotSpikeRaster(active_processing.spikes.time,'PlotType','vertline','rasterWindowOffset', 0.01,'XLimForCell',[0 0.201]);
% xlabel('Time [seconds]')
% ylabel('Unit Index')
% title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
% set(gca,'XTick',[]);


% slider_controller = fnBuildCallbackInteractiveSliderController(iscInfo, @(extantFigH, curr_i) (pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i)) );
extantFigH = figure(3);

slider_controller = fnBuildCallbackInteractiveSliderController(iscInfo, @(curr_i) (pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_i)) );

%% Plot function called as a callback on update
function plotted_figH = pho_plot_spikeRaster(active_processing, plotting_options, extantFigH, curr_windowIndex)
    if exist('extantFigH','var')
        plotted_figH = figure(extantFigH); 
    else
        plotted_figH = figure(3);
    end
    curr_rasterWindowOffset = double(curr_windowIndex - 1) * plotting_options.window_duration;
    curr_window = [0 plotting_options.window_duration];
    % RelSpikeStartTime: x-axis offset from the left edge of the plot
    % rasterWindowOffset: x-axis window start time
    plotSpikeRaster(active_processing.spikes.time,'PlotType','vertline','rasterWindowOffset', curr_rasterWindowOffset,'XLimForCell', curr_window);
    xlabel('Time [seconds]')
    ylabel('Unit Index')
    title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
    set(gca,'XTick',[]);

end
