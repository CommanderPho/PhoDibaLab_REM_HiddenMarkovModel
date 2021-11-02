
% %% From the older Pho Rotation approach:
% plotting_option.active_step_size = processing_config.step_sizes{temp.active_binsize_index};
% plotting_options.active_timesteps = temp.curr_timesteps_array{temp.active_binsize_index};

%% From running PhoDiba_BayesianDecoding2021.m:
plotting_option.active_step_size = tau;
plotting_options.active_timesteps = activeTimeBins;

plotting_options.total_duration = plotting_options.active_timesteps(end) - plotting_options.active_timesteps(1);

plotting_options.active_timesteps = plotting_options.active_timesteps(1:end-1);
% plotting_options.active_timesteps = plotting_options.active_timesteps + (plotting_option.active_step_size / 2); % This was to center the bins, but it isn't needed


extantFigH = figure(13);

%% Scrollplot mode:
curr_i = 1;
[extantFigH, currHeatmapPlotHandles, currStateMapHandle, plot_outputs] = pho_plot_firingRate(curr_flattenedOverUnits_binned_spike_firingRates, temp.curr_active_processing, filter_config, plotting_options, extantFigH, curr_i);
%currHeatmapPlotHandles.mainHeatmap.XData % comes out [1 353510]
% currHeatmapPlotHandles.mainHeatmap.XData = [0 plotting_options.active_timesteps(end)];
% xlim([0 plotting_options.active_timesteps(end)])
% currHeatmapPlotHandles.axesHandle.XLim
currPlotHandles = currHeatmapPlotHandles;

% xlim([0 35351.0593330000])


%% Build the scrollable interaction bar that sits below the main raster plot:
scrollHandles = scrollplot(currHeatmapPlotHandles.mainHeatmap, 'WindowSizeX', plotting_options.window_duration);
% scrollHandles = scrollplot(currStateMapHandle, 'WindowSizeX', plotting_options.window_duration);
% scrollHandles = scrollplot({currStateMapHandle currPlotHandles.ripplesHandle, currPlotHandles.linesHandle}, 'WindowSizeX', plotting_options.window_duration);





function [plotted_figH, heatmapPlotHandle, stateMapHandle, plot_outputs] = pho_plot_firingRate(target_data, active_processing, filter_config, plotting_options, extantFigH, curr_windowIndex)
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
    
    %% Filtering:
    target_data = target_data(plot_outputs.filter_active_units, :);
    
    %% Re-ordering/sorting
    % Re-order the units in the spike table
    target_data = target_data(plotting_options.sorting_config.unit_sort_indicies, :);
    plotting_options.unitBackgroundColors(:, plotting_options.sorting_config.unit_sort_indicies);
    
%     xx = [plotting_options.active_timesteps(1) plotting_options.active_timesteps(end)];
%     yy = [1 length(plotting_options.sorting_config.unit_sort_indicies)];
%     heatmapPlotHandle.mainHeatmap = imagesc(target_data);
%     heatmapPlotHandle.mainHeatmap = imagesc(xx, yy, target_data);
%     heatmapPlotHandle.mainHeatmap = imagesc('XData',xx,'YData', yy,'CData', target_data);
%     heatmapPlotHandle.mainHeatmap = imagesc(xx, yy, target_data);

%     heatmapPlotHandle.mainHeatmap = fnPhoMatrixPlot(target_data);
%     heatmapPlotHandle.mainHeatmap = stackedplot(plotted_figH, target_data);
    
    heatmapPlotHandle.mainHeatmap = imagesc([0 plotting_options.active_timesteps(end)], [1 size(target_data, 1)], target_data);

    heatmapPlotHandle.mainHeatmap = pcolor(X,Y,C); % returns a surface object
    axis ij
    axis square

    % colormap(mymap)


%     hold off;
%     for i = 1:size(target_data, 1)
%         heatmapPlotHandle.mainHeatmap(i) = plot(seconds(plotting_options.active_timesteps), target_data(i, :));
%         hold on;
%     end
%     colormap('jet')

    colormap('hot')

    if plotting_options.showOnlyAlwaysStableCells
        numAlwaysStableCells = sum(plot_outputs.filter_active_units, 'all');
        ylabel('Stable Unit Index')
        title(sprintf('Firing Rate Heatmap for Always Stable Units (%d of %d total)', numAlwaysStableCells, length(active_processing.spikes.time)));
    else
        ylabel('Unit Index')
        title(sprintf('Firing Rate Heatmap for %d second window from [%d, %d]', plotting_options.window_duration, curr_window(1), curr_window(end)));
    end
    xlabel('Time [seconds]')
    
    ax = gca;
    heatmapPlotHandle.axesHandle = ax;
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

    

end

