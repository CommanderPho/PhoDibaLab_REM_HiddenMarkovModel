% phoPlotSpikeRateHeatmap.m
% Pho Hale, 02/08/2021
% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded


% Produces a heatmap ...
plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;



% phoPlotSpikeRateHeatmap_temp.plot_config.y_label = 'Behavioral State Period Index';
% phoPlotSpikeRateHeatmap_temp.activeMatrix = general_results.per_behavioral_state_period.spike_rate_per_unit;


phoPlotSpikeRateHeatmap_temp.plot_config.y_label = sprintf('Bin Index (bin size: %d sec)', active_binning_resolution);
phoPlotSpikeRateHeatmap_temp.activeMatrix = active_binned_spike_data_matrix;


extantFigH = figure(9);
clf;
if plotting_options.showOnlyAlwaysStableCells
    isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    
    phoPlotSpikeRateHeatmap_temp.activeMatrix = phoPlotSpikeRateHeatmap_temp.activeMatrix(:, isAlwaysStable);
       
%     plotHandle = heatmap(phoPlotSpikeRateHeatmap_temp.activeMatrix,'GridVisible', false);
%     plotHandle = heatmap(phoPlotSpikeRateHeatmap_temp.activeMatrix,'GridVisible', false);
    plotHandle = imagesc(phoPlotSpikeRateHeatmap_temp.activeMatrix);

    xlabel('Stable Unit Index')
    title(sprintf('%s vs. Unit Spike Rate (for Always Stable Units (%d of %d total))', phoPlotSpikeRateHeatmap_temp.plot_config.y_label, ...
        numAlwaysStableCells, length(active_processing.spikes.time)));
else
    plotHandle = heatmap(phoPlotSpikeRateHeatmap_temp.activeMatrix);
    xlabel('Unit Index')
    title(sprintf('%s vs. Unit Spike Rate', phoPlotSpikeRateHeatmap_temp.plot_config.y_label));
end
ylabel(phoPlotSpikeRateHeatmap_temp.plot_config.y_label);

% xticks([])

clear phoPlotSpikeRateHeatmap_temp;
