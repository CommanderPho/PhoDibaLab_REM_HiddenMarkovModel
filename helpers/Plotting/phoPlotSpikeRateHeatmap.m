% phoPlotSpikeRateHeatmap.m
% Pho Hale, 02/08/2021
% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded


% Produces a heatmap ...
plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;

extantFigH = figure(9);
if plotting_options.showOnlyAlwaysStableCells
    isAlwaysStable = (active_processing.spikes.stability_count == 3);
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    plotHandle = heatmap(general_results.per_behavioral_state_period.spike_rate_per_unit);

    xlabel('Stable Unit Index')
    title(sprintf('Behavioral State Period vs. Unit Spike Rate (for Always Stable Units (%d of %d total))', numAlwaysStableCells, length(active_processing.spikes.time)));
else
    plotHandle = heatmap(general_results.per_behavioral_state_period.spike_rate_per_unit);
    xlabel('Unit Index')
    title('Behavioral State Period vs. Unit Spike Rate');
end
ylabel('Behavioral State Period Index');