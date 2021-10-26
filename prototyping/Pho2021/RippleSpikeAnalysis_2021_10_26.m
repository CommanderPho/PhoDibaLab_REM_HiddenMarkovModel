
% Get only the ripple spikes
active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
active_spike_ripple_indices = cellfun(@(spike_ripple_indicies, is_spike_ripple) spike_ripple_indicies(is_spike_ripple), active_processing.spikes.RippleIndex(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes

% [cnt_unique, unique_a] = hist(a,unique(a));

[active_spike_ripple_unique_count, active_spike_ripple_unique_indicies] = cellfun(@(spike_ripple_indicies) hist(spike_ripple_indicies, unique(spike_ripple_indicies)), active_spike_ripple_indices, 'UniformOutput', false); %% Filtered to only show the ripple spikes

[active_spike_flatSpikeTimes, active_spike_flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(active_spike_times);


% Test single cell array input:
[curr_flattened_table] = fnFlattenCellsToContents(active_processing.spikes.time, active_processing.spikes.isRippleSpike, active_processing.spikes.RippleIndex);
% Test multiple cell array inputs:
[curr_flattened_table] = fnFlattenCellsToContents(active_processing.spikes.time, active_processing.spikes.isRippleSpike, active_processing.spikes.RippleIndex);

% Test table input version:
curr_cell_table = table(active_processing.spikes.time, ...
        active_processing.spikes.isRippleSpike, ...
        active_processing.spikes.RippleIndex, ...
        'VariableNames', {'time','isRippleSpike','RippleIndex'});
[curr_flattened_table] = fnFlattenCellsToContents(curr_cell_table);


%% fnReconstructCellsFromFlattenedContents
% Test reconstruction of original test table from the flattened table:
[reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(curr_flattened_table);






