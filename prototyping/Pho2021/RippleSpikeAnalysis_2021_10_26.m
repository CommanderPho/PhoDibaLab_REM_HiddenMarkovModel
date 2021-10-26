
% Get only the ripple spikes
active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
active_spike_ripple_indices = cellfun(@(spike_ripple_indicies, is_spike_ripple) spike_ripple_indicies(is_spike_ripple), active_processing.spikes.RippleIndex(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes

% [cnt_unique, unique_a] = hist(a,unique(a));

[active_spike_ripple_unique_count, active_spike_ripple_unique_indicies] = cellfun(@(spike_ripple_indicies) hist(spike_ripple_indicies, unique(spike_ripple_indicies)), active_spike_ripple_indices, 'UniformOutput', false); %% Filtered to only show the ripple spikes

[active_spike_flatSpikeTimes, active_spike_flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(active_spike_times);


%% Flatten subset of spikes table for efficient ripple-related processing:
curr_cell_table = table(active_processing.spikes.time, ...
        active_processing.spikes.isRippleSpike, ...
        active_processing.spikes.RippleIndex, ...
        'VariableNames', {'time','isRippleSpike','RippleIndex'});
[curr_flattened_table] = fnFlattenCellsToContents(curr_cell_table);

% curr_flattened_table(curr_flattened_table.isRippleSpike)

%% Removes rows with missing values, meaning the rows with RippleIndex == NaN (meaning they aren't ripple spikes)
curr_ripplesOnly_flattened_table = rmmissing(curr_flattened_table);

%% Filter again to a specific ripple index
unique_ripple_indices = unique(curr_ripplesOnly_flattened_table.RippleIndex);
eachRipple_filtered_flattened_table = cell([length(unique_ripple_indices), 1]);
parfor ripple_idx = 1:length(unique_ripple_indices)
    eachRipple_filtered_flattened_table{ripple_idx} = curr_ripplesOnly_flattened_table(curr_ripplesOnly_flattened_table.RippleIndex == unique_ripple_indices(ripple_idx), :);
end




%% fnReconstructCellsFromFlattenedContents
% Test reconstruction of original test table from the flattened table:
[reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(curr_ripplesOnly_flattened_table);






