active_processing.processed_array{2, 1}.all.binned_spike_firingRates = active_processing.processed_array{2, 1}.all.binned_spike_counts ./ processing_config.step_sizes(1);


temp.active_binsize_index = 1;
processing_config.step_sizes{temp.active_binsize_index}
active_processing.processed_array{temp.active_binsize_index, 1}.all.binned_spike_firingRates = cellfun((@(x) x ./ processing_config.step_sizes{temp.active_binsize_index}), active_processing.processed_array{temp.active_binsize_index, 1}.all.binned_spike_counts, 'UniformOutput', false);
% active_processing.processed_array{temp.active_binsize_index, 1}.all.binned_spike_firingRates = fnCellContentsTranpose(active_processing.processed_array{temp.active_binsize_index, 1}.all.binned_spike_firingRates);
%% Flatten over the rows (which are units) subset of spikes table for efficient ripple-related processing:
[curr_flattenedOverUnits_binned_spike_firingRates] = fnFlattenCellsToContents(active_processing.processed_array{temp.active_binsize_index, 1}.all.binned_spike_firingRates); % Flatten the cells:


binning_resolution = 0.003; % 3 ms bins
% [data.spikemat] = bz_SpktToSpkmat(spikeTimes, 'dt', binning_resolution, 'win', [active_processing.behavioral_epochs.start_seconds(1), active_processing.behavioral_epochs.end_seconds(3)], 'units', 'rate'); % Can use 'win', [startTime endTime]
% [data.spikemat] = bz_SpktToSpkmat(spikeTimes', 'dt', binning_resolution, 'win', [active_processing.behavioral_epochs.start_seconds(1), active_processing.behavioral_epochs.end_seconds(3)], 'units', 'rate'); % Can use 'win', [startTime endTime]
% curr_flattenedOverUnits_table(:, {'time','flattened_UnitIDs'});

% Get just the [timestmap UID] pairs

[data.spikemat] = bz_SpktToSpkmat(curr_flattenedOverUnits_table{:, {'time','flattened_UnitIDs'}}, 'dt', binning_resolution, 'win', [active_processing.behavioral_epochs.start_seconds(1), active_processing.behavioral_epochs.end_seconds(3)], 'units', 'rate'); % Can use 'win', [startTime endTime]
data.spikemat.timestamps

data.spikemat.data % [ts x numUnits]

curr_flattenedOverUnits_binned_spike_firingRates = sparse(data.spikemat.data);
max(data.spikemat.data)
min(data.spikemat.data)
mean(data.spikemat.data)