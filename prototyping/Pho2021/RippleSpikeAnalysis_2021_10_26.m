
% Get only the ripple spikes
active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
active_spike_ripple_indices = cellfun(@(spike_ripple_indicies, is_spike_ripple) spike_ripple_indicies(is_spike_ripple), active_processing.spikes.RippleIndex(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes

% [cnt_unique, unique_a] = hist(a,unique(a));

[active_spike_ripple_unique_count, active_spike_ripple_unique_indicies] = cellfun(@(spike_ripple_indicies) hist(spike_ripple_indicies, unique(spike_ripple_indicies)), active_spike_ripple_indices, 'UniformOutput', false); %% Filtered to only show the ripple spikes

[active_spike_flatSpikeTimes, active_spike_flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(active_spike_times);

%% Filtering Options:
% filter_config.filter_included_cell_types = {};
filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

%% Get filter info for active units
[plot_outputs.filter_active_units, plot_outputs.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, true, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);

temp.num_active_units = sum(plot_outputs.filter_active_units, 'all');
fprintf('Filter: Including %d of %d total units\n', temp.num_active_units, length(plot_outputs.filter_active_units));


%% Flatten subset of spikes table for efficient ripple-related processing:
curr_cell_table = table(active_processing.spikes.time, ...
        active_processing.spikes.isRippleSpike, ...
        active_processing.spikes.RippleIndex, ...
        'VariableNames', {'time','isRippleSpike','RippleIndex'});
% Add the optional variables
curr_cell_table.behavioral_duration_indicies = fnCellContentsTranpose(active_processing.spikes.behavioral_duration_indicies);
curr_cell_table.behavioral_states = fnCellContentsTranpose(active_processing.spikes.behavioral_states);
curr_cell_table.behavioral_epoch = fnCellContentsTranpose(active_processing.spikes.behavioral_epoch);

[curr_flattened_table] = fnFlattenCellsToContents(curr_cell_table);

% curr_flattened_table(curr_flattened_table.isRippleSpike)

%% Removes rows with missing values, meaning the rows with RippleIndex == NaN (meaning they aren't ripple spikes)
curr_filtered_table = rmmissing(curr_flattened_table);

%% Exclude sleep states:
% curr_filtered_table = curr_filtered_table(('rem' == curr_filtered_table.behavioral_states), :); % Only those occuring during rem
curr_filtered_table = curr_filtered_table(('nrem' == curr_filtered_table.behavioral_states), :);
% curr_filtered_table = groupfilter(curr_filtered_table, 'groupID', @(x) all(x == 'rem'), 'behavioral_states');

curr_filtered_table = curr_filtered_table((('rem' == curr_filtered_table.behavioral_states) | ('nrem' == curr_filtered_table.behavioral_states)), :);


%% Filter again to a specific ripple index
unique_ripple_indices = unique(curr_filtered_table.RippleIndex);
eachRipple_filtered_flattened_table = cell([length(unique_ripple_indices), 1]);
parfor ripple_idx = 1:length(unique_ripple_indices)
    eachRipple_filtered_flattened_table{ripple_idx} = curr_filtered_table(curr_filtered_table.RippleIndex == unique_ripple_indices(ripple_idx), :);
end


%% Visualize 

eachRipple_activeSet_Matrix = sparse(zeros([height(curr_cell_table), length(unique_ripple_indices)]));
for ripple_idx = 1:length(unique_ripple_indices)
    eachRipple_activeSet_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = 1;
end

fnPhoMatrixPlot(eachRipple_activeSet_Matrix(plot_outputs.filter_active_units, :))

%% Clustering exploration
X = eachRipple_activeSet_Matrix(plot_outputs.filter_active_units, :);
[idx,V,D] = spectralcluster(X, 3);
gscatter(X(:,1),X(:,2),idx);

% plottingOptions = struct();
% [h, temp.plot_info] = fnPhoMatrixPlotDetailed(eachRipple_activeSet_Matrix(plot_outputs.filter_active_units, :), plottingOptions);

%% fnReconstructCellsFromFlattenedContents
% Test reconstruction of original test table from the flattened table:
[reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(curr_filtered_table);






