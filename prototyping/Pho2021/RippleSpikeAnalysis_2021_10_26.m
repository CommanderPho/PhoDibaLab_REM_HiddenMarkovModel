
% Get only the ripple spikes
active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
active_spike_ripple_indices = cellfun(@(spike_ripple_indicies, is_spike_ripple) spike_ripple_indicies(is_spike_ripple), active_processing.spikes.RippleIndex(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes

% [cnt_unique, unique_a] = hist(a,unique(a));

[active_spike_ripple_unique_count, active_spike_ripple_unique_indicies] = cellfun(@(spike_ripple_indicies) hist(spike_ripple_indicies, unique(spike_ripple_indicies)), active_spike_ripple_indices, 'UniformOutput', false); %% Filtered to only show the ripple spikes

[active_spike_flatSpikeTimes, active_spike_flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(active_spike_times);


% Get uniques across all of them:

unique(active_spike_ripple_unique_indicies



%% accepts a cell array of length N, with each cell containing items with the same datatype but potentially different sizes.
% curr_cell_table = table(active_processing.spikes.time, ...
%     active_processing.spikes.isRippleSpike, ...
%     active_processing.spikes.RippleIndex, ...
% 'VariableNames', {'time','isRippleSpike','RippleIndex'});
% 
% % first_cell_array = active_processing.spikes.time;
% first_cell_array = curr_cell_table{:,1};
% num_cell_arrays_to_flatten = width(curr_cell_table);
% 
% num_cells = length(first_cell_array);
% cell_content_counts = cellfun(@length, first_cell_array);
% cell_indicies = 1:num_cells;
% 
% % The total number of output items for the flattened array
% flattened_total_item_count = sum(cell_content_counts, 'all');
% 
% % flattened_UnitIDs: the unit the original entry belongs to:
% flattened_UnitIDs = repelem(cell_indicies, cell_content_counts); % the second argument specifies how many times to repeat each item in the first array
% 
% % Now flatten each table variable:
% 
% % curr_cell_table{:,1}
% 
% % Closed-form vectorized flattening for non-equal sized cells
%    
% flatTableColumns = cell([num_cell_arrays_to_flatten 1]);
% % curr_flattened_table = table();
% curr_flattened_table = table('Size', [flattened_total_item_count, (1 + num_cell_arrays_to_flatten)], 'VariableNames', ['flattened_UnitIDs', curr_cell_table.Properties.VariableNames], 'VariableTypes', ['double', varfun(@class, curr_cell_table, 'OutputFormat', 'cell')]);
% % Add the flattened unit id's as the first column
% curr_flattened_table.flattened_UnitIDs = flattened_UnitIDs';
% 
% for var_idx = 1:num_cell_arrays_to_flatten
%     flatTableColumns{var_idx} = [curr_cell_table{:,var_idx}{:}];
%     curr_variable_name = curr_cell_table.Properties.VariableNames{var_idx};
%     curr_flattened_table.(curr_variable_name) = flatTableColumns{var_idx}';
% %     curr_flattened_table.(curr_variable_name) = [curr_cell_table{:,var_idx}{:}]';
% end
% 
% % Sort the output table by the time column:
% curr_flattened_table = sortrows(curr_flattened_table,'time','ascend');
% 
% % curr_flattened_table = cell2table(flatTableCells);
% % curr_flattened_table = cell2table(flatTableColumns, ...
% %     'VariableNames', curr_cell_table.Properties.VariableNames);


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


% fnReconstructCellsFromFlattenedContents

if ~isVariable(curr_flattened_table, 'flattened_UnitIDs')
	error('flattened_table lacks the column/Variable "flattened_UnitIDs", meaning it was not produced by fnReconstructCellsFromFlattenedContents. Aborting.');
end

[includedCellIDs, unitSpikeCells, unitFlatIndicies] = fnFlatSpikesToUnitCells(curr_flattened_table{:,2}, curr_flattened_table.flattened_UnitIDs, true);
num_reconstructed_cells = length(includedCellIDs);
num_reconstructed_columns = width(curr_flattened_table) - 1; % minus the flattened_unitID column


% reconstructedCellTable = table;
reconstructed_column_types = varfun(@class, curr_flattened_table, 'OutputFormat', 'cell');
reconstructed_column_types = reconstructed_column_types(2:end);
reconstructed_column_names = curr_flattened_table.Properties.VariableNames(2:end);

reconstructedCellTable = table('Size', [num_reconstructed_cells, num_reconstructed_columns], 'VariableNames', reconstructed_column_names, 'VariableTypes', reconstructed_column_types);

var_idx = 1;
curr_variable_name = reconstructed_column_names{var_idx};
reconstructedCellTable.(curr_variable_name) = unitSpikeCells;

for var_idx = 1:num_reconstructed_columns
    curr_variable_name = reconstructed_column_names{var_idx};
%     curr_selected_data = cellfun(@(column_data, included_indicies) column_data(included_indicies), curr_flattened_table{:, (var_idx + 1)}, unitFlatIndicies, 'UniformOutput', false); %% Filtered to only show the ripple spikes

    column_data = curr_flattened_table{:, (var_idx + 1)};
    curr_selected_data = cellfun(@(included_indicies) column_data(included_indicies), unitFlatIndicies, 'UniformOutput', false); %% Filtered to only show the ripple spikes
    reconstructedCellTable.(curr_variable_name) = curr_selected_data;
    
%     reconstructedCellTable.(curr_variable_name) = [curr_flattened_table{,var_idx}{:}]';
end





