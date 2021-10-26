function [curr_flattened_table] = fnFlattenCellsToContents(varargin)
%FNFLATTENCELLSTOCONTENTS Summary of this function goes here
%   Detailed explanation goes here
%% accepts a cell array of length N, with each cell containing items with the same datatype but potentially different sizes.

%% Example:
% curr_cell_table = table(active_processing.spikes.time, ...
%         active_processing.spikes.isRippleSpike, ...
%         active_processing.spikes.RippleIndex, ...
%         'VariableNames', {'time','isRippleSpike','RippleIndex'});

numFnInputs = nargin;

if numFnInputs == 0
    error('fnFlattenCellsToContents(...) called with no arguments!');
elseif ((numFnInputs == 1) & istabular(varargin{1}))
%     if istabular(varargin{1})
        % A table object
        curr_cell_table = varargin{1};
        first_cell_array = curr_cell_table{:,1};
        num_cell_arrays_to_flatten = width(curr_cell_table);

%     elseif iscell(varargin{1})
% 		framesList = varargin{1};
% 	else
% 		error('fnFlattenCellsToContents(...) called with only one argument and it is not a list of frames! Why would you take the average of only one item?');
% 	end
else
    % Otherwise more than one input hopefully
    curr_cell_table = table(varargin{:});
    first_cell_array = curr_cell_table{:,1};
    num_cell_arrays_to_flatten = width(curr_cell_table);

%     for arg_idx = 1:numFnInputs
%         if iscell(varargin{arg_idx})
%             first_cell_array = curr_cell_table{:,1};
%             num_cell_arrays_to_flatten = width(curr_cell_table);
% 	    else
% 		    error('FnComputeAverageFrame(...) called with only one argument and it is not a list of frames! Why would you take the average of only one item?');
%         end
%     end

end

num_cells = length(first_cell_array);
cell_content_counts = cellfun(@length, first_cell_array);
cell_indicies = 1:num_cells;

% The total number of output items for the flattened array
flattened_total_item_count = sum(cell_content_counts, 'all');

% flattened_UnitIDs: the unit the original entry belongs to:
flattened_UnitIDs = repelem(cell_indicies, cell_content_counts); % the second argument specifies how many times to repeat each item in the first array

% Now flatten each table variable:


% Closed-form vectorized flattening for non-equal sized cells
flatTableColumns = cell([num_cell_arrays_to_flatten 1]);
% curr_flattened_table = table();
curr_flattened_table = table('Size', [flattened_total_item_count, (1 + num_cell_arrays_to_flatten)], 'VariableNames', ['flattened_UnitIDs', curr_cell_table.Properties.VariableNames], 'VariableTypes', ['double', varfun(@class, curr_cell_table, 'OutputFormat', 'cell')]);
% Add the flattened unit id's as the first column
curr_flattened_table.flattened_UnitIDs = flattened_UnitIDs';

for var_idx = 1:num_cell_arrays_to_flatten
    flatTableColumns{var_idx} = [curr_cell_table{:,var_idx}{:}];
    curr_variable_name = curr_cell_table.Properties.VariableNames{var_idx};
    curr_flattened_table.(curr_variable_name) = flatTableColumns{var_idx}';
%     curr_flattened_table.(curr_variable_name) = [curr_cell_table{:,var_idx}{:}]';
end

% Sort the output table by the time column:
curr_flattened_table = sortrows(curr_flattened_table,'time','ascend');

% curr_flattened_table = cell2table(flatTableCells);
% curr_flattened_table = cell2table(flatTableColumns, ...
%     'VariableNames', curr_cell_table.Properties.VariableNames);



end

