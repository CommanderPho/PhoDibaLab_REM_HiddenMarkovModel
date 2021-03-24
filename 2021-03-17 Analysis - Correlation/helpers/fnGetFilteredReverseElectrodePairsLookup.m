function [filtered_reverse_lookup_unique_electrode_pairs] = fnGetFilteredReverseElectrodePairsLookup(general_results, filtered_original_unit_indicies, filtered_other_dim_original_unit_indicies)
%fnGetFilteredReverseElectrodePairsLookup Trivial right now

% filtered_original_unit_indicies: the set of unit indicies to filter by. 
% filtered_other_dim_original_unit_indicies: an optional set of other original_unit_indicies that will be used to filter for the other second dimenion of the output. If not provided, filtered_original_unit_indicies will be used.

if ~exist('other_dim_original_unit_indicies','var')
    filtered_other_dim_original_unit_indicies = filtered_original_unit_indicies;
end

% Pre allocate new lookup info
filtered_reverse_lookup_unique_electrode_pairs = general_results.indicies.reverse_lookup_unique_electrode_pairs(filtered_original_unit_indicies, filtered_other_dim_original_unit_indicies);

% Output will be a [length(filtered_original_unit_indicies) x length(filtered_other_dim_original_unit_indicies)] matrix of electrode_pair_indicies that can be used as a lookup.


dim_1.num_valid_units = length(filtered_original_unit_indicies);
dim_2.num_valid_units = length(filtered_other_dim_original_unit_indicies);


%% Want a 7875 x 1 matrix (one for each pair) of indicies in the new plot
% Nieve strategy:
% for active_unit_A_index = 1:dim_1.num_valid_units
%     
%     for active_unit_B_index = 1:dim_2.num_valid_units
%         if plotting_options.new_xcorr_plot.plot_single_row & (active_unit_A_index > 1)
%             % do nothing, just continue
%         elseif ~ismember(active_unit_A_index, plotting_options.new_xcorr_plot.included_unit_A_indicies)
%             % do nothing also
%         else
%             active_pair_index = across_experiment_results{active_expt_index}.general_results.indicies.reverse_lookup_unique_electrode_pairs(active_unit_A_index, active_unit_B_index);
%             % Convert subplot index to incement along a column first (each row for the column) and then move to the next column instead of its default (row -> column) mode
%             
%             [curr_row, curr_col] = ind2subplot(xcorr_all_plots.num_subplot_rows, xcorr_all_plots.num_subplot_columns, temp.linear_subplot_accumulator);
%             
%             temp.linear_subplot_accumulator = temp.linear_subplot_accumulator + 1;
%         end
%         temp.linear_accumulator = temp.linear_accumulator + 1;
%     end
% end




% general_results.indicies.reverse_lookup_unique_electrode_pairs: 126 x 126
% Want to get the correct reverse_loookup_electrode_pair index for only the filtered values.


% so let's say you've filtered the cells down to the ones you want.
% You have a filter_config.original_unit_index
% [filter_config.filter_active_units, filter_config.original_unit_index] = fnFilterUnitsWithCriteria(across_experiment_results{active_expt_index}.active_processing, across_experiment_results{active_expt_index}.processing_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
%   filter_config.filter_maximum_included_contamination_level);
    
% filter_config.filter_active_units: [totalNumAllUnits x 1]
% filter_config.original_unit_index: [filteredNumValidUnits x 1]

% num_valid_units = length(filter_config.original_unit_index);

% Given the

%   Should return a numCell x numCell matrix 


% [is_pair_included, original_pair_index] = fnFilterPairsWithCriteria(general_results, included_unit_indicies)
% 
% [temp.is_pair_included, temp.original_pair_index] = fnFilterPairsWithCriteria(general_results, filter_config.original_unit_index);
    % original_pair_index: 126x86 double
end

