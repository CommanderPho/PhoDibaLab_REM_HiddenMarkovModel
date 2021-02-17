function [is_unit_included] = fnFilterUnitsWithCriteria(active_processing, showOnlyAlwaysStableCells, included_cell_types, maximum_included_contamination_level)
%fnFilterUnitsWithCriteria Summary of this function goes here
% showOnlyAlwaysStableCells: if true, returns only the cells that are stable across all three behavioral epochs.
% included_cell_types: a cell array of cell types to include.
    % e.g. {'pyramidal', 'interneurons'}
% maximum_included_contamination_level: a cell array of the maximum permitted contamination level of a unit to include.
    % e.g. {2, 2}
    
    
    num_units = height(active_processing.spikes);
    is_unit_included = logical(ones([num_units 1]));

    if exist('showOnlyAlwaysStableCells','var') & showOnlyAlwaysStableCells
        isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
        is_unit_included = is_unit_included & isAlwaysStable;
    end


    if exist('included_cell_types','var') & ~isempty(included_cell_types)
        curr_conditions.is_unit_type = (@(compare_type) (active_processing.spikes.speculated_unit_type == compare_type));
        temp.curr_is_included = cellfun(curr_conditions.is_unit_type, included_cell_types, 'UniformOutput', false);
        for i = 1:length(temp.curr_is_included)
            is_unit_included = is_unit_included & temp.curr_is_included{i};
        end
    end

    if exist('maximum_included_contamination_level','var') & ~isempty(maximum_included_contamination_level)
        curr_conditions.meets_unit_contamination_requirements = (@(compare_value) (active_processing.spikes.maximum_included_contamination_level >= compare_value));
        temp.curr_is_included = cellfun(curr_conditions.meets_unit_contamination_requirements, maximum_included_contamination_level, 'UniformOutput', false);
        for i = 1:length(temp.curr_is_included)
            is_unit_included = is_unit_included & temp.curr_is_included{i};
        end
    end


end

