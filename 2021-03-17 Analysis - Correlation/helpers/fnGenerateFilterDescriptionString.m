function [filter_description_string] = fnGenerateFilterDescriptionString(filter_config)
%FNGENERATEFILTERDESCRIPTIONSTRING Summary of this function goes here
%   Detailed explanation goes here
    if isempty(filter_config.filter_included_cell_types)
        filter_description_string.cells = 'All Cells';
    else
        num_cell_types = length(filter_config.filter_included_cell_types);
        filter_description_string.cells = join(filter_config.filter_included_cell_types,'|');
        filter_description_string.cells = filter_description_string.cells{1};
%         filter_description_string.cells = [filter_description_string.cells ' Cells'];
        if num_cell_types == 1
            % Append ' Only' to the end of the string if there's only one cell type (like 'Pyramidal Only'
            filter_description_string.cells = [filter_description_string.cells ' Only'];
        end
    end
    
    filter_description_string.cells = strcat(filter_description_string.cells);

end

