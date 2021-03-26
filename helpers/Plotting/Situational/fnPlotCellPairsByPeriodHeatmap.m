function [xcorr_fig] = fnPlotCellPairsByPeriodHeatmap(active_processing, general_results, results, filter_config, plottingOptions, reference_unit_index)
% fnPlotCellPairsByPeriodHeatmap - plots a heatmap for a given unit with cell pair on its x-axis and period on the y-axis
% Detailed explanation goes here
% 
% Syntax:  
%     [xcorr_fig] = fnPlotCellPairsByPeriodHeatmap(active_processing, general_results, results, filter_config, plottingOptions, reference_unit_index)
% 
% Inputs:
%    reference_unit_index: a valid unit index in the set of pairs
%    input2 - Description
%    input3 - Description
% 
% Outputs:
%    output1 - Description
%    output2 - Description
% 

% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 26-Mar-2021 ; Last revision: 26-Mar-2021 

% ------------- BEGIN CODE --------------

%     plotResults.figures.xcorr = figure(expt_info.index);
   
    temp.curr_ref_unit_index_string = num2str(reference_unit_index);
    
    temp.curr_active_unit_index = reference_unit_index;
%     temp.curr_active_unit_index = filter_config.filter_active_pair_values(reference_unit_index, 1);
%     temp.curr_active_pair_indicies = (temp.curr_active_unit_index == filter_config.filter_active_pair_values(:,1));
%     temp.curr_active_pair_values = filter_config.filter_active_pair_values(temp.curr_active_pair_indicies,:);
%     temp.curr_active_pair_labels = num2str(temp.curr_active_pair_values(:,2));
    
    %% New
    [temp.is_pair_included, temp.original_pair_index] = fnFilterPairsWithCriteria(general_results, temp.curr_active_unit_index);
    % original_pair_index: 126x1 double
    filtered.is_pair_included = temp.is_pair_included(filter_config.filter_active_pairs, 1); % 3655x1 logical
    filtered.original_pair_index = temp.original_pair_index(filter_config.filter_active_units, 1); % 86x1 double
    
    temp.curr_active_pair_indicies = find(filtered.is_pair_included);
    temp.curr_active_pair_values = filter_config.filter_active_pair_values(temp.curr_active_pair_indicies, :);
    temp.curr_active_pair_other_unit_index = ones([size(temp.curr_active_pair_values, 1), 1]);
    temp.curr_active_pair_other_unit_values = zeros([size(temp.curr_active_pair_values, 1), 1]);
    
    % Find the value in the pair that DOESN'T correspond to the active index (to get the other relevant unit)
    for i = 1:size(temp.curr_active_pair_values, 1)
        if (temp.curr_active_pair_values(i,1) == temp.curr_active_unit_index)
            temp.curr_active_pair_other_unit_index(i) = 2; % Use the other index for the value;
            temp.curr_active_pair_other_unit_values(i) = temp.curr_active_pair_values(i,2); 
        else
            temp.curr_active_pair_other_unit_index(i) = 1; % Use this index for the value;
            temp.curr_active_pair_other_unit_values(i) = temp.curr_active_pair_values(i,1);
        end
    end
    
    temp.curr_active_pair_labels = num2str(temp.curr_active_pair_other_unit_values);
    
    [xcorr_fig, heatmap_handle] = fnPlotAcrossREMXcorrHeatmap(results.all.per_period.xcorr_all_lags(:, temp.curr_active_pair_indicies)); % 668x3655 double
    heatmap_axis = gca;    
    xlabel('Cell Pairs')
    
%     xticklabels('manual');
    xticks(1:length(temp.curr_active_pair_labels));
    xticklabels(temp.curr_active_pair_labels);
    xticklabels('manual')
    xtickangle(90)
    [state_ax, epoch_ax] = fnPlotHelper_AddStateMapSubplot(active_processing, heatmap_axis);
    temp.curr_ref_unit_index_string = sprintf('unit[%d]', reference_unit_index);
    sgtitle([plottingOptions.curr_expt_string ' : XCorr for ' temp.curr_ref_unit_index_string ' and all units - All Periods - Good Units'])
    
end


% ------------- END OF CODE --------------
