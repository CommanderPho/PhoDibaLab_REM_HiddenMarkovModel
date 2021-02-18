    %% active_processing: pre-processed results
        %% behavioral_state_names: TODO
        %% earliest_start_timestamp: TODO
        %% behavioral_epochs: TODO
        %% behavioral_periods_table: TODO
        %% spikes: TODO
        %% processed_array: TODO 
		%% definitions: TODO 
			%% behavioral_state: {'rem','rem','quiet','active'}
			%% behavioral_epoch: {'pre_sleep','track','post_sleep'}

    %% active_results: items dependant on the bin size
        %% all: TODO
        %% by_epoch: TODO
        %% aggregates: TODO
        %% by_state: TODO 

    %% general_results: results that are independent of binning resolution
        %% indicies: TODO
        %% GroupedByState: TODO
        %% flattened_across_all_units: TODO
        %% edges: TODO
        %% per_behavioral_state_period: TODO
        %% counts: TODO 
        
        
        
       
        
        
%% TO REFACTOR:
% active_processing.processed_array.* -> active_results.all.*
% 
% temp.migrations.numElements = length(active_processing.processed_array);
% 
% for i = 1:temp.migrations.numElements
%     active_results = results_array{i};
%     
%     % Does not work because mergestruct isn't recurrsive:
% %     active_results = mergestruct(active_results, ...
% %                active_processing.processed_array{i});
%            
%            
%     temp.migrations.activeSourceFieldNames = fieldnames(active_processing.processed_array{i});    
%     temp.migrations.activeDestFieldNames = fieldnames(results_array{i});
%     
%     % Note that active_processing.processed_array{i} has the fields {'all', 'by_epoch', 'by_state'}
%         % Each of these contain {'spike_data', 'binned_spike_counts'}
%     
%     for fieldIndex = 1:length(temp.migrations.activeSourceFieldNames)
%         temp.migrations.currFieldName = temp.migrations.activeSourceFieldNames{fieldIndex};
%         
%         if ~isfield(temp.migrations.currFieldName, results_array{i})
%             % Field doesn't exist in destination structure, just add it
%             results_array{i}.(temp.migrations.currFieldName) = active_processing.processed_array{i}.(temp.migrations.currFieldName);
%         
%         else
%            % Field exists, must merge:
%            results_array{i}.(temp.migrations.currFieldName) = mergestruct(results_array{i}.(temp.migrations.currFieldName), ...
%                active_processing.processed_array{i}.(temp.migrations.currFieldName));
%         end
%         
%         % TODO: need to *merge* the fields:
%             % for example active_processing.processed_array.all exists both in 
%         
% %         active_processing.processed_array{i}.all
%     
%     end
%     
% %     results_array{i} = active_results;
% end