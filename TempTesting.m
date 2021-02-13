clc
% 
% % docGeneration.variableNames = {'active_processing.spikes'};
% % docGeneration.variableValues = cellfun((@(varName) evalin('base',varName)), docGeneration.variableNames);
% 
% 
% % docGeneration.variableMetadata.names = {'active_processing', 'active_results', 'general_results'};
% % docGeneration.variableMetadata.descriptions = {'pre-processed data','depend on the bin size','results that are independent of binning resolution'};
% 
% 
% 
% docGeneration.variableMetadata.names = {'active_results.all'};
% docGeneration.variableMetadata.descriptions = {'TODO','depend on the bin size','results that are independent of binning resolution'};
% 
% 
% 
% 
% 
% docGeneration.variableValues = cellfun((@(varName) evalin('base',varName)), docGeneration.variableMetadata.names, 'UniformOutput', false);
% 
% % variableMetadata.names = {'active_processing', 'active_processing.spikes'};
% % variableMetadata.descriptions = {'idk','idk'};
% 
% [copyableCodeCommentString, copyableCodeCommentBlocks] = fnGenerateDocumentationSkeleton(docGeneration.variableValues, docGeneration.variableMetadata);
% 
% fprintf('%s',copyableCodeCommentString);
% % disp(copyableCodeCommentString);
% 
%     %% active_processing: pre-processed results
%         %% behavioral_state_names: TODO
%         %% earliest_start_timestamp: TODO
%         %% behavioral_epochs: TODO
%         %% behavioral_periods_table: TODO
%         %% spikes: TODO
%         %% processed_array: Important, one for each bin size 
% active_processing.processed_array  
%     %% active_results: depend on the bin size
%         %% all: TODO
%             %% autocorrelations: TODO
%             %% partial_autocorrelations: TODO
%             %% pairwise_xcorrelations: TODO 
%         %% by_epoch: TODO
%         %% aggregates: TODO
%         %% by_state: TODO 
% 
%     %% general_results: results that are independent of binning resolution
%         %% indicies: TODO
%         %% GroupedByState: TODO
%         %% flattened_across_all_units: TODO
%         %% edges: TODO
%         %% per_behavioral_state_period: TODO
%         %% counts: TODO 
%         
%         
%   
%      % active_processing.processed_array{1}.by_state and .by_epoch don't seem to make much sense, do they?
%      
% % active_results.aggregates

% PhoDibaPrepare_Stage0
PhoDibaProcess_Stage1
PhoDibaProcess_Stage2
