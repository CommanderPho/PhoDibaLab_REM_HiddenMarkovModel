function [active_results] = fnProcessCorrelationalMeasures(spike_data, active_results, process_config)
%CORR Summary of this function goes here
%   Detailed explanation goes here

%% Auto-Correlations:

% % [autocor,lags] = xcorr(spike_data{1});
% [acf,lags,bounds] = autocorr(spike_data{1},'NumMA',2);

active_results.autocorrelations = cellfun((@(spike_timestamps) xcorr(spike_timestamps)), ...
        spike_data, 'UniformOutput', false);
    
active_results.partial_autocorrelations = cellfun((@(spike_timestamps) parcorr(spike_timestamps)), ...
        spike_data, 'UniformOutput', false);    


%% Cross-Correlations:

%% Extract all timeseries in the appropriate matrix format:
temp.smoothed_spike_data_matrix = cell2mat(spike_data'); % 35351x126 double

% rho: Pairwise linear correlation coefficient
% 
% [active_results.pairwise_correlations.rho, active_results.pairwise_correlations.pval] = corr(temp.smoothed_spike_data_matrix);


active_results.pairwise_xcorrelations.lag_offsets = (-process_config.max_xcorr_lag):process_config.max_xcorr_lag;
active_results.pairwise_xcorrelations.num_lag_steps = length(active_results.pairwise_xcorrelations.lag_offsets);

% [active_results.pairwise_xcorrelations.xcorr, active_results.pairwise_xcorrelations.lags] = xcorr(temp.smoothed_spike_data_matrix, process_config.max_xcorr_lag, 'normalized'); 
% xcorr: 13x15876 double ((2 × max_xcorr_lag + 1) × N^2)
% lags: 1x13 double


% pre-allocate output

% parfor unpacking
smoothed_spike_data_matrix = temp.smoothed_spike_data_matrix;
unique_electrode_index_pairs = active_results.indicies.unique_electrode_pairs;
max_xcorr_lag = process_config.max_xcorr_lag;

pairwise_xcorrelations = zeros([active_results.indicies.num_unique_pairs active_results.pairwise_xcorrelations.num_lag_steps]);
parfor i = 1:active_results.indicies.num_unique_pairs
%    temp.curr_pair = active_result.indicies.unique_electrode_pairs(i,:);
   pairwise_xcorrelations(i,:) = xcorr(smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,1)), ...
       smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,2)), ...
       max_xcorr_lag,'normalized');
end

active_results.pairwise_xcorrelations.xcorr = pairwise_xcorrelations;

% parfor cleanup
clear smoothed_spike_data_matrix unique_electrode_index_pairs max_xcorr_lag pairwise_xcorrelations

% Find maximum correlations and when they occur
[active_results.pairwise_xcorrelations.processed.MaxXCorr, active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex] = max(active_results.pairwise_xcorrelations.xcorr');
[active_results.pairwise_xcorrelations.processed.MinXCorr, active_results.pairwise_xcorrelations.processed.MinXCorrLagIndex] = min(active_results.pairwise_xcorrelations.xcorr');
active_results.pairwise_xcorrelations.processed.meanXCorr = mean(active_results.pairwise_xcorrelations.xcorr, 2)';

% active_results.pairwise_xcorrelations.processed.MaxXCorr: contains the maximum xcorr value obtained for each unit
% active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex: contains when the maximum xcorr value occurred for each unit

[temp.sorted_values, temp.sorted_index] = sortrows(active_results.pairwise_xcorrelations.processed.MaxXCorr','descend');

% Get the sorted indicies:
active_results.pairwise_xcorrelations.processed.sorted.unique_electrode_index_pairs = active_results.indicies.unique_electrode_pairs(temp.sorted_index, :);
% Get the sorted xcorr values:
active_results.pairwise_xcorrelations.processed.sorted.MaxXCorr = active_results.pairwise_xcorrelations.processed.MaxXCorr(temp.sorted_index);
active_results.pairwise_xcorrelations.processed.sorted.MaxXCorrLagIndex = active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex(temp.sorted_index);
active_results.pairwise_xcorrelations.processed.sorted.xcorr = active_results.pairwise_xcorrelations.xcorr(temp.sorted_index, :);


end

