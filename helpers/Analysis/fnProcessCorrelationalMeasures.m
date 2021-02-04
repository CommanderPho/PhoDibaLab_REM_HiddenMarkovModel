function [active_spike_data_matrix, autocorrelations, partial_autocorrelations, pairwise_xcorrelations] = fnProcessCorrelationalMeasures(spike_data, active_results_indicies, processing_config)
%CORR Summary of this function goes here
%   Detailed explanation goes here

%% Auto-Correlations:

% % [autocor,lags] = xcorr(spike_data{1});
% [acf,lags,bounds] = autocorr(spike_data{1},'NumMA',2);

autocorrelations = cellfun((@(spike_timestamps) xcorr(spike_timestamps)), ...
        spike_data, 'UniformOutput', false);
    
partial_autocorrelations = cellfun((@(spike_timestamps) parcorr(spike_timestamps)), ...
        spike_data, 'UniformOutput', false);    


%% Cross-Correlations:

%% Extract all timeseries in the appropriate matrix format:
% 1 x 126 cell should not be transposed.
% 126 x 1 cell should
if (size(spike_data, 1) > size(spike_data, 2))
   active_spike_data_matrix = cell2mat(spike_data'); % 35351x126 double
else
   active_spike_data_matrix = cell2mat(spike_data); % 35351x126 double
end

% rho: Pairwise linear correlation coefficient
% 
% [pairwise_correlations.rho, pairwise_correlations.pval] = corr(temp.smoothed_spike_data_matrix);

pairwise_xcorrelations.lag_offsets = (-processing_config.max_xcorr_lag):processing_config.max_xcorr_lag;
pairwise_xcorrelations.num_lag_steps = length(pairwise_xcorrelations.lag_offsets);

% [pairwise_xcorrelations.xcorr, pairwise_xcorrelations.lags] = xcorr(temp.smoothed_spike_data_matrix, processing_config.max_xcorr_lag, 'normalized'); 
% xcorr: 13x15876 double ((2 × max_xcorr_lag + 1) × N^2)
% lags: 1x13 double


% pre-allocate output

% parfor unpacking
smoothed_spike_data_matrix = active_spike_data_matrix;
unique_electrode_index_pairs = active_results_indicies.unique_electrode_pairs;
max_xcorr_lag = processing_config.max_xcorr_lag;

output_pairwise_xcorrelations = zeros([active_results_indicies.num_unique_pairs pairwise_xcorrelations.num_lag_steps]);
parfor i = 1:active_results_indicies.num_unique_pairs
%    temp.curr_pair = active_result.indicies.unique_electrode_pairs(i,:);
   output_pairwise_xcorrelations(i,:) = xcorr(smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,1)), ...
       smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,2)), ...
       max_xcorr_lag,'normalized');
end

pairwise_xcorrelations.xcorr = output_pairwise_xcorrelations;

% parfor cleanup
clear smoothed_spike_data_matrix unique_electrode_index_pairs max_xcorr_lag output_pairwise_xcorrelations

% Find maximum correlations and when they occur
[pairwise_xcorrelations.processed.MaxXCorr, pairwise_xcorrelations.processed.MaxXCorrLagIndex] = max(pairwise_xcorrelations.xcorr');
[pairwise_xcorrelations.processed.MinXCorr, pairwise_xcorrelations.processed.MinXCorrLagIndex] = min(pairwise_xcorrelations.xcorr');
pairwise_xcorrelations.processed.meanXCorr = mean(pairwise_xcorrelations.xcorr, 2)';

% pairwise_xcorrelations.processed.MaxXCorr: contains the maximum xcorr value obtained for each unit
% pairwise_xcorrelations.processed.MaxXCorrLagIndex: contains when the maximum xcorr value occurred for each unit

[temp.sorted_values, temp.sorted_index] = sortrows(pairwise_xcorrelations.processed.MaxXCorr','descend');

% Get the sorted indicies:
pairwise_xcorrelations.processed.sorted.unique_electrode_index_pairs = active_results_indicies.unique_electrode_pairs(temp.sorted_index, :);
% Get the sorted xcorr values:
pairwise_xcorrelations.processed.sorted.MaxXCorr = pairwise_xcorrelations.processed.MaxXCorr(temp.sorted_index);
pairwise_xcorrelations.processed.sorted.MaxXCorrLagIndex = pairwise_xcorrelations.processed.MaxXCorrLagIndex(temp.sorted_index);
pairwise_xcorrelations.processed.sorted.xcorr = pairwise_xcorrelations.xcorr(temp.sorted_index, :);


end

