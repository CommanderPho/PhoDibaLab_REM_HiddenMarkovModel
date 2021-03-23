function [output] = fnNormalizeAllPairsXCorr(xcorr_full, xcorr_lag_zero_offset)
%FNNORMALIZEALLPAIRSXCORR Normalizes an xcorr_full input matrix in several ways and returns it in output (along with the data computed related to the normalization)
%   Pho Hale 2021-03-22

% xcorr_full: a [num_unique_pairs x num_lag_steps] array
% xcorr_lag_zero_offset: an optional scalar index that represents where the zero offset for num_lag_step is. If missing, zero_offset-based normalizations are skipped

%% Output:
% output.xcorr.raw: duplicate of input matrix
% output.xcorr.globally_normalized: input matrix normalized by the global maximum


% Compute max and normalization factors:
output.global_xcorr_pair_max = max(abs(xcorr_full), [], 2); % Max for each pair, should be [num_unique_pairs 1] array
output.global_xcorr_max = max(output.global_xcorr_pair_max, [], 'all'); % Absolute global max across all pairs

if exist('xcorr_lag_zero_offset','var')
    output.global_xcorr_zero_offset_values = xcorr_full(:, xcorr_lag_zero_offset); % The zero offset value for each pair, [num_unique_pairs 1] array
    output.global_xcorr_zero_offset_max = max(abs(output.global_xcorr_zero_offset_values), [], 'all');
    % %normalize so that zero-lag has a height of 1
    % 
    % TODO
end


output.xcorr.raw = xcorr_full;
% Normalize by global maximum
output.xcorr.globally_normalized = xcorr_full ./ output.global_xcorr_max; 
% % Normalize for each pair:
% output.by_behavioral_period.curr_xcorr_allPairs.normalized = output.by_behavioral_period.curr_xcorr_allPairs.raw ./ output.global_xcorr_max; % Normalize by global maximum
% output.plotVals = 
% output.plotVals = output.plotVals ./ output.plotVals(output.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
end

