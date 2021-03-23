function [B, integration_info] = fnTemporalBias_SkaggsMcNaughton(lag_offsets, xcorr_input, integration_range)
%FNTEMPORALBIAS_SKAGGSMCNAUGHTON Temporal Bias "B" as defined in SkaggsMcNaughtonScience1996.pdf
% 
% measures the difference between the number of events in which a spike from cell i was followed within 200 ms by a spike from cell j and the
% number of events in which a spike from cell j was followed within 200 ms by a spike from cell i, possibly with other spikes of either cell in
% between.


% lag_offsets: the mapping of the lag indicies to real values.
% xcorr_input: a [num_units x num_units x lag_times] a matrix of cross-correlations X_ij(t).
% integration_range: a 2x1 array where the first element is the time before zero to integrate over, and the second is the time after 0.


%% Given a matrix of cross-correlations X_ij(t) provided as xcorr_input
if ~exist('integration_range','var')
%    integration_range = [-200, 200]; % specified in [ms] because I'm using active_results.all.pairwise_xcorrelations.lag_offsets_ms
   integration_range = [-0.2, 0.2];
end




% Xcorr lag offsets are given by -9
% active_results.all.pairwise_xcorrelations.lag_offsets
% Find the range 200ms prior and after the 0 point.
integration_info.integration_range.lowerIndicies = find((0 >= lag_offsets) & (lag_offsets >= integration_range(1)));
integration_info.integration_range.upperIndicies = find((0 <= lag_offsets) & (lag_offsets <= integration_range(2)));
integration_info.lag_zero_offset = find(lag_offsets == 0);


% Get before period:
integration_info.before_period = squeeze(sum(xcorr_input(:,:,integration_info.integration_range.lowerIndicies), 3)); % Sum over all indicies (integrate)

% Get after period:
integration_info.after_period = squeeze(sum(xcorr_input(:,:,integration_info.integration_range.upperIndicies), 3));

% Returns a scalar quantity for each pair of units [numUnits x numUnits]
B = integration_info.after_period - integration_info.before_period;

% One for each behavioral period:
%B = sum(temp.by_behavioral_period.curr_xcorr_forPair(:, temp.curr_xcorr_integration_range.upperIndicies),2) - sum(temp.by_behavioral_period.curr_xcorr_forPair(:, temp.curr_xcorr_integration_range.lowerIndicies), 2);

% For all: a scalar quantity:
%B = sum(temp.all.curr_xcorr_forPair(temp.curr_xcorr_integration_range.upperIndicies)) - sum(temp.all.curr_xcorr_forPair(temp.curr_xcorr_integration_range.lowerIndicies));




end

