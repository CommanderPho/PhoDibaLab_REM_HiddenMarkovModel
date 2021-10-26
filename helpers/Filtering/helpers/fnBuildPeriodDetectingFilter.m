function [outputFilterFunction] = fnBuildPeriodDetectingFilter(startTimes, endTimes)
%FNBUILDPERIODDETECTINGFILTER Should return a function that takes a single event time as input and returns the period index it falls within if it does fall within one and nan otherwise
%   Detailed explanation goes here

%% Testing: 
ripples_time_mat = source_data.ripple.RoyMaze1.time;
ripples_time_mat = (ripples_time_mat - source_data.behavior.RoyMaze1.time(1,1)) ./ 1e6; % Convert to relative timestamps since start
startTimes = ripples_time_mat(:,1);
endTimes = ripples_time_mat(:,2);

% Pre-allocate the cell arrays that will be used as teh table columns
is_ripple_spike = cell([height(active_processing.spikes) 1]);
spike_ripple_index = cell([height(active_processing.spikes) 1]);

for unit_idx = 1:height(active_processing.spikes)
%     unit_idx = 2;
    % functionCells{i}(active_processing.spikes.time{unit_idx})
    % find(functionCells{i}(active_processing.spikes.time{unit_idx})~=0, 1, 'first')
    % find(functionCells(active_processing.spikes.time{unit_idx})~=0, 1, 'first')
    unit_matches = cellfun(@(period_comparison_fcn) find(period_comparison_fcn(active_processing.spikes.time{unit_idx})~=0, 1, 'first'), functionCells, 'UniformOutput', false);
    unit_matched_periods = find(~cellfun(@isempty, unit_matches)); % These are the indicies of the periods that each spike falls within
    unit_matched_spike_indicies = [unit_matches{unit_matched_periods}]'; % These are the indicies of spikes for this unit that fall within a period
    
    % Finally we can go through and build the required row entry for each of the new columns for the spikes table:
    unit_is_ripple_spike = false(size(active_processing.spikes.time{unit_idx}));
    unit_is_ripple_spike(unit_matched_spike_indicies) = true;
    % Set the accumulator column:
    is_ripple_spike{unit_idx} = unit_is_ripple_spike;


    unit_spike_ripple_indicies = NaN(size(active_processing.spikes.time{unit_idx})); % Most spikes fall outside a period, and have a NaN value for the matched period index
    unit_spike_ripple_indicies(unit_matched_spike_indicies) = unit_matched_periods; % For the spikes that did match a period, set their value to the matched period index
    spike_ripple_index{unit_idx} = unit_spike_ripple_indicies;

end

% Finally, add the new columns to the table
active_processing.spikes.isRippleSpike = is_ripple_spike;
active_processing.spikes.RippleIndex = spike_ripple_index;



%% Function Start:
num_periods = length(startTimes);
functionCells = cell(num_periods,1);

for i = 1:num_periods
    currentEventFallsWithinPeriodFunction = @(t) ((startTimes(i) <= t) & (t <= endTimes(i)));
    functionCells{i} = currentEventFallsWithinPeriodFunction;
end

outputFilterFunction = inputArg1;

startTimes


ripples_time_mat


end

