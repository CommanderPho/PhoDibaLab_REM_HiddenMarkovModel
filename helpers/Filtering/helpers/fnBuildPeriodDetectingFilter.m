function [outputFilterFunction] = fnBuildPeriodDetectingFilter(startTimes, endTimes)
%FNBUILDPERIODDETECTINGFILTER Should return a function that takes a single event time as input and returns the period index it falls within if it does fall within one and nan otherwise
%   Detailed explanation goes here

%% Testing: 
ripples_time_mat = source_data.ripple.RoyMaze1.time;
ripples_time_mat = (ripples_time_mat - source_data.behavior.RoyMaze1.time(1,1)) ./ 1e6; % Convert to relative timestamps since start
startTimes = ripples_time_mat(:,1);
endTimes = ripples_time_mat(:,2);

unit_idx = 2;
i = 50;
% functionCells{i}(active_processing.spikes.time{unit_idx})
% find(functionCells{i}(active_processing.spikes.time{unit_idx})~=0, 1, 'first')
% find(functionCells(active_processing.spikes.time{unit_idx})~=0, 1, 'first')
unit_matches = cellfun(@(period_comparison_fcn) find(period_comparison_fcn(active_processing.spikes.time{unit_idx})~=0, 1, 'first'), functionCells, 'UniformOutput', false);



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

