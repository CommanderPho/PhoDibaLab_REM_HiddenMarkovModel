function [outputFilterFunction] = fnBuildPeriodDetectingFilter(startTimes, endTimes)
%FNBUILDPERIODDETECTINGFILTER Should return a function that takes a single event time as input and returns the period index it falls within if it does fall within one and nan otherwise
%   Detailed explanation goes here

    %% Function Start:
    num_periods = length(startTimes);
    functionCells = cell(num_periods,1);
    
    for i = 1:num_periods
        currentEventFallsWithinPeriodFunction = @(t) ((startTimes(i) <= t) & (t <= endTimes(i)));
        functionCells{i} = currentEventFallsWithinPeriodFunction;
    end

end

