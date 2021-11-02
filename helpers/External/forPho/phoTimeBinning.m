function [Binnedfiring, noFiringUnits] = phoTimeBinning(pbEvents, spike, qclus, binDur, fileinfo)
%% phoTimeBinning: Aims to improve the computational efficiency of timeBinning.m

% This function is intended to bin the spikes trains of pyramidal units
% (idnetifies by qclus) within each event. 

% INPUTS
% binOverlapRatio: the overlap between the successive bins

eventBeg = pbEvents(:,1);
eventEnd = pbEvents(:,2);

noEvents = size(pbEvents, 1);
binSizes = [0.001*fileinfo.Fs binDur*fileinfo.Fs]; % binning for two time resolutions (1 ms and 20 ms durations)
setofUnits = cell(1, noEvents);
totalNumofUnits = max(spike.unit);
Binnedfiring = cell(noEvents, length(binSizes));
noFiringUnits = zeros(noEvents, 1);

static_isSpikeMember = ismember(spike.qclu , qclus); % 2875634 x 1


% [spikesTable] = fnAddRippleIdentifiersToSpikesTable(spikesTable, rippleStartTimes, rippleEndTimes);
% [outputFilterFunction] = fnBuildPeriodDetectingFilter(eventBeg, eventEnd);
[outputFilterFunction, functionCells] = fnBuildPeriodDetectingFilter(eventBeg, eventEnd);
event_matches = outputFilterFunction(spike.t(static_isSpikeMember));
event_matched_periods = find(~cellfun(@isempty, event_matches)); % These are the indicies of the periods that each spike falls within
event_matched_spike_indicies = [event_matches{event_matched_periods}]'; % These are the indicies of spikes for this unit that fall within a period


% 
% % Pre-allocate the cell arrays that will be used as teh table columns
% is_ripple_spike = cell([height(spikesTable) 1]);
% spike_ripple_index = cell([height(spikesTable) 1]);
% 
% unit_matches = cellfun(@(period_comparison_fcn) find(period_comparison_fcn(spikesTable.time{unit_idx})~=0, 1, 'first'), functionCells, 'UniformOutput', false);
% event_matched_periods = find(~cellfun(@isempty, unit_matches)); % These are the indicies of the periods that each spike falls within
% event_matched_spike_indicies = [unit_matches{event_matched_periods}]'; % These are the indicies of spikes for this unit that fall within a period

spikeTimes = spike.t(static_isSpikeMember);
spikeUnit  = spike.unit(static_isSpikeMember);

spikeTimes = spikeTimes(event_matched_spike_indicies);
spikeUnit  = spikeUnit(event_matched_spike_indicies);

isbetween(t,tlower,tupper)

isbetween(t,tlower,tupper)

% filtered_spike_times = spike.t(static_isSpikeMember);
% filtered_spike_times = seconds((spike.t(static_isSpikeMember) ./ 1e6)); 
% further_filtered_spike_times = filtered_spike_times((filtered_spike_times >= eventBeg) && (filtered_spike_times < eventEnd));

% [outputFilterFunction] = fnBuildPeriodDetectingFilter(eventBeg, eventEnd);
% [outputFilterFunction] = fnBuildPeriodDetectingFilter(rippleStartTimes, rippleEndTimes);
% unit_matches = cellfun(@(period_comparison_fcn) find(period_comparison_fcn(spikesTable.time{unit_idx})~=0, 1, 'first'), functionCells, 'UniformOutput', false);
% unit_matched_periods = find(~cellfun(@isempty, unit_matches)); % These are the indicies of the periods that each spike falls within
% unit_matched_spike_indicies = [unit_matches{unit_matched_periods}]'; % These are the indicies of spikes for this unit that fall within a period

curr_timetable = timetable(seconds((spike.t(static_isSpikeMember) ./ 1e6)), spike.t(static_isSpikeMember), spike.unit(static_isSpikeMember), ...
			'VariableNames',{'t', 'unit'});
% eventTimeRanges = timerange(seconds(eventBeg ./ 1e6), seconds(eventEnd ./ 1e6), 'closed');


% noEvents = 89118
for evt = 1:noEvents
    
    %% This line is very inefficient, using 91.7% of the computation time at 652seconds and with 89118 calls
%     spikeInd   = find(spike.t >= eventBeg(evt) & spike.t < eventEnd(evt) & ismember(spike.qclu , qclus));
%     spikeInd   = find(spike.t >= eventBeg(evt) & spike.t < eventEnd(evt) & static_isSpikeMember);
%     spikeTimes = spike.t(spikeInd);
%     spikeUnit  = spike.unit(spikeInd);
%     eventTimeRanges = timerange(seconds(eventBeg ./ 1e6), seconds(eventEnd ./ 1e6), 'closed');
    curr_filtered_TT = curr_timetable(timerange(seconds(eventBeg(evt) ./ 1e6), seconds(eventEnd(evt) ./ 1e6), 'closed'), :);
%     spikeTimes = seconds(curr_filtered_TT.t) .* 1e6; % convert back to double seconds
    spikeTimes = curr_filtered_TT.t;
    spikeUnit  = curr_filtered_TT.unit;
    
    setofUnits{evt}    = unique(spikeUnit);
    noFiringUnits(evt) = numel(setofUnits{evt});

    for iBin = 1:2
        binSize = binSizes(iBin);
        binEdges = eventBeg(evt): binSize : eventEnd(evt);
        noBins   = length(binEdges) - 1;
        
        Binnedfiring{evt, iBin} = zeros(totalNumofUnits, noBins);
        for ii = 1:length(setofUnits{evt})
            unit           = setofUnits{evt}(ii);
            unitSpikeTimes = spikeTimes(spikeUnit == unit);
            temp           = histc(unitSpikeTimes, binEdges);
            temp(end)      = [];
            if ~isempty(temp) %% why should it be empty?!
                Binnedfiring{evt, iBin}(unit, :) = temp; 
            else
                Binnedfiring{evt, iBin}(unit, :) = zeros(1, noBins); 
            end
            
        end
    end
end
    

end