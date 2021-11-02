% Pho Test TimeBinning:
%% phoTimeBinning: Aims to improve the computational efficiency of timeBinning.m
%% TODO: Does not work
error('Not yet implemented')

pbEvents = primaryPBEs;
spike = spikeStruct;

eventBeg = pbEvents(:,1);
eventEnd = pbEvents(:,2);


flattened_Event_breaks = reshape(pbEvents(:,1:2), [1, 2 * noEvents])';
[flattened_Event_breaks, flattened_Event_breaks_sort_idx] = sort(flattened_Event_breaks, 1, "ascend");


% Get the indicies that span the sessions from the session delination indicies.
	[SpanningRanges, SpanningIndexMatrix] = PointsToRanges(splitIndicies' ,1, length(allData.consolidatedDataTable.interSampleDuration)); %

noEvents = size(pbEvents, 1);
binSizes = [0.001*fileinfo.Fs binDur*fileinfo.Fs]; % binning for two time resolutions (1 ms and 20 ms durations)
setofUnits = cell(1, noEvents);
totalNumofUnits = max(spike.unit);
Binnedfiring = cell(noEvents, length(binSizes));
noFiringUnits = zeros(noEvents, 1);
no_flattened_spikes = length(spike.t);

static_isSpikeMember = ismember(spike.qclu , qclus); % 2875634 x 1

TimeIntervalCombined(


% spike_included_mask = zeros([no_flattened_spikes, noEvents]);
spike_included_mask = ones([no_flattened_spikes, noEvents]);
spike_included_mask(~static_isSpikeMember, :) = 0;



spikeInd = find(static_isSpikeMember & (spike.t >= eventBeg(evt)) & (spike.t < eventEnd(evt)));

for i = 1:no_flattened_spikes
    

end

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
    
