function [Binnedfiring, noFiringUnits] = timeBinning(pbEvents, spike, qclus, binDur, fileinfo)

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

static_isSpikeMember = ismember(spike.qclu , qclus);

for evt = 1:noEvents
    
%     spikeInd   = find(spike.t >= eventBeg(evt) & spike.t < eventEnd(evt) & ismember(spike.qclu , qclus));
    spikeInd   = find(spike.t >= eventBeg(evt) & spike.t < eventEnd(evt) & static_isSpikeMember);
    
    
    
    spikeTimes = spike.t(spikeInd);
    spikeUnit  = spike.unit(spikeInd);
    
    
    setofUnits{evt}    = unique(spikeUnit);
    noFiringUnits(evt) = numel(setofUnits{evt});

    
    for iBin = 1:2

        binSize = binSizes(iBin);
        
        binEdges = eventBeg(evt): binSize : eventEnd(evt);
        noBins   = length(binEdges) - 1;
        

        Binnedfiring{evt, iBin} = zeros(totalNumofUnits, noBins);

        for ii = 1:length(setofUnits{evt})
            
            unit           = setofUnits{evt}(ii);
            unitSpikeTimes = spikeTimes(find(spikeUnit == unit));
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