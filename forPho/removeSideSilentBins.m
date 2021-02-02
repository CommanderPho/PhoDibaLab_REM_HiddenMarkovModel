function [eventsBinnedfiring2, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring, eventBeg2, binDur, fileinfo)

% This function removes the silent (non-firing) bins from the beginnig and
% end of each event.

nEvents  = size(eventsBinnedfiring, 1); %% original events

eventsBinnedfiring2 = cell(size(eventsBinnedfiring));
eventLen = zeros(nEvents, 1);

eventBeg = zeros(nEvents, 1);
eventEnd = zeros(nEvents, 1);


% Remove the surrounding silent bins from each event

for evt = 1:nEvents
    
    eventData = eventsBinnedfiring{evt, 2}; % working with the coarse time bins (e.g., 20 ms for PBEs)
    
    binnSpikes = sum(eventData, 1); % number of spikes within each bin
    
    firstBin = find(binnSpikes, 1, 'first');
    lastBin  = find(binnSpikes, 1, 'last');
    
    if ~isempty(firstBin) && ~isempty(lastBin) 
        
        eventLen(evt) = lastBin-firstBin+1; %% length of the truncated event
        
        
        % Create the new trunctaed events (for 20ms-binning and also 1ms-binnig accordingly)
        eventsBinnedfiring2{evt, 2} = eventData(:, firstBin:lastBin);
        eventsBinnedfiring2{evt, 1} = eventsBinnedfiring{evt, 1}(:, floor((firstBin-1)*binDur*1000+1) : floor(lastBin*binDur*1000));
        
        
        
        % update the start and end times of the event
        eventBeg(evt) = eventBeg2(evt)+(firstBin-1)*binDur*fileinfo.Fs;
        eventEnd(evt) = eventBeg2(evt)+lastBin*binDur*fileinfo.Fs;

    end
    
end


end
