function [eventsBinnedfiring, secondaryPBEs, noFiringUnits, eventLen] = finalBinningResult(PBEs, spikeStruct, qclus, fileinfo, binDur)

    %% binnig the spikes within an event and qualifying the event based on number of active units and length
    
    
    %%% binning the spikes within each event

    if ~exist('binDur','var')
        binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) 
        % qclus = [1 2 3]; % Only the pyramidal neurons were included
    end


    [eventsBinnedfiring2, noFiringUnits] = timeBinning(PBEs, spikeStruct, qclus, binDur, fileinfo);
%     [eventsBinnedfiring2, noFiringUnits] = phoTimeBinning(PBEs, spikeStruct, qclus, binDur, fileinfo);
    
    % 
    % binOverlapRatio = 0.5;
    % [eventsBinnedfiring3, noFiringUnits2] = timeBinning_withOverlap(PBEs, spikeStruct, qclus, binDur, binOverlapRatio, fileinfo);
    
    % Remove the flanking zero-firing periods(silent bins) from the beginnig
    % and end of each event. Since the interneurons were also involved in calculating
    % the PBEs' boundaries and here we are binning just pyramidal units' spike
    % trains, silent bins are expected. 
    
    [eventsBinnedfiring, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring2, PBEs(:, 1), binDur, fileinfo);
    
    % Qualify the events based on the minumum active units and length criteria
    % Considering only the event with number of active units firing above 10 percent of the total or 5 (whichever is greater) and
    % duration of at least 4 20ms-time bins.
    
    
    % activeUnits = unique(spikeStruct.unit);
    activeUnits = 1:size(eventsBinnedfiring2{1,1}, 1);
    
    idx_acceptedEvents = find(noFiringUnits >= 5 & eventLen >= 4); % maximum length 500 ms (???)
    
    % idx_acceptedEvents = find(noFiringUnits >= max(5, floor(0.1 * length(activeUnits))) & (eventLen >= 3 & eventLen <= 15)); % maximum length 500 ms (???)
    % idx_acceptedEvents =  1:numel(eventBeg);
    
    secondaryPBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) PBEs(idx_acceptedEvents, 3:4)]; 
    
    eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses
    noFiringUnits = noFiringUnits(idx_acceptedEvents);
    
    numEvents = size(eventsBinnedfiring, 1);
    eventsBinnedfiring(:, 3) = mat2cell((1:numEvents)', ones(1, numEvents));

end
