function spikeCounts = computeSpikeCounts(spikeTimes, tstart, tend, tau)

    nCells = length(spikeTimes);
    
    % Spike counts during each estimation time bin
    nTimeBins = ceil((tend-tstart)/tau);
    windows = linspace(tstart,tend,nTimeBins);
    
    spikeCounts = zeros(nTimeBins-1,nCells);
    
    for i = 1:nCells
        spikeCounts(:,i) = histcounts(spikeTimes{i},windows);
    end

end



