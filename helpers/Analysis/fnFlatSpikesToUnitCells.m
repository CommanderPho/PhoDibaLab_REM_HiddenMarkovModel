function [includedCellIDs, unitSpikeCells] = fnFlatSpikesToUnitCells(flatSpikeTimes, flatSpikeUnitIDs)
%fnFlatSpikesToUnitCells Converts a flat spikes representation (where each spike is represented by a time in flatSpikeTimes and a corresponding cell ID in flatSpikeUnitIDs) to a cell array of the spikes (grouped by unit).
%   Detailed explanation goes here

    if size(flatSpikeTimes) ~= size(flatSpikeUnitIDs) 
        error('flatSpikeTimes is not the same size as flatSpikeUnitIDs!')
    end
    
    includedCellIDs = unique(flatSpikeUnitIDs);
    numOutputCellIDs = length(includedCellIDs);
    unitSpikeCells = cell(numOutputCellIDs, 1);
    
    % unitSpikeCells = cellfun(
    for i = 1:numOutputCellIDs 
        unitSpikeCells{i} = flatSpikeTimes(includedCellIDs(i) == flatSpikeUnitIDs);
    end

end

