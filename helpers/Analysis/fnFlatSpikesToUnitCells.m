function [includedCellIDs, unitSpikeCells, unitFlatIndicies] = fnFlatSpikesToUnitCells(flatSpikeTimes, flatSpikeUnitIDs, includeIndicies)
%fnFlatSpikesToUnitCells Converts a flat spikes representation (where each spike is represented by a time in flatSpikeTimes and a corresponding cell ID in flatSpikeUnitIDs) to a cell array of the spikes (grouped by unit).
%   Detailed explanation goes here

    if size(flatSpikeTimes) ~= size(flatSpikeUnitIDs) 
        error('flatSpikeTimes is not the same size as flatSpikeUnitIDs!')
    end

    if ~exist('includeIndicies','var')
        includeIndicies = false;
    end
    
    includedCellIDs = unique(flatSpikeUnitIDs);
    numOutputCellIDs = length(includedCellIDs);
    unitSpikeCells = cell(numOutputCellIDs, 1);

    if includeIndicies
        unitFlatIndicies = cell(numOutputCellIDs, 1);
    else
        unitFlatIndicies = {};
    end
    
    % unitSpikeCells = cellfun(
    for i = 1:numOutputCellIDs
         if includeIndicies
            unitFlatIndicies{i} = find(includedCellIDs(i) == flatSpikeUnitIDs);
            unitSpikeCells{i} = flatSpikeTimes(unitFlatIndicies{i});
         else
            unitSpikeCells{i} = flatSpikeTimes(includedCellIDs(i) == flatSpikeUnitIDs);
         end
    end

end

