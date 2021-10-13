function [flatSpikeTimes, flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(unitSpikeCells, includedCellIDs)
%fnUnitSpikeCells2FlatSpikes Converts a cell array of the spikes (grouped by unit) to a flat spikes representation (where each spike is represented by a time in flatSpikeTimes and a corresponding cell ID in flatSpikeUnitIDs)
% Inverse of fnFlatSpikesToUnitCells

    num_unit_cells = length(unitSpikeCells);
    
    if ~exist('includedCellIDs','var')
        % If no IDs provided, generate IDs 1:N for each unit
        includedCellIDs = 1:num_unit_cells;
    end

    if size(unitSpikeCells) ~= size(includedCellIDs) 
        error('unitSpikeCells is not the same size as includedCellIDs!')
    end
    
    flatSpikeTimes = [];
    flatSpikeUnitIDs = [];
    for i = 1:num_unit_cells
        flatSpikeTimes = [flatSpikeTimes; unitSpikeCells{i}];
        flatSpikeUnitIDs = [flatSpikeUnitIDs; repmat(includedCellIDs(i),[length(unitSpikeCells{i}) 1])];
    end

    [flatSpikeTimes, sort_indicies] = sort(flatSpikeTimes);
    flatSpikeUnitIDs = flatSpikeUnitIDs(sort_indicies);
end

