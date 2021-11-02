function [sortedPeakPlaces, sortedTuningCurveIndicies] = fnSortPlaceCellSpatialTuningCurves(placeFieldTuningCurves, PositionBins, originalUnitIDs)
% fnSortPlaceCellSpatialTuningCurves - One line description of what the function performs
% Detailed explanation goes here
% 
% Syntax:  
%     [sortedPeakPlaces, sortedTuningCurveIndicies] = fnSortPlaceCellSpatialTuningCurves(placeFieldTuningCurves, PositionBins, originalUnitIDs)
% 
% Input:
%    placeFieldTuningCurves - Description
%    PositionBins - Description
%    originalUnitIDs - Description
% 
% Output:
%     - Description
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 02-Nov-2021 ; Last revision: 02-Nov-2021 

% ------------- BEGIN CODE --------------

    [~, peakIndex] = max(placeFieldTuningCurves, [], 2);
    % convert the indices stored in maxIndex to position on the track
    peakPlace = PositionBins(peakIndex);
%     [sortedPeakPlaces, sortedTuningCurveIndicies] = sort(peakPlace, 2, "ascend");
    [sortedPeakPlaces, sortedTuningCurveIndicies] = sort(peakPlace, 1, "ascend");
end


% ------------- END OF CODE --------------
