function [sortedPeakPlaces, sortedTuningCurveSortIndicies, sortedTuningCurveUnitIDs, inverseSortingIndicies] = fnSortPlaceCellSpatialTuningCurves(placeFieldTuningCurves, PositionBins, originalUnitIDs)
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
%     - sortedTuningCurveSortIndicies: the sort indices (in terms of units)
%     to sort by the spatial tuned positions. NOT the cellIDs
%     - sortedTuningCurveUnitIDs: the actual unitIDs corresponding to the
%     sorted curves. Only can be produced if you pass in the optional
%     originalUnitIDs argument
%     - inverseSortingIndicies: the reciprocal of
%     sortedTuningCurveSortIndicies. Things sorted by
%     sortedTuningCurveSortIndicies can be returned to the original unit
%     sort order by applying inverseSortingIndicies
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
    [sortedPeakPlaces, sortedTuningCurveSortIndicies] = sort(peakPlace, 1, "ascend");

    if exist('originalUnitIDs','var')
        sortedTuningCurveUnitIDs = originalUnitIDs(sortedTuningCurveSortIndicies);
    end

    % Inverse sort indicies:
    inverseSortingIndicies = zeros(size(sortedTuningCurveSortIndicies));
    for i = 1:length(sortedTuningCurveSortIndicies)
        inverseSortingIndicies(i) = find(sortedTuningCurveSortIndicies == i);
    end

end


% ------------- END OF CODE --------------
