function [linearPos2, nanIdx] = fnFixPositionNaNs(linearPos, tpos)
% fnFixPositionNaNs: if the pos includes nans interpolate to replace them
%     linearPos2 = fnFixPositionNaNs(linearPos, tpos);
    nanIdx = find(isnan(linearPos));
    
    linearPos_woNans = linearPos;
    linearPos_woNans(isnan(linearPos)) = [];
    tpos_woNans = tpos;
    tpos_woNans(isnan(linearPos)) = [];
    
    num_nans_to_interpolate = length(tpos) - length(tpos_woNans);

    temp = interp1(tpos_woNans, linearPos_woNans, tpos(nanIdx)); % interpolate over the nans
    
     
    successfully_interpolated_nans = find(~isnan(temp));
    num_nans_successfully_interpolated = length(successfully_interpolated_nans);

    linearPos2 = linearPos;
    linearPos2(nanIdx) = temp;

%     unresolved_nanIdx = find(isnan(linearPos2));        
%     tpos(~isnan(linearPos2)), linearPos2(~isnan(linearPos2))


end