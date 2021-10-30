% TrackLapTesting - One line description of what the script performs
% Detailed explanation goes here
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 29-Oct-2021 ; Last revision: 29-Oct-2021 


lapsInfo = smartload('/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/TrackLaps/trackLaps.mat', 'laps', 'lapsStruct');
%% laps are represented in absolute timestamps, convert to experiment relative timestamps
    % outputs will be aligned with the timestamps in the active_processing.position_table's timestamp column
lapsInfo.laps(:, 1:2) = ((lapsInfo.laps(:, 1:2) - fileinfo.tbegin) ./ 1e6); % Convert to relative timestamps since start
lapsInfo.numTotalLaps = size(lapsInfo.laps, 1); 

[~, ~] = fnJumpToLap(lapsInfo.laps, 5);


% ------------- END OF CODE --------------
