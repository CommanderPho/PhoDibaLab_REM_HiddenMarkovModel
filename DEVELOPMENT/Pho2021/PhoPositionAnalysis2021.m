%% Conversion Info:
% PhoPositionAnalysis2021 - This replaces Kourosh's place cell code
% Detailed explanation goes here
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 02-Nov-2021 ; Last revision: 02-Nov-2021 

% fileinfo.xyt2(:, 3) -> positionTable.lap_index
% fileinfo.xyt(:, 3) -> positionTable.t
% fileinfo.xyt2(:, 1) -> positionTable.linearPos
currDir = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/';
currSrcDir = fullfile(currDir, 'src');
animalID = 3;
sessionNumber = 1;
%% loading session data
rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3],[1 2 3], [1]}; 
mazeShape.Roy = {'linear';'linear';'linear'};
mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
mazeShape.Kevin = {'linear'};
currRat = rats{animalID};
% sessionNumber
%% session info
%%% load the data
sessionName = [currRat 'Maze' num2str(sessionNumber)];
sessionName2 = [currRat '-maze' num2str(sessionNumber)];
VarList = {'spikes','behavior','position','speed','basics','ripple'};
for var = 1 : length(VarList)
    load([currSrcDir '/wake-' VarList{var} '.mat'])
end
spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['position.' sessionName]);
speed    = eval(['speed.' sessionName]);
basics   = eval(['basics.' sessionName]);
ripple       = eval(['ripple.' sessionName]);
rippleEvents = ripple.time;
%%% initiate a unified structure for the current session data and make a
%%% master folder for storing the results from different analyses
fileinfo = struct('name', sessionName2, 'animal', currRat, 'xyt', [], 'xyt2', [], ...
    'tbegin', [], 'tend', [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'pix2cm', 0.3861); 
mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
if ~exist(mainDir,'dir')
   mkdir(mainDir) % Make the output dir 
end
%%% behavior 
if strcmp(currRat, 'Kevin') % there are two wake periods, we need to concatenate them
   behavior.time(2,2) = behavior.time(3,2);
   behavior.time(3,1) = behavior.time(4,1);
   behavior.time(3,2) = behavior.time(4,2);
   behavior.time(4,:) = [];
end
behavior.time(2,1) = behavior.time(2,1) + 1e6; % slight modification to make sure this time period is only limited to the maze (1e6 = 1 sec)
behavior.time(2,2) = behavior.time(2,2) - 1e6;
bvrTimeList = behavior.list(:,[1 2]); % start and end of each behavioral state
bvrState    = behavior.list(:,3); % the label of corresponding states
%%% recording configurations 
% Fs = basics.SampleRate;
fileinfo.Fs = 1e6; % Since in Hiro's data the spike timestamps are already in microsecond we don't need to use the original recording sampling rate. 
            % The behavioral timepoints are also in microsecond therefore easier to handle them.
% lfpSampleRate = basics.lfpSampleRate;
fileinfo.lfpSampleRate = 1e6; % agian we use the time unit after conversion
fileinfo.nCh = basics.nChannels;
%%% position info and preprocessings
fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 
% linearized positions (we need this specially for analyzing the L-shape, U-shape, and circular mazes)
% the information will be loaded into the xyt2 field of fileinfo structure
animalMazeShapes = eval(sprintf('mazeShape.%s', currRat));
currMazeShape    = animalMazeShapes{sessionNumber};
specified_maze_platform_centers = [43.4908, 36.2828; 223.7903, 32.3178];
%% Linearize Position for the current maze shape:
[linearPos, userSelectedCenters] = linearizePosition2(fileinfo.xyt, behavior.time(2,1), behavior.time(2,2), currMazeShape, specified_maze_platform_centers);
% linearPos: 970086x1 double
fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);
%% Build a smarter output structure here:
positionTable = table(fileinfo.xyt(:, 3), fileinfo.xyt(:, 1), fileinfo.xyt(:, 2), linearPos, ...
   'VariableNames', {'t', 'x', 'y', 'linearPos'});
%% Lap Timings:
direction = 'bi';
[lapsStruct, turningPeriods, occupancyInfo, trackInfo] = calculateLapTimings(positionTable.t, positionTable.linearPos, fileinfo.Fs, speed, direction, mainDir);
if length(lapsStruct.RL) > length(lapsStruct.LR)
   lapsStruct.RL(1,:) = [];
   behavior.MazeEpoch(1,:) = lapsStruct.LR(1,:);
end
totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
laps = zeros(totNumLaps, 2);
laps(1:2:totNumLaps, :)  = lapsStruct.LR; % odd lap number entries are set from the LR values
laps(2:2:totNumLaps, :)  = lapsStruct.RL; % even lap number entries are set from the RL values
laps(:, 3) = 1:size(laps, 1); % the "lapnumber" itself is not monotonically increasing when sorted by startTime or endTime
lapsStruct.lapDirectionLabels = cell([totNumLaps 1]);
[lapsStruct.lapDirectionLabels{1:2:totNumLaps}] = deal('LR');
[lapsStruct.lapDirectionLabels{2:2:totNumLaps}] = deal('RL');
lapsTable = array2table(laps,...
    'VariableNames',{'lapStartTime','lapEndTime', 'lapNumber'});
lapsTable.direction = lapsStruct.lapDirectionLabels;
lapsTable.duration_sec = ((lapsTable.lapEndTime - lapsTable.lapStartTime) ./ 1e6);
 % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)
fileinfo.xyt2(:, 3) = zeros(size(fileinfo.xyt2(:, 1)));
for ii = 1: length(laps)
   idx = find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);         
end
positionTable.lap_index = fileinfo.xyt2(:, 3);
% Save out the laps info:
subfolder = fullfile(mainDir, 'TrackLaps');
mkdir(subfolder)
% save(fullfile(subfolder, 'trackLaps.mat'), 'lapsStruct', 'turningPeriods', 'laps', 'totNumLaps', 'lapsTable', 'xyt', 'xyt2', 'currMazeShape', 'occupancyInfo', 'trackInfo')
save(fullfile(subfolder, 'trackLaps.mat'), 'lapsStruct', 'turningPeriods', 'laps', 'totNumLaps', 'lapsTable', 'positionTable', 'currMazeShape', 'occupancyInfo', 'trackInfo')
% clear xyt2
%% formating the spike info
unitTypes = 'all';
[spikeStruct, okUnits] = spikeBehaviorAnalysis(spikes, laps, rippleEvents, speed, unitTypes, fileinfo);
temp = [spikes.id];
shanks = temp(2*okUnits - 1);
% save(fullfile(mainDir, 'allVariables.mat'), 'spikeStruct', 'okUnits', 'shanks')
save(fullfile(mainDir, 'spikesVariables.mat'), 'spikeStruct', 'okUnits', 'shanks')
%% Tuning/PlaceFields Curves:
close all
subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)

runSpeedThresh = 10; % 10cm/sec
%%% 1D spatial tuning: using linearized position:
[spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR, positionBinningInfo_LR] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, behavior.time(2,1), behavior.time(2,2), [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);
[spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL, positionBinningInfo_RL] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, behavior.time(2,1), behavior.time(2,2), [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);
[spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir, positionBinningInfo_biDir] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, behavior.time(2,1), behavior.time(2,2), [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);
save(fullfile(subfolder, 'biDirectional.mat'), 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'positionBinningInfo_biDir')
%% TODO: make a good output format:
% save(fullfile(mainDir, 'allVariables.mat'), 'spikeStruct', 'okUnits', 'shanks')
%% PBEs during different behavioral periods
time_resolution = 0.001; % in second
threshZ         = 3; % sdf with 3 std deviation above the mean
qclus = [1 2 3];
velocityFilter = 1;
%% all time
fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend   = behavior.time(3,2);
subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)
exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);
save(fullfile(subfolder, 'PBEvariables.mat'), 'primaryPBEs', 'sdat', 'exclude', 'velocityFilter')
save(fullfile(mainDir, 'toAddVariables.mat'), 'fileinfo', 'behavior');

binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) 
%% Just Testing:
%% Implementation
[binnedPBEs, secondaryPBEs] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileinfo, binDur); % Very slow function
PBErippleIdx = ifContainRipples(secondaryPBEs, rippleEvents);
secondaryPBEs(:, 5) = PBErippleIdx;
nPBEs = size(binnedPBEs, 1);
secondaryPBEs = [secondaryPBEs zeros(nPBEs, 4)];
for ii= 1:nPBEs
    pbeCenter = secondaryPBEs(ii, 3);
    boutInd        = find(bvrTimeList(:,1) < pbeCenter & bvrTimeList(:,2) > pbeCenter, 1, 'first');
    boutBrainState = bvrState(boutInd);
    secondaryPBEs(ii, 5 + boutBrainState) = 1;
end
save(fullfile(subfolder, 'binnedPBEvariables.mat'), 'binnedPBEs', 'secondaryPBEs', 'qclus', 'rippleEvents','nPBEs','PBErippleIdx')
function [lapsStruct, turningPeriods] = subfn_calculateLaps()
   direction = 'bi';
   % [lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
   [lapsStruct, turningPeriods] = calculateLapTimings(fileinfo.xyt2(:, 2), fileinfo.xyt2(:, 1), fileinfo.Fs, speed, direction, mainDir);
   if length(lapsStruct.RL) > length(lapsStruct.LR)
      lapsStruct.RL(1,:) = [];
      behavior.MazeEpoch(1,:) = lapsStruct.LR(1,:);
   end
   totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
   laps = zeros(totNumLaps, 2);
   laps(1:2:totNumLaps, :)  = lapsStruct.LR;
   laps(2:2:totNumLaps, :)  = lapsStruct.RL;
   laps(:, 3) = 1:size(laps, 1); 
   % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)
   fileinfo.xyt2(:, 3) = zeros(size(fileinfo.xyt2(:, 1)));
   for ii = 1: length(laps)
      idx = find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
      fileinfo.xyt2(idx, 3) = laps(ii, 3);         
   end
end
function [occupancy, posbin] = subfn_computeOccupancyWithSpeedCutoff(positionTable, velocityTable, low_speed_cutoff)
   %  low_speed_cutoff = 10; % minimum speed cm/s
   tpos = positionTable.t;
   linearPos = positionTable.linearPos;
   tvelocity = velocityTable.t;
   velocity = velocityTable.v;
    speed_at_pos = interp1(tvelocity, velocity, tpos)';
    loSpeedPositions = linearPos2(find(speed_at_pos < low_speed_cutoff)); % positions with speed less than 10 cm/s
    %%% finding the positions with highest occupancy 
    [occupancy, posbin] = hist(loSpeedPositions, 100);
    occupancy = occupancy/sum(occupancy) * 100;
    sigma = 5; %% do a little smoothing
    halfwidth = 15;
    win = gausswindow(sigma, halfwidth);
    occupancy = conv(occupancy, win, 'same'); % smoothed occupancy
end


% ------------- END OF CODE --------------
