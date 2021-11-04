%% Conversion Info:
% PhoTuningCurveComputation2021 - This replaces Kourosh's place cell code, specifically "basicProcedure.m"
% Historical: Previously called 'PhoPositionAnalysis2021.m'
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 02-Nov-2021 ; Last revision: 02-Nov-2021 

% fileinfo.xyt2(:, 3) -> positionTable.lap_index
% fileinfo.xyt(:, 3) -> positionTable.t
% fileinfo.xyt2(:, 1) -> positionTable.linearPos
close all;
currDir = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/';
currSrcDir = fullfile(currDir, 'src');
animalID = 1;
sessionNumber = 1;
runSpeedThresh = 10; % 10cm/sec
unitTypes = 'all';
% Runs with the above settings
[sessionInfo] = PhoBayesianDecoder.getHiroExperimentName(animalID, sessionNumber);
fnPerformKouroshTrackProcessingComputations(sessionInfo, currDir, currSrcDir, runSpeedThresh, unitTypes);


function [] = fnPerformKouroshTrackProcessingComputations(sessionInfo, currDir, currSrcDir, runSpeedThresh, unitTypes)
   %% loading session data
%    %% session info
%    sessionName = [currRat 'Maze' num2str(sessionNumber)];
%    sessionName2 = [currRat '-maze' num2str(sessionNumber)];
%    animalMazeShapes = eval(sprintf('mazeShape.%s', currRat));
%    currMazeShape    = animalMazeShapes{sessionNumber};
%    rats = {'Roy','Ted', 'Kevin'};
%    allsessionNumbers = {[1 2 3],[1 2 3], [1]}; 
%    mazeShape.Roy = {'linear';'linear';'linear'};
%    mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
%    mazeShape.Kevin = {'linear'};
%    currRat = rats{animalID};

    %% OutputVariableLists:
    outputVariableNames.trackLaps = {'lapsStruct', 'turningPeriods', 'laps', 'totNumLaps', 'lapsTable', 'positionTable', 'currMazeShape', 'occupancyInfo', 'trackInfo'};
    outputVariableNames.toAddVariables = {'fileinfo', 'behavior'};
    outputVariableNames.spikesVariables = {'spikeStruct', 'okUnits', 'shanks', 'spikesTable'};
    outputVariableNames.biDirectional = {'okUnits', 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'positionBinningInfo_biDir'};
    outputVariableNames.PBEvariables = {'primaryPBEs', 'sdat', 'exclude', 'velocityFilter'};
    outputVariableNames.binnedPBEvariables = {'binnedPBEs', 'secondaryPBEs', 'qclus', 'rippleEvents','nPBEs','PBErippleIdx'};


    runOptions.shouldIncludeAnyPBE = false;
    runOptions.shouldIncludeExtendedPBE = false;
    

   %% Load from the sessionInfo object produced by PhoBayesianDecoder.getHiroExperimentName(animalID, sessionNumber):
   currRat = sessionInfo.ratName;
   currMazeShape = sessionInfo.mazeShape;
   sessionName = sessionInfo.sessionNameCamelCase;
   sessionName2 = sessionInfo.sessionNameHypenated;

   %%% load the data:
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
   fileinfo.lfpSampleRate = 1e6; % agian we use the time unit after conversion %TODO: HARDCODED_PARAMETER
   fileinfo.nCh = basics.nChannels;
   %%% position info and preprocessings
   fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 
   % linearized positions (we need this specially for analyzing the L-shape, U-shape, and circular mazes)
   % the information will be loaded into the xyt2 field of fileinfo structure
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
   %% all time
   fileinfo.tbegin = behavior.time(1,1); 
   fileinfo.tend   = behavior.time(3,2);
   save(fullfile(mainDir, 'toAddVariables.mat'), 'fileinfo', 'behavior');

   % Save out the laps info:
   subfolder = fullfile(mainDir, 'TrackLaps');
   mkdir(subfolder)
   % save(fullfile(subfolder, 'trackLaps.mat'), 'lapsStruct', 'turningPeriods', 'laps', 'totNumLaps', 'lapsTable', 'xyt', 'xyt2', 'currMazeShape', 'occupancyInfo', 'trackInfo')
   save(fullfile(subfolder, 'trackLaps.mat'), outputVariableNames.trackLaps{:})
   %% formating the spike info
   [spikeStruct, okUnits, spikesTable] = spikeBehaviorAnalysis(spikes, laps, rippleEvents, speed, unitTypes, fileinfo);
   %% TODO: to match the implementation produced by loadData.m, we need to convert the time column to relative times.
   	% Convert the first column (of timestamp offsets) to relative offsets into the experiment
% 	spikesTable.time = cellfun((@(timestamps) (timestamps ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp), ...
% 			spikesTable.time, ...
% 			'UniformOutput', false);

   temp = [spikes.id];
   shanks = temp(2*okUnits - 1);
   % save(fullfile(mainDir, 'allVariables.mat'), 'spikeStruct', 'okUnits', 'shanks')
   save(fullfile(mainDir, 'spikesVariables.mat'), outputVariableNames.spikesVariables{:})
   %% Tuning/PlaceFields Curves:
   close all
   subfolder = fullfile(mainDir, 'PlaceFields');
   mkdir(subfolder)


   tuningCurveStartTime = behavior.time(2,1);
   tuningCurveEndTime = behavior.time(2,2);
   
   %%% 1D spatial tuning: using linearized position:
%    [spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR, positionBinningInfo_LR] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, tuningCurveStartTime, tuningCurveEndTime, [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);
%    [spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL, positionBinningInfo_RL] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, tuningCurveStartTime, tuningCurveEndTime, [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);
   [spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir, positionBinningInfo_biDir] = spatialTuning_1D_Pho(spikeStruct, [1 2 3], positionTable, tuningCurveStartTime, tuningCurveEndTime, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);

   %% Per-interval mode:
   tuningCurveStartTimes = lapsTable.lapStartTime;
   tuningCurveEndTimes = lapsTable.lapEndTime;
   [spatialTunings_biDirCells, PF_sorted_biDirCells, runTemplate_biDirCells,  spatialInfo_biDirCells, conslapsRatio_biDirCells, diffWithAvg_biDirCells, positionBinningInfo_biDirCells] = computeTuningCurvesOverTimeRanges(spikeStruct, [1 2 3], positionTable, tuningCurveStartTimes, tuningCurveEndTimes, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, fileinfo.name);

   %% These "runTemplate"s are really important, they're the cellIDs for the active spatial maps!
   % spatialTunings_biDir: 60 x 105 %% The spatialTunings_* has the same number of units as 'okUnits' from the struct, but 10 of the rows are nothing but zeros, explaining why we only see 50 x 105 for the PF_sorted_biDir and the runTemplate_*
   % runTemplate_biDir: 50 x 1
   % PF_sorted_biDir: 50 x 105
   save(fullfile(subfolder, 'biDirectional.mat'), outputVariableNames.biDirectional{:})
 
   %% TODO: make a good output format:  
   
   
   %% PBEs during different behavioral periods
   if runOptions.shouldIncludeAnyPBE
       time_resolution = 0.001; % in second %TODO: HARDCODED_PARAMETER
       threshZ         = 3; % sdf with 3 std deviation above the mean %TODO: HARDCODED_PARAMETER
       qclus = [1 2 3]; %TODO: HARDCODED_PARAMETER
       velocityFilter = 1; %TODO: HARDCODED_PARAMETER

       subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
       mkdir(subfolder)
       exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4
       [primaryPBEs, sdat] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);
       save(fullfile(subfolder, 'PBEvariables.mat'), outputVariableNames.PBEvariables{:})
      
       if runOptions.shouldIncludeExtendedPBE
           binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) %TODO: HARDCODED_PARAMETER
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
           save(fullfile(subfolder, 'binnedPBEvariables.mat'), outputVariableNames.binnedPBEvariables{:})
        
%            %% Final form of output:
%            baseStruct = struct('data', [], 'p', [], 'ts', [], 'pts', []);
%         
%            PREbinnedPBEs       = baseStruct;
%            PREidx              = find(secondaryPBEs(:, 1) > behavior.time(1,1) & secondaryPBEs(:, 2) < behavior.time(1,2));
%            PREbinnedPBEs.data  = binnedPBEs(PREidx, :);
%            secondaryPBEs_PRE   = secondaryPBEs(PREidx, :);
%            % PREbinnedPBEs       = genSurrogates(PREbinnedPBEs);
%         
%            RUNbinnedPBEs       = baseStruct;
%            RUNidx              = find(secondaryPBEs(:, 1) > behavior.time(2,1) & secondaryPBEs(:, 2) < behavior.time(2,2));
%            RUNbinnedPBEs.data  = binnedPBEs(RUNidx, :);
%            secondaryPBEs_RUN   = secondaryPBEs(RUNidx, :);
%            % RUNbinnedPBEs       = genSurrogates(RUNbinnedPBEs);
%         
%            POSTbinnedPBEs       = baseStruct;
%            POSTidx              = find(secondaryPBEs(:, 1) > behavior.time(3,1) & secondaryPBEs(:, 2) < behavior.time(3,2));
%            POSTbinnedPBEs.data  = binnedPBEs(POSTidx, :);
%            secondaryPBEs_POST   = secondaryPBEs(POSTidx, :);
%            % POSTbinnedPBEs       = genSurrogates(POSTbinnedPBEs);
%            save(fullfile(mainDir, 'toAddVariables.mat'), 'fileinfo', 'behavior', 'secondaryPBEs')
       else
           fprintf('Skipping extended PBE calculation because runOptions.shouldIncludeExtendedPBE is false. binnedPBEvariables.mat will not be updated. \n');
       end
   else
           fprintf('Skipping any PBE calculations because runOptions.shouldIncludeAnyPBE is false. PBEvariables.mat, and binnedPBEvariables.mat will not be updated. \n');
   end
   fprintf('done.\n');

end % end function

function [spatialTunings, PF_sorted, template, spatialInfo, conslapsRatio, diffWithAvg, positionBinningInfo] = computeTuningCurvesOverTimeRanges(spikeStruct, qclus, positionTable, tuningCurveStartTimes, tuningCurveEndTimes, thetaPeriods, turningPeriods, speed, direction, posBinSize, runSpeedThresh, combinedUnits, timeUnit, FileBase, FileName)
    
%     spatialTunings, PF_sorted, template,  spatialInfo, conslapsRatio, diffWithAvg, positionBinningInfo = deal(cell(size(tuningCurveStartTimes)));
    
    for i = 1:length(tuningCurveStartTimes)
        tuningCurveStartTime = tuningCurveStartTimes(i);
        tuningCurveEndTime = tuningCurveEndTimes(i);
        [spatialTunings{i}, PF_sorted{i}, template{i},  spatialInfo{i}, conslapsRatio{i}, diffWithAvg{i}, positionBinningInfo{i}] = spatialTuning_1D_Pho(spikeStruct, qclus, positionTable, tuningCurveStartTime, tuningCurveEndTime, thetaPeriods, turningPeriods, speed, direction, posBinSize, runSpeedThresh, combinedUnits, timeUnit, FileBase, FileName);
    end
end

% ------------- END OF CODE --------------
