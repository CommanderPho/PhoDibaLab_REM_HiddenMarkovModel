function [spatialTunings, PF_sorted, template, spatialInfo, conslapsRatio, diffWithAvg, positionBinningInfo] = spatialTuning_1D_Pho(spikeStruct, qclus, positionTable, maze_start_time, maze_end_time, thetaPeriods, turningPeriods, speed, direction, posBinSize, runSpeedThresh, combinedUnits, timeUnit, FileBase, FileName)
% spatialTuning_1D_Pho - One line description of what the function or script performs
% Detailed explanation goes here
% 
% Syntax:  
%     [spatialTunings, PF_sorted, template, spatialInfo, conslapsRatio, diffWithAvg] = spatialTuning_1D_Pho(spikeStruct, qclus, fileinfo, behavior, thetaPeriods, turningPeriods, speed, direction, posBinSize, runSpeedThresh, combinedUnits, timeUnit, FileBase)
% 
% Inputs:
%    spikeStruct - Description
%    qclus - Description
%    fileinfo - Description
%    behavior - Description
%    thetaPeriods - Description
%    turningPeriods - Description
%    speed - Description
%    direction - Description
%    posBinSize - Description
%    runSpeedThresh - Description
%    combinedUnits - Description
%    timeUnit - Description
%    FileBase - Description
%
% 
% Outputs:
%    template - for each spatial tuning, the set of cells active in the tuning form what's referred to in these files as the 'template'. These are the IDs of the original units provided in spikeStruct that are active given the filter critera.
%    output2 - Description
% 

% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 28-Sep-2021 ; Last revision: 28-Sep-2021 

% ------------- BEGIN CODE --------------

features.shouldPlot = false; % Whether to render the tuning curve plots or not
features.shouldSave = false; % Whether to save out the tuning curve plots to disk by default
features.shouldWriteTextFile = false; % Whether to save out .clu and .txt files

% EXPLAIN THE FUNCTION
totalNumofUnits = max(spikeStruct.unit);
% including only the spikes meeting certain criteria
if strcmp(direction, 'LR')
    desiredMod = 0; % even laps
    directionNum = 1;
    unitSortingMode = 'ascend';
elseif strcmp(direction, 'RL')
    desiredMod = 1; % odd laps
    directionNum = 2;
    unitSortingMode = 'descend';
elseif strcmp(direction, 'uni')
    desiredMod = [0 1];
    directionNum = 1;
    unitSortingMode = 'ascend';
end
spikeInd = find(spikeStruct.t >= maze_start_time & spikeStruct.t < maze_end_time ... % within the RUN period
                & ismember(spikeStruct.qclu, qclus) ... % only stable pyramidal units are included
                & spikeStruct.speed > runSpeedThresh ... % only the spikes happening when the velocity of the animal is higher than the threshold (usually 10 cm/sec)
                & spikeStruct.lap > 0 & ismember(mod(spikeStruct.lap, 2), desiredMod));  % only spike during high theta power
%                 & spikeStruct.theta == 1); % limiting to the travels in specific direction
% % excluding the spikes that occurr when the animal stopped or turned around in the middle of the track
turningSpikes = [];
for ii = 1: size(turningPeriods, 1); turningSpikes = [turningSpikes; find(spikeStruct.t >= turningPeriods(ii, 1) & spikeStruct.t <= turningPeriods(ii, 2))]; end     
spikeInd = setdiff(spikeInd, turningSpikes);
               
spikePositions = spikeStruct.linearPos(spikeInd);
spikeUnit      = spikeStruct.unit(spikeInd);
spikeLap       = spikeStruct.lap(spikeInd); %% to investigate lap-by-lap consistence of unit firings with the average tuning across the laps
activeUnits    = unique(spikeUnit); 
% finding position samples pertaining to the maze period and certain direction of travel 
speedatPos = interp1(speed.t, speed.v, positionTable.t);
positionIdx    = find(positionTable.t > maze_start_time & positionTable.t < maze_end_time ... % within the run period
                    & speedatPos > runSpeedThresh ... % position bins when animal's velocity is higher than threshold
                    & positionTable.lap_index > 0 & ismember(mod(positionTable.lap_index, 2), desiredMod)); % limiting to the travels in specific direction
                
                
% % include only the position bins within the theta peridos
if ~isempty(thetaPeriods)
    thetaPositionInd = []; 
    for ii = 1: size(thetaPeriods, 1); thetaPositionInd = [thetaPositionInd; find(positionTable.t >= thetaPeriods(ii, 1) & positionTable.t <= thetaPeriods(ii, 2))]; end     
    positionIdx = intersect(positionIdx, thetaPositionInd);
end
if ~isempty(turningPeriods)
    turningPositionInd = []; % should be excluded
    for ii = 1: size(turningPeriods, 1); turningPositionInd = [turningPositionInd; find(positionTable.t >= turningPeriods(ii, 1) & positionTable.t <= turningPeriods(ii, 2))]; end     
    positionIdx = setdiff(positionIdx, turningPositionInd);
end

linearPos = [positionTable.linearPos, positionTable.lap_index];
directionLinearPos = nan(size(linearPos));          
directionLinearPos(positionIdx, :) = linearPos(positionIdx, :);

% defining the position bins
[posBinEdges, linearPoscenters, nPosBins] = PhoFallAnalysis2021.subfn_buildPositionBinInformation(positionTable.linearPos, posBinSize);
% nPosBins = floor((max(linearPos(:, 1)) - min(linearPos(:, 1)))/posBinSize);
% posBinEdges = min(linearPos(:, 1)): posBinSize: max(linearPos(:, 1)); % center of the position bins
% linearPoscenters = posBinEdges(1:end-1) + posBinSize/2;
posSamplingPeriod = median(diff(positionTable.t))/timeUnit; % 1/sampling frequency - Hiro's dataset: timeunit is microsecond

%% Build positionBinningInfo output struct
positionBinningInfo.directionLinearPos = directionLinearPos;
positionBinningInfo.posBinSize = posBinSize;
positionBinningInfo.posBinEdges = posBinEdges;
positionBinningInfo.linearPoscenters = linearPoscenters;
positionBinningInfo.nPosBins = nPosBins;
positionBinningInfo.timeUnit = timeUnit;
positionBinningInfo.posSamplingPeriod = posSamplingPeriod;

%% Preinitialize:
spatialTunings = zeros(totalNumofUnits, nPosBins);
peakPosBin     = zeros(totalNumofUnits, 1);
currDirLaps   = unique(spikeLap);
combinedflag  = zeros(totalNumofUnits, 1);
diffWithAvg   = zeros(numel(currDirLaps), totalNumofUnits);
conslapsRatio = zeros(totalNumofUnits, 1); % number of laps with consistent peak position for each unit
spatialInfo   = zeros(totalNumofUnits, 1); % spatial information of each unit
spatial_smoothing_window = gausswindow(3,5); % smoothing with a standard deviation of 6 cm (assuming that each positon bins is 2cm long)

for ii = 1: length(activeUnits)    
    unit = activeUnits(ii);
     
    posBinnSpikeCnts_laps = zeros(nPosBins, numel(currDirLaps));
    posBinDwelltime_laps  = zeros(nPosBins, numel(currDirLaps));
        
    placeFields_laps = zeros(nPosBins, numel(currDirLaps));
    peakPosBin_laps = zeros(numel(currDirLaps), 1);
    for jj = 1:numel(currDirLaps)
        unitSpikePositions = spikePositions(spikeUnit == unit & spikeLap == currDirLaps(jj));
        positions = directionLinearPos(directionLinearPos(:, 2) == currDirLaps(jj), 1);
        temp   = histc(unitSpikePositions, posBinEdges);
        posBinnSpikeCnts_laps(:, jj) = temp(1:end-1);
        temp    = histc(positions, posBinEdges) * posSamplingPeriod;
        posBinDwelltime_laps(:, jj)  = temp(1:end-1);
        
        unsmoothed_tunings = posBinnSpikeCnts_laps(:, jj) ./ posBinDwelltime_laps(:, jj);
        unsmoothed_tunings(isnan(unsmoothed_tunings)) = 0;
        unsmoothed_tunings(isinf(unsmoothed_tunings)) = 0;
        placeFields_laps(:, jj) = conv(unsmoothed_tunings, spatial_smoothing_window, 'same');
        
        if sum(placeFields_laps(:, jj)) > 0
            [~, peakPosBin_laps(jj)] = max(placeFields_laps(:, jj));
        else
            peakPosBin_laps(jj) = nan;
        end
    
    end % end for jj
    
    % average across the laps:
    posBinnSpikeCnts = sum(posBinnSpikeCnts_laps, 2);
    posBinDwelltime  = sum(posBinDwelltime_laps, 2);
    
    unsmoothed_tunings = posBinnSpikeCnts ./ posBinDwelltime;
    unsmoothed_tunings(isnan(unsmoothed_tunings)) = 0;
    unsmoothed_tunings(isinf(unsmoothed_tunings)) = 0;
    
	% Smooth the tunings:
    spatialTunings(unit , :) = conv(unsmoothed_tunings, spatial_smoothing_window, 'same');
    [~, peakPosBin(unit)] = max(spatialTunings(unit , :));
    
    % investigate consistence of firings across the laps:
    tolerationLimit = 30; % in cm
    tolerationLimit = floor(tolerationLimit/posBinSize/2);
    
    diffWithAvg(:, unit) = peakPosBin_laps - peakPosBin(unit);
    consLaps = find(diffWithAvg(:, unit) > -tolerationLimit &  diffWithAvg(:, unit) < tolerationLimit); % consistent laps 
    conslapsRatio(unit) = numel(consLaps)/numel(currDirLaps);
    
    % spatial inforamtion
    pi = posBinDwelltime/sum(posBinDwelltime);
    tempPlaceMap = spatialTunings(unit , :) + 1e-4; % add a small value to remove zero firing rates
    
    FR = tempPlaceMap * pi; % average firing rate
    spatialInfo(unit) = (tempPlaceMap/FR .* log2(tempPlaceMap/FR)) * pi;
    
    if ismember(unit, combinedUnits)
        combinedflag(unit) = 1; 
    end 
end % end for ii 
placeCells = find(spatialInfo > 0 & conslapsRatio > 0); % here we are not using the consistence and spatial info to filter cells
[~, sortIdx] = sort(peakPosBin(placeCells), unitSortingMode);
PF_sorted    = spatialTunings(placeCells(sortIdx), :); % sort the place fields based on the peaks
combinedflag = combinedflag(placeCells(sortIdx), :);
 
%  activeUnits_sorted = activeUnits(sortIdx); % keep track of the unit labels
template = placeCells(sortIdx);
% peakRates_sorted = max(PF_sorted, [], 2); % maximum firing rate
PF_sorted_norm = PF_sorted ./ repmat(max(PF_sorted, [], 2), [1 size(PF_sorted, 2)]); % Normalize the peaks to one for visulaization
% sum_PF = mean(PF_sorted_norm, 1);
% plot the place fields
if features.shouldPlot
    % units with peak firing rates below the threshold (2 Hz) will be shown in black
    figure;
    x0=0;
    y0=0;
    width=400;
    height=400* size(PF_sorted_norm, 1)/20;
    difpos = linearPoscenters(2)- linearPoscenters(1);
    set(gcf,'units','points','position',[x0,y0,width,height])
    tt = 0;
    for jj = 1 : size(PF_sorted_norm, 1)
        tt = tt + 1;
        if combinedflag(jj) == 0
        cl = 'r';
        elseif combinedflag(jj) == 1
            cl = [238,130,238]/255; % violet
        end
        fill([linearPoscenters fliplr(linearPoscenters)], [0.06*tt+PF_sorted_norm(jj, :)/20 fliplr(0.06*tt*ones(size(PF_sorted_norm(jj, :))))], cl,'LineStyle','none')
        hold on
        plot(linearPoscenters, 0.06*tt+PF_sorted_norm(jj, :)/20,'color', 'k','linewidth', 0.5);
        alpha(0.5)
        
        set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
    end
    xlim([linearPoscenters(1)-difpos/2 linearPoscenters(end)+difpos/2])
    xlabel('Position on track (cm)', 'fontsize', 10)
    h = text(linearPoscenters(1)-5*difpos, 0.06*(tt)/2, 'Unit', 'fontsize', 10, 'HorizontalAlignment', 'center');
    set(h, 'rotation', 90)
    if ~strcmp(direction, 'uni')
        ha = annotation('arrow');  
        ha.Parent = gca;  
        if strcmp(direction, 'LR')
            ha.X = [linearPoscenters(1) linearPoscenters(15)]; 
        elseif strcmp(direction, 'RL')
            ha.X = [linearPoscenters(end) linearPoscenters(end-15)]; 
        end
        ha.Y = [0.06*(tt+2) 0.06*(tt+2)];   
        ha.LineWidth  = 2;          % make the arrow bolder for the picture
        ha.HeadWidth  = 10;
        ha.HeadLength = 10;
    end

    if ~isempty(FileBase) & features.shouldSave
        if ~strcmp(direction, 'uni')
            filename = [FileName '_placeFields1D_' direction];
        else
            filename = [FileName '_placeFields1D_uni'];
        end
        savefig(gcf, fullfile(FileBase, filename))
    %     savepdf(gcf, fullfile(FileBase, filename), '-dpng')
    %     exportgraphics(gcf,[fullfile(FileBase, filename) '.pdf'],'BackgroundColor','none','ContentType','vector')
    end 
end % end if shouldPlot

%%% making clu and res files
if features.shouldWriteTextFile
    rest = [];
    cluorder = [];
    shank = zeros(length(template), 1);
    cluster = zeros(length(template), 1);
    for ii = 1 : length(template)
        spkInd = find(spikeStruct.unit == template(ii));
        currSpiketimes = floor(spikeStruct.t(spkInd)*20000);
        rest = [rest; currSpiketimes];
        cluorder = [cluorder; ii*ones(length(currSpiketimes),1)];
        shank(ii) = 1; %unique(spikeStruct.shank(spkInd));
        cluster(ii) = unique(spikeStruct.unit(spkInd));
    end
    [sortedRest, ind] = sort(rest);
    cluorder = cluorder(ind)-1;
    numclu = [shank cluster];
    numclu = numclu';
    % Saveres([FileBase '/' FileName direction 'Active.10' num2str(directionNum) '.res.' num2str(directionNum)],sortedRest);
    % SaveClu([FileBase '/' FileName direction 'Active.10' num2str(directionNum) '.clu.' num2str(directionNum)],cluorder);
    % dlmwrite([FileBase '/' FileName direction 'Active.10' num2str(directionNum) '.numclu.' num2str(directionNum)],numclu);
    writematrix(sortedRest, [FileBase '/' FileName direction 'Active.10' num2str(directionNum) '.res.' num2str(directionNum) '.txt'])
    writematrix(cluorder, [FileBase '/' FileName direction 'Active.10' num2str(directionNum) '.clu.' num2str(directionNum) '.txt'])
end

end


% ------------- END OF CODE --------------
