function [laps, turningPeriods] = calculateLapTimings(tpos, linearPos, Fs, speed, direction, FileBase)
% function [laps, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, FileBase)

    % speed: struct
        %% t: 3243057x1
        %% v: 3243057x1
    % direction: string, specificies the directionality
    % values: {'bi', 'uni'}
    % FileBase: string, file base path
    % 
        
    %% Outputs:
    % laps
    % turningPeriods
    %% position and speed data

    % Old Call style:
    % [laps, turningPeriods] = calculateLapTimings(fileinfo.xyt2(:, 2), fileinfo.xyt2(:, 1), fileinfo.Fs, speed, direction, FileBase)

    % tpos = fileinfo.xyt2(:, 2);
    % linearPos = fileinfo.xyt2(:, 1);

    tvelocity = speed.t;
    velocity = speed.v;

    % if the pos includes nans interpolate to replace them
    linearPos2 = subfn_fixPositionNaNs(linearPos, tpos);

    %% find the dimensions of the track

    %%% to do this I tried to find the positions where the animal spent most of
    %%% the times; which assumed to be the resting box and either of the
    %%% platforms

    low_speed_cutoff = 10; % minimum speed cm/s
    % [occupancy, posbin] = subfn_computeOccupancyWithSpeedCutoff(speed, tpos, linearPos2, low_speed_cutoff);
    [occupancy, posbin, localMaxima] = subfn_computeOccupancyWithSpeedCutoff(tpos, linearPos2, tvelocity, velocity, low_speed_cutoff);

    %% Plot occupancy figure:
    figure
    hold on

    plot(posbin, occupancy, 'k', 'linewidth', 2)
    plot(posbin(localMaxima), occupancy(localMaxima), '.r', 'markersize', 15)

    hold off
    xlabel('x', 'fontsize', 20)
    ylabel('occupancy(%)', 'fontsize', 20)
    set(gca, 'fontsize', 16)

    % saveas(gcf, [FileBase '/restOccupancyHist.fig'])


    %%% define the boundaries of the track
    hiOccPositions = posbin(localMaxima); %% high occupancy positions; the first 
                                        %%% and last could correspond to either ends of the track
                                        %%% should be confirmed visually

    trackBounds = zeros(1,2);
    trackBounds(1) = hiOccPositions(1) + 0.1*(hiOccPositions(end) - hiOccPositions(1)); %% positions behind this point would be considered within the left platform
    trackBounds(2) = hiOccPositions(end) - 0.1*(hiOccPositions(end) - hiOccPositions(1)); %% positions beyond this point would be considered within the right platform


    %% calculate the laps timings

    % to define the laps we will consider periods (could be comprised of 
    % little quiescence bouts as well) in whcih animal travels from one end to 
    % the other (travels with return don't count)


    trackPart = zeros(length(tpos), 1); %% assigning an index based on the location of the position points relative to the track ends

    trackPart(find(linearPos2 < trackBounds(1))) = 1; %% if within track_end 1 then the index is 1 
    trackPart(find(linearPos2 > trackBounds(2))) = 2; 
    %%% if neither of the ends (could be on the track or rest box assign 0)

                                                
    boundPass = [0; diff(trackPart)];   % type of transition at the positions
                                                        % end1 outbound --> -1
                                                        % end1 inbound  --> +1
                                                        % end2 outbound --> -2
                                                        % end2 inbound  --> +2
                                                        % within each part
                                                        % --> 0

    passPnts = find(boundPass); %% points with non-zero boundPass
    passTypes = boundPass(passPnts); %% the type of pass at the corresponding points


    %%% extracting the primary beginning and end points of the laps

    % LR laps (left to right)

    tempInd = find(passTypes(1:end-1) == -1 & passTypes(2:end) == +2); 

    LRbegins = passPnts(tempInd); %% the start positions of Left(end1) to right(end2) traversals 
    LRends = passPnts(tempInd + 1); %% the just next passing points

    LRlaps_prim = [LRbegins LRends]; 


    % RL laps (right to left)

    tempInd = find(passTypes(1:end-1) == -2 & passTypes(2:end) == +1); 

    RLbegins = passPnts(tempInd); %% the start positions of right(end2) to left(end1) traversals 
    RLends = passPnts(tempInd + 1); %% the just next passing points

    RLlaps_prim = [RLbegins RLends];


    %%% extending the currently found laps to the point the animal stops or
    %%% returns or there is nan

    % smoothing the position to compensate for small position recording errors 

    sigma = 5; %% do a little smoothing
    halfwidth = 15;
    win = gausswindow(sigma, halfwidth);
    linearPos2 = conv(linearPos2, win, 'same');


    LRlaps_sec = zeros(size(LRlaps_prim));
    RLlaps_sec = zeros(size(RLlaps_prim));

    maxExtension = 1000; %% the number of points for extension of the primary lap timings (1000 pnts := 1000/30 ~ 33 sec)

    turningPnts = [];

    rmIdx = zeros(size(LRlaps_prim, 1), 1);
    for lap = 1 : size(LRlaps_prim, 1)
        
        
        startpoint = LRlaps_prim(lap, 1);
        endpoint = LRlaps_prim(lap, 2);
        
        
        precedpoints = max(startpoint - maxExtension, 1) : startpoint;
        posdiff = [0; diff(linearPos2(precedpoints))];
        
        nanPnt = find(isnan(linearPos(precedpoints)), 1, 'last');
        peakPnt = find(posdiff .* [posdiff(2:end); 0] < 0, 1, 'last');
        
        peakPnt = max([nanPnt, peakPnt]);
        
        if ~isempty(peakPnt)
            LRlaps_sec(lap, 1) = precedpoints(peakPnt);
        else
            LRlaps_sec(lap, 1) = precedpoints(1);
        end
        
        
        postpoints = endpoint : min(endpoint + maxExtension, length(linearPos2));
        posdiff = [0; diff(linearPos2(postpoints))];
        
        nanPnt = find(isnan(linearPos(postpoints)), 1, 'first');
        peakPnt = find(posdiff .* [posdiff(2:end); 0] < 0, 1, 'first');
        
        peakPnt = min([nanPnt, peakPnt]);
        
        if ~isempty(peakPnt)
            LRlaps_sec(lap, 2) = postpoints(peakPnt);
        else
            LRlaps_sec(lap, 2) = postpoints(end);
        end
        
        
        if LRlaps_sec(lap,:) == LRlaps_prim(lap, :)
            rmIdx(lap) = 1;
        end
        
        % detect the periods the animal turns or explore the around
        
        posdiff = [0; diff(linearPos2(LRlaps_sec(lap, 1): LRlaps_sec(lap, 2)))];
        
        turnpoints = find(posdiff .* [posdiff(2:end); 0] < 0);
        
        turnPnt1 = turnpoints(1:2:end);
        turnPnt2 = turnpoints(2:2:end);
        
        
        try
            if numel(turnPnt1) > 1
                posatTurnPnt1 = linearPos2(LRlaps_sec(lap, 1)+turnPnt1);
                posatTurnPnt2 = linearPos2(LRlaps_sec(lap, 1)+turnPnt2);

                temp = [nan; abs(posatTurnPnt1(2:end)-posatTurnPnt2(1:end-1))];
                tempInices = find(temp < 20); % in cm

                turnPnt1(tempInices)   = [];
                turnPnt2(tempInices-1) = [];
            end
        
            if ~isempty(turnPnt1) && ~isempty(turnPnt2) 
                turningPnts = [turningPnts; [turnPnt1 turnPnt2]+LRlaps_sec(lap, 1)];
            end
        catch
            lap
        end
    end
    LRlaps_sec(find(rmIdx), :) = [];



    rmIdx = zeros(size(RLlaps_prim, 1), 1);
    for lap = 1 : size(RLlaps_prim, 1)
        
        startpoint = RLlaps_prim(lap, 1);
        endpoint = RLlaps_prim(lap, 2);
        
        precedpoints = max(startpoint - maxExtension, 1) : startpoint;
        posdiff = [0; diff(linearPos2(precedpoints))];
        
        peakPnt = find(posdiff .* [posdiff(2:end); 0] < 0, 1, 'last');
        nanPnt = find(isnan(linearPos(precedpoints)), 1, 'last');
        
        peakPnt = max([nanPnt, peakPnt]);
        
        if ~isempty(peakPnt)
            RLlaps_sec(lap, 1) = precedpoints(peakPnt);
        else
            RLlaps_sec(lap, 1) = precedpoints(1);
        end

        postpoints = endpoint : min(endpoint + maxExtension, length(linearPos2));
        posdiff = [0; diff(linearPos2(postpoints))];
        
        nanPnt = find(isnan(linearPos(postpoints)), 1, 'first');
        peakPnt = find(posdiff .* [posdiff(2:end); 0] < 0, 1, 'first');
        
        peakPnt = min([nanPnt, peakPnt]);
        
        if ~isempty(peakPnt)
            RLlaps_sec(lap, 2) = postpoints(peakPnt);
        else
            RLlaps_sec(lap, 2) = postpoints(end);
        end
        
        
        % detect the periods the animal turns or explore the around
        
        posdiff = [0; diff(linearPos2(RLlaps_sec(lap, 1): RLlaps_sec(lap, 2)))];
        
        turnpoints = find(posdiff .* [posdiff(2:end); 0] < 0);
        
        turnPnt1 = turnpoints(1:2:end);
        
        turnPnt2 = turnpoints(2:2:end);
        
        try
            if numel(turnPnt1) > 1

                posatTurnPnt1 = linearPos2(RLlaps_sec(lap, 1)+turnPnt1);
                posatTurnPnt2 = linearPos2(RLlaps_sec(lap, 1)+turnPnt2);

                temp = [nan; abs(posatTurnPnt1(2:end)-posatTurnPnt2(1:end-1))];
                tempInices = find(temp < 20); % in cm

                turnPnt1(tempInices)   = [];
                turnPnt2(tempInices-1) = [];
            end

        
            if ~isempty(turnPnt1) && ~isempty(turnPnt2) 
                turningPnts = [turningPnts; [turnPnt1 turnPnt2]+RLlaps_sec(lap, 1)];
            end
        catch
            lap
        end

        if RLlaps_sec(lap,:) == RLlaps_prim(lap, :)
            rmIdx(lap) = 1;
        end
        
    end
    RLlaps_sec(find(rmIdx), :) = [];



    if strcmp(direction, 'bi')
    laps.RL = tpos(RLlaps_sec); lapDur = diff(laps.RL')'; laps.RL(lapDur < 1, :) = []; RLlaps_sec(lapDur < 1, :) = [];
    laps.LR = tpos(LRlaps_sec); lapDur = diff(laps.LR')'; laps.LR(lapDur < 1, :) = []; LRlaps_sec(lapDur < 1, :) = [];

    elseif strcmp(direction, 'uni')
    
    if length(RLlaps_sec) > length(LRlaps_sec)
        laps     = tpos(RLlaps_sec); 
        lapsPnts = RLlaps_sec;
    else
        laps     = tpos(LRlaps_sec);
        lapsPnts = LRlaps_sec;
    end   
    lapDur = diff(laps')'; laps(lapDur < 1, :) = []; lapsPnts(lapDur < 1, :) = []; 
    end

    turningPeriods = tpos(turningPnts);


    figure;
    hold on
    plot(tpos/Fs, linearPos, '.k', 'markersize', 2)

    if strcmp(direction, 'bi')
        for ii = 1:size(RLlaps_sec, 1); plot(tpos(RLlaps_sec(ii, 1):RLlaps_sec(ii, 2))/Fs, linearPos2(RLlaps_sec(ii, 1):RLlaps_sec(ii, 2)), 'r', 'linewidth', 3); end
        for ii = 1:size(LRlaps_sec, 1); plot(tpos(LRlaps_sec(ii, 1):LRlaps_sec(ii, 2))/Fs, linearPos2(LRlaps_sec(ii, 1):LRlaps_sec(ii, 2)), 'b', 'linewidth', 3); end
        
    elseif strcmp(direction, 'uni')
        for ii = 1:size(lapsPnts, 1); plot(tpos(lapsPnts(ii, 1):lapsPnts(ii, 2))/Fs, linearPos2(lapsPnts(ii, 1):lapsPnts(ii, 2)), 'b', 'linewidth', 3); end

    end 


    for ii = 1:size(turningPeriods, 1); plot(tpos(turningPnts(ii, 1):turningPnts(ii, 2))/Fs, linearPos2(turningPnts(ii, 1):turningPnts(ii, 2)), 'color', [0.7 0.7 0.7], 'linewidth', 3, 'linestyle', '-'); end

    % plot(tpos/Fs, linearPos2+1, 'color', 'g', 'linewidth', 2, 'linestyle', '-')

    hold off
    set(gca, 'fontsize', 12)
    xlabel('time(sec)', 'fontsize', 14)
    ylabel('position', 'fontsize', 14)

    % Extracted out into separate file:
    mkdir([FileBase '/lapInfo'])
    saveas(gcf, [FileBase '/lapInfo/laps.fig'])



end

function linearPos2 = subfn_fixPositionNaNs(linearPos, tpos)
    % if the pos includes nans interpolate to replace them
%     linearPos2 = subfn_fixPositionNaNs(linearPos, tpos);

    nanIdx = find(isnan(linearPos));
    
    linearPos_woNans = linearPos;
    linearPos_woNans(isnan(linearPos)) = [];
    tpos_woNans = tpos;
    tpos_woNans(isnan(linearPos)) = [];
    
    temp = interp1(tpos_woNans, linearPos_woNans, tpos(nanIdx));
    
    linearPos2 = linearPos;
    linearPos2(nanIdx) = temp;
end

function [occupancy, posbin, localMaxima] = subfn_computeOccupancyWithSpeedCutoff(tpos, linearPos2, tvelocity, velocity, low_speed_cutoff)
    speed_at_pos = interp1(tvelocity, velocity, tpos)';
    loSpeedPositions = linearPos2(find(speed_at_pos < low_speed_cutoff)); % positions with speed less than 10 cm/s
    
    %%% finding the positions with highest occupancy 
    [occupancy, posbin] = hist(loSpeedPositions, 100);
    occupancy = occupancy/sum(occupancy) * 100;
    
    
    sigma = 5; %% do a little smoothing
    halfwidth = 15;
    win = gausswindow(sigma, halfwidth);
    occupancy = conv(occupancy, win, 'same');

    %%% finding the local maxima of the occupancy (three maxima are expected coresponding to restbox and either ends of the track)
    diffOcc = [0; diff(occupancy)']; %% descending or ascending to the current point
    diffOcc_sb = [diffOcc(2:end); 0]; %% ascending or descending to the next point
    localMaxima = find(diffOcc > 0 & diffOcc_sb < 0);
    if numel(localMaxima) == 1
        localMaxima = [localMaxima length(diffOcc)];
    end

end