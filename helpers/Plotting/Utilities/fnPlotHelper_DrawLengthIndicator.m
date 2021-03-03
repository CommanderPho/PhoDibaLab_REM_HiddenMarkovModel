function [plotHandle] = fnPlotHelper_DrawLengthIndicator(centerPoint, plottingOptions)
%FNPLOTHELPER_DRAWLENGTHINDICATOR Draws a measurement/length indicator line with an optional label
%   Detailed explanation goes here
%% Example:
    % figure
    % clf
    % 
    % centerPoint = [0 0];
    % plottingOptions.majorAxisStartEnd = [0 2];
    % [plotHandle] = fnPlotHelper_DrawLengthIndicator(centerPoint, plottingOptions);
    % 
    % xlim([-1 5]);
    % ylim([-1 5]);

    ax = gca;
    hold on;
    if ~exist('centerPoint','var')
        centerPoint = [0.0, 0.0]; % temporary, allow position.
    end

    if ~exist('plottingOptions','var')
         plottingOptions = struct();
    end
    plottingOptions = fnAddDefaultOptionalArgs({'orientation', 'majorAxisStartEnd', 'majorAxisSize', 'minorAxisSize', 'barbDirections', 'label', 'lineFormat'}, ...
        {'vertical', ...
        [0, 1.0], ...
        10.0, ...
        0.25, ...
        [true, true], ... % first element indicates to include barb in positive direction of minor axis, second in the negative dir.
        'trial1', ...
        struct();}, ...
        plottingOptions);
    
%     testPlottingOptions.plotMode = 'plot';
    testPlottingOptions.plotMode = 'figure';
    
    
    if ~isfield(plottingOptions, 'barbLineFormat')
       plottingOptions.barbLineFormat = plottingOptions.lineFormat;
    end
    
    % Specify [start end] along the major axis
    minorAxisMidpointOffset = (plottingOptions.minorAxisSize ./ 2.0);
    
    numIndicators = size(plottingOptions.majorAxisStartEnd, 1);
    
    includePositiveBarb = plottingOptions.barbDirections(1);
    includeNegativeBarb = plottingOptions.barbDirections(2);
    
    if strcmpi(plottingOptions.orientation, 'vertical')
        barbMinMaxMinorAxisPoints = [centerPoint(1), centerPoint(1)];
                
    else
        barbMinMaxMinorAxisPoints = [centerPoint(2), centerPoint(2)];
    end
    
    % barbMinMaxMinorAxisPoints: for horizontal, barbBottomTopPoints
    if includePositiveBarb
        barbMinMaxMinorAxisPoints = barbMinMaxMinorAxisPoints + [0, (minorAxisMidpointOffset)];
    end
    if includeNegativeBarb
        barbMinMaxMinorAxisPoints = barbMinMaxMinorAxisPoints + [(-minorAxisMidpointOffset), 0];
    end
    
     % Convert to argument lists:
    lineFormat = struct2argsList(plottingOptions.lineFormat);
    barLineFormat = struct2argsList(plottingOptions.barbLineFormat);
    
        
    numPreallocPoints = numIndicators * 2 * 3; % 2 points for each line, 3 lines total, times numIndicators
    %% Pre-allocate points:
    xPoints = [];
    yPoints = [];
    
    barbs.xPoints = [];
    barbs.yPoints = [];
    
    
    for i = 1:numIndicators
%         majorAxisStartEnd
 
        if strcmpi(plottingOptions.orientation, 'vertical')
            % If vertically oriented, flip the x and y values
            
            error('not yet implemented!')
                
        else
            % Horizontal (x-axis) is major axis
            % Main/Base Line:
            xPoints(end+1:end+3) = [plottingOptions.majorAxisStartEnd(i,:), nan];
            yPoints(end+1:end+3) = [centerPoint(2), centerPoint(2), nan];
            
            %% Barbs:
            % Left Barb
            barbs.xPoints(end+1:end+3) = [plottingOptions.majorAxisStartEnd(i,1), plottingOptions.majorAxisStartEnd(i,1), nan]; 
            barbs.yPoints(end+1:end+3) = [barbMinMaxMinorAxisPoints(1), barbMinMaxMinorAxisPoints(2), nan];
            
            
            % Right Barb:
            barbs.xPoints(end+1:end+3) = [plottingOptions.majorAxisStartEnd(i,2), plottingOptions.majorAxisStartEnd(i,2), nan];
            barbs.yPoints(end+1:end+3) = [barbMinMaxMinorAxisPoints(1), barbMinMaxMinorAxisPoints(2), nan];
            
        end
        
        % Draw Major Axis Line:
        
        
        % Draw Barbs (Along minor axis)
    end
    
   
    
    if strcmpi(testPlottingOptions.plotMode, 'plot')
        plotHandle.MainLine = plot(xPoints, yPoints, 'k', lineFormat{:});
        plotHandle.Barbs = plot(barbs.xPoints, barbs.yPoints, 'b', barLineFormat{:});

    else
        set(gcf,'Units','normalized');
        
        %% Convert to Figure Coordinates:
        figureXPoints = reshape(xPoints,3, [])';
        figureYPoints = reshape(yPoints,3, [])';
       
        
        barbs.figureXPoints = reshape(barbs.xPoints,3, [])';
        barbs.figureYPoints = reshape(barbs.yPoints,3, [])';
        
         % Drop the NaN Columns:
        figureXPoints(:,end) = [];
        figureYPoints(:,end) = [];
        
        barbs.figureXPoints(:,end) = [];
        barbs.figureYPoints(:,end) = [];
        
%         
%         figureXPoints = reshape(figureXPoints,[],3);
%         figureYPoints = reshape(figureYPoints,[],3);
                
        
        [figureXPoints, figureYPoints] = axescoord2figurecoord(figureXPoints, figureYPoints);
        [barbs.figureXPoints, barbs.figureYPoints] = axescoord2figurecoord(barbs.figureXPoints, barbs.figureYPoints);

        for i = 1:numIndicators
    %     h = annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', titles{s_idx});

            plotHandle.MainLine = annotation('line', figureXPoints(i,:), figureYPoints(i,:), lineFormat{:});
            
            minBarbIndex = 2*i - 1;
            plotHandle.Barbs = annotation('line', barbs.figureXPoints(minBarbIndex,:), barbs.figureYPoints(minBarbIndex,:), barLineFormat{:});
            maxBarbIndex = 2*i;
            plotHandle.Barbs = annotation('line', barbs.figureXPoints(maxBarbIndex,:), barbs.figureYPoints(maxBarbIndex,:), barLineFormat{:});
            
            %% TODO: Draw the other barb too, as there are two barbs added for each main line:
            

        end
        
    end
end

