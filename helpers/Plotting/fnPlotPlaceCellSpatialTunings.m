function [fig, h] = fnPlotPlaceCellSpatialTunings(spatialTunings, varargin)
% function [fig, h] = fnPlotPlaceCellSpatialTunings(spatialTunings, plotting_options)
%FNPLOTPLACECELLSPATIALTUNINGS plot the place fields
%   Detailed explanation goes here

%% History: Extracted from Kourosh's spatialTuning_1D_tempModifications.m file by PhoHale on 2021-10-27

%% Example:
% [~, ~] = fnPlotPlaceCellSpatialTunings(PF_sorted_biDir,'linearPoscenters', linearPoscenters, 'unitLabels', num2cellstr(plot_outputs.original_unit_index));

    numSpatialTuningCurves = size(spatialTunings, 1);

    p = inputParser;
    addParameter(p,'linearPoscenters', 1:size(spatialTunings,2), @isnumeric)
    addParameter(p,'unitLabels', num2cellstr(1:numSpatialTuningCurves), @iscellstr)
    addParameter(p,'unitColors', colormap(jet(numSpatialTuningCurves)), @isnumeric) % should be [numUnits x 3] RGB triplets
    addParameter(p,'sortOrder', 1:numSpatialTuningCurves, @isnumeric) % should be [numUnits 1] set of indicies to sort the spatialTunings with
    addParameter(p,'colorSortOrder', 1:numSpatialTuningCurves, @isnumeric) % should be [numUnits 1] set of indicies to sort the colors with with
    addParameter(p,'peaks', [], @isnumeric)
    addParameter(p,'plotting_options', struct(), @isstruct)
    
    
    % addParameter(p,'minPeakRate',3,@isnumeric)
    parse(p, varargin{:})
    
    linearPoscenters = p.Results.linearPoscenters;
    unitLabels = p.Results.unitLabels;
    dynamic_colors = p.Results.unitColors;
    sortOrder = p.Results.sortOrder;
    colorSortOrder = p.Results.colorSortOrder;
    peaks = p.Results.peaks;
    plotting_options = p.Results.plotting_options;
    
%     if ~isempty(sortOrder)
%         % Sort the spatialTunings if that option is specified
%         spatialTunings = spatialTunings(sortOrder , :);
%         unitLabels = unitLabels(sortOrder); % resort the labels too, so when plotted they correctly correspond to the plot
%     end
% 
%     if ~isempty(colorSortOrder)
%         dynamic_colors = dynamic_colors(colorSortOrder, :);
%     end

    if isempty(sortOrder)
        sortOrder = 1:numSpatialTuningCurves;
    end

    if isempty(colorSortOrder)
        colorSortOrder = 1:numSpatialTuningCurves;
    end

    [spatialTunings_unnormalizedPeakValues, spatialTunings_peakIndices] = max(spatialTunings, [], 2);
    spatialTunings_norm = spatialTunings ./ repmat(spatialTunings_unnormalizedPeakValues, [1 size(spatialTunings, 2)]); % Normalize the peaks to one for visulaization

    if ~isempty(peaks)
        spatialTunings_normalizedPeakValues = zeros([numSpatialTuningCurves 1]);
    end

    %% An unknown magic constant that was initially set to 20
    magic_unit = 20;
    %magic_unit = numSpatialTuningCurves;

    % plot the place fields
    % units with peak firing rates below the threshold (2 Hz) will be shown in black
    fig = figure;
    x0=0;
    y0=0;
    width=400;
    height=400* numSpatialTuningCurves/magic_unit; % why 20?

    difpos = linearPoscenters(2)- linearPoscenters(1);
    set(gcf,'units','points','position',[x0,y0,width,height])
    tt = 0;

    %% Loop over the units    
    for jj = 1 : numSpatialTuningCurves
        
        tt = tt + 1;
        % Get the appropriate color for the curve based on its status in the combinedflag array
%             if combinedflag(jj) == 0
%                 cl = 'r';
%             elseif combinedflag(jj) == 1
%                 cl = [238,130,238]/255; % violet
%             end

        cl = dynamic_colors(colorSortOrder(jj), :); % get the specific color from the map
        
        %% Config the y-spacing between each subplot
        curr_y_offset_factor = 0.06*tt; % Normal y-offset
%             curr_y_offset_factor = 0; % stacked-up graph mode
%             curr_y_offset_factor = 0.01*tt; % small offset only graph mode
        curr_y_offset_factor_midpoint = (0.06*(tt) + 0.03);
        next_y_offset_factor = 0.06*(tt+1);

        curr_y_offset = curr_y_offset_factor + spatialTunings_norm(sortOrder(jj), :)/magic_unit;

        % first is numUnits x 2
        % second term is 1 x (2.0 * numUnits)
        h.fill = fill([linearPoscenters fliplr(linearPoscenters)], [curr_y_offset fliplr(curr_y_offset_factor*ones(size(spatialTunings_norm(sortOrder(jj), :))))], cl,'LineStyle','none');
        alpha(0.5)
        hold on
        % first is 92x1
        % second is 1x108 double
        h.line = plot(linearPoscenters, curr_y_offset, 'color', 'k','linewidth', 0.5);
        alpha(0.5)

        if ~isempty(peaks)
            spatialTunings_normalizedPeakValues(sortOrder(jj)) = spatialTunings_norm(sortOrder(jj), spatialTunings_peakIndices(sortOrder(jj)));
            curr_y_peak_height = curr_y_offset_factor + spatialTunings_normalizedPeakValues(sortOrder(jj))/magic_unit;
            h.peakLine = line([peaks(sortOrder(jj)), peaks(sortOrder(jj))], [curr_y_offset_factor, curr_y_peak_height],'Color','black','LineWidth', 1);
            alpha(0.5);
        end
            
        set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
        if exist('unitLabels','var')
            
            text(0.0, curr_y_offset_factor_midpoint, unitLabels{sortOrder(jj)}, 'fontsize', 11, 'HorizontalAlignment', 'center', 'Color', 'k'); % black background text for legibility
            text(0.0, curr_y_offset_factor_midpoint, unitLabels{sortOrder(jj)}, 'fontsize', 10, 'HorizontalAlignment', 'center', 'Color', cl);
        end
    end
    xlim([linearPoscenters(1)-difpos/2 linearPoscenters(end)+difpos/2])
    xlabel('Position on track (cm)', 'fontsize', 10)

end

