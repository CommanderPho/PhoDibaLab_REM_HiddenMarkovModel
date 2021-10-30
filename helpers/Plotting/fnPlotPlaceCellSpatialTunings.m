function [fig, h] = fnPlotPlaceCellSpatialTunings(spatialTunings, varargin)
% function [fig, h] = fnPlotPlaceCellSpatialTunings(spatialTunings, plotting_options)
%FNPLOTPLACECELLSPATIALTUNINGS plot the place fields
%   Detailed explanation goes here

%% History: Extracted from Kourosh's spatialTuning_1D_tempModifications.m file by PhoHale on 2021-10-27

%% Example:
% [~, ~] = fnPlotPlaceCellSpatialTunings(PF_sorted_biDir,'linearPoscenters', linearPoscenters, 'unitLabels', num2cellstr(plot_outputs.original_unit_index));


    p = inputParser;
    addParameter(p,'linearPoscenters', 1:size(spatialTunings,2), @isnumeric)
    addParameter(p,'unitLabels', num2str(1:size(spatialTunings,1)), @iscellstr)
    addParameter(p,'plotting_options', struct(), @isstruct)
    
    % addParameter(p,'minPeakRate',3,@isnumeric)
    parse(p, varargin{:})
    
    linearPoscenters = p.Results.linearPoscenters;
    unitLabels = p.Results.unitLabels;
    plotting_options = p.Results.plotting_options; 
    
    PF_sorted_norm = spatialTunings ./ repmat(max(spatialTunings, [], 2), [1 size(spatialTunings, 2)]); % Normalize the peaks to one for visulaization

    % plot the place fields
    % units with peak firing rates below the threshold (2 Hz) will be shown in black
    fig = figure;
    x0=0;
    y0=0;
    width=400;
    height=400* size(PF_sorted_norm, 1)/20;

    %% TODO: stuff that would need to be brought in from the original file to work completely:
%     linearPos = fileinfo.xyt2(:, [1 3]); % the third row indicates lap indices of the positions 
%     % defining the position bins
% %     nPosBins = floor((max(linearPos(:, 1)) - min(linearPos(:, 1)))/posBinSize);
%     posBinEdges = min(linearPos(:, 1)): posBinSize: max(linearPos(:, 1)); % center of the position bins
%     linearPoscenters = posBinEdges(1:end-1) + posBinSize/2;


    difpos = linearPoscenters(2)- linearPoscenters(1);
    set(gcf,'units','points','position',[x0,y0,width,height])
    tt = 0;

    dynamic_colors = colormap(hsv(size(PF_sorted_norm, 1))); % should be the RGB triplets

    %% Loop over the units
    for jj = 1 : size(PF_sorted_norm, 1)
        
    %     if peakRates_sorted(jj) > 1
            
            tt = tt + 1;
            % Get the appropriate color for the curve based on its status in the combinedflag array
%             if combinedflag(jj) == 0
%                 cl = 'r';
%             elseif combinedflag(jj) == 1
%                 cl = [238,130,238]/255; % violet
%             end

            cl = dynamic_colors(jj, :); % get the specific color from the map
            
            %% Config the y-spacing between each subplot
            curr_y_offset_factor = 0.06*tt; % Normal y-offset
%             curr_y_offset_factor = 0; % stacked-up graph mode
%             curr_y_offset_factor = 0.01*tt; % small offset only graph mode

            curr_y_offset = curr_y_offset_factor + PF_sorted_norm(jj, :)/20;

            % first is 92x2
            % second term is 1x216 for some reason
            h.fill = fill([linearPoscenters fliplr(linearPoscenters)], [curr_y_offset fliplr(curr_y_offset_factor*ones(size(PF_sorted_norm(jj, :))))], cl,'LineStyle','none');
            alpha(0.5)
            hold on
            % first is 92x1
            % second is 1x108 double
            h.line = plot(linearPoscenters, curr_y_offset,'color', 'k','linewidth', 0.5);
            alpha(0.5)
            
    %         sparsity(tt) = mean(PF_sorted_norm(jj, :))^2 / mean(PF_sorted_norm(jj, :).^2);
    %         peakBins(tt) = peakPosBin_sorted(jj) * posBinSize;
            
    %     end
        
        
        set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
        text(linearPoscenters(1)-5*difpos, curr_y_offset_factor, unitLabels{jj}, 'fontsize', 10, 'HorizontalAlignment', 'center');
%         text(linearPoscenters(end)+5*difpos, 0.06*jj, sprintf('%.1f Hz', peakRates_sorted(jj)), 'fontsize', 7, 'HorizontalAlignment', 'center');
    end
    % tt = tt + 3;
    % cl = 'k';
    % 
    % fill([linearPoscenters fliplr(linearPoscenters)], [0.06*tt+sum_PF/2 fliplr(0.06*tt*ones(size(sum_PF)))], cl,'LineStyle','none')
    % 
    % hold on
    % plot(linearPoscenters, 0.06*tt+sum_PF/2,'color', 'k','linewidth', 0.5);
    % 
    % alpha(0.5)
    xlim([linearPoscenters(1)-difpos/2 linearPoscenters(end)+difpos/2])
    xlabel('Position on track (cm)', 'fontsize', 10)
%     h = text(linearPoscenters(1)-5*difpos, 0.06*(tt)/2, 'Unit', 'fontsize', 10, 'HorizontalAlignment', 'center');
%     set(h, 'rotation', 90)
    % 
    % h = text(linearPoscenters(1)-5*difpos, 0.06*tt, 'Norm pooled', 'fontsize', 10, 'HorizontalAlignment', 'center');
    % set(h, 'rotation', 90)

%% Add direction arrow indicator:    
%     if ~strcmp(direction, 'uni')
%         ha = annotation('arrow');  
%         ha.Parent = gca;  
%         if strcmp(direction, 'LR')
%             ha.X = [linearPoscenters(1) linearPoscenters(15)]; 
%         elseif strcmp(direction, 'RL')
%             ha.X = [linearPoscenters(end) linearPoscenters(end-15)]; 
%         end
%         ha.Y = [0.06*(tt+2) 0.06*(tt+2)];   
%         ha.LineWidth  = 2;          % make the arrow bolder for the picture
%         ha.HeadWidth  = 10;
%         ha.HeadLength = 10;
%     end


end

