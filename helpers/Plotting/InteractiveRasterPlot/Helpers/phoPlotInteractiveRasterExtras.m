classdef phoPlotInteractiveRasterExtras
    %PHOPLOTINTERACTIVERASTEREXTRAS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end % end regular methods

    methods(Static)

        function [annotations] = testAnnotations_phoTrialRectangles()
            %% testAnnotations_phoTrialRectangles: build test annotations
            % annotations: cell array of RasterplotAnnotation objects.
            annotationLists.unitIDs = [{1, 5, 9}, {2, 3, 4, 7}]; % The absolute unit ID (original/unfiltered ID) of each unit included
            annotationLists.referenceTimes = [1080, nan; 1440, 1467]; % The index timestamps of which to make a line-annotation
            % the scalar is an example of a point annotation (to highlight an event) and the second a window annotation (to highlight a range of timestamps)
            annotationLists.notes = {'', ''};
            
            for i = 1:length(annotationLists.notes)
                startTimestamp = annotationLists.referenceTimes(i, 1);
                endTimestamp = annotationLists.referenceTimes(i, 2);
                comment = annotationLists.notes{i};
                typeName = 'Temp';
                annotations{i} = RasterplotAnnotation(typeName, startTimestamp, endTimestamp, comment, annotationLists.unitIDs{i});
            end

        end

        function obj = phoPlotInteractiveRasterExtras(inputArg1, inputArg2)
            %PHOPLOTINTERACTIVERASTEREXTRAS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end


        function [linearPos, unitSpikeLaps, unitIsRippleSpike] = processSpikeStructExtendedExtras(spikeStruct)
            if ~exist('spikeStruct','var')
                error('spikeStruct not defined!')
                load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/allVariables.mat')
            end
              [includedCellIDs, unitSpikeCells, unitFlatIndicies] = fnFlatSpikesToUnitCells(spikeStruct.t, spikeStruct.unit, true);
              linearPos = cellfun(@(x) spikeStruct.linearPos(find(x)), unitFlatIndicies, 'UniformOutput', false);
              unitIsRippleSpike = cellfun(@(x) spikeStruct.ripple(find(x)), unitFlatIndicies, 'UniformOutput', false);
              unitSpikeLaps = cellfun(@(x) spikeStruct.ripple(find(x)), unitFlatIndicies, 'UniformOutput', false);
        end

        function [t, t_rel, x, y, linearPos] = processFileInfoPositionExtendedExtras(fileinfo)
            % Loads the position info out of a 'fileinfo' structure
            if ~exist('fileinfo','var')
                error('fileinfo not defined!')
                load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/toAddVariables.mat')
            end
              x = fileinfo.xyt(:,1);
              y = fileinfo.xyt(:,2);
              t = fileinfo.xyt(:,3);
              t_rel = ((t-fileinfo.tbegin) ./ 1e6); % Convert to relative timestamps since start
              % t_rel is aligned with the timestamps in the active_processing.position_table's timestamp column
              linearPos = fileinfo.xyt2(:,1);
        end



        function [extraControlsPanel] = addExtraControls(plotted_figH, mainRasterPlotAx)
%             f = figure( 'Position', 200*ones(1,4) );
%             vbox = uix.VBox( 'Parent', f );
%             axes( 'Parent', vbox );
%             hbox = uix.HButtonBox( 'Parent', vbox, 'Padding', 5 );
%             uicontrol( 'Parent', hbox, ...
%                 'String', 'Button 1' );
%             uicontrol( 'Parent', hbox, ...
%                 'String', 'Button 2' );
%             set( vbox, 'Heights', [-1 35] )
%             extraControlsPanel = f;

            % Add the contents
%             mainLayout = uix.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
            mainLayout = uix.VBoxFlex('Parent', plotted_figH, 'Spacing', 3 );
            
            % Create the panels
            gui.controlPanel = uix.BoxPanel( ...
               'Parent', mainLayout, ...
               'Title', 'Select a demo:' );
            hbox = uix.HButtonBox( 'Parent', gui.controlPanel, 'Padding', 5 );
            uicontrol( 'Parent', hbox, ...
                'String', 'Button 1' );
            uicontrol( 'Parent', hbox, ...
                'String', 'Button 2' );



            gui.ViewPanel = uix.BoxPanel( ...
               'Parent', mainLayout, ...
               'Title', 'Viewing: ???', ...
               'HelpFcn', @onDemoHelp );
            
            set(mainRasterPlotAx, 'Parent', gui.ViewPanel);

            % Adjust the main layout
            set( mainLayout, 'Widths', [-1,-2]  );
    
            extraControlsPanel = gui;
        end
        

        function [subplots_rect] = reallocateForAddingSubplots(mainRasterPlotAx, subplots_normalized_height)
            %% Resizes the main raster plot figure to add subplots above or below it
            % Get the position of the main raster plot axes:
            temp.updated_main_axes_pos = mainRasterPlotAx.Position;    
            if ~exist('subplots_normalized_height','var')
                subplots_normalized_height = 0.05;
            end
            main_plot_reduction_needed = subplots_normalized_height - 0.05; % 0.05 (5%) fits within the normal raster plots top margin, but anything greater than that will require reducing the rasterplot by that height
            if (main_plot_reduction_needed > 0.001)
                mainRasterPlotAx.Position(4) = temp.updated_main_axes_pos(4) - main_plot_reduction_needed; 
                mainRasterPlotAx.Position(2) = temp.updated_main_axes_pos(2) - main_plot_reduction_needed; % Also need to shift the raster plot downward by the reduction amount (I think)
                % Get the updated size for positioning the subplot:
                temp.updated_main_axes_pos = mainRasterPlotAx.Position; 
            end
            %% Puts Above the main raster plot box:
            subplots_rect = temp.updated_main_axes_pos;
            subplots_rect(2) = temp.updated_main_axes_pos(2) + temp.updated_main_axes_pos(4);
            subplots_rect(4) = subplots_normalized_height; % Manually constrain height, while keeping width the same as the main raster figure
        end

        function [stateMapHandle] = addStateMapSubplot(active_processing, mainRasterPlotAx, positionRectangle)
            % positionRectangle: the frame where the subplot will be rendered.

            %% Add the state_map:
            state_statemapPlottingOptions.orientation = 'horizontal';
            state_statemapPlottingOptions.plot_variable = 'behavioral_state';
            state_statemapPlottingOptions.vertical_state_mode = 'combined';
            state_statemapPlottingOptions.x_axis = 'timestamp'; % Timestamp-appropriate relative bins
            
            % Get the position of the main raster plot axes:
            if ~exist('positionRectangle','var')
                temp.updated_main_axes_pos = mainRasterPlotAx.Position;    
                
                %% Puts Above the main raster plot box:
                temp.statemap_pos = temp.updated_main_axes_pos;
                temp.statemap_pos(2) = temp.updated_main_axes_pos(2) + temp.updated_main_axes_pos(4);
                temp.statemap_pos(4) = 0.05;
            else
                % If provided, installs the state map in the position rectangle
                temp.statemap_pos = positionRectangle;
            end
            
            stateMapHandle = [];
            subplot('Position', temp.statemap_pos);
            [stateMapHandle] = fnPlotStateDiagram(active_processing, state_statemapPlottingOptions);
            
            % Link the state map to the main raster plot. Important for when the raster plot is scrolled after the scrollHandles are added. 
            linkaxes([mainRasterPlotAx stateMapHandle],'x'); 
        end

        function [ax, h] = addPositionSubplot(timestamps_array, timeseries_cells, mainRasterPlotAx, currPlotHandles, positionRectangle)
            %% Add the positions:
            % timeseries_cells: normall a cell array containing {x, y} with each cell having the same length as timestamps_array
            curr_num_plots = length(timeseries_cells);
            
            % Subdivide its height into equal rectangles
            temp.currNumOfHeightSubdivisions = curr_num_plots;

            if ~exist('positionRectangle','var')
            % Get the position of the main raster plot axes:
                temp.updated_main_axes_pos = mainRasterPlotAx.Position;    
                temp.subplotHeight = (temp.updated_main_axes_pos(3)-0.12) ./ temp.currNumOfHeightSubdivisions;
            else
                temp.updated_main_axes_pos = positionRectangle;
                temp.subplotHeight = positionRectangle(4) ./ temp.currNumOfHeightSubdivisions;
            end
            
            % Get the subplot's y-offset for each subplot:
            hold on;

            for i = 1:curr_num_plots
            %     ax(i) = subplot(numBlurredSpikeOutputs,1,i);
                % Need to convert to parent-space:
                currSubplotPositionRect = temp.updated_main_axes_pos;
                currSubplotPositionRect(4) = temp.subplotHeight; % Set to the common height

                if ~exist('positionRectangle','var')
                    currSubplotPositionRect(2) = temp.updated_main_axes_pos(2) + ((i-1) * temp.subplotHeight);
                else
%                     currSubplotPositionRect(2) = (temp.updated_main_axes_pos(2) - positionRectangle(4)) + ((i-1) * temp.subplotHeight); % Need to get bottom of the subplot rect for offsets to be right
                    currSubplotPositionRect(2) = temp.updated_main_axes_pos(2) + ((i-1) * temp.subplotHeight); % Need to get bottom of the subplot rect for offsets to be right
                end

                %% Reversed Y-axis:
            %     currSubplotPositionRect(2) = temp.currRasterAxisPosition(3) - currSubplotPositionRect(2);
            %     [figureXPoints, figureYPoints] = axescoord2figurecoord(figureXPoints, figureYPoints);
                ax(i) = axes('Position', currSubplotPositionRect,'Color','none');
                % Normalize the blurredSpikeOutputs down to unit height for plotting:
                h(i) = plot(ax(i), seconds(timestamps_array), timeseries_cells{i});
                ax(i).Color = 'none';
                xlabel(ax(i), [])
                xticks(ax(i), [])
                box off
                yticks(ax(i), [])
                ylabel(ax(i),'')
            end
        
            %% Set the current window to the specified range:
            % xlim(ax, xlim(scrollHandles.ParentAxesHandle))
            linkaxes([currPlotHandles.axesHandle ax],'x'); % Link all blurred axes to the main rasterplot axes
        end



        function [ax, h] = addBlurredSpikeOverlays(unitStatistics, scrollHandles, currPlotHandles, plotting_options)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            numBlurredSpikeOutputs = length(unitStatistics.blurredSpikeOutputs);
            unitStatistics.blurredStats.max = cellfun(@max, unitStatistics.blurredSpikeOutputs);
            unitStatistics.blurredStats.min = cellfun(@min, unitStatistics.blurredSpikeOutputs);
            unitStatistics.blurredStats.range = unitStatistics.blurredStats.max - unitStatistics.blurredStats.min;
        
            %% Axes:
        
            % temp.currRasterAxisPosition = currPlotHandles.axesHandle.Position;
            temp.currRasterAxisPosition = scrollHandles.ParentAxesHandle.Position;
            % Subdivide its height into equal rectangles
            temp.currNumOfHeightSubdivisions = length(plotting_options.trialSelection.TrialBackgroundRects.pos);
        
            temp.subplotHeight = (temp.currRasterAxisPosition(3)-0.12) ./ temp.currNumOfHeightSubdivisions;
            % Get the subplot's y-offset for each subplot:
            hold on;
            for i = 1:numBlurredSpikeOutputs
            %     ax(i) = subplot(numBlurredSpikeOutputs,1,i);
                % Need to convert to parent-space:
                currSubplotPositionRect = temp.currRasterAxisPosition;
                currSubplotPositionRect(4) = temp.subplotHeight; % Set to the common height
                currSubplotPositionRect(2) = temp.currRasterAxisPosition(2) + ((i-1) * temp.subplotHeight);
        
                %% Reversed Y-axis:
            %     currSubplotPositionRect(2) = temp.currRasterAxisPosition(3) - currSubplotPositionRect(2);
            %     [figureXPoints, figureYPoints] = axescoord2figurecoord(figureXPoints, figureYPoints);
                ax(i) = axes('Position', currSubplotPositionRect,'Color','none');
                % Normalize the blurredSpikeOutputs down to unit height for plotting:
                h(i) = plot(ax(i), seconds(temp.curr_timesteps_array{2}), (unitStatistics.blurredSpikeOutputs{i} ./ unitStatistics.blurredStats.range(i)));
                ax(i).Color = 'none';
                xlabel(ax(i), [])
                xticks(ax(i), [])
                box off
                yticks(ax(i), [])
                ylabel(ax(i),'')
            end
        
            %% Set the current window to the specified range:
            % xlim(ax, xlim(scrollHandles.ParentAxesHandle))
            linkaxes([currPlotHandles.axesHandle ax],'x'); % Link all blurred axes to the main rasterplot axes
        end


        function [ax, h] = addPeriodOverlays(startTimes, endTimes, color)
            % Plots periods along an axis
            % startTimes: T x 1 list of period start timestamps
            % endTimes: T x 1 list of period end timestamps
            % color: the color to display for each period:

            %% Example:
            % source_data.ripple.RoyMaze1.time


        end

    end % end static methods block

    %% Static Testing Method Block:
    methods(Static)
    
        function [outputs] = test_addPeriodOverlays()
            ripples_time_mat = evalin("caller", 'source_data.ripple.RoyMaze1.time');
            start_times = ripples_time_mat(:,1);
            end_times = ripples_time_mat(:,2);
            color = 'b';
            [outputs.ax, outputs.h] = phoPlotInteractiveRasterExtras.addPeriodOverlays(start_times, end_times, color);
        end


    end % end static Testing methods block

end

