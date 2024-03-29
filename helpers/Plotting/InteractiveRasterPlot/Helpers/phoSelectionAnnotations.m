% +------------------------------------------------------+
% |      Interactive Section Selection on 2D Plot        |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        09/18/16 | 
% +------------------------------------------------------+
% 
% function: disableDraggingHandles(hPlot)
%
% Input:
% hPlot - plot handle where the interactive section selection must be
%         applied
% 
% Output:
% N/A
%
% Instructions for use:
% - Double click the left mouse button to place the first mark line to a 
%   new possition;
% - Double click the right mouse button to place the second mark line to a 
%   new possition;
% - Press "Enter" key in order to save the selected 2D plot section as a
%   mat. data file;
% - Press "End" key in order to close the interactive section selection
%   regime.

function phoSelectionAnnotations(hPlot, plottingOptions)

    %% plottingOptions: a structure with the following properties:
        %% Rectangle Selection Properties:
        % plottingOptions.trialSelection.TrialBackgroundRects.handles
        % plottingOptions.trialSelection.TrialBackgroundRects.pos % [nTrials x 4]
        % plottingOptions.trialSelection.TrialBackgroundRects.pos(:, 2) % y-offsets [0.5, 1.5, 2.5, ...]


    % determine the axes limits                 
    axlms = axis(get(hPlot, 'Parent'));                                                                
    xlen = axlms(2) - axlms(1);	

    %% Creates Invisible Lines to begin with ('Color', 'none') which are repositioned when they're made visible.
    % section start-line creation
    hLine1 = line([axlms(1)+xlen/3, axlms(1)+xlen/3], [axlms(3), axlms(4)],...
                  'LineStyle', '-.', 'LineWidth', 2, 'Color', 'none'); 

    % section end-line creation      
    hLine2 = line([axlms(1)+2*xlen/3, axlms(1)+2*xlen/3], [axlms(3), axlms(4)],...
                  'LineStyle', '-.', 'LineWidth', 2, 'Color', 'none'); 

    % merge the mark lines' handles          
    hLines = [hLine1, hLine2]; 

    % set the callback functions
    set(gcf, 'WindowButtonDownFcn', @mouseDownCallback_phoSelectionAnnotations)
    set(gcf, 'KeyPressFcn', @keyPressCallback_phoSelectionAnnotations)

    % Probably need to use
%     'WindowButtonMotionFcn'
    

    % update the Handles Structure
    handles.hPlot = hPlot;
    handles.hLines = hLines;
    
    
    % Pack Gui Data:
    dataPack.handles = handles;
    dataPack.plottingOptions = plottingOptions;
    
    %% Build the required variables for interaction
    dataPack.plottingOptions.trialSelection.TrialBackgroundRects.numTrials = length(dataPack.plottingOptions.trialSelection.TrialBackgroundRects.handles);
    dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected = false([dataPack.plottingOptions.trialSelection.TrialBackgroundRects.numTrials 1]);
    
    
    guidata(gcf, dataPack)
    
    
    % Render any passed-in annotations:
    performRenderAnnotations_phoSelectionAnnotations(gcf, dataPack.plottingOptions.annotations);
    

end


function mouseDownCallback_phoSelectionAnnotations(src, eventdata)

    % Note: Bellow, a special code is developed in order to determine which  
    % of the two mouse button is double clicked. The Matlab do not support  
    % such feature as an build-in functionality.

    % retrive the Gui data Structure
    dataPack = guidata(gcf);
    % Un-Pack Gui Data:
    handles = dataPack.handles;
    
    % make the variable "left" persistent - its value must 
    % retained in memory between calls of the function
    persistent left

    % determine which button is double clicked - the left or the right
    if strcmp(get(src, 'SelectionType'), 'normal')
        left = true;
        pause(0.5) 
    elseif strcmp(get(src, 'SelectionType'), 'alt')
        left = false;
        pause(0.5)   
    elseif strcmp(get(src, 'SelectionType'), 'open')
        % get the current pointer possition
        ppos = get(gca, 'CurrentPoint');
        
        %% Test Trial Rectangle Functionality
        temp.rectYOffsets = dataPack.plottingOptions.trialSelection.TrialBackgroundRects.pos(:, 2); % y-offsets [0.5, 1.5, 2.5, ...]
        % Need height of elements by subtracting y-offsets (any two, they should all be the same
        temp.rectHeight = temp.rectYOffsets(2) - temp.rectYOffsets(1);
        temp.rectHalfHeight = temp.rectHeight ./ 2;
        
%         temp.rectYOffsetStarts = temp.rectYOffsets - temp.rectHalfHeight;
%         temp.rectYOffsetEnds = temp.rectYOffsets + temp.rectHalfHeight;
%         
        inRectIndex = max(1, floor(ppos(1,2) ./ temp.rectHeight)); % don't 0 index
        fprintf('Click was in rect[%d]\n', inRectIndex);
        
        % Toggle the isSelected property for this rectangle/trial:
%         dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected(inRectIndex) = ~dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected(inRectIndex);
        if (dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected(inRectIndex))
            % Deselect it:
            paramCell = struct2argsList(dataPack.plottingOptions.trialSelection.RectangleProperties.normal);
            dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected(inRectIndex) = false;
  
        else
             % Select:
            paramCell = struct2argsList(dataPack.plottingOptions.trialSelection.RectangleProperties.selected);
            dataPack.plottingOptions.trialSelection.TrialBackgroundRects.isSelected(inRectIndex) = true;
            
        end
        % Set the property on the rectangle to visually update its selection:
        set(dataPack.plottingOptions.trialSelection.TrialBackgroundRects.handles(inRectIndex), paramCell{:});


        
        %% Normal Functionality
        if left
            % the left mouse button is double clicked, so place the first mark
            % line to the new pointer possition 
            set(handles.hLines(1), 'XData', [ppos(1,1), ppos(1,1)], 'Color', dataPack.plottingOptions.lineVisibleColor)   
        else
            % the right mouse button is double clicked, so place the first mark
            % line to the new pointer possition 
            set(handles.hLines(2), 'XData', [ppos(1,1), ppos(1,1)], 'Color', dataPack.plottingOptions.lineVisibleColor) 
        end

        % updata the figure
        drawnow
    end

    % Pack Gui Data:
    dataPack.handles = handles;
    % update the Gui Data Structure
    guidata(gcf, dataPack)

end

%% keyPressCallback
function keyPressCallback_phoSelectionAnnotations(src, eventdata)

    switch eventdata.Key

        case 'l'
            %% Load from workspace
            loadAnnotations_phoSelectionAnnotations(gcf);
            
        case 's'
            %% Save to workspace:
            performSaveAnnotations_phoSelectionAnnotations(gcf);
                    
%         case 'return'
            
            
            
%             secstart = get(handles.hLines(1), 'XData');
%             secend = get(handles.hLines(2), 'XData');
%             lpos = sort([secstart(1), secend(1)]);
% 
%             % determine the selected data
%             xdata = get(handles.hPlot, 'XData');
%             ydata = get(handles.hPlot, 'YData');
%             xsel = xdata(xdata >= lpos(1) & xdata <= lpos(2));
%             ysel = ydata(xdata >= lpos(1) & xdata <= lpos(2));
%             data = [xsel(:), ysel(:)]; %#ok<NASGU>

        case 'end'
            % end of the selsec function operation
            set(gcf, 'WindowButtonDownFcn', '') % TODO: Do I need to restore the default one, or the one set by scrollplot(...)?
            set(gcf, 'KeyPressFcn', '')
            % retrive the Gui data Structure
            dataPack = guidata(gcf);
            % Un-Pack Gui Data:
            handles = dataPack.handles;

            delete(handles.hLines)
            % TODO: set the figure guidata to clear it properly?
            

    end
    
end




function addAnnotations_phoSelectionAnnotations(hFig, annotations)
    %% addAnnotations_phoSelectionAnnotations: add one or more annotation objects to the current annotations list. Should be called by the double click function and such.
    % annotations: cell array of RasterplotAnnotation objects.
%     obj = RasterplotAnnotation(typeName, startTimestamp, endTimestamp, comment, unitIDs);
    %% TODO: Develop the data model that keeps track of the annotations. Maybe look at //Matlab-Pho-Helper-Tools/Interactive/UserAnnotationsManager/UserAnnotationsManager.m for inspiration.
    

end


function performRenderAnnotations_phoSelectionAnnotations(hFig, desiredAnnotations)
    %% performRenderAnnotations_phoSelectionAnnotations: render the current set of annotations out to the figure as lines and windows
    %% TODO: Implement so that it doesn't redraw pointlessly every time a single annotation is added, but it can render the annotations easily.
    
    % retrive the Gui data Structure
    dataPack = guidata(hFig);
    
    %% Need Trial Rectangle Positions to Set y-positions of created objects:
    temp.rectYOffsets = dataPack.plottingOptions.trialSelection.TrialBackgroundRects.pos(:, 2); % y-offsets [0.5, 1.5, 2.5, ...]
    % Need height of elements by subtracting y-offsets (any two, they should all be the same
    temp.rectHeight = temp.rectYOffsets(2) - temp.rectYOffsets(1);
    temp.rectHalfHeight = temp.rectHeight ./ 2;

    numUpdatedRenderedAnnotations = length(desiredAnnotations);
    lineAnnotations.xPositions = []; % Create an empty array to hold line annotation x-positions
    lineAnnotations.yPositions = []; % Create an empty array to hold line annotation y-positions
    
    lineAnnotations.associatedData = {};
    
    for i = 1:numUpdatedRenderedAnnotations
        curr_start_timestamp = desiredAnnotations{i}.timestamp;
        
        if ~isnan(desiredAnnotations{i}.endTimestamp)
            % Has an end timestamp, meaning that it's a window
            % TODO:
            warning('Skipping rectangles for now!');
            
        else
            % No end timestamp means that it's a line
%             inRectIndex = floor(ppos(1,2) ./ temp.rectHeight);
%             fprintf('Click was in rect[%d]\n', inRectIndex);
            
            lineAnnotations.xPositions(end+1) = curr_start_timestamp;
            lineAnnotations.associatedData{end+1} = struct;
        end
        
        
    end
    
    
    %% convert to the form that the build lines function can read and build the lines from them:
    buildLines_phoSelectionAnnotations(lineAnnotations.xPositions, lineAnnotations.associatedData, hFig);
    
    %% TODO: Add rectangles option too:
    
    
    
end



function loadAnnotations_phoSelectionAnnotations(hFig)
    %% loadAnnotations_phoSelectionAnnotations: Loads the annotations from the base workspace
    % hFig: the figure handle
    
    fprintf('phoSelectionAnnotations: Loading annotations from base workspace (variable named phoSelectionAnnotations_selections)\n');
    % Load the phoSelectionAnnotations_selections variable from the base workspace
    phoSelectionAnnotations_selections = evalin('base','phoSelectionAnnotations_selections');
    
    %% loadAnnotations_phoSelectionAnnotations: Loads the annotations from the workspace or a file and builds lines from them:
    buildLines_phoSelectionAnnotations(phoSelectionAnnotations_selections.xPositions, phoSelectionAnnotations_selections.associatedData, hFig);
    
    drawnow;

end

function performSaveAnnotations_phoSelectionAnnotations(hFig)
    %% performSaveAnnotations_phoSelectionAnnotations: performs the task of saving the annotations to the base workspace or saving them to a file.
    temp.persistToBaseWorkspace = true;
    temp.persistToMatFile = false;
    
    % retrive the Gui data Structure
    dataPack = guidata(gcf);
    % Un-Pack Gui Data:
    handles = dataPack.handles;
    
    % determine the mark lines' abscissa values
    numLines = length(handles.hLines);
    xdata = get(handles.hPlot, 'XData');
    ydata = get(handles.hPlot, 'YData');

    phoSelectionAnnotations_selections.xPositions = zeros([numLines 1]);
    for i = 1:numLines
        currLineXData = get(handles.hLines(i), 'XData');
        currLineXData = currLineXData(1);
%         currLineXData = ceil(currLineXData(1));
%         phoSelectionAnnotations_selections.xPositions(i) = xdata(currLineXData); % Get the actual x-axis value that corresponds to the found index.
        phoSelectionAnnotations_selections.xPositions(i) = currLineXData;
    end

    %% TODO: this would overwrite the metadata, so be sure to change how this works before entering any info. It should be merged with the other data.
    %% TODO: save some sort of color/label data too?
    phoSelectionAnnotations_selections.associatedData = cell([numLines 1]);
    
    if temp.persistToBaseWorkspace
        % saving to workspace variable instead of out to file.
        fprintf('Lines saved to base workspace in variable named "phoSelectionAnnotations_selections"\n');
        assignin('base', 'phoSelectionAnnotations_selections', phoSelectionAnnotations_selections);
    end

    % Save out to file if that option is set:
    if temp.persistToMatFile
        % organize a save dialog menu
        [filename, pathname] = uiputfile({'*.mat', 'MAT-file (*.mat)'}, 'Save as...');
        if filename == 0, return, end

        % try to save the selected data in a .mat file
        try
            export_path = [pathname filename];
            fprintf('phoSelectionAnnotations: Saving annotations to mat file %s.\n', export_path);
            save(export_path, 'phoSelectionAnnotations_selections')
        catch
            errordlg('An error was emerge!', 'Error', 'modal');
            return
        end
    end

end

 
function buildLines_phoSelectionAnnotations(xPositions, yTrials, associatedData, hFig)
    %% buildLines_phoSelectionAnnotations: creates or updates the lines in the handles of the given figure
    
    % xPositions: a vector of x-positions at which to draw the lines
    % yTrials: a list of indicies to draw the track over. (TODO: how to deal with discontinuities?
    % associatedData: a structure-array containing other associated properties (like line colors, labels, etc) which to attach to the lines.
    
    
    %% TODO: figure out better way of drawing the lines over discontinuous tracks. Ideally they wouldn't all be separate annotations, or if they would they could be readjusted/deleted/labeled etc
    % could do a system of linking annotations as related in a group or something, but this got complicated pretty fast.
    
    % the unfiltered y-lims (spanning all units) is [axlms(3), axlms(4)]
    if ~exist('yTrials')
        
    end
    
    
    if ~exist('hFig','var')
        hFig = gcf;
    end
    % retrive the Gui data Structure
    dataPack = guidata(hFig);
    % Un-Pack Gui Data:
    handles = dataPack.handles;
    
    numUpdatedLines = length(xPositions);
        
    if ~isfield(handles, 'hLines')
       handles.hLines = gobjects(numUpdatedLines, 1); % Initialize new objects array
    end
    
    numExtantLines = length(handles.hLines);
    
    numNewLinesNeeded = numUpdatedLines - numExtantLines;
    if numNewLinesNeeded > 0
        % Need to create new lines
        % determine the axes limits                 
        axlms = axis(get(handles.hPlot, 'Parent'));                                                                
        xlen = axlms(2) - axlms(1);	

        for j = numExtantLines+1:numUpdatedLines
            %% Creates Invisible Lines to begin with ('Color', 'none') which are repositioned when they're made visible.
            handles.hLines(j) = line([axlms(1)+xlen/3, axlms(1)+xlen/3], [axlms(3), axlms(4)],...
                      'LineStyle', '-.', 'LineWidth', 2, 'Color', 'none'); 
        end
        
    elseif numNewLinesNeeded < 0
        % We need to remove lines:
        for j = numExtantLines+1:numUpdatedLines
            %% Creates Invisible Lines to begin with ('Color', 'none') which are repositioned when they're made visible.
%             set(handles.hLines(j), 'Color', 'none');
            delete(handles.hLines(j));
        end
        
    end
    
    % Loop through all lines to position them, including extant ones.
    for i = 1:length(xPositions)
        set(handles.hLines(i), 'XData', [xPositions(i), xPositions(i)], 'Color', dataPack.plottingOptions.lineVisibleColor)       
    end
    
    % Pack Gui Data:
    dataPack.handles = handles;
    % update the Gui Data Structure
    guidata(gcf, dataPack)

end
    
