% +------------------------------------------------------+
% |      Interactive Section Selection on 2D Plot        |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        09/18/16 | 
% +------------------------------------------------------+
% 
% function: secsel(hPlot)
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

function secsel(hPlot)

    % determine the axes limits                 
    axlms = axis(get(hPlot, 'Parent'));                                                                
    xlen = axlms(2) - axlms(1);	

    % section start-line creation
    hLine1 = line([axlms(1)+xlen/3, axlms(1)+xlen/3], [axlms(3), axlms(4)],...
                  'LineStyle', '-.', 'LineWidth', 2, 'Color', 'none'); 

    % section end-line creation      
    hLine2 = line([axlms(1)+2*xlen/3, axlms(1)+2*xlen/3], [axlms(3), axlms(4)],...
                  'LineStyle', '-.', 'LineWidth', 2, 'Color', 'none'); 

    % merge the mark lines' handles          
    hLines = [hLine1, hLine2]; 

    % set the callback functions
    set(gcf, 'WindowButtonDownFcn', @createlines)
    set(gcf, 'KeyPressFcn', @savedata)

    % update the Handles Structure
    handles.hPlot = hPlot;
    handles.hLines = hLines;
    guidata(gcf, handles)

end


function createlines(src, eventdata)

    % Note: Bellow, a special code is developed in order to determine which  
    % of the two mouse button is double clicked. The Matlab do not support  
    % such feature as an build-in functionality.

    % retrive the Handles Structure
    handles = guidata(gcf);

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
        if left
            % the left mouse button is double clicked, so place the first mark
            % line to the new pointer possition 
            set(handles.hLines(1), 'XData', [ppos(1,1), ppos(1,1)], 'Color', 'k')   
        else
            % the right mouse button is double clicked, so place the first mark
            % line to the new pointer possition 
            set(handles.hLines(2), 'XData', [ppos(1,1), ppos(1,1)], 'Color', 'k') 
        end

        % updata the figure
        drawnow
    end

    % update the Handles Structure
    guidata(gcf, handles)

end


function savedata(src, eventdata)

    % retrive the Handles Structure
    handles = guidata(gcf);

    switch eventdata.Key

        case 'return'
            % determine the mark lines' abscissa values
            secstart = get(handles.hLines(1), 'XData');
            secend = get(handles.hLines(2), 'XData');
            lpos = sort([secstart(1), secend(1)]);

            % determine the selected data
            xdata = get(handles.hPlot, 'XData');
            ydata = get(handles.hPlot, 'YData');
            xsel = xdata(xdata >= lpos(1) & xdata <= lpos(2));
            ysel = ydata(xdata >= lpos(1) & xdata <= lpos(2));
            data = [xsel(:), ysel(:)]; %#ok<NASGU>

            %% TODO: Enable saving to workspace variable instead of out to file.
            % organize a save dialog menu
            [filename, pathname] = uiputfile({'*.mat', 'MAT-file (*.mat)'}, 'Save as...');
            if filename == 0, return, end

            % try to save the selected data in a .mat file
            try
                save([pathname filename], 'data')
            catch
                errordlg('An error was emerge!', 'Error', 'modal');
                return
            end

        case 'end'
            % end of the selsec function operation
            set(gcf, 'WindowButtonDownFcn', '') % TODO: Do I need to restore the default one, or the one set by scrollplot(...)?
            set(gcf, 'KeyPressFcn', '')
            delete(handles.hLines)

    end

end
        