function [ax] = fnPlotStateDiagram(active_processing, plottingOptions)
%fnPlotStateDiagram Plots a behavioral-state indicator as a function of time
%   Detailed explanation goes here


    if ~exist('plottingOptions','var')
         plottingOptions = struct();
    end
    
    plottingOptions = fnAddDefaultOptionalArgs({'orientation', 'vertical_state_mode', 'plot_variable', 'x_axis'}, ...
        {'vertical', 'combined', 'behavioral_state', 'timestamp'}, ...
        plottingOptions);
    
    
%         plottingOptions.orientation = 'vertical';
% %         plottingOptions.orientation = 'horizontal';
%         
%     %     plottingOptions.vertical_state_mode = 'stacked';
%         plottingOptions.vertical_state_mode = 'combined';
%         
%         plottingOptions.plot_variable = 'behavioral_epoch';
%     %     plottingOptions.plot_variable = 'behavioral_state';
%         
%         plottingOptions.x_axis = 'timestamp'; % Timestamp-appropriate relative bins
% %         plottingOptions.x_axis = 'index'; % Equal-sized bins
%         


    if strcmpi(plottingOptions.plot_variable, 'behavioral_epoch')
        % behavioral_epoch
        if ~isfield(active_processing.definitions.behavioral_epoch, 'classColors')
            active_processing.definitions.behavioral_epoch.classColors = [0.0, 0.5, 0.0
               0.2, 1.0, 0.2
               0.0, 0.2, 0.0];
        end
        
        state_names = active_processing.definitions.behavioral_epoch.classNames;
        
        if strcmpi(plottingOptions.x_axis, 'timestamp')
            states = [active_processing.behavioral_periods_table.epoch_start_seconds, ...
                active_processing.behavioral_periods_table.epoch_end_seconds, ...
                double(active_processing.behavioral_periods_table.behavioral_epoch)];
        else
%             linear_start_indicies = 0:1:(length(active_processing.behavioral_periods_table.behavioral_epoch)-1);
%             linear_end_indicies = linear_start_indicies + 1;
%                        
            linear_end_indicies = [1:length(active_processing.behavioral_periods_table.behavioral_epoch)];
            linear_start_indicies = linear_end_indicies - 1;
            states = [linear_start_indicies', ...
                linear_end_indicies', ...
                double(active_processing.behavioral_periods_table.behavioral_epoch)];
            
        end
        
        color_state = active_processing.definitions.behavioral_epoch.classColors;
           
    elseif strcmpi(plottingOptions.plot_variable, 'behavioral_state')
        % behavioral_state
        if ~isfield(active_processing.definitions.behavioral_state, 'classColors')
            active_processing.definitions.behavioral_state.classColors = [0.5, 0.5, 1.0
               0.7, 0.7, 1.0
               1.0, 0.7, 0.7
               1.0, 0.0, 0.0];
        end
        state_names = active_processing.definitions.behavioral_state.classNames;
        
        if strcmpi(plottingOptions.x_axis, 'timestamp')
            states = [active_processing.behavioral_periods_table.epoch_start_seconds, ...
                active_processing.behavioral_periods_table.epoch_end_seconds, ...
                double(active_processing.behavioral_periods_table.type)];
        else
            linear_end_indicies = [1:length(active_processing.behavioral_periods_table.behavioral_epoch)];
            linear_start_indicies = linear_end_indicies - 1;
                    
            states = [linear_start_indicies', ...
                linear_end_indicies', ...
                double(active_processing.behavioral_periods_table.type)];
        end
        color_state = active_processing.definitions.behavioral_state.classColors;
    else
        error('invalid plottingOptions.plot_variable');
    end
            
    t = 0:1:states(end, 2);
    
    for s_idx = 1:size(states,1)
        if strcmpi(plottingOptions.vertical_state_mode, 'stacked')
            curr_s_y = states(s_idx,3)-0.5;
        else
            % in combined mode they all have the same y-position.
            curr_s_y = 1.0-0.5;
        end
        
        % 
        if strcmpi(plottingOptions.orientation, 'vertical')
            % If vertically oriented, flip the x and y values
            temp.rect_pos = [curr_s_y, states(s_idx,1), 1, diff(states(s_idx,1:2))];
        else
            temp.rect_pos = [states(s_idx,1), curr_s_y, diff(states(s_idx,1:2)), 1];
        end
        
        rectangle('Position', temp.rect_pos,...
                    'LineStyle','none','facecolor',color_state(states(s_idx,3),:))
    end

    ax=gca;
    
    if strcmpi(plottingOptions.vertical_state_mode, 'stacked')
        state_stack_dim.Tick = 1:length(state_names);
        state_stack_dim.TickLabel = state_names;
        state_stack_dim.Lim = 0.5+[0,length(state_names)];
        
    else
        % combined mode:
        state_stack_dim.Tick = [];
        state_stack_dim.TickLabel = '';
        state_stack_dim.Lim = 0.5 + [0,1];
        
    end
     
    state_stack_dim.Dir='reverse';
    
    length_dim.Lim = t([1,end]);
    length_dim.Axis.Visible = 'off';

     if strcmpi(plottingOptions.orientation, 'vertical')
        % If vertically oriented, flip the x and y values
        ax.XTick = state_stack_dim.Tick;
        ax.XTickLabel = state_stack_dim.TickLabel;
        ax.XLim = state_stack_dim.Lim;
        ax.XDir = state_stack_dim.Dir;
        
        ax.YLim = length_dim.Lim;
        ax.YAxis.Visible = length_dim.Axis.Visible;
        
        
     else
        ax.YTick = state_stack_dim.Tick;
        ax.YTickLabel = state_stack_dim.TickLabel;
        ax.YLim = state_stack_dim.Lim;
        ax.YDir = state_stack_dim.Dir;
        
        ax.XLim = length_dim.Lim;
        ax.XAxis.Visible = length_dim.Axis.Visible;
     end
        
     

   
end

