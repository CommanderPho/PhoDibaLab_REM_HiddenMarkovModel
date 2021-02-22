function [ax] = fnPlotStateDiagram(active_processing, plottingOptions)
%fnPlotStateDiagram Plots a behavioral-state indicator as a function of time
%   Detailed explanation goes here


if ~exist('plottingOptions','var')
%     plottingOptions.vertical_state_mode = 'stacked';
    plottingOptions.vertical_state_mode = 'combined';
    
    
end


%     state_names={'wake','nrem','rem'}; 

    
    % Table:
%     state_changes = across_experiment_results{1, 1}.active_processing.behavioral_periods_table.epoch_start_seconds;
    
%     state_changes=[0;cumsum(rand(10,1)*10)];
%     states=[state_changes(1:end-1),state_changes(2:end),mod(cumsum(randi(2,10,1)),3)+1];


    % behavioral_epoch
%     state_names = across_experiment_results{1, 1}.active_processing.definitions.behavioral_epoch.classNames;

    % behavioral_state
    state_names = active_processing.definitions.behavioral_state.classNames;
    states = [active_processing.behavioral_periods_table.epoch_start_seconds, ...
        active_processing.behavioral_periods_table.epoch_end_seconds, ...
        double(active_processing.behavioral_periods_table.type)];
    color_state=[0.5, 0.5, 1.0
           0.7, 0.7, 1.0
           1.0, 0.7, 0.7
           1.0, 0.0, 0.0];
       
    
    t = 0:1:states(end, 2);
    
    for s_idx = 1:size(states,1)
        if strcmpi(plottingOptions.vertical_state_mode, 'stacked')
            curr_s_y = states(s_idx,3)-0.5;
        else
            % in combined mode they all have the same y-position.
            curr_s_y = 1.0-0.5;
        end
        rectangle('Position',[states(s_idx,1), curr_s_y, diff(states(s_idx,1:2)),1],...
                    'LineStyle','none','facecolor',color_state(states(s_idx,3),:))
    end

    ax=gca;
    
    if strcmpi(plottingOptions.vertical_state_mode, 'stacked')
        ax.YTick=1:length(state_names);
        ax.YTickLabel=state_names;
        ax.YLim=0.5+[0,length(state_names)];
    else
        % combined mode:
        ax.YLim=0.5+[0,1];
        ax.YTickLabel='state';
        
    end
        
    ax.YDir='reverse';
    
  
    ax.XLim=t([1,end]);
    ax.XAxis.Visible = 'off';

   
end

