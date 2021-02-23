function [state_ax, epoch_ax] = fnPlotAddStateMapSubplot(active_processing, curr_axis)
%FNPLOTADDSTATEMAPSUBPLOT adds a state map/partition subplot to indicate the state as a function of time or period index.
    % By default adds them to the right side of the plot.

    % Resize current plots:
    temp.curr_heatmap_pos = curr_axis.Position;
    temp.updated_heatmap_pos = curr_axis.Position;
    temp.updated_heatmap_pos(1) = 0.1; 
%     temp.updated_heatmap_pos(3) = temp.curr_heatmap_pos(3) * 0.9; % Set to 90% of original width
    curr_axis.Position =  temp.updated_heatmap_pos;
    
    
    %% Add the behavioral period map:
    state_statemapPlottingOptions.orientation = 'vertical';
    state_statemapPlottingOptions.vertical_state_mode = 'combined';
    state_statemapPlottingOptions.plot_variable = 'behavioral_state';
    
    temp.statemap_pos = temp.updated_heatmap_pos;
    temp.statemap_pos(1) = 0.9;
    temp.statemap_pos(1) = temp.updated_heatmap_pos(1) + temp.updated_heatmap_pos(3);
    temp.statemap_pos(3) = 0.05;
    
    subplot('Position', temp.statemap_pos);
    [state_ax] = fnPlotStateDiagram(active_processing, state_statemapPlottingOptions);

    
    %% Plot Epoch Position map
    epoch_statemapPlottingOptions.orientation = 'vertical';
    epoch_statemapPlottingOptions.vertical_state_mode = 'combined';
    epoch_statemapPlottingOptions.plot_variable = 'behavioral_epoch';
    temp.epoch_statemap_pos = temp.statemap_pos;
    temp.epoch_statemap_pos(1) = temp.statemap_pos(1) + temp.statemap_pos(3);
    temp.epoch_statemap_pos(3) = 0.025;
    
    subplot('Position', temp.epoch_statemap_pos);
    [epoch_ax] = fnPlotStateDiagram(active_processing, epoch_statemapPlottingOptions);
end

