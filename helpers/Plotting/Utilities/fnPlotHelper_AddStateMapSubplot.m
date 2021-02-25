function [state_ax, epoch_ax] = fnPlotHelper_AddStateMapSubplot(active_processing, curr_axis, plottingOptions)
%fnPlotHelper_AddStateMapSubplot adds a state map/partition subplot to indicate the state as a function of time or period index.
    % By default adds them to the right side of the plot.

    
% %     common_statemapPlottingOptions.x_axis = 'timestamp';% Timestamp-appropriate relative bins (default)
%     common_statemapPlottingOptions.x_axis = 'index'; % Equal-sized bins
%     common_statemapPlottingOptions.orientation = 'vertical';
%     common_statemapPlottingOptions.vertical_state_mode = 'combined';
    
    
     if ~exist('plottingOptions','var')
         plottingOptions = struct();
     end
    
     common_statemapPlottingOptions = fnAddDefaultOptionalArgs({'x_axis', 'orientation', 'vertical_state_mode', 'include_state', 'include_epoch'}, ...
        {'index', 'vertical', 'combined', true, true}, ...
        plottingOptions);
    
    
    
    % Resize current plots:
    % temp.curr_heatmap_pos = curr_axis.Position;

    temp.updated_heatmap_pos = curr_axis.Position;
    temp.updated_heatmap_pos(1) = 0.1; % Move left edge (x) to 0.1 in relative coords.
%     temp.updated_heatmap_pos(3) = temp.curr_heatmap_pos(3) * 0.9; % Set to 90% of original width
    curr_axis.Position =  temp.updated_heatmap_pos;
    
    
    %% Add the behavioral period map:
    if common_statemapPlottingOptions.include_state
		state_statemapPlottingOptions.orientation = common_statemapPlottingOptions.orientation;
		state_statemapPlottingOptions.vertical_state_mode = common_statemapPlottingOptions.vertical_state_mode;
		state_statemapPlottingOptions.plot_variable = 'behavioral_state';
		state_statemapPlottingOptions.x_axis = common_statemapPlottingOptions.x_axis;
		
		temp.statemap_pos = temp.updated_heatmap_pos;
		temp.statemap_pos(1) = 0.9;
		temp.statemap_pos(1) = temp.updated_heatmap_pos(1) + temp.updated_heatmap_pos(3);
		temp.statemap_pos(3) = 0.05;
		
		subplot('Position', temp.statemap_pos);
		[state_ax] = fnPlotStateDiagram(active_processing, state_statemapPlottingOptions);
	else
		state_ax = [];
	end
    
    %% Plot Epoch Position map
	 if common_statemapPlottingOptions.include_epoch
		epoch_statemapPlottingOptions.orientation = common_statemapPlottingOptions.orientation;
		epoch_statemapPlottingOptions.vertical_state_mode = common_statemapPlottingOptions.vertical_state_mode;
		epoch_statemapPlottingOptions.plot_variable = 'behavioral_epoch';
		epoch_statemapPlottingOptions.x_axis = common_statemapPlottingOptions.x_axis;
		
		temp.epoch_statemap_pos = temp.statemap_pos;
		temp.epoch_statemap_pos(1) = temp.statemap_pos(1) + temp.statemap_pos(3);
		temp.epoch_statemap_pos(3) = 0.025;
		
		subplot('Position', temp.epoch_statemap_pos);
		[epoch_ax] = fnPlotStateDiagram(active_processing, epoch_statemapPlottingOptions);
	else
		epoch_ax = [];
	end

end

