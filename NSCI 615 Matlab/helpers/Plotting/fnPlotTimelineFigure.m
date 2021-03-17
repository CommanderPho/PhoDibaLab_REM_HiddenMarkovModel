
function [ax, outputs] = fnPlotTimelineFigure(timeline, plottingOptions)
%fnPlotTimelineFigure Plots a series of blocks that represent a fixed-duration stage of an experiment as a function of time

    % timeline: a struct containing the fields:

    if ~exist('timeline','var')
        error('timeline variable must be specified!')
    %      timeline = struct();
    end
    % timeline = fnAddDefaultOptionalArgs({'titles', 'title_options', 'durations', 'colors'}, ...
    %     {{'CS-on','Reward','Intertrial Interval'}, ...
    %     [90, 90, 0], ...
    %     [seconds(8), seconds(1), seconds(60)], ...
    %     [0.0, 0.5, 0.0; 0.9, 0.1, 0.1; 0.0, 0.2, 0.0]}, ...
    %     timeline);


    if ~exist('plottingOptions','var')
         plottingOptions = struct();
    end
    plottingOptions = fnAddDefaultOptionalArgs({'orientation', 'vertical_state_mode', 'plot_variable', 'x_axis'}, ...
        {'vertical', 'combined', 'behavioral_state', 'timestamp'}, ...
        plottingOptions);


    outputs.positions.totalDuration = sum(timeline.durations, 'all');
    outputs.positions.proportionalDurations = timeline.durations ./ outputs.positions.totalDuration;

    numBlocks = length(outputs.positions.proportionalDurations);

    % Cumulative sum start timesteps
    outputs.positions.end_timesteps = cumsum(timeline.durations);
    outputs.positions.start_timesteps = outputs.positions.end_timesteps - timeline.durations;

    outputs.positions.state_y_offsets = ones([numBlocks 1]) - 0.5;
    states = [seconds(outputs.positions.start_timesteps)', ...
        seconds(outputs.positions.end_timesteps)'];


    outputs.positions.t = 0:1:seconds(outputs.positions.totalDuration);

    %% Pre-allocation of output objects:
    
    
    outputs.Rectangles.pos = zeros([numBlocks 4]);
    outputs.Rectangles.handles = gobjects([numBlocks 1]);
%     outputs.Labels.handles = gobjects([numBlocks 1]);  
            
    for s_idx = 1:numBlocks

    %     curr_s_y = 1.0-0.5;
        curr_s_y = outputs.positions.state_y_offsets(s_idx);

        % 
        if strcmpi(plottingOptions.orientation, 'vertical')
            % If vertically oriented, flip the x and y values
            outputs.Rectangles.pos(s_idx,:) = [curr_s_y, states(s_idx,1), 1, diff(states(s_idx,1:2))];
        else
            outputs.Rectangles.pos(s_idx,:) = [states(s_idx,1), curr_s_y, diff(states(s_idx,1:2)), 1];
        end

        outputs.Rectangles.handles(s_idx) = rectangle('Position', outputs.Rectangles.pos(s_idx,:),...
                    'LineStyle','none','facecolor', timeline.colors(s_idx,:));
    end

    %% Cleanup Axes:
    ax = gca;
    % combined mode:
    state_stack_dim.Tick = [];
    state_stack_dim.TickLabel = '';
    state_stack_dim.Lim = 0.5 + [0,1];
    % state_stack_dim.Lim = [0,1];
    state_stack_dim.Dir='reverse';

    length_dim.Lim = outputs.positions.t([1,end]);
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

    %% Add Labels/Text:
    midpoint_timestep_offsets = (outputs.positions.end_timesteps - outputs.positions.start_timesteps) ./ 2.0;
    outputs.positions.midpoint_timesteps = outputs.positions.start_timesteps + midpoint_timestep_offsets;

    y_offset = (state_stack_dim.Lim(2) - state_stack_dim.Lim(1)) / 2;
    y_midpoint = min(state_stack_dim.Lim) + y_offset;
    y_midpoints = repelem(y_midpoint, numBlocks);

    outputs.Labels.handles = text(seconds(outputs.positions.midpoint_timesteps), y_midpoints, timeline.titles, ...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'FontSize',14,'Color','white');


    for s_idx = 1:numBlocks
    %     h = text(X(:,10),X(:,i),labels,'VerticalAlignment','top','HorizontalAlignment','right');
    %     h = annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', titles{s_idx});
        % Rotate text:
        set(outputs.Labels.handles(s_idx), 'Rotation', timeline.title_options.Rotations(s_idx));
    end

end
