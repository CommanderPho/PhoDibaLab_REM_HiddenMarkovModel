
%% General Plot Config Options:
plottingOptions.orientation = 'horizontal';
plottingOptions.plot_variable = 'behavioral_state';
plottingOptions.vertical_state_mode = 'combined';
plottingOptions.x_axis = 'timestamp'; % Timestamp-appropriate relative bins


plottingOptions.lineFormat.Color = 'white';
plottingOptions.lineFormat.LineWidth = 4.0;

plottingOptions.barbLineFormat.Color = 'cyan';
plottingOptions.barbLineFormat.LineWidth = 2.0;

%% Trial Timeline Figure:
trialTimeline.titles = {'CS-on','Reward','Intertrial Interval'};
trialTimeline.title_options.Rotations = [90, 90, 0];

trialTimeline.durations = [seconds(8), seconds(1), seconds(60)];
% colors = [0.0, 0.5, 0.0
%                0.2, 1.0, 0.2 % Lime Green
%                0.0, 0.2, 0.0];
trialTimeline.colors = [0.0, 0.5, 0.0 % Lighter Dark Green
               0.9, 0.1, 0.1 % Red
               0.4, 0.4, 0.4]; % Dark Green

currIdx = 1;
figure(currIdx)
clf;
%% Build Time Plot:
[ax(currIdx), outputs] = fnPlotTimelineFigure(trialTimeline, plottingOptions);

centerPoint = [0.75 0.75];
% plottingOptions.majorAxisStartEnd = [0 2];
% Get the line start/stop positions from the outputs of the timeline plotting function:
plottingOptions.majorAxisStartEnd = [seconds(outputs.positions.start_timesteps)', seconds(outputs.positions.end_timesteps)'];
[plotHandle] = fnPlotHelper_DrawLengthIndicator(centerPoint, plottingOptions);




% %% Trial Timeline Figure:
% 
% numTrialsPerSession = 25;
% sessionTrialRepeatDims = [1 numTrialsPerSession];
% sessionTimeline.titles = repmat({'', '', ''}, sessionTrialRepeatDims);
% sessionTimeline.title_options.Rotations = repmat([0, 0, 0], sessionTrialRepeatDims);
% 
% sessionTimeline.durations = repmat(trialTimeline.durations, sessionTrialRepeatDims);
% sessionTimeline.colors = repmat(trialTimeline.colors, [numTrialsPerSession 1]);
% 
% currIdx = 2;
% figure(currIdx)
% clf;
% %% Build Time Plot:
% [ax(currIdx)] = fnPlotTimelineFigure(sessionTimeline, plottingOptions);
% % 
