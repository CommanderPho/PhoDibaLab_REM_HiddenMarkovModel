% PhoDibaTest_PositionalAnalysis.m
% Peform analysis of animal position and state transitions in wake on the track

% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded
% Stage 3 of the processing pipeline

% Makes use of active_processing.position_table, active_processing.behavioral_epochs.start_seconds

addpath(genpath('../../helpers'));
addpath(genpath('../../libraries/buzcode/'));

enable_plotting = false;


clear temp

if ~exist('data_config','var')
    Config;
end


if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading data from %s...\n', data_config.output.intermediate_file_paths{2});
    load(data_config.output.intermediate_file_paths{2}, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end

% if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
%     fprintf('loading results from %s...\n', data_config.output.results_file_path);
%     load(data_config.output.results_file_path, 'results_array', 'general_results');
%     fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
% else
%     fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
% end

%% Begin:
fprintf('PhoDibaTest_PositionalAnalysis ready to process!\n');

%% Binning Options:
active_expt_index = 1;


% We have active_processing.position_table and active_processing.speed_table

%% Computer displacements per timestep, velocities per timestep, etc.
positionalAnalysis.displacement.dt = [diff(active_processing.position_table.timestamp)];
positionalAnalysis.displacement.dx = [diff(active_processing.position_table.x)];
positionalAnalysis.displacement.dy = [diff(active_processing.position_table.y)];

positionalAnalysis.displacement.speeds = sqrt(((positionalAnalysis.displacement.dx .^ 2)) + (positionalAnalysis.displacement.dy .^ 2)) ./ positionalAnalysis.displacement.dt;

% Add the initial zero (indicating no change on the initial point) to keep length
positionalAnalysis.displacement.dt = [0; diff(active_processing.position_table.timestamp)];
positionalAnalysis.displacement.dx = [0; diff(active_processing.position_table.x)];
positionalAnalysis.displacement.dy = [0; diff(active_processing.position_table.y)];
positionalAnalysis.displacement.speeds = [0; positionalAnalysis.displacement.speeds];


%% Find the start and end times the animal was put on the track, which is the period in which positions are relevant.
positionalAnalysis.track_epoch.start_end_delay = 0.5; % positionalAnalysis.track_epoch.start_end_delay: a buffer at the start and end of the epoch to account for the recording not being setup quite yet.
positionalAnalysis.track_epoch.begin = active_processing.behavioral_epochs.start_seconds(2) + positionalAnalysis.track_epoch.start_end_delay;
positionalAnalysis.track_epoch.end = active_processing.behavioral_epochs.end_seconds(2) - positionalAnalysis.track_epoch.start_end_delay;

%% FEATURE: 'use_focus_time': if true, displays a only a range of the data prior to the focus_time, leaving a "recency-trial" of places the animal has previously been that fades away as the focus_time advances.
positionalAnalysis.plotting.filtering.use_focus_time = false; % If true, highlights the time range preceeding the focus time.
% active_step_sizes


% Filter the timestamps where the animal was on the track:
positionalAnalysis.track_indicies = ((positionalAnalysis.track_epoch.begin < active_processing.position_table.timestamp) & (active_processing.position_table.timestamp < positionalAnalysis.track_epoch.end));

positionalAnalysis.track_position.t = active_processing.position_table.timestamp(positionalAnalysis.track_indicies);
positionalAnalysis.track_position.x = active_processing.position_table.x(positionalAnalysis.track_indicies);
positionalAnalysis.track_position.y = active_processing.position_table.y(positionalAnalysis.track_indicies);
positionalAnalysis.track_position.speeds = positionalAnalysis.displacement.speeds(positionalAnalysis.track_indicies);

% Filter the displacents in case we want those as well:
positionalAnalysis.displacement.dt = positionalAnalysis.displacement.dt(positionalAnalysis.track_indicies);
positionalAnalysis.displacement.dx = positionalAnalysis.displacement.dx(positionalAnalysis.track_indicies);
positionalAnalysis.displacement.dy = positionalAnalysis.displacement.dy(positionalAnalysis.track_indicies);

positionalAnalysis.track_position.xyv = [positionalAnalysis.track_position.x, positionalAnalysis.track_position.y, positionalAnalysis.track_position.speeds]; % a matrix with 3 columns corresponding to the x-pos, y-pos, and speed
% To have complete info, we need:
% positionalAnalysis.track_position.t
% positionalAnalysis.track_position.xyv

% Computer the new bounds restricted to the valid (track) positions:
[positionalAnalysis.plotting.bounds.x(1), positionalAnalysis.plotting.bounds.x(2)] = bounds(positionalAnalysis.track_position.x);
[positionalAnalysis.plotting.bounds.y(1), positionalAnalysis.plotting.bounds.y(2)] = bounds(positionalAnalysis.track_position.y);


if (enable_plotting)
    positionalAnalysis.plotting.speed_fig = figure(18);
    clf(positionalAnalysis.plotting.speed_fig);
    title('Animal Speed vs. time')
    % plot(positionalAnalysis.displacement.speeds(positionalAnalysis.track_indicies));
    stackedaxes(positionalAnalysis.plotting.speed_fig, positionalAnalysis.track_position.t, positionalAnalysis.track_position.xyv, ...
        'DisplayLabels', {'x','y','speed'})
    % xline([positionalAnalysis.track_epoch.begin positionalAnalysis.track_epoch.end],'--r',{'<track>', '</track>'})
    xlabel('t')
end

if positionalAnalysis.plotting.filtering.use_focus_time
    current_binning_index = 2;
    active_binning_resolution = across_experiment_results{active_expt_index}.processing_config.step_sizes{current_binning_index}; % This isn't right, it should instead be the dt

    positionalAnalysis.plotting.filtering.previous_path_duration = 100 / active_binning_resolution; % Show 100 timesteps prior to the filter time
    positionalAnalysis.plotting.filtering.focus_time = positionalAnalysis.track_epoch.begin + positionalAnalysis.plotting.filtering.previous_path_duration; % The end time to show
    positionalAnalysis.plotting.filtering.initial_focus_time = positionalAnalysis.plotting.filtering.focus_time - positionalAnalysis.plotting.filtering.previous_path_duration;

    % Build the focused indicies 
    positionalAnalysis.plotting.filtering.focused_indicies = (positionalAnalysis.plotting.filtering.initial_focus_time <= active_processing.position_table.timestamp) & (active_processing.position_table.timestamp <= positionalAnalysis.plotting.filtering.focus_time);
    % Compute the logical-AND of the two sets of indicies by performing element-wise multiplication
    
    %% Show subset only mode:
%     positionalAnalysis.plotting.active_indicies = logical(track_indicies .* positionalAnalysis.plotting.filtering.focused_indicies); % 1x971952
    
    % Opacity mode:
    positionalAnalysis.plotting.filtering.focused_indicies = logical(positionalAnalysis.track_indicies .* positionalAnalysis.plotting.filtering.focused_indicies); % 1x971952 
    positionalAnalysis.plotting.filtering.num_active_points = sum(positionalAnalysis.plotting.filtering.focused_indicies, 'all');
    positionalAnalysis.plotting.active_indicies = positionalAnalysis.track_indicies; % 1x971952
    
else
    positionalAnalysis.plotting.active_indicies = positionalAnalysis.track_indicies; % 1x971952
    positionalAnalysis.plotting.filtering.num_active_points = sum(positionalAnalysis.plotting.active_indicies, 'all');
end

if (enable_plotting)
    % Perform the plotting of the position heatmap:
    positionalAnalysis.plotting.fig = figure(19);
    clf(positionalAnalysis.plotting.fig);
end

% h = plot(active_processing.position_table.timestamp', [active_processing.position_table.x', active_processing.position_table.y']);
positionalAnalysis.plotting.num_location_points = length(active_processing.position_table.timestamp(positionalAnalysis.plotting.active_indicies)); % 971952

if positionalAnalysis.plotting.filtering.use_focus_time
    % Opacity mode:    
    positionalAnalysis.plotting.alphaMap = zeros([1 positionalAnalysis.plotting.num_location_points]);
    positionalAnalysis.plotting.alphaMap(positionalAnalysis.plotting.filtering.focused_indicies) = linspace(0.1, 0.9, positionalAnalysis.plotting.filtering.num_active_points);
    positionalAnalysis.plotting.colorMap = zeros([1 positionalAnalysis.plotting.num_location_points]);    
    positionalAnalysis.plotting.colorMap(positionalAnalysis.plotting.filtering.focused_indicies) = linspace(1, 10, positionalAnalysis.plotting.filtering.num_active_points); % Define colormap relative to active points  
else
    positionalAnalysis.plotting.alphaMap = linspace(0.1,0.8,positionalAnalysis.plotting.num_location_points);
    positionalAnalysis.plotting.colorMap = linspace(1, 10, positionalAnalysis.plotting.num_location_points);
    
end

if (enable_plotting)
    positionalAnalysis.plotting.h = scatter(active_processing.position_table.x(positionalAnalysis.plotting.active_indicies)', active_processing.position_table.y(positionalAnalysis.plotting.active_indicies)',...
        10,...
        positionalAnalysis.plotting.colorMap,...
        'filled'); % to color based on z-axis values
    positionalAnalysis.plotting.h.AlphaData = positionalAnalysis.plotting.alphaMap;
    positionalAnalysis.plotting.h.MarkerFaceAlpha = 'flat';
    axis equal
    grid on;
    % xlabel('timestamp')
    xlabel('x-position')
    ylabel('y-position')
    title('Position and Speed vs. Time');
    xlim(positionalAnalysis.plotting.bounds.x);
    ylim(positionalAnalysis.plotting.bounds.y);
    % hold on
    % xline(track_epoch.begin,'-.','TRACK Epoch');
    % xline(track_epoch.end,'-.','End TRACK');

    % xlim([track_epoch.begin, track_epoch.end]);
end


% %% Main Procedure:
% PhoDibaTest_PositionalAnalysis_config.K = 5; % K: Number of factors

% grab_fig.figure_indicies = [1 15 19];
% grab_fig.num_figs = length(grab_fig.figure_indicies);
% grab_fig.figure_handles = gobjects();
% for fig_i = 1:grab_fig.num_figs
%    grab_fig.figure_handles(fig_i) = figure(grab_fig.figure_indicies(fig_i)); % Make the current figure active
% end


% Save out positionalAnalysis data for Python:
export_root_path = '/Users/pho/repo/Python Projects/PhoNeuronGillespie2021CodeRepo/PhoMatlabDataScripting/ExportedData';

if ~exist('active_experiment_name','var')
    active_experiment_name = active_experiment_names{active_expt_index};
end
active_experiment_export_root_path = fullfile(export_root_path, active_experiment_name);
mkdir(active_experiment_export_root_path);

if (enable_plotting)
    % Remove figure handles before saving
    positionalAnalysis = rmfield(positionalAnalysis, 'plotting.fig', 'plotting.h');
end


fprintf('Saving positional analysis data to %s...\n', fullfile(active_experiment_export_root_path, 'positionAnalysis.mat'));
save(fullfile(active_experiment_export_root_path, 'positionAnalysis.mat'), 'positionalAnalysis');
fprintf('done!\n');

fprintf('PhoDibaTest_PositionalAnalysis complete!\n');
