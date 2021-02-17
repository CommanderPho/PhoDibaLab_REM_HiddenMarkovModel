% PhoDibaTest_PositionalAnalysis.m
% Peform analysis of animal position and state transitions in wake on the track

% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded
% Stage 3 of the processing pipeline

addpath(genpath('../../helpers'));
addpath(genpath('../../libraries/buzcode/'));

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

if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'results_array', 'general_results');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('PhoDibaTest_PositionalAnalysis ready to process!\n');

% We have active_processing.position_table and active_processing.speed_table

positionalAnalysis.track_epoch.begin = active_processing.behavioral_epochs.start_seconds(2);
positionalAnalysis.track_epoch.end = active_processing.behavioral_epochs.end_seconds(2);


positionalAnalysis.plotting.filtering.use_focus_time = true; % If true, highlights the time range preceeding the focus time.

positionalAnalysis.plotting.filtering.previous_path_duration = 30 / active_binning_resolution; % Show 30 timesteps prior to the filter time
positionalAnalysis.plotting.filtering.focus_time = positionalAnalysis.track_epoch.begin + positionalAnalysis.plotting.filtering.previous_path_duration; % The end time to show
positionalAnalysis.plotting.filtering.initial_focus_time = positionalAnalysis.plotting.filtering.focus_time - positionalAnalysis.plotting.filtering.previous_path_duration;


% Filter the timestamps where the animal was on the track:
% track_indicies = find((track_epoch.begin < active_processing.position_table.timestamp) & (active_processing.position_table.timestamp < track_epoch.end));
track_indicies = ((positionalAnalysis.track_epoch.begin < active_processing.position_table.timestamp) & (active_processing.position_table.timestamp < positionalAnalysis.track_epoch.end));

[positionalAnalysis.plotting.bounds.x(1), positionalAnalysis.plotting.bounds.x(2)] = bounds(active_processing.position_table.x);
[positionalAnalysis.plotting.bounds.y(1), positionalAnalysis.plotting.bounds.y(2)] = bounds(active_processing.position_table.y);




if positionalAnalysis.plotting.filtering.use_focus_time
    % Build the focused indicies 
    positionalAnalysis.plotting.filtering.focused_indicies = (positionalAnalysis.plotting.filtering.initial_focus_time <= active_processing.position_table.timestamp) & (active_processing.position_table.timestamp <= positionalAnalysis.plotting.filtering.focus_time);
    % Compute the logical-AND of the two sets of indicies by performing element-wise multiplication
    
    %% Show only mode:
%     positionalAnalysis.plotting.active_indicies = logical(track_indicies .* positionalAnalysis.plotting.filtering.focused_indicies); % 1x971952
    
    % Opacity mode:
    positionalAnalysis.plotting.filtering.focused_indicies = logical(track_indicies .* positionalAnalysis.plotting.filtering.focused_indicies); % 1x971952 
    positionalAnalysis.plotting.filtering.num_active_points = sum(positionalAnalysis.plotting.filtering.focused_indicies, 'all');
    positionalAnalysis.plotting.active_indicies = track_indicies; % 1x971952
    
else
    positionalAnalysis.plotting.active_indicies = track_indicies; % 1x971952
    positionalAnalysis.plotting.filtering.num_active_points = sum(positionalAnalysis.plotting.active_indicies, 'all');
end


positionalAnalysis.plotting.fig = figure(19);
clf(positionalAnalysis.plotting.fig);
% h = plot(active_processing.position_table.timestamp', [active_processing.position_table.x', active_processing.position_table.y']);
positionalAnalysis.plotting.num_location_points = length(active_processing.position_table.timestamp(positionalAnalysis.plotting.active_indicies)); % 971952
% h = plot(active_processing.position_table.x', active_processing.position_table.y','r', 'LineWidth',5);

% % plot(x, y, 'bo-', 'LineWidth', 2, 'MarkerSize', 12);
% % modified jet-colormap
% cd = [uint8(jet(num_location_points)*255) uint8(ones(num_location_points,1))].';
% drawnow
% set(h.Edge, 'ColorBinding','interpolated', 'ColorData',cd)



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
hold on


% xline(track_epoch.begin,'-.','TRACK Epoch');
% xline(track_epoch.end,'-.','End TRACK');

% xlim([track_epoch.begin, track_epoch.end]);



% plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;
% 
% 
% current_binning_index = 1;
% active_binning_resolution = processing_config.step_sizes{current_binning_index};
% temp.curr_timestamps = timesteps_array{current_binning_index};
% temp.curr_processed = active_processing.processed_array{current_binning_index};
% active_results = results_array{current_binning_index};
% 
% fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);
% 
% 
% % Binned Spike rates per time:
% [PhoDibaTest_PositionalAnalysis_temp.activeMatrix] = fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts);  % 35351x126 double
% 
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds = [11512.6074066973, 12000];
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds ./ active_binning_resolution;
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = [floor(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1)), ...
%     ceil(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2))];
% 
% if processing_config.showOnlyAlwaysStableCells
%     isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
%     numAlwaysStableCells = sum(isAlwaysStable, 'all');
%     PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, isAlwaysStable);
% end
% 
% PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix'; % should be units x time
% % Extract only the portion specified as the training portion
% PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1):PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2));
% 
% 
% %% Main Procedure:
% PhoDibaTest_PositionalAnalysis_config.K = 5; % K: Number of factors

% grab_fig.figure_indicies = [1 15 19];
% grab_fig.num_figs = length(grab_fig.figure_indicies);
% grab_fig.figure_handles = gobjects();
% for fig_i = 1:grab_fig.num_figs
%    grab_fig.figure_handles(fig_i) = figure(grab_fig.figure_indicies(fig_i)); % Make the current figure active
% end




fprintf('PhoDibaTest_PositionalAnalysis complete!\n');
