% Requires the datastructures from "PhoDibaProcess_Stage1.m" to be loaded
% Stage 2 of the processing pipeline

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

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
    load(data_config.output.results_file_path, 'results_array');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('PhoDibaProcess_Stage2 ready to process!\n');


%% Find units that are stable across all sessions:
active_processing.spikes.stability_count = sum(active_processing.spikes.isStable, 2);
isAlwaysStable = (active_processing.spikes.stability_count == 3);
numAlwaysStableCells = sum(isAlwaysStable, 'all');




%% Compute ISIs:
active_processing.spikes.ISIs = cellfun((@(spikes_timestamps) diff(spikes_timestamps)), ...
 active_processing.spikes.time, 'UniformOutput', false);

active_processing.spikes.meanISI = cellfun((@(spikes_ISIs) mean(spikes_ISIs)), ...
 active_processing.spikes.ISIs, 'UniformOutput', false);

active_processing.spikes.ISIVariance = cellfun((@(spikes_ISIs) var(spikes_ISIs)), ...
 active_processing.spikes.ISIs, 'UniformOutput', false);

% active_processing.spikes.ISIs = cellfun((@(spikes_timestamps) diff(spikes_timestamps)), ...
%  active_processing.spikes.time, 'UniformOutput', false);


% Cluster ISIs:
% kmeanscluster(


% Get the duration of the epochs {'pre_sleep','track','post_sleep'}
% active_processing.curr_activity_table.behavioral_epochs.duration

% Get the duration of the states 
temp.GroupedByState.groups = findgroups(active_processing.curr_activity_table.type);
temp.GroupedByState.durations = splitapply(@sum, active_processing.curr_activity_table.duration, temp.GroupedByState.groups);

% active_processing.spikes.behavioral_duration_indicies

%% Across all cells:
temp.flattened_across_all_units.spike_timestamp = cat(2, active_processing.spikes.time{:});
temp.flattened_across_all_units.spike_state_index = cat(1, active_processing.spikes.behavioral_duration_indicies{:});
temp.flattened_across_all_units.spike_state = cat(1, active_processing.spikes.behavioral_states{:});
temp.flattened_across_all_units.spike_epoch = cat(1, active_processing.spikes.behavioral_epoch{:});


%% For each behavioral state period, every unit fires a given number of spikes.
%	The average firing rate for each unit within that period is given by this number of spikes divided by the duration of that period.
    
num_of_behavioral_state_periods = height(active_processing.curr_activity_table);
temp.edges = 1:num_of_behavioral_state_periods;

temp.per_behavioral_state_period.num_spikes_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double
temp.per_behavioral_state_period.spike_rate_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double


for unit_index = 1:num_of_electrodes

    % temp.edges = unique(active_processing.spikes.behavioral_duration_indicies{unit_index});
    temp.counts = histc(active_processing.spikes.behavioral_duration_indicies{unit_index}(:), temp.edges);
    temp.per_behavioral_state_period.num_spikes_per_unit(:, unit_index) = temp.counts;
    temp.per_behavioral_state_period.spike_rate_per_unit(:, unit_index) = temp.counts ./ active_processing.curr_activity_table.duration;
    % active_processing.spikes.behavioral_duration_indicies{unit_index}

    % curr_activity_table.duration

end


phoPlotSpikeRateHeatmap;


% 
% 
% % Take the above figure and collapse the y-axis by behavioral state.
% num_of_behavioral_state_types = length(active_processing.behavioral_state_names);
% temp.per_behavioral_state_type.num_spikes_per_unit = zeros([num_of_behavioral_state_types num_of_electrodes]); %% num_results: 4x126 double
% temp.per_behavioral_state_type.spike_rate_per_unit = zeros([num_of_behavioral_state_types num_of_electrodes]); %% num_results: 4x126 double
% 
% 
% for i = 1:length(active_processing.behavioral_state_names)
%     temp.curr_state_name =  active_processing.behavioral_state_names{i};
%     temp.per_behavioral_state_type.num_spikes_per_unit(i,:) = mean(temp.per_behavioral_state_period.num_spikes_per_unit,1);
% 
% end
% 
% 
% phoPlotInteractiveSpikesScrollableRaster;



% for current_binning_index = 1:length(processing_config.step_sizes)
% 	temp.curr_timestamps = timesteps_array{current_binning_index};
% 	temp.curr_processed = active_processing.processed_array{current_binning_index};
% 	temp.active_results = results_array{current_binning_index};
% 
% 
% 
% 
% 
% 	% results_array{current_binning_index} = active_results;
% 	
% end % end for processing_config.step_sizes loop

% fprintf('writing out results to %s...\n', data_config.output.results_file_path);
% save(data_config.output.results_file_path, 'results_array');
% fprintf('done.\n');

fprintf('PhoDibaProcess_Stage2 complete!\n');



