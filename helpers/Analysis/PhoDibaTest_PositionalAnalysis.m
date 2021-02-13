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



plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;


current_binning_index = 1;
active_binning_resolution = processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = timesteps_array{current_binning_index};
temp.curr_processed = active_processing.processed_array{current_binning_index};
active_results = results_array{current_binning_index};

fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);


% Binned Spike rates per time:
[PhoDibaTest_seqNMF_temp.activeMatrix] = fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts);  % 35351x126 double

PhoDibaTest_seqNMF_config.training_subset_start_stop_seconds = [11512.6074066973, 12000];
PhoDibaTest_seqNMF_config.training_subset_start_stop_bins = PhoDibaTest_seqNMF_config.training_subset_start_stop_seconds ./ active_binning_resolution;
PhoDibaTest_seqNMF_config.training_subset_start_stop_bins = [floor(PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(1)), ...
    ceil(PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(2))];

if processing_config.showOnlyAlwaysStableCells
    isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix(:, isAlwaysStable);
end

PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix'; % should be units x time
% Extract only the portion specified as the training portion
PhoDibaTest_seqNMF_temp.activeMatrix = PhoDibaTest_seqNMF_temp.activeMatrix(:, PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(1):PhoDibaTest_seqNMF_config.training_subset_start_stop_bins(2));


%% Main Procedure:
PhoDibaTest_PositionalAnalysis_config.K = 5; % K: Number of factors



fprintf('PhoDibaTest_PositionalAnalysis complete!\n');
