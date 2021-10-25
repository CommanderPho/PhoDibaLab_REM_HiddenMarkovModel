% PhoDibaConvert_SpikesAnalysis.m

% Makes use of active_processing.position_table, active_processing.behavioral_epochs.start_seconds

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

% if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
%     fprintf('loading results from %s...\n', data_config.output.results_file_path);
%     load(data_config.output.results_file_path, 'results_array', 'general_results');
%     fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
% else
%     fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
% end

%% Begin:
fprintf('PhoDibaConvert_SpikesAnalysis ready to process!\n');

%% Binning Options:
active_expt_index = 1;

plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;

current_binning_index = 1;
active_binning_resolution = processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = timesteps_array{current_binning_index};
temp.curr_processed = active_processing.processed_array{current_binning_index};
% active_results = results_array{current_binning_index};

fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);

% [ones(size(PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_i}')), PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_i}']
% cellfun(@(cell_i) [ones(size(PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_i}')), PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_i}'], 

% cellfun(@(cell_arr) [ones(size(cell_arr')), cell_arr'], PhoDibaTest_PositionalAnalysis_temp.active_spikes,'UniformOutput',false)
% num2cell([1:numAlwaysStableCells]')
% {1:numAlwaysStableCells}

% arrayfun(
% active_processing.spikes % Table

% Only get the spike data during the track period:
% active_processing.processed_array{2, 1}.by_epoch.track.spike_data  

% [out_times out_groups out_spikesmat] = spikes2sorted(active_processing.spikes.time);
% [spikemat] = bz_SpktToSpkmat(active_processing.spikes, varargin)
% [firingMaps] = bz_firingMap1D(varargin)
% [ fields ] = bz_getPlaceFields(varargin)
% [ fields ] = bz_getPlaceFields1D(varargin)


% [positionDecodingBayesian] = bz_positionDecodingBayesian(varargin)
% [Pr, prMax] = placeBayes(Cr, rateMap, binLength)
% [positionDecodingGLM] = bz_positionDecodingGLM(varargin)


% Binned Spike rates per time:
[PhoDibaTest_PositionalAnalysis_temp.activeMatrix] = fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts);  % 35351x126 double

PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds = [11512.6074066973, 12000];
PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds ./ active_binning_resolution;
PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = [floor(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1)), ...
    ceil(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2))];

if processing_config.showOnlyAlwaysStableCells
    isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, isAlwaysStable);
    PhoDibaTest_PositionalAnalysis_temp.active_spikes = active_processing.spikes.time(isAlwaysStable);
    numActiveCells = numAlwaysStableCells;
else
    PhoDibaTest_PositionalAnalysis_temp.active_spikes = active_processing.spikes.time;
    numActiveCells = length(active_processing.spikes.isAlwaysStable);
end

PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix'; % should be units x time
% Extract only the portion specified as the training portion
PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1):PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2));

% Extract the cells as an array
PhoDibaTest_PositionalAnalysis_temp.cell_indicies = num2cell([1:numActiveCells]);
PhoDibaTest_PositionalAnalysis_temp.spike_cells = cellfun(@(cell_idx) [(cell_idx .* ones(size(PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_idx}'))), PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_idx}'], ...
    PhoDibaTest_PositionalAnalysis_temp.cell_indicies, ...
    'UniformOutput', false);

% [includedCellIDs, unitSpikeCells, unitFlatIndicies] = fnFlatSpikesToUnitCells(spikeStruct.t, spikeStruct.unit, true);
% cellfun(@(x) spikeStruct.ripple(find(x)), unitFlatIndicies, 'UniformOutput', false)


% Save out positionalAnalysis data for Python:
export_root_path = '/Users/pho/repo/Python Projects/PhoNeuronGillespie2021CodeRepo/PhoMatlabDataScripting/ExportedData';
active_experiment_export_root_path = fullfile(export_root_path, active_experiment_names{active_expt_index});
mkdir(active_experiment_export_root_path);


fprintf('Saving spikes analysis data to %s...\n', fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'));
spike_cells = PhoDibaTest_PositionalAnalysis_temp.spike_cells;
spike_matrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix;
save(fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'), 'spike_matrix', 'spike_cells')
fprintf('done!\n');

fprintf('PhoDibaConvert_SpikesAnalysis complete!\n');



