% PhoNeuroPyConvert_ExportAllToPython_MAIN.m
% PhoNeuroPyConvert_ExportAllToPython_MAIN - Aims to export the .mat-format
% data to NeuroPy 2022-07-08 Session format.
% Makes use of active_processing.position_table, active_processing.behavioral_epochs.start_seconds
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 08-07-2022 ; Last revision: 08-07-2022 
% History: was based off of "PhoDibaConvert_SpikesToPython.m"

%%%%%%%%% CONFIG
% export_root_path = '/Users/pho/repo/Python Projects/PhoNeuronGillespie2021CodeRepo/PhoMatlabDataScripting/ExportedData';
export_root_path = 'R:\rMBP Python Repos 2022-07-07\PhoNeuronGillespie2021CodeRepo\PhoMatlabDataScripting\ExportedData';
%% Filtering Options:
filter_config.filter_included_cell_types = {};
% filter_config.filter_maximum_included_contamination_level = {2};
filter_config.filter_maximum_included_contamination_level = {};
% filter_config.showOnlyAlwaysStableCells = true;
filter_config.showOnlyAlwaysStableCells = false;



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
%% Begin:
fprintf('PhoNeuroPyConvert_ExportAllToPython_MAIN ready to process!\n');


%% Binning Options:
active_expt_index = 1;
if exist('active_experiment_names','var')
    active_experiment_name = active_experiment_names{active_expt_index};
else
    active_experiment_name = processing_config.active_expt.name; % get from processing config
end


% plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;
% current_binning_index = 1;
% active_binning_resolution = processing_config.step_sizes{current_binning_index};
% temp.curr_timestamps = timesteps_array{current_binning_index};
% temp.curr_processed = active_processing.processed_array{current_binning_index};
% % active_results = results_array{current_binning_index};
% fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);

% Binned Spike rates per time:
% [PhoDibaTest_PositionalAnalysis_temp.activeMatrix] = fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts);  % 35351x126 double
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds = [11512.6074066973, 12000];
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_seconds ./ active_binning_resolution;
% PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins = [floor(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1)), ...
    % ceil(PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2))];


%% Get filter info for active units
[filtered_outputs.filter_active_units, filtered_outputs.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, filter_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);
temp.num_active_units = sum(filtered_outputs.filter_active_units, 'all');
fprintf('Filter: Including %d of %d total units\n', temp.num_active_units, length(filtered_outputs.filter_active_units));

%% Apply the filters:
PhoDibaTest_PositionalAnalysis_temp.active_spikes = active_processing.spikes.time(filtered_outputs.filter_active_units);

%% This seems to be unused:
% numActiveCells = temp.num_active_units;
% PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, filtered_outputs.filter_active_units);
% PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix'; % should be units x time
% % Extract only the portion specified as the training portion
% PhoDibaTest_PositionalAnalysis_temp.activeMatrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix(:, PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(1):PhoDibaTest_PositionalAnalysis_config.training_subset_start_stop_bins(2));
% % Extract the cells as an array
% % PhoDibaTest_PositionalAnalysis_temp.cell_indicies = num2cell([1:numActiveCells]);


%% This is the essential loading function: relies on .active_spikes
PhoDibaTest_PositionalAnalysis_temp.cell_indicies = filtered_outputs.original_unit_index; % Get the indicies of the remaining cells:
PhoDibaTest_PositionalAnalysis_temp.spike_cells = cellfun(@(cell_idx) [(cell_idx .* ones(size(PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_idx}'))), PhoDibaTest_PositionalAnalysis_temp.active_spikes{cell_idx}'], ...
    num2cell([1:numActiveCells]), ...
    'UniformOutput', false);


%% Pho Spikes Cell Info 2022-07-08

PhoDibaTest_PositionalAnalysis_temp.speculated_unit_type = cellstr(active_processing.spikes.speculated_unit_type(filtered_outputs.filter_active_units));
% behavioral_epochs = [[0:(height(active_processing.behavioral_epochs)-1)]', table2array(active_processing.behavioral_epochs)];

%% Properties are ['shank', 'cluster', 'aclu', 'qclu']
PhoDibaTest_PositionalAnalysis_temp.shank = active_processing.spikes.id(filtered_outputs.filter_active_units, 1); % shank
PhoDibaTest_PositionalAnalysis_temp.cluster = active_processing.spikes.id(filtered_outputs.filter_active_units, 2); % cluster
PhoDibaTest_PositionalAnalysis_temp.aclu = PhoDibaTest_PositionalAnalysis_temp.cell_indicies; % (aclu)
PhoDibaTest_PositionalAnalysis_temp.qclu = active_processing.spikes.quality(filtered_outputs.filter_active_units); % quality (qclu)
PhoDibaTest_PositionalAnalysis_temp.speculated_unit_contamination_level = active_processing.spikes.speculated_unit_contamination_level(filtered_outputs.filter_active_units); % speculated_unit_contamination_level


% Make the path for the active experiment:
active_experiment_export_root_path = fullfile(export_root_path, active_experiment_name, 'ExportedData');
mkdir(active_experiment_export_root_path);



%% Spikes:
% The critical properties are: 'spike_cells'

fprintf('Saving spikes analysis data to %s...\n', fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'));
spike_cells_ids = PhoDibaTest_PositionalAnalysis_temp.cell_indicies;
spike_cells = PhoDibaTest_PositionalAnalysis_temp.spike_cells;
% spike_matrix = PhoDibaTest_PositionalAnalysis_temp.activeMatrix;
shank = PhoDibaTest_PositionalAnalysis_temp.shank;
cluster = PhoDibaTest_PositionalAnalysis_temp.cluster;
aclu = PhoDibaTest_PositionalAnalysis_temp.aclu;
qclu = PhoDibaTest_PositionalAnalysis_temp.qclu;
speculated_unit_contamination_level = PhoDibaTest_PositionalAnalysis_temp.speculated_unit_contamination_level;
speculated_unit_type = PhoDibaTest_PositionalAnalysis_temp.speculated_unit_type;

% save(fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'), 'spike_matrix', 'spike_cells', 'spike_cells_ids')
% save(fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'), 'spike_matrix', 'spike_cells', 'spike_cells_ids', 'shank', 'cluster', 'aclu', 'qclu', 'speculated_unit_contamination_level', 'speculated_unit_type')
save(fullfile(active_experiment_export_root_path, 'spikesAnalysis.mat'), 'spike_cells', 'spike_cells_ids', 'shank', 'cluster', 'aclu', 'qclu', 'speculated_unit_contamination_level', 'speculated_unit_type')

fprintf('done!\n');
fprintf('PhoDibaConvert_SpikesAnalysis complete!\n');


% Save out positionalAnalysis data for Python:




% ------------- END OF CODE --------------