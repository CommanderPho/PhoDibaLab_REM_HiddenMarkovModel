% Saves out the behavioral_epochs table data and the behavioral_periods
% table for use in python.
if ~exist('timesteps_array','var')
    load('PhoResults_Expt1_RoyMaze1.mat', 'timesteps_array')
    load('PhoResults_Expt1_RoyMaze1.mat', 'active_processing')
end

% Save out positionalAnalysis data for Python:
% export_root_path = '/Users/pho/repo/Python Projects/PhoNeuronGillespie2021CodeRepo/PhoMatlabDataScripting/ExportedData';
export_root_path = 'R:\rMBP Python Repos 2022-07-07\PhoNeuronGillespie2021CodeRepo\PhoMatlabDataScripting\ExportedData';
active_experiment_export_root_path = fullfile(export_root_path, active_experiment_name, 'ExportedData');
mkdir(active_experiment_export_root_path);
fprintf('Saving extras analysis data to %s...\n', fullfile(active_experiment_export_root_path, 'extrasAnalysis.mat'));

% Epoch names from table Row header:
behavioral_epoch_names = active_processing.behavioral_epochs.Row;

% Numerical table version:
behavioral_epochs = [[0:(height(active_processing.behavioral_epochs)-1)]', table2array(active_processing.behavioral_epochs)];
behavioral_periods = [[0:(height(active_processing.behavioral_periods_table)-1)]', double(active_processing.behavioral_periods_table.epoch_start_seconds), double(active_processing.behavioral_periods_table.epoch_end_seconds), double(active_processing.behavioral_periods_table.duration), double(active_processing.behavioral_periods_table.behavioral_epoch), double(active_processing.behavioral_periods_table.type)];

save(fullfile(active_experiment_export_root_path, 'extrasAnalysis.mat'), 'behavioral_epochs', 'behavioral_periods', 'behavioral_epoch_names')
fprintf('done!\n');
fprintf('Extras export complete!\n');

