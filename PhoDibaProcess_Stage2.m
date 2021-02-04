% Requires the datastructures from "PhoDibaProcess_Stage1.m" to be loaded
% Stage 2 of the processing pipeline

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

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

% Add the results file path:
data_config.output.results_file_name = 'PhoResults.mat';
data_config.output.results_file_path = fullfile(data_config.source_root_path, data_config.output.results_file_name);

if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'results_array');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('PhoDibaProcess_Stage2 ready to process!\n');
