% Requires the datastructures from "PhoDibaProcess_Stage1.m" to be loaded
% Stage 2 of the processing pipeline

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

%% Configure Graphics and Plotting:
process_config.show_graphics = false;
% Options for tightening up the subplots:
plotting_options.should_use_custom_subplots = true;

if plotting_options.should_use_custom_subplots
    plotting_options.subtightplot.gap = [0.01 0.01]; % [intra_graph_vertical_spacing, intra_graph_horizontal_spacing]
    plotting_options.subtightplot.width_h = [0.01 0.05]; % Looks like [padding_bottom, padding_top]
    plotting_options.subtightplot.width_w = [0.025 0.01];
    plotting_options.opt = {plotting_options.subtightplot.gap, plotting_options.subtightplot.width_h, plotting_options.subtightplot.width_w}; % {gap, width_h, width_w}
    subplot_cmd = @(m,n,p) subtightplot(m, n, p, plotting_options.opt{:});
else
    subplot_cmd = @(m,n,p) subplot(m, n, p);
end

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.output.intermediate_file_name = 'PhoIntermediate.mat';
data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);
data_config.output.intermediate_file_path = fullfile(data_config.source_root_path, data_config.output.intermediate_file_name);



if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading data from %s...\n', data_config.output.intermediate_file_path);
    load(data_config.output.intermediate_file_path, 'data_config', 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end

% Add the results file path:
data_config.output.results_file_name = 'PhoResults.mat';
data_config.output.results_file_path = fullfile(data_config.source_root_path, data_config.output.results_file_name);

if ~exist('active_results','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'active_results');
    fprintf('done.\n');
else
    fprintf('active_results already exists in workspace. Using extant data.\n');
end

%% Begin:

