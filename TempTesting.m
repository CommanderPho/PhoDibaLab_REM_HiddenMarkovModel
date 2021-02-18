clc
addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));


% clear all;

temp.possibleRootPaths = {'/Volumes/iNeo/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets', ...
	'/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets'};

temp.possibleExperiments = {'RoyMaze1',	'RoyMaze2',	'RoyMaze3',	'TedMaze1',	'TedMaze2',	'TedMaze3',	'KevinMaze1'};

% active_expt_name = 'TedMaze3';
active_expt_name = 'KevinMaze1';
% auxillary pink rMBP:
active_root_path = temp.possibleRootPaths{1};
% % main rMBP:
% active_root_path = temp.possibleRootPaths{2};

active_step_sizes = {0.1, 1.0};
fprintf('fnRunPipeline starting for active_expt_name: %s...\n', active_expt_name);
[data_config, processing_config, plotting_options] = fnDataConfigForExptName(active_expt_name, active_step_sizes, active_root_path);
PhoDibaPrepare_Stage0 % active_processing, source_data, timesteps_array
PhoDibaProcess_Stage1 % general_results, results_array
PhoDibaProcess_Stage2 % general_results

fprintf('fnRunPipeline completed for active_expt_name: %s\n', active_expt_name);
