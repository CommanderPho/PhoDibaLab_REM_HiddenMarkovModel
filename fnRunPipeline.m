function [source_data, timesteps_array, active_processing, general_results, results_array] = fnRunPipeline(active_expt_name, active_step_sizes)
%FNRUNPIPELINE Run all stages of the pipeline for a given animal name
%   Detailed explanation goes here

fprintf('fnRunPipeline starting for active_expt_name: %s...\n', active_expt_name);
[data_config, processing_config, plotting_options] = fnDataConfigForExptName(active_expt_name, active_step_sizes);
% run(PhoDibaPrepare_Stage0); % it will load all the variables inside the function workspace
% run(PhoDibaProcess_Stage1); % it will load all the variables inside the function workspace
% run(PhoDibaProcess_Stage2); % it will load all the variables inside the function workspace

PhoDibaPrepare_Stage0 % active_processing, source_data, timesteps_array
PhoDibaProcess_Stage1 % general_results, results_array
PhoDibaProcess_Stage2 % general_results

fprintf('fnRunPipeline completed for active_expt_name: %s\n', active_expt_name);
end

