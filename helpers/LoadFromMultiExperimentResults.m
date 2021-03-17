%% LoadFromMultiExperimentResults.m
% Loads active_processing and results_array from the specified experiment index for compatibility with single-experiment scripts.

% multi_experiment_config.active_experiment_id: the index of the active experiment to load
multi_experiment_config.active_experiment_id = 1;
% current_binning_index = 2;
% active_binning_resolution = processing_config.step_sizes{current_binning_index};
% temp.curr_timestamps = timesteps_array{current_binning_index};
% temp.curr_processed = active_processing.processed_array{current_binning_index};
% active_results = results_array{current_binning_index};


if ~exist('data_config','var')
    Config;
end

if ~exist('across_experiment_results','var')
   error('ERROR: across_experiment_results must be loaded!') 
   
else
    if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
        temp.variables_to_load = {'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array','results_array','general_results'};
        for i = 1:length(temp.variables_to_load)
            assignin('base',temp.variables_to_load{i}, across_experiment_results{multi_experiment_config.active_experiment_id}.(temp.variables_to_load{i}));
        end
    else
        error('active_processing will not be overwritten by this function!');
        fprintf('active_processing already exists in workspace. Using extant data.\n');
    end

end


