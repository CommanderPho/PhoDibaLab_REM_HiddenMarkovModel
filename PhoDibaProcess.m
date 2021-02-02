% Requires the datastructures from "PhoDibaAnalyze.m" to be loaded

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));
% clear all;

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.output.intermediate_file_name = 'PhoIntermediate.mat';

data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);
data_config.output.intermediate_file_path = fullfile(data_config.source_root_path, data_config.output.intermediate_file_name);


if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading from %s...\n', data_config.output.intermediate_file_path);
    load(data_config.output.intermediate_file_path, 'active_processing', 'data_config', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps');
    fprintf('done.\n');

else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end





%% Split based on experiment epoch:
for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes = cell2mat(active_processing.processed.by_epoch.(temp.curr_epoch_name).smoothed_spike_data);
    
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).across_all_cells.count = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 2);
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 'all');
    
    fprintf('epoch: %s\n total_counts: %d\n', temp.curr_epoch_name, active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts);
end


%timesteps


figure(1);

%% Split based on behavioral state:
for i = 1:length(active_processing.behavioral_state_names)
    temp.curr_state_name =  active_processing.behavioral_state_names{i};
    active_results.aggregates.by_state.(temp.curr_state_name).spikes = cell2mat(active_processing.processed.by_state.(temp.curr_state_name).smoothed_spike_data);
    
    
    active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 2);
    active_results.aggregates.by_state.(temp.curr_state_name).total_counts = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 'all');
    
    
    fprintf('state: %s\n total_counts: %d\n', temp.curr_state_name, active_results.aggregates.by_state.(temp.curr_state_name).total_counts);
    
    subplot(4,1,i);
    plot(active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count);
    ylabel(temp.curr_state_name);
    xlabel('');
end

xlim([timesteps(1), timesteps(end)]);
title('behavioral state spike counts')
