% Requires the datastructures from "PhoDibaAnalyze.m" to be loaded

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));
% clear all;

process_config.show_graphics = false;


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



%% Pairwise Indexing:
% Generates all unique pairs of indicies for pairwise comparisons (without replacement or repetition)
active_result.indicies.unique_electrode_pairs = nchoose2([1:num_of_electrodes]);
active_results.indicies.num_unique_pairs = length(active_result.indicies.unique_electrode_pairs);



%% Auto-Correlations:

% % [autocor,lags] = xcorr(active_processing.processed.smoothed_spike_data{1});
% [acf,lags,bounds] = autocorr(active_processing.processed.smoothed_spike_data{1},'NumMA',2);

active_results.autocorrelations = cellfun((@(spike_timestamps) xcorr(spike_timestamps)), ...
        active_processing.processed.smoothed_spike_data, 'UniformOutput', false);
    
active_results.partial_autocorrelations = cellfun((@(spike_timestamps) parcorr(spike_timestamps)), ...
        active_processing.processed.smoothed_spike_data, 'UniformOutput', false);    

% figure(2)
% stem(lags,autocor)
%corr2

%% Cross-Correlations:
% active_processing.processed.smoothed_spike_data

%% Extract all timeseries in the appropriate matrix format:
temp.smoothed_spike_data_matrix = cell2mat(active_processing.processed.smoothed_spike_data'); % 35351x126 double

% rho: Pairwise linear correlation coefficient
% 
% [active_results.pairwise_correlations.rho, active_results.pairwise_correlations.pval] = corr(temp.smoothed_spike_data_matrix);
process_config.max_xcorr_lag = 6; % Specified the maximum pairwise cross-correlation lag in seconds, the output ranges from -maxlag to maxlag
active_results.pairwise_xcorrelations.lag_offsets = (-process_config.max_xcorr_lag):process_config.max_xcorr_lag;
active_results.pairwise_xcorrelations.num_lag_steps = length(active_results.pairwise_xcorrelations.lag_offsets);

% [active_results.pairwise_xcorrelations.xcorr, active_results.pairwise_xcorrelations.lags] = xcorr(temp.smoothed_spike_data_matrix, process_config.max_xcorr_lag, 'normalized'); 
% xcorr: 13x15876 double ((2 × max_xcorr_lag + 1) × N^2)
% lags: 1x13 double


% pre-allocate output
active_results.pairwise_xcorrelations.xcorr = zeros([active_result.indicies.num_unique_pairs active_results.pairwise_xcorrelations.num_lag_steps]);
parfor i = 1:active_result.indicies.num_unique_pairs
   temp.curr_pair = active_result.indicies.unique_electrode_pairs(i,:);
   active_results.pairwise_xcorrelations.xcorr(i,:) = xcorr(temp.smoothed_spike_data_matrix(:,temp.curr_pair(1)), temp.smoothed_spike_data_matrix(:,temp.curr_pair(2)), process_config.max_xcorr_lag);
end



temp.active_idx = 1:5;

stem(active_results.pairwise_xcorrelations.lags, active_results.pairwise_xcorrelations.xcorr(:, temp.active_idx));






%% Aggregate Measures:

%% Split based on experiment epoch:
for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes = cell2mat(active_processing.processed.by_epoch.(temp.curr_epoch_name).smoothed_spike_data);
    
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).across_all_cells.count = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 2);
    active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 'all');
    
    fprintf('epoch: %s\n total_counts: %d\n', temp.curr_epoch_name, active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts);
end


%timesteps

if process_config.show_graphics
    figure(1);
end

%% Split based on behavioral state:
for i = 1:length(active_processing.behavioral_state_names)
    temp.curr_state_name =  active_processing.behavioral_state_names{i};
    active_results.aggregates.by_state.(temp.curr_state_name).spikes = cell2mat(active_processing.processed.by_state.(temp.curr_state_name).smoothed_spike_data);
    
    
    active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 2);
    active_results.aggregates.by_state.(temp.curr_state_name).total_counts = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 'all');
    
    
    fprintf('state: %s\n total_counts: %d\n', temp.curr_state_name, active_results.aggregates.by_state.(temp.curr_state_name).total_counts);
    
    if process_config.show_graphics
        subplot(4,1,i);
        plot(active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count);
        ylabel(temp.curr_state_name);
        xlabel('');
    end
end

if process_config.show_graphics
    xlim([timesteps(1), timesteps(end)]);
    title('behavioral state spike counts')
end