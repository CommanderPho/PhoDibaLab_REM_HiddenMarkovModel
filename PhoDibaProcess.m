% Requires the datastructures from "PhoDibaAnalyze.m" to be loaded

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));
% clear all;

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
data_config.output.results_file_name = 'PhoResults.mat';


data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);
data_config.output.intermediate_file_path = fullfile(data_config.source_root_path, data_config.output.intermediate_file_name);
data_config.output.results_file_path = fullfile(data_config.source_root_path, data_config.output.results_file_name);

if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading from %s...\n', data_config.output.intermediate_file_path);
    load(data_config.output.intermediate_file_path, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end



%% Pairwise Indexing:
% Generates all unique pairs of indicies for pairwise comparisons (without replacement or repetition)
active_results.indicies.unique_electrode_pairs = nchoose2([1:num_of_electrodes]);
active_results.indicies.num_unique_pairs = length(active_results.indicies.unique_electrode_pairs);

% Build a reverse lookup matrix
active_results.indicies.reverse_lookup_unique_electrode_pairs = zeros(num_of_electrodes);
for linear_pair_idx = 1:active_results.indicies.num_unique_pairs
    curr_pair = active_results.indicies.unique_electrode_pairs(linear_pair_idx,:);
    active_results.indicies.reverse_lookup_unique_electrode_pairs(curr_pair(1), curr_pair(2)) = linear_pair_idx;
    active_results.indicies.reverse_lookup_unique_electrode_pairs(curr_pair(2), curr_pair(1)) = linear_pair_idx;
end


%% Build the Correlational Results:
process_config.max_xcorr_lag = 9; % Specified the maximum pairwise cross-correlation lag in seconds, the output ranges from -maxlag to maxlag
[active_results] = fnProcessCorrelationalMeasures(active_processing.processed.smoothed_spike_data, active_results, process_config);

%% Display the Correlational Results:
[temp.fig, temp.h] = fnPhoPlotCorrelationalResults(active_results, plottingOptions);






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

fprintf('writing out results to %s...\n', data_config.output.results_file_path);
save(data_config.output.results_file_path, 'active_results');
fprintf('done.\n');
