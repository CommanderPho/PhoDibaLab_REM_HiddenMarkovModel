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
    load(data_config.output.intermediate_file_path, 'active_processing', 'data_config', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps');
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
process_config.max_xcorr_lag = 9; % Specified the maximum pairwise cross-correlation lag in seconds, the output ranges from -maxlag to maxlag
active_results.pairwise_xcorrelations.lag_offsets = (-process_config.max_xcorr_lag):process_config.max_xcorr_lag;
active_results.pairwise_xcorrelations.num_lag_steps = length(active_results.pairwise_xcorrelations.lag_offsets);

% [active_results.pairwise_xcorrelations.xcorr, active_results.pairwise_xcorrelations.lags] = xcorr(temp.smoothed_spike_data_matrix, process_config.max_xcorr_lag, 'normalized'); 
% xcorr: 13x15876 double ((2 × max_xcorr_lag + 1) × N^2)
% lags: 1x13 double


% pre-allocate output

% parfor unpacking
smoothed_spike_data_matrix = temp.smoothed_spike_data_matrix;
unique_electrode_index_pairs = active_results.indicies.unique_electrode_pairs;
max_xcorr_lag = process_config.max_xcorr_lag;

pairwise_xcorrelations = zeros([active_results.indicies.num_unique_pairs active_results.pairwise_xcorrelations.num_lag_steps]);
parfor i = 1:active_results.indicies.num_unique_pairs
%    temp.curr_pair = active_result.indicies.unique_electrode_pairs(i,:);
   pairwise_xcorrelations(i,:) = xcorr(smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,1)), ...
       smoothed_spike_data_matrix(:, unique_electrode_index_pairs(i,2)), ...
       max_xcorr_lag,'normalized');
end

active_results.pairwise_xcorrelations.xcorr = pairwise_xcorrelations;

% parfor cleanup
clear smoothed_spike_data_matrix unique_electrode_index_pairs max_xcorr_lag pairwise_xcorrelations

% Find maximum correlations and when they occur
[active_results.pairwise_xcorrelations.processed.MaxXCorr, active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex] = max(active_results.pairwise_xcorrelations.xcorr');
[active_results.pairwise_xcorrelations.processed.MinXCorr, active_results.pairwise_xcorrelations.processed.MinXCorrLagIndex] = min(active_results.pairwise_xcorrelations.xcorr');
active_results.pairwise_xcorrelations.processed.meanXCorr = mean(active_results.pairwise_xcorrelations.xcorr, 2)';

% active_results.pairwise_xcorrelations.processed.MaxXCorr: contains the maximum xcorr value obtained for each unit
% active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex: contains when the maximum xcorr value occurred for each unit

[temp.sorted_values, temp.sorted_index] = sortrows(active_results.pairwise_xcorrelations.processed.MaxXCorr','descend');

% Get the sorted indicies:
active_results.pairwise_xcorrelations.processed.sorted.unique_electrode_index_pairs = active_results.indicies.unique_electrode_pairs(temp.sorted_index, :);
% Get the sorted xcorr values:
active_results.pairwise_xcorrelations.processed.sorted.MaxXCorr = active_results.pairwise_xcorrelations.processed.MaxXCorr(temp.sorted_index);
active_results.pairwise_xcorrelations.processed.sorted.MaxXCorrLagIndex = active_results.pairwise_xcorrelations.processed.MaxXCorrLagIndex(temp.sorted_index);
active_results.pairwise_xcorrelations.processed.sorted.xcorr = active_results.pairwise_xcorrelations.xcorr(temp.sorted_index, :);



%% Display the Correlational Results:
%%%%%%%%%%%%%%%%%%%%%
figure(3);
clf
% temp.active_idx = 1:5;
% temp.active_idx = num_of_electrodes-5:num_of_electrodes;
temp.active_idx = 10;

% temp.found_lin_idx = find((active_results.indicies.unique_electrode_pairs(:,1) == temp.active_idx) | (active_results.indicies.unique_electrode_pairs(:,2) == temp.active_idx));
temp.found_lin_idx = active_results.indicies.reverse_lookup_unique_electrode_pairs(temp.active_idx, :); % 1x126 double

% temp.preview_subset = temp.active_idx+1:temp.active_idx+80;
temp.preview_subset = 1:num_of_electrodes;
temp.preview_subset_size = length(temp.preview_subset);
temp.active_found_lin_idx = temp.found_lin_idx(temp.preview_subset);

% Overlay version:
% stem(active_results.pairwise_xcorrelations.lags, active_results.pairwise_xcorrelations.processed.sorted.xcorr(temp.active_found_lin_idx, :)','filled');

% % Subplot version:
% % tiledlayout(temp.preview_subset_size, 1)
% for i = 1:temp.preview_subset_size
%     % Subplot
%     temp.isLastSession = (i == temp.preview_subset_size);
%     ax(i) = subplot_cmd(temp.preview_subset_size, 1, i);
%     stem(ax(i), active_results.pairwise_xcorrelations.lags, active_results.pairwise_xcorrelations.processed.sorted.xcorr(temp.active_found_lin_idx(i), :)')
%     
%     set(gca,'ytick',[],'yticklabel',[])
%     
%     if temp.isLastSession
%        xlabel('t');
%     else
% %         title('');
% %         xlabel('');
% %         ylabel('');
%         set(gca,'xtick',[],'xticklabel',[])
% %         set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
%     end
%     
% end


%% Image Matrix Version:
% fnPhoMatrixPlot(active_results.pairwise_xcorrelations.xcorr');

temp.excluded_indicies = find(temp.active_found_lin_idx < 1);
temp.active_found_lin_idx(temp.excluded_indicies) = 1; % temporarily fill it in with a valid index (referring to row 1), but then replace it afterwards with NaN
temp.plot_matrix = active_results.pairwise_xcorrelations.xcorr(temp.active_found_lin_idx, :);
temp.plot_matrix(temp.excluded_indicies, :) = NaN;

% % Cell formatter: "[unitID - timeOffset]: xcorr_value"
% plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('[%d - %d]: %d', row_index, active_results.pairwise_xcorrelations.lag_offsets(column_index), temp.plot_matrix(row_index, column_index)); 

% Cell formatter: "unitID[@t=timeOffset]"
plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('%d[@t=%d]', row_index, active_results.pairwise_xcorrelations.lag_offsets(column_index)); 

[h, temp.plot_info] = fnPhoMatrixPlotDetailed(temp.plot_matrix, plottingOptions);
xlabel('[seconds]');
ylabel('xcorr');
title(sprintf('Pairwise XCorr for Unit %d',  temp.active_idx));






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
