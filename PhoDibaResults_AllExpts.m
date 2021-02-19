% PhoDibaResults_AllExpts.m
% 02/18/2021 by Pho Hale
% Displays results across all experiments (animals x days)
% Requires: across_experiment_results, active_experiment_names, active_step_sizes

%%% Note: This is a very inefficient implementation that is to serve prototyping purposes, but should be refactored for deployment scenarios.

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

if ~exist('across_experiment_results','var')
    %% Set the path to the combined across_experiment file:
    % data_config.output.all_expts_combined_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets';
    data_config.output.all_expts_combined_parent_path = '/Volumes/iNeo/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results';
    data_config.output.all_expts_combined_results_file_name = sprintf('PhoResults_AllExperiments.mat');

    data_config.output.all_expts_combined_results_file_path = fullfile(data_config.output.all_expts_combined_parent_path, data_config.output.all_expts_combined_results_file_name);

    %% Load Combined experiment results:
    fprintf('\t loading combined experiment results from %s... (this will take quite a while, up to 10 minutes)...\n', data_config.output.all_expts_combined_results_file_path);
    temp.variablesList = {'across_experiment_results', 'active_experiment_names', 'active_step_sizes'};
    load(data_config.output.all_expts_combined_results_file_path, temp.variablesList{:});

else
    warning('across_experiment_results already exists, using extant value.');
end

    
%% Binning Options:
current_binning_index = 1;

%% Filtering Options:
filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

temp.filter_states = {'rem'};
temp.filter_epochs = {'pre_sleep', 'post_sleep'};

active_binning_resolution = active_step_sizes{current_binning_index};
active_num_experiments = length(active_experiment_names);

% Loop through each experiment:
for expt_index = 1:active_num_experiments
    expt_info.index = expt_index;
    expt_info.name = active_experiment_names{expt_index};
    
    temp.curr_timestamps = across_experiment_results{expt_index}.timesteps_array{current_binning_index};
    temp.curr_processed = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index};
    active_results = across_experiment_results{expt_index}.results_array{current_binning_index};

    [xcorr_fig, filtered, results] = fnPlotAcrossREMTesting(across_experiment_results{expt_index}.active_processing, ...
        across_experiment_results{expt_index}.general_results, ...
        across_experiment_results{expt_index}.results_array{current_binning_index}, ...
        across_experiment_results{expt_index}.processing_config, filter_config, expt_info);
   
end





function [xcorr_fig, filtered, results] = fnPlotAcrossREMTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info)
    % Run main analysis


    [filter_config.filter_active_units] = fnFilterUnitsWithCriteria(active_processing, processing_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);
    fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));


    [filtered.pre_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'pre_sleep'}, {'rem'}); % 668x1
    [filtered.post_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'post_sleep'}, {'rem'}); % 668x1

    [filtered.any_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, [], {'rem'}); % 668x1
    filtered.all_except_REM_indicies = ~filtered.any_REM_indicies; % 668x1, the complement of the REM indicies


    results.any_REM.num_behavioral_periods = sum(filtered.any_REM_indicies,'all');
    results.pre_sleep_REM.num_behavioral_periods = sum(filtered.pre_sleep_REM_indicies,'all');
    results.post_sleep_REM.num_behavioral_periods = sum(filtered.post_sleep_REM_indicies,'all');
    fprintf('any_REM: %d periods\n pre_sleep_REM: %d periods\n post_sleep_REM: %d periods\n', results.any_REM.num_behavioral_periods, results.pre_sleep_REM.num_behavioral_periods, results.post_sleep_REM.num_behavioral_periods);

    % The number of spikes per unit
    % general_results.per_behavioral_state_period.num_spikes_per_unit(filter_config.filter_active_units, filtered.pre_sleep_REM_indicies); % 668x126
    % general_results.per_behavioral_state_period.spike_rate_per_unit(filter_config.filter_active_units, filtered.pre_sleep_REM_indicies); % 668x126


    results.pre_sleep_REM.per_period.durations = active_processing.behavioral_periods_table.duration(filtered.pre_sleep_REM_indicies);
    results.post_sleep_REM.per_period.durations = active_processing.behavioral_periods_table.duration(filtered.post_sleep_REM_indicies);

    results.pre_sleep_REM.per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(filtered.pre_sleep_REM_indicies);
    results.post_sleep_REM.per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(filtered.post_sleep_REM_indicies);

    results.pre_sleep_REM.per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(filtered.pre_sleep_REM_indicies);
    results.post_sleep_REM.per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(filtered.post_sleep_REM_indicies);

    % Compute the center of the epochs to plot the firing rates along an appropriately scaled x-axis:
    results.pre_sleep_REM.per_period.epoch_center_seconds = (results.pre_sleep_REM.per_period.epoch_start_seconds + floor(results.pre_sleep_REM.per_period.durations ./ 2.0));
    results.post_sleep_REM.per_period.epoch_center_seconds = (results.post_sleep_REM.per_period.epoch_start_seconds + floor(results.post_sleep_REM.per_period.durations ./ 2.0));


    % Leave in terms of the spike rates per unit (14x92 double):
    results.pre_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(filtered.pre_sleep_REM_indicies, filter_config.filter_active_units);
    results.post_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(filtered.post_sleep_REM_indicies, filter_config.filter_active_units);

    % results.pre_sleep_REM.spike_rate_per_unit: (14x92 double)
    % results.post_sleep_REM.spike_rate_per_unit: (9x92 double)

    % Average Across all of the units
    results.pre_sleep_REM.spike_rate_all_units.mean = mean(results.pre_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
    results.post_sleep_REM.spike_rate_all_units.mean = mean(results.post_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
    results.pre_sleep_REM.spike_rate_all_units.stdDev = std(results.pre_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double
    results.post_sleep_REM.spike_rate_all_units.stdDev = std(results.post_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double
    
    % Compute the average across the REM sessions in each epoch
    results.pre_sleep_REM.baseline_spike_rate_across_all.mean = mean(results.pre_sleep_REM.spike_rate_all_units.mean);
    results.post_sleep_REM.baseline_spike_rate_across_all.mean = mean(results.post_sleep_REM.spike_rate_all_units.mean);
    results.pre_sleep_REM.baseline_spike_rate_across_all.stdDev = std(results.pre_sleep_REM.spike_rate_all_units.mean);
    results.post_sleep_REM.baseline_spike_rate_across_all.stdDev = std(results.post_sleep_REM.spike_rate_all_units.mean);

    %% xcorr_all_pairs
    % Get REM only pairs:
    results.pre_sleep_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(filtered.pre_sleep_REM_indicies, :); % 14x19
    results.post_sleep_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(filtered.post_sleep_REM_indicies, :); % 9x19
    % Get other pairs:
    results.all_except_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(filtered.all_except_REM_indicies, :); % (remainder)x19

    %% xcorr_all_pairs_AND_lags:
    % Get REM only pairs:
    results.pre_sleep_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(filtered.pre_sleep_REM_indicies, :); % 14x1
    results.post_sleep_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(filtered.post_sleep_REM_indicies, :); % 9x1
    % Get other pairs:
    results.all_except_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(filtered.all_except_REM_indicies, :); % (remainder)x1

    %% TODO: Plot all others (states that aren't REM at all) as a separate series

    %% Display the Correlational Results:
    xcorr_fig = figure(expt_info.index);
    [xcorr_fig, h1, h2] = fnPlotAcrossREMXcorrHeatmap(results.pre_sleep_REM.per_period.xcorr_all_pairs, ...
        results.post_sleep_REM.per_period.xcorr_all_pairs);
    
    temp.curr_expt_string = sprintf('experiment[%d]: %s', expt_info.index, expt_info.name);
    sgtitle([temp.curr_expt_string ' : XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units'])
    
end



function [xcorr_fig, h1, h2] = fnPlotAcrossREMXcorrHeatmap(v1, v2)
    % Plot a heatmap of the xcorr    
%     xcorr_fig = figure(15);
    xcorr_fig = gcf;
    clf
    subplot(2,1,1);

    h1 = heatmap(v1);
    ylabel('Filtered Trial Index')
    xlabel('Time Lag')
    title(sprintf('PRE sleep REM periods: %d', size(v1,1)));


    subplot(2,1,2);
    h2 = heatmap(v2);
    ylabel('Filtered Trial Index')
    xlabel('Time Lag')    
    title(sprintf('POST sleep REM periods: %d', size(v2,1)));

%     sgtitle('XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units')

end

