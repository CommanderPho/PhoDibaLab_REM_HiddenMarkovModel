% PhoDibaResults_AllExpts.m
% 02/18/2021 by Pho Hale
% Displays results across all experiments (animals x days)
% Requires: across_experiment_results, active_experiment_names, active_step_sizes

%%% Note: This is a very inefficient implementation that is to serve prototyping purposes, but should be refactored for deployment scenarios.

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));


%% Binning Options:
current_binning_index = 1;

%% Filtering Options:
temp.filter_included_cell_types = {'pyramidal'};
% temp.filter_included_cell_types = {'interneurons'};
temp.filter_maximum_included_contamination_level = {2};


active_binning_resolution = active_step_sizes{current_binning_index};
active_num_experiments = length(active_experiment_names);

% Loop through each experiment:
for expt_index = 1:active_num_experiments
    temp.active_expt_name = active_experiment_names{expt_index};

    temp.curr_timestamps = across_experiment_results{expt_index}.timesteps_array{current_binning_index};
    temp.curr_processed = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index};
    active_results = across_experiment_results{expt_index}.results_array{current_binning_index};

    [temp.filter_active_units] = fnFilterUnitsWithCriteria(across_experiment_results{expt_index}.active_processing, ...
        processing_config.showOnlyAlwaysStableCells, temp.filter_included_cell_types, ...
    temp.filter_maximum_included_contamination_level);


end

















function [h] = fnPlotAcrossREMTesting(active_results, v2, v3, v4)

    %% xcorr_all_pairs
    % Get REM only pairs:
    temp.results.pre_sleep_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(temp.filtered.pre_sleep_REM_indicies, :); % 14x19
    temp.results.post_sleep_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(temp.filtered.post_sleep_REM_indicies, :); % 9x19
    % Get other pairs:
    temp.results.all_except_REM.per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(temp.filtered.all_except_REM_indicies, :); % (remainder)x19

    %% xcorr_all_pairs_AND_lags:
    % Get REM only pairs:
    temp.results.pre_sleep_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(temp.filtered.pre_sleep_REM_indicies, :); % 14x1
    temp.results.post_sleep_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(temp.filtered.post_sleep_REM_indicies, :); % 9x1
    % Get other pairs:
    temp.results.all_except_REM.per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(temp.filtered.all_except_REM_indicies, :); % (remainder)x1

    %% TODO: Plot all others (states that aren't REM at all) as a separate series

    %% Display the Correlational Results:
    xcorr_fig = figure(15);
    [xcorr_fig, h1, h2] = fnPlotAcrossREMXcorrHeatmap(temp.results.pre_sleep_REM.per_period.xcorr_all_pairs, ...
        temp.results.post_sleep_REM.per_period.xcorr_all_pairs);
    sgtitle('XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units')
    
    
end



function [xcorr_fig, h1, h2] = fnPlotAcrossREMXcorrHeatmap(v1, v2)
    % Plot a heatmap of the xcorr    
%     xcorr_fig = figure(15);
    clf
    subplot(2,1,1);

    h1 = heatmap(v1);
    ylabel('Filtered Trial Index')
    xlabel('Time Lag')
    title(sprintf('PRE sleep REM periods: %d', temp.results.pre_sleep_REM.num_behavioral_periods));


    subplot(2,1,2);
    h2 = heatmap(v2);
    ylabel('Filtered Trial Index')
    xlabel('Time Lag')    
    title(sprintf('POST sleep REM periods: %d', temp.results.post_sleep_REM.num_behavioral_periods));

%     sgtitle('XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units')

end

