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

%% Plotting Options:
plottingOptions.plottingXAxis = 'index';
% plottingOptions.plottingXAxis = 'timestamp';
plottingOptions.plottingYlim = [];
% plottingOptions.plottingYlim = [2 4.25];
% plottingOptions.plottingYlim = [0.2 1.4];

plottingOptions.plottingMode = 'scatter';
% plottingOptions.plottingMode = 'errorbar';
% plottingOptions.plottingMode = 'distributionPlot'; % distributionPlot should display the variance across neurons

plottingOptions.outputs.rootPath = '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/Figures';


active_binning_resolution = active_step_sizes{current_binning_index};
active_num_experiments = length(active_experiment_names);

% Loop through each experiment:
for expt_index = 1:active_num_experiments
    expt_info.index = expt_index;
    expt_info.name = active_experiment_names{expt_index};
    
    temp.curr_timestamps = across_experiment_results{expt_index}.timesteps_array{current_binning_index};
    temp.curr_processed = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index};
    active_results = across_experiment_results{expt_index}.results_array{current_binning_index};

    [xcorr_fig, filtered, results] = fnPerformAcrossREMTesting(across_experiment_results{expt_index}.active_processing, ...
        across_experiment_results{expt_index}.general_results, ...
        across_experiment_results{expt_index}.results_array{current_binning_index}, ...
        across_experiment_results{expt_index}.processing_config, filter_config, expt_info, plottingOptions);
   
end





function [plotResults, filtered, results] = fnPerformAcrossREMTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions)
    %%% fnPerformAcrossREMTesting: Run main analysis
    %%%
    
    plotResults.exports = {};
    plotResults.configs = {};
    
    %% Get filter info for active units
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

    
    %%% Plotting Results:
    temp.curr_expt_string = sprintf('experiment[%d]: %s', expt_info.index, expt_info.name);
    plottingOptions.outputs.curr_expt_filename_string = sprintf('%s_', expt_info.name);
    
    
    %% TODO: Plot all others (states that aren't REM at all) as a separate series

    %% Error bars are across units:
    plotResults.figures.meanSpikeRateFig = figure(9+expt_info.index);
    clf
    subplot(2,1,1);

    if strcmpi(plottingOptions.plottingXAxis, 'index')
        plottingOptions.x = [1:results.pre_sleep_REM.num_behavioral_periods];
    else
        plottingOptions.x = results.pre_sleep_REM.per_period.epoch_center_seconds;
    end

    [h1] = fnPlotAcrossREMTesting(plottingOptions.plottingMode, plottingOptions.x, ...
        results.pre_sleep_REM.spike_rate_all_units.mean, ...
        results.pre_sleep_REM.spike_rate_all_units.stdDev, ...
        results.pre_sleep_REM.spike_rate_per_unit);

    title(sprintf('PRE sleep REM periods: %d', results.pre_sleep_REM.num_behavioral_periods));
    if strcmpi(plottingOptions.plottingXAxis, 'index')
        xlabel('Filtered Period Index')
    else
        xlabel('Period Timestamp Offset (Seconds)')
    end
    ylabel('mean spike rate')
    if ~isempty(plottingOptions.plottingYlim)
        ylim(plottingOptions.plottingYlim)
    end

    subplot(2,1,2);

    if strcmpi(plottingOptions.plottingXAxis, 'index')
        plottingOptions.x = [1:results.post_sleep_REM.num_behavioral_periods];
    else
        plottingOptions.x = results.post_sleep_REM.per_period.epoch_center_seconds;
    end
    [h2] = fnPlotAcrossREMTesting(plottingOptions.plottingMode, plottingOptions.x, ...
        results.post_sleep_REM.spike_rate_all_units.mean, ...
        results.post_sleep_REM.spike_rate_all_units.stdDev, ...
        results.post_sleep_REM.spike_rate_per_unit);

    title(sprintf('POST sleep REM periods: %d', results.post_sleep_REM.num_behavioral_periods));
    if strcmpi(plottingOptions.plottingXAxis, 'index')
        currPlotConfig.xlabel = 'Filtered Period Index';
        
    else
        currPlotConfig.xlabel = 'Period Timestamp Offset (Seconds)';

    end
    currPlotConfig.ylabel = 'mean spike rate';
    
    xlabel(currPlotConfig.xlabel)
    ylabel(currPlotConfig.ylabel)
    if ~isempty(plottingOptions.plottingYlim)
        ylim(plottingOptions.plottingYlim)
    end
    sgtitle([temp.curr_expt_string ' : Spike Rates - PRE vs Post Sleep REM Periods - Period Index - Pyramidal Only'])
    % Figure Name:
    %'Spike Rates - PRE vs Post Sleep REM Periods - Period Index';
    %'Spike Rates - PRE vs Post Sleep REM Periods - Timestamp Offset';

%     'Spike Rates - PRE vs Post Sleep REM Periods - Period Index - Pyramidal Only'
    % MainPlotVariableBeingCompared - PurposeOfComparison - IndependentVariable - FiltersAndConstraints
    
    % Build Figure Export File path:
    currPlotConfig.curr_expt_filename_string = sprintf('%s - %s - %s - %s - %s', ...
        plottingOptions.outputs.curr_expt_filename_string, ...
        'Spike Rates', ...
        'PRE vs Post Sleep REM Periods', ...
        currPlotConfig.xlabel, ...
        'Spike Rates'...
        );

    currPlotConfig.curr_expt_parentPath = plottingOptions.outputs.rootPath;
    currPlotConfig.curr_expt_path = fullfile(currPlotConfig.curr_expt_parentPath, currPlotConfig.curr_expt_filename_string);

    % Perform the export:
    [plotResults.exports{end+1}.export_result] = fnSaveFigureForExport(plotResults.figures.meanSpikeRateFig, currPlotConfig.curr_expt_path, true, false, false, true);
    plotResults.configs{end+1} = currPlotConfig;
    
    %% Display the Correlational Results:
    plotResults.figures.xcorr = figure(expt_info.index);
    [plotResults.figures.xcorr, h1, h2] = fnPlotAcrossREMXcorrHeatmap(results.pre_sleep_REM.per_period.xcorr_all_pairs, ...
        results.post_sleep_REM.per_period.xcorr_all_pairs);
    sgtitle([temp.curr_expt_string ' : XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units'])
    
    % Build Figure Export File path:
    currPlotConfig.curr_expt_filename_string = sprintf('%s - %s - %s - %s - %s', ...
        plottingOptions.outputs.curr_expt_filename_string, ...
        'XCorr for all pairs', ...
        'PRE vs Post Sleep REM Periods', ...
        currPlotConfig.xlabel, ...
        'All Units'...
        );

    currPlotConfig.curr_expt_parentPath = plottingOptions.outputs.rootPath;
    currPlotConfig.curr_expt_path = fullfile(currPlotConfig.curr_expt_parentPath, currPlotConfig.curr_expt_filename_string);
    % Perform the export:
    [plotResults.exports{end+1}.export_result] = fnSaveFigureForExport(plotResults.figures.xcorr, currPlotConfig.curr_expt_path, true, false, false, true);
    plotResults.configs{end+1} = currPlotConfig;

    
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


function [h] = fnPlotAcrossREMTesting(mode, v1, v2, v3, v4)

    if strcmpi(mode, 'errorbar')
        h = errorbar(v1, ...
            v2, ...
            v3);

    elseif strcmpi(mode, 'scatter')
        h(1) = scatter(v1, ...
            v2);
        
        hold on;
        h(2) = errorbar(v1, v2, v3, 'LineStyle','none');
        
        
    elseif strcmpi(mode, 'distributionPlot')
        h = distributionPlot(v4'); % defaults 

    elseif strcmpi(mode, 'bar')
        h = bar(v1, v2);
        
    else
       error('Invalid mode input!') 
    end
end

