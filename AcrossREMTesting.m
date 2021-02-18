% AcrossREMTesting.m
% 02-17-2020 by Pho Hale

% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

clear temp

if ~exist('data_config','var')
    Config;
end


if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading data from %s...\n', data_config.output.intermediate_file_paths{2});
    load(data_config.output.intermediate_file_paths{2}, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end

if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'results_array', 'general_results');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('AcrossREMTesting ready to process!\n');

clear temp;

num_of_behavioral_state_periods = height(active_processing.behavioral_periods_table);
% 
% general_results.per_behavioral_state_period.spike_rate_per_unit

% All units with qualities from 1-4 are pyramidal.
% The higher the number in this range, the higher is the contamination, so 1 and 2 are well-separated pyramidal units and if your analysis is not much sensitive to contaminations you can consider 3 and 4 as well. For my analysis, I considered 1 to 3. 8 and 9 are interneurons.
% 

% general_results.per_behavioral_state_period.num_spikes_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double
% general_results.per_behavioral_state_period.spike_rate_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double


temp.filter_included_cell_types = {'pyramidal'};
% temp.filter_included_cell_types = {'interneurons'};
temp.filter_maximum_included_contamination_level = {2};
[temp.filter_active_units] = fnFilterUnitsWithCriteria(active_processing, processing_config.showOnlyAlwaysStableCells, temp.filter_included_cell_types, ...
    temp.filter_maximum_included_contamination_level);


fprintf('Filter: Including %d of %d total units\n', sum(temp.filter_active_units, 'all'), length(temp.filter_active_units));


%% Filter by Epoch:



%% What I was looking for before. Can filter by the epoch and state indicies and interest and collapse across trials
% general_results.GroupedByState.groups = findgroups(active_processing.behavioral_periods_table.type);
% general_results.GroupedByEpoch.groups = findgroups(active_processing.behavioral_periods_table.behavioral_epoch);

temp.filter_states = {'rem'};
temp.filter_epochs = {'pre_sleep', 'post_sleep'};

[temp.filtered.pre_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'pre_sleep'}, {'rem'}); % 668x1
[temp.filtered.post_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'post_sleep'}, {'rem'}); % 668x1

temp.results.pre_sleep_REM.num_behavioral_periods = sum(temp.filtered.pre_sleep_REM_indicies,'all');
temp.results.post_sleep_REM.num_behavioral_periods = sum(temp.filtered.post_sleep_REM_indicies,'all');
fprintf('pre_sleep_REM: %d periods\n post_sleep_REM: %d periods\n', temp.results.pre_sleep_REM.num_behavioral_periods, temp.results.post_sleep_REM.num_behavioral_periods);


% Alternative: splitapply workflow:
% splitapply(@mean,Height,G)

% The number of spikes per unit
% general_results.per_behavioral_state_period.num_spikes_per_unit(temp.filter_active_units, temp.filtered.pre_sleep_REM_indicies); % 668x126
% general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filter_active_units, temp.filtered.pre_sleep_REM_indicies); % 668x126


temp.results.pre_sleep_REM.per_period.durations = active_processing.behavioral_periods_table.duration(temp.filtered.pre_sleep_REM_indicies);
temp.results.post_sleep_REM.per_period.durations = active_processing.behavioral_periods_table.duration(temp.filtered.post_sleep_REM_indicies);

temp.results.pre_sleep_REM.per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(temp.filtered.pre_sleep_REM_indicies);
temp.results.post_sleep_REM.per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(temp.filtered.post_sleep_REM_indicies);

temp.results.pre_sleep_REM.per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(temp.filtered.pre_sleep_REM_indicies);
temp.results.post_sleep_REM.per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(temp.filtered.post_sleep_REM_indicies);

% Compute the center of the epochs to plot the firing rates along an appropriately scaled x-axis:
temp.results.pre_sleep_REM.per_period.epoch_center_seconds = (temp.results.pre_sleep_REM.per_period.epoch_start_seconds + floor(temp.results.pre_sleep_REM.per_period.durations ./ 2.0));
temp.results.post_sleep_REM.per_period.epoch_center_seconds = (temp.results.post_sleep_REM.per_period.epoch_start_seconds + floor(temp.results.post_sleep_REM.per_period.durations ./ 2.0));


% Leave in terms of the spike rates per unit (14x92 double):
temp.results.pre_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filtered.pre_sleep_REM_indicies, temp.filter_active_units);
temp.results.post_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filtered.post_sleep_REM_indicies, temp.filter_active_units);

% temp.results.pre_sleep_REM.spike_rate_per_unit: (14x92 double)
% temp.results.post_sleep_REM.spike_rate_per_unit: (9x92 double)


% Average Across all of the units
temp.results.pre_sleep_REM.spike_rate_all_units.mean = mean(temp.results.pre_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
temp.results.post_sleep_REM.spike_rate_all_units.mean = mean(temp.results.post_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
temp.results.pre_sleep_REM.spike_rate_all_units.stdDev = std(temp.results.pre_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double
temp.results.post_sleep_REM.spike_rate_all_units.stdDev = std(temp.results.post_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double

% temp.results.pre_sleep_REM.num_behavioral_periods = length(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
% temp.results.post_sleep_REM.num_behavioral_periods = length(temp.results.post_sleep_REM.spike_rate_all_units.mean);


% Compute the average across the REM sessions in each epoch
temp.results.pre_sleep_REM.baseline_spike_rate_across_all.mean = mean(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
temp.results.post_sleep_REM.baseline_spike_rate_across_all.mean = mean(temp.results.post_sleep_REM.spike_rate_all_units.mean);
temp.results.pre_sleep_REM.baseline_spike_rate_across_all.stdDev = std(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
temp.results.post_sleep_REM.baseline_spike_rate_across_all.stdDev = std(temp.results.post_sleep_REM.spike_rate_all_units.mean);


temp.plottingOptions.plottingXAxis = 'index';
% temp.plottingOptions.plottingXAxis = 'timestamp';
temp.plottingOptions.plottingYlim = [];
% temp.plottingOptions.plottingYlim = [2 4.25];
% temp.plottingOptions.plottingYlim = [0.2 1.4];

temp.plottingOptions.plottingMode = 'scatter';
% temp.plottingOptions.plottingMode = 'errorbar';
% temp.plottingOptions.plottingMode = 'distributionPlot'; % distributionPlot should display the variance across neurons


%% Error bars are across units:
figure(9);
clf
subplot(2,1,1);

if strcmpi(temp.plottingOptions.plottingXAxis, 'index')
    temp.plottingOptions.x = [1:temp.results.pre_sleep_REM.num_behavioral_periods];
else
    temp.plottingOptions.x = temp.results.pre_sleep_REM.per_period.epoch_center_seconds;
end

[h1] = fnPlotAcrossREMTesting(temp.plottingOptions.plottingMode, temp.plottingOptions.x, ...
    temp.results.pre_sleep_REM.spike_rate_all_units.mean, ...
    temp.results.pre_sleep_REM.spike_rate_all_units.stdDev, ...
    temp.results.pre_sleep_REM.spike_rate_per_unit);

title(sprintf('PRE sleep REM periods: %d', temp.results.pre_sleep_REM.num_behavioral_periods));
if strcmpi(temp.plottingOptions.plottingXAxis, 'index')
    xlabel('Filtered Trial Index')
else
    xlabel('Trial Timestamp Offset (Seconds)')
end
ylabel('mean spike rate')
if ~isempty(temp.plottingOptions.plottingYlim)
    ylim(temp.plottingOptions.plottingYlim)
end

subplot(2,1,2);

if strcmpi(temp.plottingOptions.plottingXAxis, 'index')
    temp.plottingOptions.x = [1:temp.results.post_sleep_REM.num_behavioral_periods];
else
    temp.plottingOptions.x = temp.results.post_sleep_REM.per_period.epoch_center_seconds;
end
[h2] = fnPlotAcrossREMTesting(temp.plottingOptions.plottingMode, temp.plottingOptions.x, ...
    temp.results.post_sleep_REM.spike_rate_all_units.mean, ...
    temp.results.post_sleep_REM.spike_rate_all_units.stdDev, ...
    temp.results.post_sleep_REM.spike_rate_per_unit);


title(sprintf('POST sleep REM periods: %d', temp.results.post_sleep_REM.num_behavioral_periods));
if strcmpi(temp.plottingOptions.plottingXAxis, 'index')
    xlabel('Filtered Trial Index')
else
    xlabel('Trial Timestamp Offset (Seconds)')
end
ylabel('mean spike rate')
if ~isempty(temp.plottingOptions.plottingYlim)
    ylim(temp.plottingOptions.plottingYlim)
end
sgtitle('Spike Rates - PRE vs Post Sleep REM Periods - Period Index - Pyramidal Only')
% Figure Name:
%'Spike Rates - PRE vs Post Sleep REM Periods - Period Index';
%'Spike Rates - PRE vs Post Sleep REM Periods - Timestamp Offset';


% 
% figure(10);
% clf
% subplot(2,1,1);
% [h1] = fnPlotAcrossREMTesting('bar', [1:temp.results.pre_sleep_REM.num_behavioral_periods], ...
%     temp.results.pre_sleep_REM.per_period.durations);
% 
% title(sprintf('Period Durations - PRE sleep REM periods: %d', temp.results.pre_sleep_REM.num_behavioral_periods));
% xlabel('Filtered Trial Index')
% ylabel('period duration')
% 
% subplot(2,1,2);
% [h2] = fnPlotAcrossREMTesting('bar', [1:temp.results.post_sleep_REM.num_behavioral_periods], ...
%     temp.results.post_sleep_REM.per_period.durations);
% 
% title(sprintf('Period Durations - POST sleep REM periods: %d', temp.results.post_sleep_REM.num_behavioral_periods));
% xlabel('Filtered Trial Index')
% ylabel('period duration')

% Figure Name:
%'Period Duration - PRE vs Post Sleep REM Periods';
% Conclusion: Duration is not sufficient to explain periodic behavior in REM periods for either subplot



%%%%%%% SECTION: Computer Per-Period correlation coefficients to compare them across REM periods:
%%%%%

current_binning_index = 2;
active_binning_resolution = processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = timesteps_array{current_binning_index};
temp.curr_processed = active_processing.processed_array{current_binning_index};
active_results = results_array{current_binning_index};

% Binned Spike rates per time:
% [active_binned_spike_data_matrix] = sparse(fnUnitDataCells2mat(temp.curr_processed.all.binned_spike_counts));  % 35351x126 double
[active_binned_spike_data_matrix] = fnUnitDataCells2mat(temp.curr_processed.all.binned_spike_counts);  % 35351x126 double

% active_binned_spike_data_matrix = active_binned_spike_data_matrix(:, temp.filter_active_units); % Filter down to only the active units.

temp.curr_timestamp_period_map = zeros(size(temp.curr_timestamps));

if ~isfield('by_behavioral_period', active_results) | ~isfield('pairwise_xcorrelations', active_results.by_behavioral_period) | ~isfield('xcorr_full', active_results.by_behavioral_period.pairwise_xcorrelations)
    active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_full = zeros([num_of_behavioral_state_periods, general_results.indicies.num_unique_pairs, pairwise_xcorrelations.num_lag_steps]);   
    active_results.by_behavioral_period.pairwise_xcorrelations.xcorr = zeros([num_of_behavioral_state_periods, general_results.indicies.num_unique_pairs]);   

    pairwise_xcorrelations.lag_offsets = (-processing_config.max_xcorr_lag):active_binning_resolution:processing_config.max_xcorr_lag;
    pairwise_xcorrelations.num_lag_steps = length(pairwise_xcorrelations.lag_offsets);
    %% max_xcorr_lag must be specified in terms of samples (num unit timesteps), not seconds, so we must convert by dividing by the currStepSize
    max_xcorr_lag_unit_time = processing_config.max_xcorr_lag / active_binning_resolution;

    % Loop over behavioral periods
    for behavioral_period_index = 1:num_of_behavioral_state_periods
        temp.curr_state_start = active_processing.behavioral_periods_table.epoch_start_seconds(behavioral_period_index);
        temp.curr_state_end = active_processing.behavioral_periods_table.epoch_end_seconds(behavioral_period_index);
        temp.curr_state_type = active_processing.behavioral_periods_table.type(behavioral_period_index);
        temp.curr_epoch_type = active_processing.behavioral_periods_table.behavioral_epoch(behavioral_period_index);

        % general_results.per_behavioral_state_period
        % active_processing.processed_array{1, 1}.by_state.rem.binned_spike_counts

        temp.curr_state_timesteps_start = ceil(temp.curr_state_start / active_binning_resolution) + 1;
        temp.curr_state_timesteps_end = floor(temp.curr_state_end / active_binning_resolution) + 1;

        temp.curr_timestamp_period_map(temp.curr_state_timesteps_start:temp.curr_state_timesteps_end) = behavioral_period_index;
    %     active_binned_spike_data_matrix(temp.curr_state_timesteps_start:temp.curr_state_timesteps_end, :); % Filter down to only the portion of the matrix that applies to the behavioral period.

        % Filter down to only the portion of the matrix that applies to the behavioral period.
    %     [~, ~, ~, active_results.by_behavioral_period.pairwise_xcorrelations.xcorr(behavioral_period_index, :)] = fnProcessCorrelationalMeasures(active_binned_spike_data_matrix(temp.curr_state_timesteps_start:temp.curr_state_timesteps_end, :), ...
    %         general_results.indicies, processing_config, active_binning_resolution);

        temp.output_pairwise_xcorrelations = zeros([general_results.indicies.num_unique_pairs pairwise_xcorrelations.num_lag_steps]); % 7875x181 double
        for j = 1:general_results.indicies.num_unique_pairs
           temp.output_pairwise_xcorrelations(j,:) = xcorr(active_binned_spike_data_matrix(temp.curr_state_timesteps_start:temp.curr_state_timesteps_end, general_results.indicies.unique_electrode_pairs(j,1)), ...
               active_binned_spike_data_matrix(temp.curr_state_timesteps_start:temp.curr_state_timesteps_end, general_results.indicies.unique_electrode_pairs(j,2)), ...
               max_xcorr_lag_unit_time,'normalized'); % 181x1 double
        end

        % temp.output_pairwise_xcorrelations is 7875x19 (one for each timestep)
        active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_full(behavioral_period_index, :, :) = 
        active_results.by_behavioral_period.pairwise_xcorrelations.xcorr(behavioral_period_index, :) = mean(temp.output_pairwise_xcorrelations, 2);
    %     active_processing.processed_array{1}.all.binned_spike_counts


    %     [~, active_results.by_epoch.(temp.curr_epoch_name).autocorrelations, active_results.by_epoch.(temp.curr_epoch_name).partial_autocorrelations, active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations] = fnProcessCorrelationalMeasures(temp.curr_processed.by_epoch.(temp.curr_epoch_name).binned_spike_counts, ...
    % 			general_results.indicies, processing_config, active_binning_resolution);

        fprintf('behavioral state progress: %d/%d\n', behavioral_period_index, num_of_behavioral_state_periods);

    end
    
    results_array{current_binning_index} = active_results;
    

else
    fprintf('active_results.by_behavioral_period.pairwise_xcorrelations.xcorr already exists. Not overwritting.\n');
end % end if ~isfield pairwise_xcorrelations

fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);
active_results.by_behavioral_period.pairwise_xcorrelations.xcorr


%% Display the Correlational Results:
xcorr_fig = figure(14);
plotting_options.plotMode = 'xcorr';
[temp.fig, temp.h] = fnPhoPlotCorrelationalResults(active_processing, general_results, active_results, plotting_options, xcorr_fig);
    



% active_processing.spikes.behavioral_states
fprintf('AcrossREMTesting done.\n');

function [h] = fnPlotAcrossREMTesting(mode, v1, v2, v3, v4)

    if strcmpi(mode, 'errorbar')
        h = errorbar(v1, ...
            v2, ...
            v3);

    elseif strcmpi(mode, 'scatter')
        h = scatter(v1, ...
            v2);

    elseif strcmpi(mode, 'distributionPlot')
        h = distributionPlot(v4'); % defaults 

    elseif strcmpi(mode, 'bar')
        h = bar(v1, v2);
        
    else
       error('Invalid mode input!') 
    end
end