% Iterate across REM states

clear temp;

num_of_behavioral_state_periods = height(active_processing.behavioral_periods_table);
% 
% general_results.per_behavioral_state_period.spike_rate_per_unit

% All units with qualities from 1-4 are pyramidal.
% The higher the number in this range, the higher is the contamination, so 1 and 2 are well-separated pyramidal units and if your analysis is not much sensitive to contaminations you can consider 3 and 4 as well. For my analysis, I considered 1 to 3. 8 and 9 are interneurons.
% 

% general_results.per_behavioral_state_period.num_spikes_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double
% general_results.per_behavioral_state_period.spike_rate_per_unit = zeros([num_of_behavioral_state_periods num_of_electrodes]); %% num_results: 668x126 double


% for unit_index = 1:num_of_electrodes
% 
%     % temp.edges = unique(active_processing.spikes.behavioral_duration_indicies{unit_index});
%     general_results.counts = histc(active_processing.spikes.behavioral_duration_indicies{unit_index}(:), general_results.edges);
%     general_results.per_behavioral_state_period.num_spikes_per_unit(:, unit_index) = general_results.counts;
%     general_results.per_behavioral_state_period.spike_rate_per_unit(:, unit_index) = general_results.counts ./ active_processing.behavioral_periods_table.duration;
%     % active_processing.spikes.behavioral_duration_indicies{unit_index}
% 
%     % behavioral_periods_table.duration
% 
% end


% if processing_config.showOnlyAlwaysStableCells
%     isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
%     numAlwaysStableCells = sum(isAlwaysStable, 'all');
%     temp.filter_active_units = isAlwaysStable;
%     
% else
%     temp.filter_active_units = logical(ones([num_of_electrodes 1]));
% end

temp.filter_included_cell_types = {};
temp.filter_maximum_included_contamination_level = {};

[temp.testing.filter_active_units] = fnFilterUnitsWithCriteria(active_processing, processing_config.showOnlyAlwaysStableCells, temp.filter_included_cell_types, ...
    temp.filter_maximum_included_contamination_level);

all(temp.testing.filter_active_units == temp.filter_active_units)

%% Filter by Epoch:
active_processing.spikes


%% What I was looking for before. Can filter by the epoch and state indicies and interest and collapse across trials
% general_results.GroupedByState.groups = findgroups(active_processing.behavioral_periods_table.type);
% general_results.GroupedByEpoch.groups = findgroups(active_processing.behavioral_periods_table.behavioral_epoch);

temp.filter_states = {'rem'};
temp.filter_epochs = {'pre_sleep', 'post_sleep'};
temp.indicies.states.REM = (active_processing.behavioral_periods_table.type == 'rem');
temp.indicies.epochs.pre_sleep = (active_processing.behavioral_periods_table.behavioral_epoch == 'pre_sleep');
temp.indicies.epochs.post_sleep = (active_processing.behavioral_periods_table.behavioral_epoch == 'post_sleep');

% Element-wise multiplication acting as a logical AND
temp.filtered.pre_sleep_REM_indicies = logical(temp.indicies.states.REM .* temp.indicies.epochs.pre_sleep); % 668x1
temp.filtered.post_sleep_REM_indicies = logical(temp.indicies.states.REM .* temp.indicies.epochs.post_sleep); % 668x1




[temp.testing.filtered.pre_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'pre_sleep'}, {'rem'});
[temp.testing.filtered.post_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, {'post_sleep'}, {'rem'});

% Testing new functions:
all(temp.testing.filtered.pre_sleep_REM_indicies == temp.filtered.pre_sleep_REM_indicies)
all(temp.testing.filtered.post_sleep_REM_indicies == temp.filtered.post_sleep_REM_indicies)


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
temp.plottingOptions.plottingXAxis = 'timestamp';
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
ylim([2 4.25])
xlim([2000 

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
ylim([2 4.25])
% Figure Name:
%'Spike Rates - PRE vs Post Sleep REM Periods - Period Index';
%'Spike Rates - PRE vs Post Sleep REM Periods - Timestamp Offset';



figure(10);
clf
subplot(2,1,1);
[h1] = fnPlotAcrossREMTesting('bar', [1:temp.results.pre_sleep_REM.num_behavioral_periods], ...
    temp.results.pre_sleep_REM.per_period.durations);

title(sprintf('Period Durations - PRE sleep REM periods: %d', temp.results.pre_sleep_REM.num_behavioral_periods));
xlabel('Filtered Trial Index')
ylabel('period duration')

subplot(2,1,2);
[h2] = fnPlotAcrossREMTesting('bar', [1:temp.results.post_sleep_REM.num_behavioral_periods], ...
    temp.results.post_sleep_REM.per_period.durations);

title(sprintf('Period Durations - POST sleep REM periods: %d', temp.results.post_sleep_REM.num_behavioral_periods));
xlabel('Filtered Trial Index')
ylabel('period duration')

% Figure Name:
%'Period Duration - PRE vs Post Sleep REM Periods';
% Conclusion: Duration is not sufficient to explain periodic behavior in REM periods for either subplot




    
% 
% % Loop over behavioral periods
% for state_index = 1:num_of_behavioral_state_periods
%     temp.curr_state_start = active_processing.behavioral_periods_table.epoch_start_seconds(state_index);
%     temp.curr_state_end = active_processing.behavioral_periods_table.epoch_end_seconds(state_index);
%     temp.curr_state_type = active_processing.behavioral_periods_table.type(state_index);
%     temp.curr_epoch_type = active_processing.behavioral_periods_table.behavioral_epoch(state_index);
%     
%     if strcmpi(temp.curr_state_type, 'REM')
%         % It's a REM state
%         temp.curr_state_spikes = cell(num_of_electrodes, 1);
%         % Extract the spike train for each electrode
%         for electrode_index = 1:num_of_electrodes
%             % Convert spike times to relative to expt start and scale to seconds.
%             temp.curr_electrode_spikes = active_processing.spikes.time{electrode_index};
%             % Get the spike times that belong to this particular state.
% %             temp.curr_state_spikes_idx{electrode_index} = find((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
%             temp.curr_state_spikes{electrode_index} = temp.curr_electrode_spikes((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
% 
%             
%         end
% 
%     
%         
%         % Compute Period Average Firing Rate
%         
%         
%         % Compute Period Std Dev/ Variation
%         
%         
%     end
%     fprintf('behavioral state progress: %d/%d\n', state_index, temp.curr_num_of_behavioral_states);
%     
% 
% end

% active_processing.spikes.behavioral_states


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