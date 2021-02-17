% Iterate across REM states

clear temp;

num_of_behavioral_state_periods = height(active_processing.behavioral_periods_table);
% 
% general_results.per_behavioral_state_period.spike_rate_per_unit


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


if processing_config.showOnlyAlwaysStableCells
    isAlwaysStable = active_processing.spikes.isAlwaysStable; % 126x1
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    temp.filter_active_units = isAlwaysStable;
    
else
    temp.filter_active_units = logical(ones([num_of_electrodes 1]));
end



%% Filter by Epoch:



%% What I was looking for before. Can filter by the epoch and state indicies and interest and collapse across trials
general_results.GroupedByState.groups = findgroups(active_processing.behavioral_periods_table.type);
general_results.GroupedByEpoch.groups = findgroups(active_processing.behavioral_periods_table.behavioral_epoch);

temp.filter_states = {'rem'};
temp.filter_epochs = {'pre_sleep', 'post_sleep'};
temp.indicies.states.REM = (active_processing.behavioral_periods_table.type == 'rem');
temp.indicies.epochs.pre_sleep = (active_processing.behavioral_periods_table.behavioral_epoch == 'pre_sleep');
temp.indicies.epochs.post_sleep = (active_processing.behavioral_periods_table.behavioral_epoch == 'post_sleep');

% Element-wise multiplication acting as a logical AND
temp.filtered.pre_sleep_REM_indicies = logical(temp.indicies.states.REM .* temp.indicies.epochs.pre_sleep); % 668x1
temp.filtered.post_sleep_REM_indicies = logical(temp.indicies.states.REM .* temp.indicies.epochs.post_sleep); % 668x1


sum(temp.filtered.pre_sleep_REM_indicies,'all')
sum(temp.filtered.post_sleep_REM_indicies,'all')



% Alternative: splitapply workflow:
% splitapply(@mean,Height,G)

% The number of spikes per unit
% general_results.per_behavioral_state_period.num_spikes_per_unit(temp.filter_active_units, temp.filtered.pre_sleep_REM_indicies); % 668x126
% general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filter_active_units, temp.filtered.pre_sleep_REM_indicies); % 668x126


% Leave in terms of the spike rates per unit (14x92 double):
temp.results.pre_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filtered.pre_sleep_REM_indicies, temp.filter_active_units);
temp.results.post_sleep_REM.spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(temp.filtered.post_sleep_REM_indicies, temp.filter_active_units);

% Average Across all of the units
temp.results.pre_sleep_REM.spike_rate_all_units.mean = mean(temp.results.pre_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
temp.results.post_sleep_REM.spike_rate_all_units.mean = mean(temp.results.post_sleep_REM.spike_rate_per_unit, 2); % 14x1 double
temp.results.pre_sleep_REM.spike_rate_all_units.stdDev = std(temp.results.pre_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double
temp.results.post_sleep_REM.spike_rate_all_units.stdDev = std(temp.results.post_sleep_REM.spike_rate_per_unit, 0, 2); % 14x1 double

temp.results.pre_sleep_REM.num_behavioral_periods = length(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
temp.results.post_sleep_REM.num_behavioral_periods = length(temp.results.post_sleep_REM.spike_rate_all_units.mean);


% Compute the average across the REM sessions in each epoch
temp.results.pre_sleep_REM.baseline_spike_rate_across_all.mean = mean(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
temp.results.post_sleep_REM.baseline_spike_rate_across_all.mean = mean(temp.results.post_sleep_REM.spike_rate_all_units.mean);
temp.results.pre_sleep_REM.baseline_spike_rate_across_all.stdDev = std(temp.results.pre_sleep_REM.spike_rate_all_units.mean);
temp.results.post_sleep_REM.baseline_spike_rate_across_all.stdDev = std(temp.results.post_sleep_REM.spike_rate_all_units.mean);





%% Error bars are across units:
figure(9);
subplot(2,1,1);
% h1 = scatter([1:temp.results.pre_sleep_REM.num_behavioral_periods], ...
%         temp.results.pre_sleep_REM.spike_rate_all_units.mean);
h1 = errorbar([1:temp.results.pre_sleep_REM.num_behavioral_periods], ...
        temp.results.pre_sleep_REM.spike_rate_all_units.mean, ...
        temp.results.pre_sleep_REM.spike_rate_all_units.stdDev);
title(sprintf('PRE sleep REM periods: %d', temp.results.pre_sleep_REM.num_behavioral_periods));
xlabel('Filtered Trial Index')
ylabel('mean spike rate')

subplot(2,1,2);
% h2 = scatter([1:temp.results.post_sleep_REM.num_behavioral_periods], ...
%         temp.results.post_sleep_REM.spike_rate_all_units.mean);
    
    
h2 = errorbar([1:temp.results.post_sleep_REM.num_behavioral_periods], ...
        temp.results.post_sleep_REM.spike_rate_all_units.mean, ...
        temp.results.post_sleep_REM.spike_rate_all_units.stdDev);
    
    
% [h2] = fnPlotAcrossREMTesting('errorbar', [1:temp.results.post_sleep_REM.num_behavioral_periods], ...
%     temp.results.post_sleep_REM.spike_rate_all_units.mean, ...
%     temp.results.post_sleep_REM.spike_rate_all_units.stdDev);


title(sprintf('POST sleep REM periods: %d', temp.results.post_sleep_REM.num_behavioral_periods));
xlabel('Filtered Trial Index')
ylabel('mean spike rate')






    
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


% function [h] = fnPlotAcrossREMTesting(mode, x, y, z)
% 
%     if strcmpi(mode, 'errorbar')
%         h = errorbar(x, ...
%             y, ...
%             z);
% 
%     elseif strcmpi(mode, 'scatter')
%         h = scatter(x, ...
%             y);
% 
%     elseif strcmpi(mode, 'distributionPlot')
%         h = distributionPlot(x); % defaults 
% 
%     else
%        error('Invalid mode input!') 
%     end
% end