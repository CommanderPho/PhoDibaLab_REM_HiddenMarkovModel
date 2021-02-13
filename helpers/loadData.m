function [active_processing, source_data] = loadData(data_config, processing_config)
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here
% ,'y'
% source_data.spikes = load(fullfile(data_config.source_root_path,'wake-spikes.mat'), 'spikes').spikes;
% source_data.behavior = load(fullfile(data_config.source_root_path, 'wake-behavior.mat'), 'behavior').behavior;


source_data.variableList = {'spikes','behavior','position','speed','basics','ripple'};
for variableIndex = 1:length(source_data.variableList)
    curr_variable_name = source_data.variableList{variableIndex};
    source_data.(curr_variable_name) = load(fullfile(data_config.source_root_path, sprintf('wake-%s.mat', curr_variable_name)), ...
        curr_variable_name).(curr_variable_name);
end

% Begin Pre-processing
processing_config.active_expt.spikes_list = source_data.spikes.(processing_config.active_expt.name);
processing_config.active_expt.behavior_list = source_data.behavior.(processing_config.active_expt.name).list;
processing_config.active_expt.position_list = source_data.position.(processing_config.active_expt.name); % x, y, t
processing_config.active_expt.speed_list = source_data.speed.(processing_config.active_expt.name); % t, v;

active_processing.behavioral_state_names = source_data.behavior.(processing_config.active_expt.name).name;


%% Determine the experiment start timestamp to convert the timestamps into experiment-start-relative durations.
active_processing.earliest_start_timestamp = Inf; % Set the initial timestamp to be infinity, so that any timestamp it's compared to will be less than it

for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    temp.curr_epoch_start_stop_absolute = source_data.behavior.(processing_config.active_expt.name).time(i,:) ./ data_config.conversion_factor; % Absolute
    
    active_processing.earliest_start_timestamp = min(active_processing.earliest_start_timestamp, temp.curr_epoch_start_stop_absolute(1));
    temp.curr_epoch_duration = (temp.curr_epoch_start_stop_absolute(2) - temp.curr_epoch_start_stop_absolute(1));
    active_processing.behavioral_epochs.(temp.curr_epoch_name) = [temp.curr_epoch_start_stop_absolute temp.curr_epoch_duration];
end

fprintf('earliest recording timestamp is %d\n', active_processing.earliest_start_timestamp)

temp.behavioral_epochs_flattened = [active_processing.behavioral_epochs.pre_sleep; active_processing.behavioral_epochs.track; active_processing.behavioral_epochs.post_sleep];
temp.behavioral_epochs_flattened = [temp.behavioral_epochs_flattened(:,1:2), (temp.behavioral_epochs_flattened(:,1:2) - active_processing.earliest_start_timestamp), temp.behavioral_epochs_flattened(:,3)];

active_processing.behavioral_epochs = array2table(temp.behavioral_epochs_flattened, ...
    'RowNames', {'pre_sleep','track','post_sleep'}, ...
    'VariableNames',{'start_seconds_absolute', 'end_seconds_absolute', 'start_seconds', 'end_seconds', 'duration'});

% We have active_processing.position_table and active_processing.speed_table
active_processing.position_table = table(((processing_config.active_expt.position_list.t ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    processing_config.active_expt.position_list.x, ...
    processing_config.active_expt.position_list.y, ...
    'VariableNames',{'timestamp', 'x', 'y'});

active_processing.speed_table =  table(((processing_config.active_expt.speed_list.t ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    processing_config.active_expt.speed_list.v, ...
    'VariableNames',{'timestamp', 'speed'});


temp.durations = ((processing_config.active_expt.behavior_list(:,2) - processing_config.active_expt.behavior_list(:,1)) ./ data_config.conversion_factor);

% Seconds relative to start of recording based:
active_processing.behavioral_periods_table = table(((processing_config.active_expt.behavior_list(:,1) ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    ((processing_config.active_expt.behavior_list(:,2) ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    temp.durations, ...
    categorical(processing_config.active_expt.behavior_list(:,3), [1:length(active_processing.behavioral_state_names)], active_processing.behavioral_state_names),  ...
    'VariableNames',{'epoch_start_seconds', 'epoch_end_seconds', 'duration', 'type'});


%% For each behavioral period in behavioral_periods_table:
% we want to be able to extract:
%% Any spikes that occur within that period
%% the experimental_phase it belongs in {pre_sleep, track, post_sleep}

temp.edges = [active_processing.behavioral_epochs.start_seconds(1), active_processing.behavioral_epochs.start_seconds(2), active_processing.behavioral_epochs.start_seconds(3), active_processing.behavioral_epochs.end_seconds(3)];
temp.behavior_types = discretize(active_processing.behavioral_periods_table.epoch_start_seconds, temp.edges);
% Add the categorical data to the table
active_processing.behavioral_periods_table.behavioral_epoch = categorical(temp.behavior_types, [1:length(data_config.behavioral_epoch_names)], data_config.behavioral_epoch_names);

% Build Spikes table:
active_processing.spikes = struct2table(processing_config.active_expt.spikes_list);

% Convert the first column (of timestamp offsets) to relative offsets into the experiment
active_processing.spikes.time = cellfun((@(timestamps) (timestamps ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp), ...
        active_processing.spikes.time, ...
        'UniformOutput', false);


    
% Count the spikes for each unit
active_processing.spikes.num_spikes = cellfun((@(timestamps) length(timestamps)), ...
        active_processing.spikes.time);
    
   

    
%% Preprocessing:
% num_of_electrodes = height(active_processing.spikes);







%% Loading complete.


end

