function [active_processing, source_data] = loadData(data_config, processing_config)
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here
% ,'y'
source_data.spikes = load(fullfile(data_config.source_root_path,'wake-spikes.mat'), 'spikes').spikes;
source_data.behavior = load(fullfile(data_config.source_root_path, 'wake-behavior.mat'), 'behavior').behavior;

% Begin Pre-processing
processing_config.active_expt.spikes_list = source_data.spikes.(processing_config.active_expt.name);
processing_config.active_expt.behavior_list = source_data.behavior.(processing_config.active_expt.name).list;

active_processing.behavioral_state_names = source_data.behavior.(processing_config.active_expt.name).name;

active_processing.earliest_start_timestamp = Inf; % Set the initial timestamp to be infinity, so that any timestamp it's compared to will be less than it

for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    temp.curr_epoch_start_stop = source_data.behavior.(processing_config.active_expt.name).time(i,:);
    active_processing.earliest_start_timestamp = min(active_processing.earliest_start_timestamp, (temp.curr_epoch_start_stop(1) / data_config.conversion_factor));
    temp.curr_epoch_duration = ((temp.curr_epoch_start_stop(2) - temp.curr_epoch_start_stop(1)) / data_config.conversion_factor);
    active_processing.behavioral_epochs.(temp.curr_epoch_name) = [temp.curr_epoch_start_stop temp.curr_epoch_duration];
end


fprintf('earliest recording timestamp is %d\n', active_processing.earliest_start_timestamp)

temp.durations = ((processing_config.active_expt.behavior_list(:,2) - processing_config.active_expt.behavior_list(:,1)) ./ data_config.conversion_factor);

% Seconds relative to start of recording based:
active_processing.curr_activity_table = table(((processing_config.active_expt.behavior_list(:,1) ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    ((processing_config.active_expt.behavior_list(:,2) ./ data_config.conversion_factor) - active_processing.earliest_start_timestamp),  ...
    temp.durations, ...
    categorical(processing_config.active_expt.behavior_list(:,3), [1:length(active_processing.behavioral_state_names)], active_processing.behavioral_state_names),  ...
    'VariableNames',{'epoch_start_seconds', 'epoch_end_seconds', 'duration', 'type'});

active_processing.spikes = struct2table(processing_config.active_expt.spikes_list);
%% Loading complete.


end

