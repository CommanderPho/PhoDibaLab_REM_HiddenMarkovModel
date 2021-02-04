function [output] = fnBinSpikeData(active_processing, timesteps)
%FNBINSPIKEDATA Summary of this function goes here
%   Detailed explanation goes here
output.all.smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.all.spike_data, 'UniformOutput', false);

%% Split based on experiment epoch:
for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    output.by_epoch.(temp.curr_epoch_name).smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.by_epoch.(temp.curr_epoch_name).spike_data, 'UniformOutput', false);
end

%% Split based on behavioral state:
for i = 1:length(active_processing.behavioral_state_names)
    temp.curr_state_name =  active_processing.behavioral_state_names{i};
    output.by_state.(temp.curr_state_name).smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.by_state.(temp.curr_state_name).spike_data, 'UniformOutput', false);
end

% output.all.smoothed_spike_data = cell([num_of_electrodes, 1]);
% for electrode_index = 1:num_of_electrodes
%     % Gaussian smoothing of spike data:
%     fprintf('progress: %d/%d\n', electrode_index, num_of_electrodes);
%     % This retimed version is pretty slow: > 30 seconds execution time.
%     output.normalized_spike_data{electrode_index} = retime(active_processing.processed.all.spike_data{electrode_index},'regular','count','TimeStep', seconds(1));    
%     output.all.smoothed_spike_data{electrode_index} = smoothdata(output.normalized_spike_data{electrode_index},'gaussian', seconds(1.5));
% 
%     histcounts(active_processing.processed.all.spike_data{electrode_index}, timesteps);
%     
%     
%     %% Split based on experiment epoch:
%     for i = 1:length(data_config.behavioral_epoch_names)
%         temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
%         output.by_epoch.(temp.curr_epoch_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_epoch == temp.curr_epoch_name), :);
%     end
%     
%     %% Split based on behavioral state:
%     for i = 1:length(active_processing.behavioral_state_names)
%         temp.curr_state_name =  active_processing.behavioral_state_names{i};
%         output.by_state.(temp.curr_state_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_state == temp.curr_state_name), :);
%     end
%     
%     
% end


end

