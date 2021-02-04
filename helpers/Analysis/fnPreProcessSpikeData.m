function [output] = fnPreProcessSpikeData(active_processing, data_config, num_of_electrodes)
%fnPreProcessSpikeData Summary of this function goes here
%   Detailed explanation goes here

%% Processing:
output.all.spike_data = cell([1, num_of_electrodes]);

for electrode_index = 1:num_of_electrodes
    % Convert spike times to relative to expt start and scale to seconds. 
    fprintf('electrode progress: %d/%d\n', electrode_index, num_of_electrodes);
    
    temp.curr_timetable = timetable(seconds(active_processing.spikes.time{electrode_index}'), active_processing.spikes.behavioral_epoch{electrode_index}, active_processing.spikes.behavioral_states{electrode_index}, active_processing.spikes.behavioral_duration_indicies{electrode_index}, ...
        'VariableNames',{'behavioral_epoch','behavioral_state','behavioral_period_index'});
    
    output.all.spike_data{electrode_index} = temp.curr_timetable;
   
    %% Split based on experiment epoch:
    for i = 1:length(data_config.behavioral_epoch_names)
        temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
        output.by_epoch.(temp.curr_epoch_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_epoch == temp.curr_epoch_name), :);
    end
    
    %% Split based on behavioral state:
    for i = 1:length(active_processing.behavioral_state_names)
        temp.curr_state_name =  active_processing.behavioral_state_names{i};
        output.by_state.(temp.curr_state_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_state == temp.curr_state_name), :);
    end

end % end for


[output] = fnBinSpikeData(active_processing, data_config, timesteps)




end

