% READ
addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));
clear all;

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.output.intermediate_file_name = 'PhoIntermediate.mat';


data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);
data_config.output.intermediate_file_path = fullfile(data_config.source_root_path, data_config.output.intermediate_file_name);

% microseconds (10^6): 1000000
% nanoseconds (10^9): 1000000000
data_config.conversion_factor = (10^6);
data_config.behavioral_epoch_names = {'pre_sleep', 'track', 'post_sleep'};

% Process one of the experiments: 
processing_config.active_expt.name = 'RoyMaze1';

if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    [active_processing, source_data] = loadData(data_config, processing_config);
end

%% Preprocessing:

% Each entry in active_processing.spikes has a variable number of double entries, indicating the relative offset (in seconds) the spike occured for each unit.
num_of_electrodes = height(active_processing.spikes);

%% For each behavioral period in curr_activity_table:
% we want to be able to extract:
%% Any spikes that occur within that period
%% the experimental_phase it belongs in {pre_sleep, track, post_sleep}

%% Partition spikes based on behavioral state:
temp.curr_num_of_behavioral_states = height(active_processing.curr_activity_table);

temp.spikes_behavioral_states = cell([num_of_electrodes, 1]);
temp.spikes_behavioral_epoch = cell([num_of_electrodes, 1]);

% Pre-allocation: Loop over electrodes
for electrode_index = 1:num_of_electrodes
    % Convert spike times to relative to expt start and scale to seconds. 
    temp.spikes_behavioral_states{electrode_index} = categorical(ones([active_processing.spikes.num_spikes(electrode_index), 1]), [1:length(active_processing.behavioral_state_names)], active_processing.behavioral_state_names);
    temp.spikes_behavioral_epoch{electrode_index} = categorical(ones([active_processing.spikes.num_spikes(electrode_index), 1]), [1:length(data_config.behavioral_epoch_names)], data_config.behavioral_epoch_names);
    active_processing.spikes.behavioral_duration_indicies{electrode_index} = zeros([active_processing.spikes.num_spikes(electrode_index), 1]); % to store the index of the corresponding behavioral state the spike belongs to
end


% Loop over behavioral activities
for state_index = 1:temp.curr_num_of_behavioral_states
    temp.curr_state_start = active_processing.curr_activity_table.epoch_start_seconds(state_index);
    temp.curr_state_end = active_processing.curr_activity_table.epoch_end_seconds(state_index);
    temp.curr_state_type = active_processing.curr_activity_table.type(state_index);
    temp.curr_epoch_type = active_processing.curr_activity_table.behavioral_epoch(state_index);
    
    fprintf('progress: %d/%d\n', state_index, temp.curr_num_of_behavioral_states);
    
    temp.curr_state_spikes = cell(num_of_electrodes, 1);
    % Extract the spike train for each electrode
    for electrode_index = 1:num_of_electrodes
        % Convert spike times to relative to expt start and scale to seconds.
        temp.curr_electrode_spikes = active_processing.spikes.time{electrode_index};
        % Get the spike times that belong to this particular state.
        temp.curr_state_spikes_idx{electrode_index} = find((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
        
        temp.curr_state_spikes{electrode_index} = temp.curr_electrode_spikes((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
        
        temp.spikes_behavioral_states{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = temp.curr_state_type;
        temp.spikes_behavioral_epoch{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = temp.curr_epoch_type;
        
        active_processing.spikes.behavioral_duration_indicies{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = state_index;
    end

end

active_processing.spikes.behavioral_states = temp.spikes_behavioral_states;
active_processing.spikes.behavioral_epoch = temp.spikes_behavioral_epoch;


%% Processing:
active_processing.processed.spike_data = cell([num_of_electrodes, 1]);

for electrode_index = 1:num_of_electrodes
    % Convert spike times to relative to expt start and scale to seconds. 
    fprintf('electrode progress: %d/%d\n', electrode_index, num_of_electrodes);
    
    temp.curr_timetable = timetable(seconds(active_processing.spikes.time{electrode_index}'), active_processing.spikes.behavioral_epoch{electrode_index}, active_processing.spikes.behavioral_states{electrode_index}, active_processing.spikes.behavioral_duration_indicies{electrode_index}, ...
        'VariableNames',{'behavioral_epoch','behavioral_state','behavioral_period_index'});
    
    active_processing.processed.spike_data{electrode_index} = temp.curr_timetable;
   
    %% Split based on experiment epoch:
    for i = 1:length(data_config.behavioral_epoch_names)
        temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
        active_processing.processed.by_epoch.(temp.curr_epoch_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_epoch == temp.curr_epoch_name), :);
    end
    
    %% Split based on behavioral state:
    for i = 1:length(active_processing.behavioral_state_names)
        temp.curr_state_name =  active_processing.behavioral_state_names{i};
        active_processing.processed.by_state.(temp.curr_state_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_state == temp.curr_state_name), :);
    end

end

fprintf('writing out to %s...\n', data_config.output.intermediate_file_path);
save(data_config.output.intermediate_file_path, 'active_processing', 'data_config', 'processing_config', 'num_of_electrodes', 'source_data');
fprintf('done.\n');


% Smooth with a gaussian window
% Perform Gaussian Binning: (sigma = 1.5 sec)


timesteps = seconds(active_processing.behavioral_epochs.start_seconds(1):active_processing.behavioral_epochs.end_seconds(end));


active_processing.processed.smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.spike_data, 'UniformOutput', false);

%% Split based on experiment epoch:
for i = 1:length(data_config.behavioral_epoch_names)
    temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
    active_processing.processed.by_epoch.(temp.curr_epoch_name).smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.by_epoch.(temp.curr_epoch_name).spike_data, 'UniformOutput', false);
end

%% Split based on behavioral state:
for i = 1:length(active_processing.behavioral_state_names)
    temp.curr_state_name =  active_processing.behavioral_state_names{i};
    active_processing.processed.by_state.(temp.curr_state_name).smoothed_spike_data = cellfun((@(ttable) histcounts(ttable.Time, timesteps)'), active_processing.processed.by_state.(temp.curr_state_name).spike_data, 'UniformOutput', false);
end

fprintf('writing out to %s...\n', data_config.output.intermediate_file_path);
save(data_config.output.intermediate_file_path, 'active_processing', 'data_config', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps');
fprintf('done.\n');

% active_processing.processed.smoothed_spike_data = cell([num_of_electrodes, 1]);
% for electrode_index = 1:num_of_electrodes
%     % Gaussian smoothing of spike data:
%     fprintf('progress: %d/%d\n', electrode_index, num_of_electrodes);
%     % This retimed version is pretty slow: > 30 seconds execution time.
%     active_processing.processed.normalized_spike_data{electrode_index} = retime(active_processing.processed.spike_data{electrode_index},'regular','count','TimeStep', seconds(1));    
%     active_processing.processed.smoothed_spike_data{electrode_index} = smoothdata(active_processing.processed.normalized_spike_data{electrode_index},'gaussian', seconds(1.5));
% 
%     histcounts(active_processing.processed.spike_data{electrode_index}, timesteps);
%     
%     
%     %% Split based on experiment epoch:
%     for i = 1:length(data_config.behavioral_epoch_names)
%         temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
%         active_processing.processed.by_epoch.(temp.curr_epoch_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_epoch == temp.curr_epoch_name), :);
%     end
%     
%     %% Split based on behavioral state:
%     for i = 1:length(active_processing.behavioral_state_names)
%         temp.curr_state_name =  active_processing.behavioral_state_names{i};
%         active_processing.processed.by_state.(temp.curr_state_name).spike_data{electrode_index} = temp.curr_timetable((temp.curr_timetable.behavioral_state == temp.curr_state_name), :);
%     end
%     
%     
% end


