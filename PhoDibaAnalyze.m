% READ
addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));
clear all;

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);

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

spikes_behavioral_states = cell([num_of_electrodes, 1]);
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
active_processing.spikes.behavioral_epoch = temp.spikes_behavioral_epoch';








% [xPoints, yPoints] = plotSpikeRaster({active_processing.spikes(1).time});
% plot(xPoints, yPoints)

    
% TimeWindows
%     FileLoaderOpenEphys


%     table(indexArray', final_data_explorer_obj.uniqueComps, final_data_explorer_obj.cellROIIndex_mapper.compIDsArray, final_data_explorer_obj.multiSessionCellRoi_CompListIndicies,...
% 		manualRoiFilteringResults.final_is_Excluded, manualRoiFilteringResults.final_quality_of_tuning, ...
% 'datetime'

% 
% curr_activity_timetable = timetable(all_actigraphy_files_output_data{i}.concatenatedTimestamps, ...
%             all_actigraphy_files_output_data{i}.concatenatedResults.num_changed_pixels, all_actigraphy_files_output_data{i}.concatenatedResults.total_sum_changed_pixel_value, ...
%             'VariableNames',{'NumChangedPixels', 'TotalSumChangedPixelValues'});
%         
%         
%         
% processing_config.active_expt.spikes_list

