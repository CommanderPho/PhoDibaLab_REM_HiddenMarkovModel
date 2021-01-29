% READ
clear all;

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);


% microseconds (10^6): 1000000
% nanoseconds (10^9): 1000000000
conversion_factor = (10^6);

% ,'y'
source_data.spikes = load(fullfile(data_config.source_root_path,'wake-spikes.mat'), 'spikes').spikes;
source_data.behavior = load(fullfile(data_config.source_root_path, 'wake-behavior.mat'), 'behavior').behavior;

% Process one of the experiments: 
processing_config.active_expt.name = 'RoyMaze1';

processing_config.active_expt.spikes_list = source_data.spikes.(processing_config.active_expt.name);
processing_config.active_expt.behavior_list = source_data.behavior.(processing_config.active_expt.name).list;


behavioral_state_names = source_data.behavior.(processing_config.active_expt.name).name;

behavioral_epoch_names = {'pre_sleep', 'track', 'post_sleep'};

active_processing.earliest_start_timestamp = Inf; % Set the initial timestamp to be infinity, so that any timestamp it's compared to will be less than it

for i = 1:length(behavioral_epoch_names)
    temp.curr_epoch_name = behavioral_epoch_names{i};
    temp.curr_epoch_start_stop = source_data.behavior.(processing_config.active_expt.name).time(i,:);
    
    % microseconds (10^6): 1000000
    % nanoseconds (10^9): 1000000000
    active_processing.earliest_start_timestamp = min(active_processing.earliest_start_timestamp, (temp.curr_epoch_start_stop(1) / conversion_factor));
    temp.curr_epoch_duration = ((temp.curr_epoch_start_stop(2) - temp.curr_epoch_start_stop(1)) / conversion_factor);
    disp(temp.curr_epoch_duration)
    active_processing.behavioral_epochs.(temp.curr_epoch_name) = [temp.curr_epoch_start_stop temp.curr_epoch_duration];
end


fprintf('earliest recording timestamp is %d\n', active_processing.earliest_start_timestamp)
% curr_activity_table = array2table(processing_config.active_expt.behavior_list, ...
%         'VariableNames',{'epoch_start_time', 'epoch_end_time', 'type'});
    

temp.durations = ((processing_config.active_expt.behavior_list(:,2) - processing_config.active_expt.behavior_list(:,1)) ./ conversion_factor);


% Datetime based:
curr_activity_table = table(((processing_config.active_expt.behavior_list(:,1) ./ conversion_factor) - active_processing.earliest_start_timestamp),  ...
    ((processing_config.active_expt.behavior_list(:,2) ./ conversion_factor) - active_processing.earliest_start_timestamp),  ...
    temp.durations, ...
    categorical(processing_config.active_expt.behavior_list(:,3), [1:length(behavioral_state_names)], behavioral_state_names),  ...
    'VariableNames',{'epoch_start_seconds', 'epoch_end_seconds', 'duration', 'type'});


% % Datetime based:
% curr_activity_table = table(convert_date(processing_config.active_expt.behavior_list(:,1), conversion_factor),  ...
%     convert_date(processing_config.active_expt.behavior_list(:,2), conversion_factor),  ...
%     temp.durations, ...
%     categorical(processing_config.active_expt.behavior_list(:,3), [1:length(behavioral_state_names)], behavioral_state_names),  ...
%     'VariableNames',{'epoch_start_time', 'epoch_end_time', 'duration', 'type'});
    
% TimeWindows
%     FileLoaderOpenEphys
function [date_out] = convert_date(aDateTimestamp, conversion_factor)
%% Raw Timestamps Mode:
return_raw_timestamps = false;
    timestamps_seconds = aDateTimestamp ./ conversion_factor;
if return_raw_timestamps
   date_out = timestamps_seconds
else
    %     date_out = datetime(aDateTimestamp,'ConvertFrom','epochtime','TicksPerSecond',1000);
%     date_out = datetime(timestamps_seconds,'ConvertFrom','epochtime','Epoch','2000-01-01','TicksPerSecond', conversion_factor); % Prety close
    date_out = datetime(timestamps_seconds,'ConvertFrom','epochtime');
    %     date_out = datetime(aDateTimestamp,'ConvertFrom','datenum');
end % end if
end % end function convert_date


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