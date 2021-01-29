% READ
clear all;

data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);


% ,'y'
source_data.spikes = load(fullfile(data_config.source_root_path,'wake-spikes.mat'), 'spikes').spikes;
source_data.behavior = load(fullfile(data_config.source_root_path, 'wake-behavior.mat'), 'behavior').behavior;

% Process one of the experiments: 
processing_config.active_expt.name = 'RoyMaze1';

processing_config.active_expt.spikes_list = source_data.spikes.(processing_config.active_expt.name);
processing_config.active_expt.behavior_list = source_data.behavior.(processing_config.active_expt.name).list;


behavioral_state_names = source_data.behavior.(processing_config.active_expt.name).name;

behavioral_epoch_names = {'pre_sleep', 'track', 'post_sleep'};

for i = 1:length(behavioral_epoch_names)
    curr_epoch_name = behavioral_epoch_names{i};
    curr_epoch_start_stop = source_data.behavior.(processing_config.active_expt.name).time(i,:);
    
    % microseconds (10^6): 1000000
    % nanoseconds (10^9): 1000000000
    
    curr_epoch_duration = (curr_epoch_start_stop(2) - curr_epoch_start_stop(1)) / 1000000;
    disp(curr_epoch_duration)
    behavioral_epochs.(curr_epoch_name) = [curr_epoch_start_stop curr_epoch_duration];
end

% curr_activity_table = array2table(processing_config.active_expt.behavior_list, ...
%         'VariableNames',{'epoch_start_time', 'epoch_end_time', 'type'});
    
curr_activity_table = table(convert_date(processing_config.active_expt.behavior_list(:,1)),  ...
    convert_date(processing_config.active_expt.behavior_list(:,2)),  ...
    categorical(processing_config.active_expt.behavior_list(:,3), [1:length(behavioral_state_names)], behavioral_state_names),  ...
    'VariableNames',{'epoch_start_time', 'epoch_end_time', 'type'});
    
TimeWindows
    FileLoaderOpenEphys
function [date_out] = convert_date(aDateTimestamp)
%% Raw Timestamps Mode:
return_raw_timestamps = true;

if return_raw_timestamps
   date_out = aDateTimestamp;
else
    %     date_out = datetime(aDateTimestamp,'ConvertFrom','epochtime','TicksPerSecond',1000);
    date_out = datetime(aDateTimestamp,'ConvertFrom','epochtime','Epoch','2000-01-01','TicksPerSecond', 1000000); % Prety close
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