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

% curr_activity_table = array2table(processing_config.active_expt.behavior_list, ...
%         'VariableNames',{'epoch_start_time', 'epoch_end_time', 'type'});
    
curr_activity_table = table(datetime(processing_config.active_expt.behavior_list(:,1),'ConvertFrom','epochtime','TicksPerSecond',1000),  ...
    datetime(processing_config.active_expt.behavior_list(:,2),'ConvertFrom','epochtime','TicksPerSecond',1000),  ...
    processing_config.active_expt.behavior_list(:,3),  ...
    'VariableNames',{'epoch_start_time', 'epoch_end_time', 'type'});
    

    
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