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

% active_processing.behavioral_state_names




num_of_electrodes = height(active_processing.spikes);

% Loop over electrodes
for electrode_index = 1:num_of_electrodes
    temp.curr_spikes = active_processing.spikes.time{electrode_index};
    subplot(num_of_electrodes, 1, electrode_index);
    histogram(temp.curr_spikes)
end

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