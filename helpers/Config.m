data_config.root_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject';
data_config.source_data_prefix = 'Hiro_Datasets';

data_config.output.intermediate_file_names = {'PhoIntermediate_Stage0_0.mat', 'PhoIntermediate_Stage0_1.mat'};


data_config.source_root_path = fullfile(data_config.root_parent_path, data_config.source_data_prefix);
data_config.output.intermediate_file_paths = cellfun((@(filename) fullfile(data_config.source_root_path, filename)), data_config.output.intermediate_file_names, 'UniformOutput', false);


% microseconds (10^6): 1000000
% nanoseconds (10^9): 1000000000
data_config.conversion_factor = (10^6);
data_config.behavioral_epoch_names = {'pre_sleep', 'track', 'post_sleep'};

% Process one of the experiments: 
processing_config.active_expt.name = 'RoyMaze1';
processing_config.step_sizes = {0.1, 1.0}; % Step Sizes in seconds 
processing_config.num_step_sizes = length(processing_config.step_sizes);