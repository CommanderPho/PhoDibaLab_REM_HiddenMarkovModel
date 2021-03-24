% 'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run'
data_config.output.March24OutputFilepath = '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/ResultsTemp/March 24/smoothed.mat';

temp.smoothedVariablesList = {'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run'};

if ~exist('gau_sdf','var')
    %% Load from the file:
    load(data_config.output.March24OutputFilepath, 'gau_sdf', 'tbin_edges', 'tbin_centers', 'spk_count', 'run');
end

plot(tbin_centers,gau_sdf,'g');