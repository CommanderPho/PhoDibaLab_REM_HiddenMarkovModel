% PhoDiba_BayesianDecoding2021 - Uses some example code to compute the maximum likelihood of the data
% Detailed explanation goes here
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 28-Oct-2021 ; Last revision: 29-Oct-2021 

import PhoFallAnalysis2021.*

override_should_recompute = false;

% smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1_v7_3.mat');
% smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/toAddVariables.mat');
% smartload('/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/PlaceFields/biDirectional.mat', 'PF_sorted_biDir', 'spatialTunings_biDir');
% smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1_v7_3.mat', ...
%     'active_processing', 'general_results', 'num_of_electrodes', 'processing_config', 'results_array', 'source_data', 'timesteps_array');
% smartload('/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_30-Oct-2021/Roy-maze1/toAddVariables.mat', ...
%     'behavior', 'fileinfo', 'secondaryPBEs');
% smartload('/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_30-Oct-2021/Roy-maze1/PlaceFields/biDirectional.mat', ...
%     'PF_sorted_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'spatialTunings_biDir');

% experimentName = 'RoyMaze1';
% experimentName = 'Roy-maze1';
% experimentName = 'KevinMaze1';
experimentName = 'Kevin-maze1';
% '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/Kevin-maze1'
% smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1_v7_3.mat', ...
%     'active_processing', 'general_results', 'num_of_electrodes', 'processing_config', 'results_array', 'source_data', 'timesteps_array');

parentFolder = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/';
smartload([parentFolder experimentName '/toAddVariables.mat'], ...
    'behavior', 'fileinfo', '-f');
smartload([parentFolder experimentName '/PlaceFields/biDirectional.mat'], ...
    'PF_sorted_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'spatialTunings_biDir', 'positionBinningInfo_biDir');
smartload([parentFolder experimentName '/TrackLaps/trackLaps.mat']);

smartload([parentFolder experimentName '/spikesVariables.mat']);



% smartload('C:\Share\data\RoyMaze1\analysesResults_13-Oct-2021\PlaceFields\biDirectional.mat', 'PF_sorted_biDir')


addpath(genpath('helpers'));
% [maxL] = 
% 
% PF_sorted_biDir,
% 
% % PF_sorted_biDir
% spikeTimes = active_processing.processed_array{1, 1}.by_epoch.track.spike_data;
spikeTimes = active_processing.spikes.time(plot_outputs.filter_active_units);

%% Get filter info for active units
% filter_config.filter_included_cell_types = {};
% filter_config.filter_maximum_included_contamination_level = {2};
% [pf_filter.filter_active_units, pf_filter.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, true, filter_config.filter_included_cell_types, filter_config.filter_maximum_included_contamination_level);
% temp.num_active_units = sum(pf_filter.filter_active_units, 'all');
% fprintf('PF Filter: Including %d of %d total units\n', temp.num_active_units, length(pf_filter.filter_active_units));


% PhoFallAnalysis2021.showPlaceCellSpatialTuningsPlot();


[PositionBins, placeFieldTuningCurves, sortedTuningCurveIndicies] = PhoFallAnalysis2021.showPlaceCellSpatialTuningPlots(spatialTunings_biDir, fileinfo, plot_outputs.original_unit_index);



%% Previously using: 
% placeFieldTuningCurves = PF_sorted_biDir;



% Actual timesteps if we want those
% timesteps_array{1, 1}
tau = 0.1; % bin size (seconds)
% tau = 0.25; % bin size (seconds)
% tau = 1.0; % bin size (seconds)
% tau = 10.0;

%% Define the active time range:

% Whole period:
TrialStart = active_processing.behavioral_epochs.start_seconds(1);
TrialEnd = active_processing.behavioral_epochs.end_seconds(3);

% % Track only:
% TrialStart = active_processing.behavioral_epochs.start_seconds(2);
% TrialEnd = active_processing.behavioral_epochs.end_seconds(2);


%% Compute the actual likelihoods:
if ~exist('likelihood','var') || (exist('override_should_recompute','var') && override_should_recompute)
    disp('recomputing maximum likelihood...')
    [maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, placeFieldTuningCurves, PositionBins, TrialStart, TrialEnd, tau);
    disp('done.')
else
    warning('variable "likelihood" already exists, and since it is expensive to compute the existing value will be used.');
end

%% Plotting and Validation:
%%%%%%%%
%% decode the animal's position using the maximum likelihood estimate
nTimeBins = ceil((TrialEnd-TrialStart)/tau);
activeTimeBins = linspace(TrialStart, TrialEnd, (nTimeBins-1)); % 1x353520
% t_rel is aligned with the timestamps in the active_processing.position_table's timestamp column
t_rel = ((fileinfo.xyt2(:, 2)-fileinfo.tbegin) ./ 1e6); % Convert to relative timestamps since start
% TimestampPosition = t_rel;
AnimalPosition = fileinfo.xyt2(:,1);
% Main Visualization of output: compares the animal's actual recorded position to the maximum likelihood predected position at each timepoint
subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, AnimalPosition);

% t_rel, fileinfo.xyt2(:, 1)

% % determine the index of the preferred location for each place cell from its place field
% [maxFiringRate, maxIndex] = max(placeFieldTuningCurves,[],2);
% % convert the indices stored in maxIndex to position on the track
% maxPlace = PositionBins(maxIndex);
% %% Firing Rate Plot
% % compute firing rates for all place cells during Trial 1
% firingRates = computeFiringRates(spikeTimes, TrialStart, TrialEnd, tau);
% % display the firing rates of four place cells during one traversal 
% % of the animal along the linear track
% %% Winner-take-all maximum place cell detection:
% % determine which place cell has the maximum firing rate at each moment, 
% % and then determine the preferred location for this place cell
% [maxRate, maxIndex] = max(firingRates);
% posTrial = maxPlace(maxIndex);



% Vizualization of likelihoods (plots samples)
% subfn_plotSampleTrajectories(PositionBins, likelihood);


% %% Visualize Firing Rate and Winner-take-all results:
% %
% %% Testing: Plot some sample winner-take-all place fields:
% % cell number (row index) for a subset of 12 place cells
% cellNum = [1, 2, 3, 6, 8, 9, 11, 13, 15, 17, 19, 20];
% figure('Position', [100 100 1000 600])
% for i = 1:12
%     subplot(3,4,i)
%     bar(PositionBins, placeFieldTuningCurves(cellNum(i), :));
%     hold on
%     pos = maxPlace(cellNum(i));
%     plot([pos pos], [0 maxFiringRate(cellNum(i))], 'r', 'LineWidth', 2);
% %     axis([0 4 0 30])
%     xlabel('Position (m)')
%     ylabel('Firing Rate (Hz)')
%     title(['Cell #' num2str(cellNum(i))])
% end
% %% Firing Rate plot:
% figure('Position', [100 100 1000 600])
% cellNum = [48 44 11 51];
% for i = 1:4
%     subplot(4,1,i)
%     rates = firingRates(cellNum(i),:);
%     times = linspace(TrialStart, TrialEnd, length(rates));
%     bar(times, rates);
%     xlim([TrialStart TrialEnd]);
%     xlabel('Time (sec)')
%     ylabel('Firing Rate (Hz)')
%     title(['Cell #' num2str(cellNum(i))])
% end
% % display the computed trajectory of the animal for Trial 1
% figure
% % timesTrial = linspace(TrialStart, TrialEnd, 75); % 1x75 double
% times = linspace(TrialStart, TrialEnd, length(posTrial));
% scatter(times, posTrial, 'r*')
% hold on
% plot(TimestampPosition, AnimalPosition, '-b', 'LineWidth', 2);
% xlim([TrialStart, TrialEnd])
% xlabel('Time (sec)')
% ylabel('Position (m)')

%% end main function body


% ------------- END OF CODE --------------
