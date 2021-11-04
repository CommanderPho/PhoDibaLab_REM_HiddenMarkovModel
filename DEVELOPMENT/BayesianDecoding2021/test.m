%% This script tests PhoBayesianDecoder.m's ability to load data saved to file by PhoPositionAnalaysis2021.m and use it to independent compute place curves

% function [] = test()
%TEST Summary of this function goes here
%   Detailed explanation goes here

%     experimentName = 'RoyMaze1';
%     experimentIdentifier = 'RoyMaze1';
%     experimentName = 'Roy-maze1';
    
    experimentIdentifier = 'TedMaze2';
    experimentName = 'Ted-maze2';

    % experimentName = 'KevinMaze1';
    %experimentName = 'Kevin-maze1';
    
    %experimentIdentifier = 'KevinMaze1';
    experimentName = experimentIdentifier; %override on windows to get by naming problems

    % '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/Kevin-maze1'
    % smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1_v7_3.mat', ...
    %     'active_processing', 'general_results', 'num_of_electrodes', 'processing_config', 'results_array', 'source_data', 'timesteps_array');
    
    %parentFolder = 'C:\Share\data\analysesResults\KevinMaze1';
    parentFolder = 'C:\Share\data\analysesResults';

    outputFiguresFolder = fullfile(parentFolder, experimentName, 'Figures');
    if ~exist(outputFiguresFolder, 'dir')
        mkdir(outputFiguresFolder)
    end
%     parentFolder = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/';
    obj = PhoBayesianDecoder();
    obj.performLoadDataHiroFormat(parentFolder, experimentName);
    obj.build(); % build the parameters from the loaded data
    % Config 1:
    sigma = [6];
    bin_size = [3]; % spatial bin size (cm)
    f_base = 2; % base firing rate (Hz)
    min_t_occ = 0.5;
    % [] = buildTuningCurves(obj, spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma, f_base, min_t_occ);
    obj.buildTuningCurves(bin_size, sigma, f_base, min_t_occ);
    %% Should update obj.TuningCurves when done!
    [outFilePath] = obj.performSaveComputedData(parentFolder, experimentName, experimentIdentifier);

    obj.plotKouroshLoadedPlaceFieldSpatialTunings(experimentName, outputFiguresFolder)
    obj.plotPlaceFieldSpatialTunings(outputFiguresFolder);

    %% Just decode over the track period for testing:
    t_0 = obj.Loaded.behavior.time(2,1);
    t_f = obj.Loaded.behavior.time(2,2);

    %% Track-specific decoding:
%     active_lap_index = 9;    
%     t_0 = obj.Loaded.lapsTable.lapStartTime(active_lap_index);
% %     t_f = obj.Loaded.lapsTable.lapEndTime(active_lap_index);
%     t_f = obj.Loaded.lapsTable.lapEndTime(active_lap_index+2); % go through the next 2 laps

    t_step_seconds = 0.25; % 250ms steps
    % Perform neural decoding for the specified time ranges:
    obj.performNeuralDecode(t_0, t_f, t_step_seconds);

    obj.plotPosteriors();
% end

