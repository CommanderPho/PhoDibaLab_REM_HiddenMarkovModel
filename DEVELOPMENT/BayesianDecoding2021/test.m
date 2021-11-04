%% This script tests PhoBayesianDecoder.m's ability to load data saved to file by PhoPositionAnalaysis2021.m and use it to independent compute place curves

% function [] = test()
%TEST Summary of this function goes here
%   Detailed explanation goes here

%     sessionNameHypenated = 'RoyMaze1';
%     sessionNameCamelCase = 'RoyMaze1';
%     sessionNameHypenated = 'Roy-maze1';
    
%     sessionNameCamelCase = 'TedMaze2';
%     sessionNameHypenated = 'Ted-maze2';

    % sessionNameHypenated = 'KevinMaze1';
    %sessionNameHypenated = 'Kevin-maze1';
    
    %sessionNameCamelCase = 'KevinMaze1';
%     sessionNameHypenated = sessionNameCamelCase; %override on windows to get by naming problems

  
    [activeSessionInfo] = PhoBayesianDecoder.getHiroExperimentName(1, 1);
     parentFolder = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analyses';
%     parentFolder = 'C:\Share\data\analysesResults';

    do_decoding = false;

    outputFiguresFolder = fullfile(parentFolder, activeSessionInfo.sessionNameHypenated, 'Figures');
    if ~exist(outputFiguresFolder, 'dir')
        mkdir(outputFiguresFolder)
    end
%     parentFolder = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/';
    obj = PhoBayesianDecoder();
    obj.performLoadDataHiroFormat(parentFolder, activeSessionInfo.sessionNameHypenated);
    
    % Config 1:
    sigma = [6];
    bin_size = [3]; % spatial bin size (cm)
    f_base = 2; % base firing rate (Hz)
    min_t_occ = 0.5;
    % [] = buildTuningCurves(obj, spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma, f_base, min_t_occ);
    obj.buildTuningCurves(bin_size, sigma, f_base, min_t_occ);
    %% Should update obj.TuningCurves when done!
    [outFilePath] = obj.performSaveComputedData(parentFolder, activeSessionInfo.sessionNameHypenated, activeSessionInfo.sessionNameCamelCase);

    %% Filtering Options:
    % filter_config.filter_included_cell_types = {};
    filter_config.filter_included_cell_types = {'pyramidal'};
    % filter_config.filter_included_cell_types = {'interneurons'};
    filter_config.filter_maximum_included_contamination_level = {2};
    obj.applyFilter(filter_config);
    obj.clearFilter();
    obj.plotKouroshLoadedPlaceFieldSpatialTunings(activeSessionInfo.sessionNameHypenated, outputFiguresFolder);
    obj.plotPlaceFieldSpatialTunings(outputFiguresFolder);

    %% DECODING:

    if do_decoding
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
    end
% end

