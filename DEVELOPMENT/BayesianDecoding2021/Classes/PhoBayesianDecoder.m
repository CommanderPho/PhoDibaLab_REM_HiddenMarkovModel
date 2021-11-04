classdef PhoBayesianDecoder < handle
    %PhoBayesianDecoder Summary of this class goes here
    %   Detailed explanation goes here

    properties
        TuningCurves
        InformationContentCurves
        Loaded
        Parameters
        Posteriors
    end

    methods
        function obj = PhoBayesianDecoder()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.TuningCurves = struct;
            obj.InformationContentCurves = struct;
            obj.Loaded = struct;
            obj.Parameters = struct;
            obj.Posteriors = struct;
            
        end


        function [] = build(obj)
            %%% New mode:
            %% Pho's comments:
            % Looks like I need some position data in X
            obj.Parameters = struct;

            % spikes = {source_data.spikes.RoyMaze1.time};
            [~, unitSpikeCells, ~] = fnFlatSpikesToUnitCells(obj.Loaded.spikeStruct.t, obj.Loaded.spikeStruct.unit, false);
%             spikes = unitSpikeCells; % the timestamps we need
             obj.Parameters.spikes = unitSpikeCells'; % for rebuilt flat spike arrays, we need just the transpose of the whole cell array, not the elements
%              obj.Parameters.spikes = fnCellContentsTranpose(unitSpikeCells)'; % for active_processing.spikes.time type objects, we need the tranpose of the contents and the whole cell array
             % obj.Parameters.spikes should be a 1x76 cell with contents:
                % obj.Parameters.spikes{1}: 80828x1 double
                % obj.Parameters.spikes{2}: 17874x1 double
                ...
            % Absolute timestamp version are available here: source_data.position.RoyMaze1.t            
            
            % t_rel is aligned with the timestamps in the active_processing.position_table's timestamp column
            % t = ((fileinfo.xyt2(:, 2)-fileinfo.tbegin) ./ 1e6)'; % Convert to relative timestamps since start
            obj.Loaded.positionTable.t_rel = ((obj.Loaded.positionTable.t - obj.Loaded.fileinfo.tbegin) ./ 1e6); % Convert to relative timestamps since start
            
            obj.Parameters.X = obj.Loaded.positionTable.linearPos';
%             obj.Parameters.t = obj.Loaded.positionTable.t_rel;
            obj.Parameters.t = obj.Loaded.positionTable.t';
            
            obj.Parameters.sample_rate = obj.Loaded.fileinfo.Fs;

            obj.Parameters.t_start = obj.Loaded.behavior.time(2,1);
            obj.Parameters.t_end = obj.Loaded.behavior.time(2,2);

%             obj.buildTuningCurves(obj.Parameters.spikes, obj.Parameters.X, obj.Parameters.t, obj.Parameters.sample_rate, ...
%                 obj.Parameters.t_start, obj.Parameters.t_end, 

        end

        function [] = buildTuningCurves(obj, bin_size, sigma, f_base, min_t_occ, t_start, t_end)
%             t_start

%             sigma = [5];
%             bin_size = [2]; % spatial bin size (cm)
%             f_base = 2; % base firing rate (Hz)
%             min_t_occ = 0.5;

            if exist('t_start', 'var')
                obj.Parameters.t_start = t_start; % use the user's value
            end
            if exist('t_end', 'var')
                obj.Parameters.t_end = t_end; % use the user's value
            end
            obj.Parameters.bin_size = bin_size;
            obj.Parameters.sigma = sigma;
            obj.Parameters.f_base = f_base;
            obj.Parameters.min_t_occ = min_t_occ;

            obj.performBuildTuningCurves(obj.Parameters.spikes, obj.Parameters.X, obj.Parameters.t, obj.Parameters.sample_rate, ...
                obj.Parameters.t_start, obj.Parameters.t_end, bin_size, sigma, f_base, min_t_occ);
    
            % Compute the obj.TuningCurves.sortedOriginalUnitIndicies that would be required to order the tuning curves by their location on the track:
            [obj.TuningCurves.sortedPeakPlaces, obj.TuningCurves.sortIndicies, obj.TuningCurves.sortedOriginalUnitIDs, obj.TuningCurves.inverseSortingIndicies] = fnSortPlaceCellSpatialTuningCurves(obj.TuningCurves.lambda, ...
                obj.TuningCurves.coords{1}, obj.Loaded.okUnits);
            obj.TuningCurves.originalUnitIDs = obj.Loaded.okUnits;
            % Unchanged, the units will be colored by ascending ID:
            obj.TuningCurves.sortDynamicColors = colormap(jet(size(obj.TuningCurves.sortedOriginalUnitIDs, 1))); % should be the RGB triplets
            % This will cause the units to be colored by their tuning position 
%             obj.TuningCurves.sortDynamicColors = obj.TuningCurves.sortDynamicColors(obj.TuningCurves.sortedOriginalUnitIDs, :);

            fprintf('Successfully built tuning curves.\n');
        end

        function [] = performBuildTuningCurves(obj, spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma, f_base, min_t_occ)
            % performBuildTuningCurves(obj, ...):
            [~, obj.TuningCurves.lambda] = build_tuning_curves(spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma);
            %{
            lambda: (N+1)-dimensional array of mean firing rate for each bin,
                    where the size of the first dimension is the number of units.
            %}
            [obj.TuningCurves.coords, obj.TuningCurves.alpha, obj.TuningCurves.beta] = build_tuning_curves(spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma, f_base, min_t_occ);
            %{
            coords: N x 1 cell array of coordinate locations, where N is the number of stimulus dimensions:
                        coords{1}: N-dimensional array of the coordinate position for
                        each bin in stimulus dimension 1.
                        coords{2}: N-dimensional array of the coordinate position for
                        each bin in stimulus dimension 2.
            
            alpha: (N+1)-dimensional array of alpha parameter for each bin,
                    where the size of the first dimension is the number of units.
                    
            beta: N-dimensional array of beta parameter for each bin.
                    (if f_base and min_t_occ are supplied)
            %}
            %% Get information content curves.
            obj.InformationContentCurves = get_IC_curves(obj.TuningCurves.alpha, obj.TuningCurves.beta, f_base, min_t_occ);
        end


        function performLoadDataHiroFormat(obj, parentFolder, experimentName)
            %performLoadDataHiroFormat Summary of this method goes here
            %   Detailed explanation goes here
            obj.Loaded = struct; % clear the loaded if it hasn't already been done.
            obj.Loaded.parentFolder = parentFolder;
            obj.Loaded.experimentName = experimentName;
            
            temp.L1 = smartload(fullfile(parentFolder, experimentName, 'toAddVariables.mat'), ...
                'behavior', 'fileinfo', '-f');
            temp.L2 = smartload(fullfile(parentFolder, experimentName, '/PlaceFields/biDirectional.mat'), ...
                'PF_sorted_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'spatialTunings_biDir', 'positionBinningInfo_biDir');
            temp.L3 = smartload(fullfile(parentFolder, experimentName, '/TrackLaps/trackLaps.mat'), ...
                'lapsStruct', 'turningPeriods', 'laps', 'totNumLaps', 'lapsTable', 'positionTable', 'currMazeShape', 'occupancyInfo', 'trackInfo');
            temp.L4 = smartload(fullfile(parentFolder, experimentName, '/spikesVariables.mat'), ...
                'spikeStruct', 'okUnits', 'shanks');
            [obj.Loaded] = fnMergeStructs(obj.Loaded, temp.L1, temp.L2, temp.L3, temp.L4);
            
            fprintf('Successfully loaded data.\n');
        end


        function [outFilePath] = performSaveComputedData(obj, parentFolder, experimentName, experimentIdentifier)
            %performSaveComputedData Summary of this method goes here
            %   [outFilePath] = obj.performSaveComputedData(parentFolder, experimentName, experimentIdentifier);
               
            if ~exist('experimentName', 'var')
                experimentName = obj.Loaded.experimentName;
            end
            if ~exist('experimentIdentifier', 'var')
                experimentIdentifier = obj.Loaded.experimentIdentifier;
            end
           if ~exist('parentFolder', 'var')
                parentFolder = obj.Loaded.parentFolder;
           end

%             input_variables_names = {'spikes', 'X', 't', 'sample_rate', 't_start', 't_end', 'bin_size', 'sigma', 'f_base', 'min_t_occ'};
            output_variable_names = {'lambda', 'coords', 'alpha', 'beta'};
%             for i = 1:length(input_variables_names)
%                 BayesocampusResults.(experimentIdentifier).Inputs.(input_variables_names{i}) = eval(input_variables_names{i}); 
%             end
            
            for i = 1:length(output_variable_names)
                BayesocampusResults.(experimentIdentifier).Outputs.(output_variable_names{i}) = obj.TuningCurves.(output_variable_names{i}); 
            end
            outFilePath = fullfile(parentFolder, experimentName, sprintf('Bayesiocampus_Results_%s_2021_11_01-6.mat', experimentName));
            fprintf('Saving data data to %s...\n', outFilePath);
            save(outFilePath, "BayesocampusResults");
            fprintf('\t done.\n');
        end


        function [] = performNeuralDecode(obj, t_0, t_f, t_step_seconds)
            %performNeuralDecode Decodes the positions at the given time steps using the computed tuning curves and the spikes matrix.
            %   Detailed explanation goes here
           obj.Posteriors = struct;

%             t_0 = 1000;
%             t_f = 1030;
%             t_step = 0.25;

            t_step = t_step_seconds * 1e6; % must convert to microseconds, since that's what t_0 and t_f are assumed to be.

            obj.Posteriors.t_start = t_0:t_step:t_f; % if t_start and t_end are 1x27 double, then posteriors will be 27x73 double
            obj.Posteriors.t_end = (t_0:t_step:t_f) + t_step;
            
            obj.Posteriors.poiss_posterior = bayesian_decode(obj.Parameters.spikes, obj.Posteriors.t_start, obj.Posteriors.t_end, obj.TuningCurves.lambda);
            obj.Posteriors.nb_posterior = bayesian_decode(obj.Parameters.spikes, obj.Posteriors.t_start, obj.Posteriors.t_end, obj.TuningCurves.alpha, obj.TuningCurves.beta);
        end

        function [] = plotPosteriors(obj)
            num_sample_timesteps = length(obj.Posteriors.t_start);
            temporal_midpoints = obj.Posteriors.t_start + ((obj.Posteriors.t_end - obj.Posteriors.t_start) ./ 2.0); % Get the centers of the temporal bins the posteriors were calculated for
            
            num_position_bins = size(obj.Posteriors.poiss_posterior, 2); 

            figure(1);
%             imagesc(obj.Posteriors.poiss_posterior);
%             title('Posterior')
%             xlabel('Position Bin Likelihood')
%             ylabel('Timestep')

            [~, maxLikelyPositionBins] = max(obj.Posteriors.poiss_posterior, [], 2); % should get a 27x1 vector of the most likely positions for each timestamp
            activeSpatialLinearPositions = obj.TuningCurves.coords{1};
            maxLikelyPositions = activeSpatialLinearPositions(maxLikelyPositionBins);

%             % Main Visualization of output: compares the animal's actual recorded position to the maximum likelihood predected position at each timepoint
%             hold off;
%             % plot the most likely trajectory
%             scatter(temporal_midpoints, maxLikelyPositions, 'r', 'MarkerEdgeAlpha', 0.4);
%         %     plot(activeTimeBins, maxL, 'r');
%             hold on
%             % Plot the animal's actual trajectory
%             plot(obj.Loaded.positionTable.t, obj.Loaded.positionTable.linearPos, 'b', 'LineWidth', 0.5); 
%             hold off;
%             xlim([temporal_midpoints(1) temporal_midpoints(end)]);

            PhoFallAnalysis2021.subfn_plotTrajectoryComparison(temporal_midpoints, maxLikelyPositions, obj.Loaded.positionTable.t, obj.Loaded.positionTable.linearPos);
            xlim([temporal_midpoints(1) temporal_midpoints(end)]);

        end



        function [figH, h] = plotKouroshLoadedPlaceFieldSpatialTunings(obj, experimentName, figureSaveParentPath)
            % Plots the tuning curves loaded from the biDirectional.mat
            % file's PF_sorted_biDir variable, which were computed using
            % Kourosh's old code.
            if ~exist('experimentName', 'var')
                experimentName = obj.Loaded.experimentName;
            end

            %% Plot tuning curves:
            [figH, h] = fnPlotPlaceCellSpatialTunings(obj.Loaded.PF_sorted_biDir, 'linearPoscenters', obj.Loaded.positionBinningInfo_biDir.linearPoscenters, 'unitLabels', num2cellstr(obj.Loaded.runTemplate_biDir));
            curr_fig_name = sprintf('Kourosh Style - %s - Sorted Position Tuning Curves', experimentName);
            title(curr_fig_name)
            if exist('figureSaveParentPath', 'var')
               curr_fig_output_basename = strrep(curr_fig_name, ' ', ''); % Remove Spaces from filename
               curr_fig_output_basename = strrep(curr_fig_output_basename, '/', '|'); % Remove Slashes from filename
               figureSavePath = fullfile(figureSaveParentPath, curr_fig_output_basename);
               savefig(figH, figureSavePath,'compact');
               fprintf('Figure saved to %s\n', figureSavePath);
            end
        end

        function [figH, h] = plotPlaceFieldSpatialTunings(obj, figureSaveParentPath)
            % plotPlaceFieldSpatialTunings: Plots the tuning curves internally computed from within this class
            %% Plot tuning curves:
%             activeSpatialTunings = obj.TuningCurves.lambda(obj.TuningCurves.sortedOriginalUnitIndicies, :);
%             activeSpatialLinearPositions = obj.TuningCurves.alpha;

            activePeakLocations = obj.TuningCurves.sortedPeakPlaces(obj.TuningCurves.inverseSortingIndicies);
            
            %% Sorted by unitID:
            activeOriginalUnitIDs = obj.TuningCurves.originalUnitIDs;
            activeSpatialTunings = obj.TuningCurves.lambda;
            activeUnitLabels = num2cellstr(activeOriginalUnitIDs);

            % Sorted by Unit ID:
            activeSortOrder = [];
            activeColorSortOrder = [];

%             % Sorted by tuned position, colored by original index:
            activeSortOrder = obj.TuningCurves.sortIndicies;
%             activeColorSortOrder = obj.TuningCurves.sortIndicies;

            % Colored by unitID:
%             activeColorSortOrder = 1:length(obj.TuningCurves.originalUnitIDs);


            %% Sorted by tuning place:
%             activeColorSortOrder = obj.TuningCurves.sortIndicies;
            %activeOriginalUnitIDs = obj.TuningCurves.sortedOriginalUnitIDs;
            %activeSpatialTunings = obj.TuningCurves.lambda(obj.TuningCurves.sortIndicies, :);

            activeSpatialLinearPositions = obj.TuningCurves.coords{1}';
%             activeSpatialLinearPositions = 1:size(activeSpatialTunings, 2);



            [figH, h] = fnPlotPlaceCellSpatialTunings(activeSpatialTunings, 'linearPoscenters', activeSpatialLinearPositions, ...
                'unitLabels', activeUnitLabels, 'unitColors', obj.TuningCurves.sortDynamicColors, 'colorSortOrder', activeColorSortOrder, 'sortOrder', activeSortOrder, 'peaks', activePeakLocations);
        
            curr_fig_name = sprintf('PhoBayesianDecoder Style - %s - Sorted Position Tuning Curves', obj.Loaded.experimentName);
            title(curr_fig_name)
            if exist('figureSaveParentPath', 'var')
               curr_fig_output_basename = strrep(curr_fig_name, ' ', ''); % Remove Spaces from filename
               curr_fig_output_basename = strrep(curr_fig_output_basename, '/', '|'); % Remove Slashes from filename
               figureSavePath = fullfile(figureSaveParentPath, curr_fig_output_basename);
               savefig(figH, figureSavePath,'compact'); 
               fprintf('Figure saved to %s\n', figureSavePath);
            end
        end


    end % end methods

    %% StaticFunction Prototypes BLock
    % function signature required to match the one in the function's .m file
    methods (Static)
        function [] = test()
        %TEST Summary of this function goes here
        %   Detailed explanation goes here
        
            % experimentName = 'RoyMaze1';
            % experimentName = 'Roy-maze1';
            % experimentName = 'KevinMaze1';
            experimentName = 'Kevin-maze1';
%             experimentIdentifier = 'KevinMaze1';
            experimentName = experimentIdentifier; %override on windows to get by naming problems
            % '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/Kevin-maze1'
            % smartload('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1_v7_3.mat', ...
            %     'active_processing', 'general_results', 'num_of_electrodes', 'processing_config', 'results_array', 'source_data', 'timesteps_array');
            
%             parentFolder = '/Volumes/iNeo/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_02-Nov-2021/';
            parentFolder = 'C:\Share\data\analysesResults';
            obj = PhoBayesianDecoder();
            obj.performLoadDataHiroFormat(parentFolder, experimentName);
            obj.build(); % build the parameters from the loaded data
            % Config 1:
            sigma = [5];
            bin_size = [3]; % spatial bin size (cm)
            f_base = 2; % base firing rate (Hz)
            min_t_occ = 0.5;

            obj.buildTuningCurves(bin_size, sigma, f_base, min_t_occ);
            %% Should update obj.TuningCurves when done!
            [outFilePath] = obj.performSaveComputedData(parentFolder, experimentName, experimentIdentifier);
            obj.plotKouroshLoadedPlaceFieldSpatialTunings(experimentName)
            obj.plotPlaceFieldSpatialTunings();

        end % end function test
        
    end % end static method block

    
end