classdef PhoBayesianDecoder < handle
    %PhoBayesianDecoder Summary of this class goes here
    %   Detailed explanation goes here

    properties
        TuningCurves
        InformationContentCurves
        Loaded
        Parameters
    end

    methods
        function obj = PhoBayesianDecoder()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.TuningCurves = struct;
            obj.InformationContentCurves = struct;
            obj.Loaded = struct;
            obj.Parameters = struct;
            
        end

        function outputArg = performNeuralDecode(obj, t_0, t_f, t_step)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
%             t_0 = 1000;
%             t_f = 1030;
%             t_step = 0.25;

            t_start = t_0:t_step:t_f;
            t_end = (t_0:t_step:t_f) + t_step;
             
            poiss_posterior = bayesian_decode(spikes, t_start, t_end, lambda);
            nb_posterior = bayesian_decode(spikes, t_start, t_end, alpha, beta);

        end


        function [] = build(obj)
            %%% New mode:
            %% Pho's comments:
            % Looks like I need some position data in X
            obj.Parameters = struct;

            % spikes = {source_data.spikes.RoyMaze1.time};
            [~, unitSpikeCells, ~] = fnFlatSpikesToUnitCells(obj.Loaded.spikeStruct.t, obj.Loaded.spikeStruct.unit, false);
%             spikes = unitSpikeCells; % the timestamps we need
             obj.Parameters.spikes = fnCellContentsTranpose(unitSpikeCells)'; % the timestamps we need
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

%             sigma = [5];
%             bin_size = [2]; % spatial bin size (cm)
%             f_base = 2; % base firing rate (Hz)
%             min_t_occ = 0.5;

            obj.performBuildTuningCurves(obj.Parameters.spikes, obj.Parameters.X, obj.Parameters.t, obj.Parameters.sample_rate, ...
                obj.Parameters.t_start, obj.Parameters.t_end, bin_size, sigma, f_base, min_t_occ);
    
        end

        function [] = performBuildTuningCurves(obj, spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma, f_base, min_t_occ)
%             [~, obj.TuningCurves.lambda] = build_tuning_curves(spikes, X, t, sample_rate, t_start, t_end, bin_size, sigma);

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
            temp.L1 = smartload([parentFolder experimentName '/toAddVariables.mat'], ...
                'behavior', 'fileinfo', '-f');
            temp.L2 = smartload([parentFolder experimentName '/PlaceFields/biDirectional.mat'], ...
                'PF_sorted_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'spatialTunings_biDir', 'positionBinningInfo_biDir');
            temp.L3 = smartload([parentFolder experimentName '/TrackLaps/trackLaps.mat']);
            temp.L4 = smartload([parentFolder experimentName '/spikesVariables.mat']);
            [obj.Loaded] = fnMergeStructs(obj.Loaded, temp.L1, temp.L2, temp.L3, temp.L4);
        end


        function performSaveComputedData(obj, parentFolder, experimentName, experimentIdentifier)
            %performSaveComputedData Summary of this method goes here
            %   Detailed explanation goes here
               
            input_variables_names = {'spikes', 'X', 't', 'sample_rate', 't_start', 't_end', 'bin_size', 'sigma', 'f_base', 'min_t_occ'};
            output_variable_names = {'lambda', 'coords', 'alpha', 'beta', 'IC_curves'};
            
%             experimentIdentifier = 'KevinMaze1';
            for i = 1:length(input_variables_names)
                BayesocampusResults.(experimentIdentifier).Inputs.(input_variables_names{i}) = eval(input_variables_names{i}); 
            end
            
            for i = 1:length(output_variable_names)
                BayesocampusResults.(experimentIdentifier).Outputs.(output_variable_names{i}) = eval(output_variable_names{i}); 
            end
            
            save('Bayesiocampus_Results_2021_11_01-6.mat', "BayesocampusResults")            
        end




    end % end methods

    
end