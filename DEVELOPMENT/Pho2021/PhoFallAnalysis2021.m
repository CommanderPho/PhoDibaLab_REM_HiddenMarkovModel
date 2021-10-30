classdef PhoFallAnalysis2021

 %________________________________ 
% PhoFallAnalysis2021 - One line description of what the function performs
% Detailed explanation goes here
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 30-Oct-2021 ; Last revision: 30-Oct-2021 

% ------------- BEGIN CODE --------------


    %PhoFallAnalysis2021 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Property1
    end
    methods
        function obj = PhoFallAnalysis2021(inputArg1,inputArg2)
            %UNTITLED8 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    methods (Static)
        
       function [filter_active_units, original_unit_index] = showPlaceCellSpatialTuningsPlot(PlaceFields, PositionBins)
           % Like spatialTunings_biDir
          
            [~, ~] = fnPlotPlaceCellSpatialTunings(PlaceFields, 'linearPoscenters', linearPoscenters, 'unitLabels', num2cellstr(plot_outputs.original_unit_index));
       end
    end % end static method block
    %% Refactored from PhoDiba_BayesianDecoding2021.m
    methods (Static)
        function [PositionBins, placeFieldTuningCurves, sortedTuningCurveIndicies] = showPlaceCellSpatialTuningPlots(spatialTunings_biDir, fileinfo, unitIndicies)
            import PhoFallAnalysis2021.*
            %% Needs fileinfo.xyz2 and spatialTunings_biDir to work:
            %% Position Bins:
            % Build the needed PositionBins value from fileinfo:
            [posBinEdges, linearPoscenters, linearPos] = subfn_buildMissingPositionBinInformation(fileinfo.xyt2);
            % Finally set the needed PositionBins value:
            PositionBins = linearPoscenters; 
            
            %% Place Fields
            valid_place_fields = (sum(spatialTunings_biDir, 2) > 0);
            placeFieldTuningCurves = spatialTunings_biDir(valid_place_fields, :); 
            
            [peakVal, peakIndex] = max(placeFieldTuningCurves, [], 2);
            % convert the indices stored in maxIndex to position on the track
            peakPlace = PositionBins(peakIndex);
            [sortedPeakPlaces, sortedTuningCurveIndicies] = sort(peakPlace, 2, "ascend");
            
            [~, ~] = fnPlotPlaceCellSpatialTunings(placeFieldTuningCurves,'linearPoscenters', PositionBins, 'unitLabels', num2cellstr(unitIndicies));
            title('Unsorted Spatial Tunings')
            
            [~, ~] = fnPlotPlaceCellSpatialTunings(placeFieldTuningCurves(sortedTuningCurveIndicies, :),'linearPoscenters', PositionBins, 'unitLabels', num2cellstr(sortedTuningCurveIndicies));
            title('Sorted Spatial Tunings')
        end
        function [maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, PlaceFields, PositionBins, TrialStart, TrialEnd, tau)
            import PhoFallAnalysis2021.*
            nTimeBins = ceil((TrialEnd-TrialStart)/tau);
            activeTimeBins = linspace(TrialStart, TrialEnd, (nTimeBins-1)); % 1x353520
            % rebuild it using their function:
            spikeCounts = computeSpikeCounts(spikeTimes, TrialStart, TrialEnd, tau); % spikeCounts: 353510x68
            if (size(spikeCounts, 1) ~= length(activeTimeBins))
                error('Oh no!')
            end
            % PlaceFields is a 68 x 108 matrix that contains the place fields of 68 place cells. 
            % Each row of the matrix stores the place field for a different cell. The place field is specified by the average firing rate of the cell as a function of the animal's position on the track.
        %     PlaceFields = placeFieldTuningCurves; % placeFieldTuningCurves: 68x108
            %% Compute Maximum Likelihoods:
            likelihood = computeLikelihood(spikeCounts, PlaceFields, tau); % 108x353510 double
            [~, index] = max(likelihood, [], 1);
            maxL = PositionBins(index);
        end % end subfn_computeMaximumLikelihood
        function [posBinEdges, linearPoscenters, linearPos] = subfn_buildMissingPositionBinInformation(xyt2)
            % xyt2: a cell array of the format that's present in fileinfo.xyt2
            posBinSize = 2; % in cm
            linearPos = xyt2(:, [1 3]); % the third row indicates lap indices of the positions        
            % defining the position bins
            nPosBins = floor((max(linearPos(:, 1)) - min(linearPos(:, 1)))/posBinSize); % x ranges from 5.3492 to 222.22
            posBinEdges = min(linearPos(:, 1)): posBinSize: max(linearPos(:, 1)); % edges of the position bins - 1x109 double
            linearPoscenters = posBinEdges(1:end-1) + posBinSize/2; % 1 x 108 double - center of the position bins
        end
        function subfn_plotTrajectoryComparison(activeTimeBins, maxL, animalPositionTimestamps, animalPosition)
            % subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, fileinfo);
            figure(1)
            clf
            hold off;
            % plot the most likely trajectory
            scatter(activeTimeBins, maxL, 'r', 'MarkerEdgeAlpha', 0.4);
        %     plot(activeTimeBins, maxL, 'r');
            hold on
            % Plot the animal's actual trajectory
            plot(animalPositionTimestamps, animalPosition, 'b', 'LineWidth', 0.5); 
            hold off;
            xlim([15030 15330]);
            
            % plot(1:length(maxL), maxL, 'r');
        %     xlim([active_processing.behavioral_epochs.start_seconds(2) 0.1 * active_processing.behavioral_epochs.end_seconds(2)]); % [15141, 15162] are good
            title('most-likely vs. observed trajectory comparsion');
        end
        function subfn_plotSampleTrajectories(PositionBins, likelihood)
            % subfn_plotSampleTrajectories: display the likely trajectory for
            % several time bins
          
            % display the likelihoods for all position bins
            % temp.start_idx = 11999;
            temp.start_idx = 1;
            timeBins = temp.start_idx:temp.start_idx+50;
            figure('Position', [100 100 800 600])
            for i = 1:5
                subplot(5,1,i)
                plot(PositionBins, likelihood(:, timeBins(i)))
                xlabel('Position (m)');
                ylabel('Likelihood');
                %     xlim([0 3.6])
                xlim([0 PositionBins(end)])
            end
            title('most-likely sample trajectories')
        end % end subfn_plotSampleTrajectories
    end % end 
end %% end static methods block


% ------------- END OF CODE --------------
