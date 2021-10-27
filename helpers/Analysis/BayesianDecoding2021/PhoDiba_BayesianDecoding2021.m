


addpath(genpath('helpers'));

% [maxL] = 
% 
% PF_sorted_biDir,
% 
% % PF_sorted_biDir


% spikeTimes = active_processing.processed_array{1, 1}.by_epoch.track.spike_data;
spikeTimes = active_processing.spikes.time(plot_outputs.filter_active_units);

% Actual timesteps if we want those
% timesteps_array{1, 1}
tau = 0.1; % bin size (seconds)

%% Define the active time range:
TrialStart = active_processing.behavioral_epochs.start_seconds(1);
TrialEnd = active_processing.behavioral_epochs.end_seconds(3);

% Build the needed PositionBins value from fileinfo:
[posBinEdges, linearPoscenters, linearPos] = subfn_buildMissingPositionBinInformation(fileinfo.xyt2);
% Finally set the needed PositionBins value:
PositionBins = linearPoscenters; 

[maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, PF_sorted_biDir, PositionBins, TrialStart, TrialEnd, tau);

%% end main function body


function [maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, PlaceFields, PositionBins, TrialStart, TrialEnd, tau)
    
    nTimeBins = ceil((TrialEnd-TrialStart)/tau);
    activeTimeBins = linspace(TrialStart, TrialEnd, (nTimeBins-1)); % 1x353520
    
    % % use what we already have
    % spikeCounts = computeSpikeCounts(active_processing.processed_array{1, 1}.all.spike_data, Trial1Start, Trial1End, tau);
    
    % rebuild it using their function:
    spikeCounts = computeSpikeCounts(spikeTimes, TrialStart, TrialEnd, tau); % spikeCounts: 353510x68
    
    if (size(spikeCounts, 1) ~= length(activeTimeBins))
        error('Oh no!')
    end
    
    % PlaceFields is a 68 x 108 matrix that contains the place fields of 68 place cells. 
    % Each row of the matrix stores the place field for a different cell. The place field is specified by the average firing rate of the cell as a function of the animal's position on the track.
%     PlaceFields = PF_sorted_biDir; % PF_sorted_biDir: 68x108
    
    %% Compute Maximum Likelihoods:
    likelihood = computeLikelihood(spikeCounts, PlaceFields, tau); % 108x353510 double
    [~, index] = max(likelihood, [], 1);
    maxL = PositionBins(index);


    % 
    % %% decode the animal's position using the maximum likelihood estimate
    % 
    % t_rel = ((fileinfo.xyt2(:, 2)-fileinfo.tbegin) ./ 1e6); % Convert to relative timestamps since start
    % % t_rel is aligned with the timestamps in the active_processing.position_table's timestamp column
    % 
    % % TimestampPosition = fileinfo.xyt2(:,2);
    % TimestampPosition = t_rel;
    % AnimalPosition = fileinfo.xyt2(:,1);
    % 
    % % % Spike counts during each estimation time bin
    % % nTimeBins = ceil((TrialEnd-TrialStart)/tau);
    % % windows = linspace(TrialStart, TrialEnd, nTimeBins);
    % 
    % 
    % nTimeBins = ceil((TrialEnd-TrialStart)/tau);
    % windows = linspace(TrialStart, TrialEnd, nTimeBins-1);
    % 
    % 
    % subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, fileinfo);
    % 
    % 
    % subfn_plotSampleTrajectories(linearPoscenters, likelihood);
    % 


end % end subfn_computeMaximumLikelihood




function [posBinEdges, linearPoscenters, linearPos] = subfn_buildMissingPositionBinInformation(xyt2)
    % xyt2: a cell array of the format that's present in fileinfo.xyt2
    posBinSize = 2; % in cm
    linearPos = xyt2(:, [1 3]); % the third row indicates lap indices of the positions        
    % defining the position bins
    nPosBins = floor((max(linearPos(:, 1)) - min(linearPos(:, 1)))/posBinSize);
    posBinEdges = min(linearPos(:, 1)): posBinSize: max(linearPos(:, 1)); % center of the position bins
    linearPoscenters = posBinEdges(1:end-1) + posBinSize/2; % 1 x 108
end


function subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, fileinfo)
    % subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, fileinfo);
    figure
    clf
    hold off;
    % plot the most likely trajectory
    % scatter(activeTimeBins, maxL, 'r', 'MarkerEdgeAlpha',0.4');
    plot(activeTimeBins, maxL, 'r');
    hold on
    
    % Plot the animal's actual trajectory
    plot(t_rel, fileinfo.xyt2(:, 1), 'b', 'LineWidth', 0.5); 
    hold off;
    xlim([15030 15330]);
    
    % plot(1:length(maxL), maxL, 'r');
    xlim([active_processing.behavioral_epochs.start_seconds(2) 0.1 * active_processing.behavioral_epochs.end_seconds(2)]); % [15141, 15162] are good
end

function subfn_plotSampleTrajectories(PositionBins, likelihood)
    % 
    
    % display the likelihoods for all position bins
    % temp.start_idx = 11999;
    temp.start_idx = 12999;
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

end % end subfn_plotSampleTrajectories