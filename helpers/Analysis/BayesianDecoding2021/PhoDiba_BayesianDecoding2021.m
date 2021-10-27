
% 
% load('C:\Share\data\RoyMaze1\analysesResults_13-Oct-2021\PlaceFields\biDirectional.mat', 'PF_sorted_biDir')

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
% tau = 1.0; % bin size (seconds)

%% Define the active time range:
TrialStart = active_processing.behavioral_epochs.start_seconds(1);
TrialEnd = active_processing.behavioral_epochs.end_seconds(3);

% Build the needed PositionBins value from fileinfo:
[posBinEdges, linearPoscenters, linearPos] = subfn_buildMissingPositionBinInformation(fileinfo.xyt2);
% Finally set the needed PositionBins value:
PositionBins = linearPoscenters; 

%% Compute the actual likelihoods:
if ~exist('likelihood','var') || (exist('override_should_recompute','var') && override_should_recompute)
    disp('recomputing maximum likelihood...')
    [maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, PF_sorted_biDir, PositionBins, TrialStart, TrialEnd, tau);
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

TimestampPosition = t_rel;
AnimalPosition = fileinfo.xyt2(:,1);


% determine the index of the preferred location for each place cell 
% from its place field
[maxFiringRate, maxIndex] = max(PF_sorted_biDir,[],2);
% convert the indices stored in maxIndex to position on the track
maxPlace = PositionBins(maxIndex);

%% Firing Rate Plot
% compute firing rates for all place cells during Trial 1
firingRates = computeFiringRates(spikeTimes, TrialStart, TrialEnd, tau);
% display the firing rates of four place cells during one traversal 
% of the animal along the linear track

%% Winner-take-all maximum place cell detection:
% determine which place cell has the maximum firing rate at each moment, 
% and then determine the preferred location for this place cell
[maxRate, maxIndex] = max(firingRates);
posTrial = maxPlace(maxIndex);



% Main Visualization of output: compares the animal's actual recorded position to the maximum likelihood predected position at each timepoint
subfn_plotTrajectoryComparison(activeTimeBins, maxL, t_rel, fileinfo);

% Vizualization of likelihoods (plots samples)
subfn_plotSampleTrajectories(PositionBins, likelihood);


%% Visualize Firing Rate and Winner-take-all results:
%

%% Testing: Plot some sample winner-take-all place fields:
% cell number (row index) for a subset of 12 place cells
cellNum = [1, 2, 3, 6, 8, 9, 11, 13, 15, 17, 19, 20];
figure('Position', [100 100 1000 600])
for i = 1:12
    subplot(3,4,i)
    bar(PositionBins, PF_sorted_biDir(cellNum(i), :));
    hold on
    pos = maxPlace(cellNum(i));
    plot([pos pos], [0 maxFiringRate(cellNum(i))], 'r', 'LineWidth', 2);
%     axis([0 4 0 30])
    xlabel('Position (m)')
    ylabel('Firing Rate (Hz)')
    title(['Cell #' num2str(cellNum(i))])
end

%% Firing Rate plot:
figure('Position', [100 100 1000 600])
cellNum = [48 44 11 51];
for i = 1:4
    subplot(4,1,i)
    rates = firingRates(cellNum(i),:);
    times = linspace(TrialStart, TrialEnd, length(rates));
    bar(times, rates);
    xlim([TrialStart TrialEnd]);
    xlabel('Time (sec)')
    ylabel('Firing Rate (Hz)')
    title(['Cell #' num2str(cellNum(i))])
end


% display the computed trajectory of the animal for Trial 1
figure
% timesTrial = linspace(TrialStart, TrialEnd, 75); % 1x75 double
times = linspace(TrialStart, TrialEnd, length(posTrial));
scatter(times, posTrial, 'r*')
hold on
plot(TimestampPosition, AnimalPosition, '-b', 'LineWidth', 2);
xlim([TrialStart, TrialEnd])
xlabel('Time (sec)')
ylabel('Position (m)')






%% end main function body


function [maxL, likelihood, activeTimeBins] = subfn_computeMaximumLikelihood(spikeTimes, PlaceFields, PositionBins, TrialStart, TrialEnd, tau)
    nTimeBins = ceil((TrialEnd-TrialStart)/tau);
    activeTimeBins = linspace(TrialStart, TrialEnd, (nTimeBins-1)); % 1x353520
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
%     xlim([active_processing.behavioral_epochs.start_seconds(2) 0.1 * active_processing.behavioral_epochs.end_seconds(2)]); % [15141, 15162] are good
    title('most-likely vs. observed trajectory comparsion');
end

function subfn_plotSampleTrajectories(PositionBins, likelihood)
    % subfn_plotSampleTrajectories: display the likely trajectory for
    % several time bins
  
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
    title('most-likely sample trajectories')
end % end subfn_plotSampleTrajectories