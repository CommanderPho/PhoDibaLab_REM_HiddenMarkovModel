% PhoDibaResults_AllExpts.m
% 02/18/2021 by Pho Hale
% Displays results across all experiments (animals x days)
% Requires: across_experiment_results, active_experiment_names, active_step_sizes

%%% Note: This is a very inefficient implementation that is to serve prototyping purposes, but should be refactored for deployment scenarios.

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

if ~exist('across_experiment_results','var')
    %% Set the path to the combined across_experiment file:
    % data_config.output.all_expts_combined_parent_path = '/Users/pho/Dropbox/Classes/Spring 2020/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets';
    data_config.output.all_expts_combined_parent_path = '/Volumes/iNeo/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results';
    data_config.output.all_expts_combined_results_file_name = sprintf('PhoResults_AllExperiments.mat');

    data_config.output.all_expts_combined_results_file_path = fullfile(data_config.output.all_expts_combined_parent_path, data_config.output.all_expts_combined_results_file_name);

    %% Load Combined experiment results:
    fprintf('\t loading combined experiment results from %s... (this will take quite a while, up to 10 minutes)...\n', data_config.output.all_expts_combined_results_file_path);
    temp.variablesList = {'across_experiment_results', 'active_experiment_names', 'active_step_sizes'};
    load(data_config.output.all_expts_combined_results_file_path, temp.variablesList{:});

else
    warning('across_experiment_results already exists, using extant value.');
end

    
%% Binning Options:
current_binning_index = 1;

%% Filtering Options:
filter_config.filter_included_cell_types = {}; % All
% filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

% temp.filter_states = {'rem'};
% temp.filter_epochs = {'pre_sleep', 'post_sleep'};

%% Plotting Options:
% plottingOptions.plottingXAxis = 'index';
plottingOptions.plottingXAxis = 'timestamp';
plottingOptions.plottingYlim = [];
% plottingOptions.plottingYlim = [2 4.25];
% plottingOptions.plottingYlim = [0.2 1.4];

plottingOptions.plottingMode = 'scatter';
% plottingOptions.plottingMode = 'errorbar';
% plottingOptions.plottingMode = 'distributionPlot'; % distributionPlot should display the variance across neurons
% plottingOptions.plottingMode = 'stem';

plottingOptions.outputs.rootPath = '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/Figures';


active_binning_resolution = active_step_sizes{current_binning_index};
temp.curr_active_experiment_names = active_experiment_names(1); % Only get the first one for testing.
active_num_experiments = length(temp.curr_active_experiment_names);

temp.outputPlottingResults = cell([active_num_experiments 1]);

% Before/After Range:
transitionAnalysis.config.beforeAfterDurations = [seconds(9), seconds(9)];

% Build the range used for plotting the offset from transition points:
transitionAnalysis.results.beforeAfterDurationRange = (-transitionAnalysis.config.beforeAfterDurations(1)):seconds(active_binning_resolution):transitionAnalysis.config.beforeAfterDurations(2);

% Loop through each experiment:
for expt_index = 1:active_num_experiments
    expt_info.index = expt_index;
    expt_info.name = temp.curr_active_experiment_names{expt_index};
    
    temp.curr_timestamps = across_experiment_results{expt_index}.timesteps_array{current_binning_index};
    temp.curr_processed = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index};
    active_results = across_experiment_results{expt_index}.results_array{current_binning_index};
    
    [temp.active_binned_spike_data_matrix] = fnUnitDataCells2mat(temp.curr_processed.all.binned_spike_counts);  % 35351x126 double;
    

    
    %% Look at Transitions from/to REM:
    % Concept: Look at the periods immediately surrounding the transitions to characterize differences between states.
    %   Comparing differences at transitions is more representitive 
    temp.currentEpochs = across_experiment_results{expt_index}.active_processing.behavioral_periods_table.type;
    temp.nextEpoch = ['<undefined>'; across_experiment_results{expt_index}.active_processing.behavioral_periods_table.type];
    % Drop the last item, which has no transition
    temp.nextEpoch(end) = [];
    
    % first entry is invalid!
    
    % Build the all-valids transition matrix (for the first N-1 times)
%     temp.currentTransitions = [temp.currentEpochs(2:end), temp.nextEpoch(2:end)];
    temp.currentTransitions = [temp.currentEpochs, temp.nextEpoch];
    % Drop first row:
    temp.currentTransitions(1,:) = [];
    
%     temp.currentTransitions(2:end,:) = 
    
    temp.indicies.fromRemIndicies = (temp.currentTransitions(:,1) == 'rem');
    
    % Drop rows that don't end up in REM:
    temp.indicies.toRemIndicies = (temp.currentTransitions(:,2) == 'rem');
    
    
    state_names = {'rem', 'nrem','quiet','active'};
    for i = 1:length(state_names)
        % From X to REM:
        temp.indicies.to.(state_names{i}) = (temp.currentTransitions(:,2) == state_names{i});
        % From REM to X
        temp.indicies.from.(state_names{i}) = (temp.currentTransitions(:,1) == state_names{i});
    end
    
    % Get All transitions from REM to anything else:
    temp.indicies.active.fromRemToAny = temp.indicies.fromRemIndicies & (temp.indicies.to.nrem | temp.indicies.to.quiet | temp.indicies.to.active);
    temp.indicies.active.toRemFromAny = temp.indicies.toRemIndicies & (temp.indicies.from.nrem | temp.indicies.from.quiet | temp.indicies.from.active);
    
    sum(temp.indicies.active.fromRemToAny, 'all')
    
    sum(temp.indicies.active.toRemFromAny, 'all')
    
    transitionAnalysis.config.transitionResultNames = {'toRem','fromRem'};
    transitionAnalysis.config.transitionIndiciesNames = {'toRemFromAny','fromRemToAny'};
    
    for transitionResultIndex = 1:length(transitionAnalysis.config.transitionResultNames)
        temp.currTransitionResultName = transitionAnalysis.config.transitionResultNames{transitionResultIndex};
		temp.currTransitionIndiciesName = transitionAnalysis.config.transitionIndiciesNames{transitionResultIndex};
	
    
        % From Period Indicies where transitions occured, get the timestamp:
        transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestamps = seconds(across_experiment_results{expt_index}.active_processing.behavioral_periods_table.epoch_start_seconds(temp.indicies.active.(temp.currTransitionIndiciesName)));

        % Get the range of start/stop timestamps surrounding the transition period.
        transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestampRange = [(transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestamps - transitionAnalysis.config.beforeAfterDurations(1)), (transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestamps + transitionAnalysis.config.beforeAfterDurations(2))];
        % [temp.filtered.any_REM_indicies] = fnFilterPeriodsWithCriteria(active_processing, [], {'rem'}); % 668x1

        %% Now I can use the above start/stop ranges to average activity around the transitions.
        % Could also get nearest bins and operate on those binned results.

        % Translate each timestamp into its nearest bin index:
        transitionAnalysis.results.(temp.currTransitionResultName).numTransitions = length(transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestamps);
        transitionAnalysis.results.(temp.currTransitionResultName).transitionBinTimesteps = round(seconds(transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestamps) ./ active_binning_resolution);
        transitionAnalysis.results.(temp.currTransitionResultName).transitionBinTimestepRange = round(seconds(transitionAnalysis.results.(temp.currTransitionResultName).transitionTimestampRange) ./ active_binning_resolution);

        transitionAnalysis.results.(temp.currTransitionResultName).transitionBinMask = zeros([length(temp.curr_timestamps)-1, 1]);

        % Get the binned spike counts:

        for transitionIndex = 1:transitionAnalysis.results.(temp.currTransitionResultName).numTransitions

            transitionAnalysis.results.(temp.currTransitionResultName).transitionBinIndicies(transitionIndex,:) = transitionAnalysis.results.(temp.currTransitionResultName).transitionBinTimestepRange(transitionIndex,1):transitionAnalysis.results.(temp.currTransitionResultName).transitionBinTimestepRange(transitionIndex,2); % Results in 23x181 array        
            transitionAnalysis.results.(temp.currTransitionResultName).transitionBinMask(transitionAnalysis.results.(temp.currTransitionResultName).transitionBinIndicies(transitionIndex,:)) = 1; % Include the regions within the bounds of the timestep range.
    %         transitionAnalysis.results.(temp.currTransitionResultName).alignedSpikes = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index}.all.binned_spike_counts(;

        end

        % Filter down to the aligned transitionBinMasks for all cells
        transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedSpikeCounts = temp.active_binned_spike_data_matrix(logical(transitionAnalysis.results.(temp.currTransitionResultName).transitionBinMask), :); % 4163x126 double

        % 8177         126
        
        % Getting only the relevant ranges means that the transition ranges are back-to-back aligned for all units. We want to shape them correctly.
        transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedSpikeCounts = reshape(transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedSpikeCounts, length(transitionAnalysis.results.beforeAfterDurationRange), [], 126); % 181x23x126 double

        transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts = squeeze(mean(transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedSpikeCounts, 2)); % Averaged over all the found ranges, preserving transition-offset and unit




         %% Get filter info for active units
        [filter_config.filter_active_units, filter_config.filter_active_unit_original_indicies] = fnFilterUnitsWithCriteria(across_experiment_results{expt_index}.active_processing,...
            across_experiment_results{expt_index}.processing_config.showOnlyAlwaysStableCells,...
            filter_config.filter_included_cell_types, ...
            filter_config.filter_maximum_included_contamination_level);
        fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));


        % Use all units:
    %     temp.activeData = transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts;
        % Use filtered units only:
        temp.activeData = transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts(:, filter_config.filter_active_units);

        % Collapse over unit data:
        transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedMeanSpikeCounts = squeeze(mean(temp.activeData, 2));
    %     transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedMeanSpikeCounts = squeeze(mean(transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts, 2)); 


    %     temp.activeData = transitionAnalysis.results.(temp.currTransitionResultName).alignedBinnedMeanSpikeCounts;
    end
    
    %% Plotting
    figure(19);
    clf;
    
    zero_timestep = find(transitionAnalysis.results.beforeAfterDurationRange == 0);

    temp.labels = sprintf('%s', transitionAnalysis.results.beforeAfterDurationRange);
    temp.handles.meanAlignedSpikeCountsHandle = imagesc(temp.activeData');

    ylabel('Unit Index')
    xlabel('Offset from Transition (Seconds)')
    
    temp.tickIndicies = [1, zero_timestep, length(transitionAnalysis.results.beforeAfterDurationRange)];
    temp.tickValues = seconds(transitionAnalysis.results.beforeAfterDurationRange(temp.tickIndicies));
%     temp.tickLabels = [temp.tickValues 'sec'];
%     temp.tickLabels = {'-18 sec', '0 sec', '18 sec'};
    temp.tickLabels = {'-9 sec', '0 sec', '9 sec'};
    
    xticks(temp.tickIndicies);
    xticklabels(temp.tickLabels)
    title(sprintf('Transitions %s from Any: mean binned spike counts', temp.currTransitionResultName));
    xline(zero_timestep,'-.','Transition', 'Color', 'white', 'LineWidth', 3.5);
    
%     ['rem'], {'nrem','quiet','active'}]
    
    
    %% Original Main Plots:
%     [temp.outputPlottingResults{expt_index}, filtered, results] = fnPerformAcrossREMTesting(across_experiment_results{expt_index}.active_processing, ...
%         across_experiment_results{expt_index}.general_results, ...
%         across_experiment_results{expt_index}.results_array{current_binning_index}, ...
%         across_experiment_results{expt_index}.processing_config, filter_config, expt_info, plottingOptions);
   

%     [temp.outputPlottingResults{expt_index}, filtered, results] = fnPerformAcrossAllTesting(across_experiment_results{expt_index}.active_processing, ...
%             across_experiment_results{expt_index}.general_results, ...
%             across_experiment_results{expt_index}.results_array{current_binning_index}, ...
%             across_experiment_results{expt_index}.processing_config, filter_config, expt_info, plottingOptions);
    
    
end

% General way:
% temp.figureTypesList = fieldnames(temp.outputPlottingResults{expt_index}.figures);
% temp.currFiguresToLayout = cellfun(@(aFigType) temp.outputPlottingResults{expt_index}.figures.(aFigType), temp.figureTypesList);

% temp.figureTypesList = fieldnames(temp.outputPlottingResults{expt_index}.figures);
% temp.currFiguresToLayout = cellfun(@(aFigResults) aFigResults.figures.meanSpikeRateAllPeriodsFig, temp.outputPlottingResults);
% 
% 
% % Align the figures:
% % align_figure(linkedFigureHandles);
% figureLayoutManager.figuresSize.width = 880;
% figureLayoutManager.figuresSize.height = 600;
% figureLayoutManager.verticalSpacing = 30;
% figureLayoutManager.horizontalSpacing = 5;
% 
% align_figure(temp.currFiguresToLayout, 1, figureLayoutManager.figuresSize.width, figureLayoutManager.figuresSize.height,...
%     100, figureLayoutManager.verticalSpacing, figureLayoutManager.horizontalSpacing, 100);





function [plotResults, filtered, results] = fnPerformAcrossREMTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions)
    % A REM specific setup for fnPerformAcrossPeriodTesting 
    plottingOptions.group_name_list = {'pre_sleep_REM', 'post_sleep_REM', 'any_REM', 'all_except_REM'};
    plottingOptions.group_indexArray_variableName_list = strcat(plottingOptions.group_name_list, '_indicies');
    plottingOptions.num_groups = length(plottingOptions.group_name_list);
    plottingOptions.group_included_epochs = {{'pre_sleep'}, {'post_sleep'}, {}, {}};
    plottingOptions.group_included_states = {{'rem'}, {'rem'}, {'rem'}, {'nrem', 'quiet', 'active'}};
    
    [plotResults, filtered, results] = fnPerformAcrossPeriodTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions);
end


function [plotResults, filtered, results] = fnPerformAcrossAllTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions)
    % A all-states specific setup for fnPerformAcrossPeriodTesting 
    plottingOptions.group_name_list = {'nrem', 'rem', 'quiet', 'active', 'all'};
    plottingOptions.group_indexArray_variableName_list = strcat(plottingOptions.group_name_list, '_indicies');
    plottingOptions.num_groups = length(plottingOptions.group_name_list);
    plottingOptions.group_included_epochs = {{}, {}, {}, {}, {}};
    plottingOptions.group_included_states = {{'nrem'}, {'rem'}, {'quiet'}, {'active'}, {}};
    plottingOptions.group_plot_types.meanfiringrate = {true, true, true, true, false};
    plottingOptions.group_plot_types.xcorrheatmap = {false, false, false, false, true};
    
    
    [plotResults, filtered, results] = fnPerformAcrossPeriodTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions);
end






