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

% Loop through each experiment:
for expt_index = 1:active_num_experiments
    expt_info.index = expt_index;
    expt_info.name = temp.curr_active_experiment_names{expt_index};
    
    temp.curr_timestamps = across_experiment_results{expt_index}.timesteps_array{current_binning_index};
    temp.curr_processed = across_experiment_results{expt_index}.active_processing.processed_array{current_binning_index};
    active_results = across_experiment_results{expt_index}.results_array{current_binning_index};

    
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
    
    
%     ['rem'], {'nrem','quiet','active'}]
    
    
    

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




function [plotResults, filtered, results] = fnPerformAcrossPeriodTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions)
    %%% fnPerformAcrossREMTesting: Run main analysis
    %%%
 
    plotResults.exports = {};
    plotResults.configs = {};
    
    %% Get filter info for active units
    [filter_config.filter_active_units, filter_config.filter_active_unit_original_indicies] = fnFilterUnitsWithCriteria(active_processing, processing_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
        filter_config.filter_maximum_included_contamination_level);
    fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));
    
    temp.curr_pairs_indicies = find(filter_config.filter_active_units);
    filter_config.filter_active_pairs =  ismember(general_results.indicies.unique_electrode_pairs(:,1), temp.curr_pairs_indicies) & ismember(general_results.indicies.unique_electrode_pairs(:,2), temp.curr_pairs_indicies);
    filter_config.filter_active_pair_values = general_results.indicies.unique_electrode_pairs(filter_config.filter_active_pairs, :);
    
    
    %% Compute fields for all groups:
    for i = 1:plottingOptions.num_groups
        temp.curr_group_name = plottingOptions.group_name_list{i};
        
        %% Group Indexing:
        [filtered.(plottingOptions.group_indexArray_variableName_list{i})] = fnFilterPeriodsWithCriteria(active_processing, plottingOptions.group_included_epochs{i}, plottingOptions.group_included_states{i}); % 668x1
        temp.curr_group_indicies = filtered.(plottingOptions.group_indexArray_variableName_list{i});
        
        % Computations:
        results.(temp.curr_group_name).num_behavioral_periods = sum(temp.curr_group_indicies,'all');
        results.(temp.curr_group_name).per_period.durations = active_processing.behavioral_periods_table.duration(temp.curr_group_indicies);
        results.(temp.curr_group_name).per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(temp.curr_group_indicies);
        results.(temp.curr_group_name).per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(temp.curr_group_indicies);

        % Compute the center of the epochs to plot the firing rates along an appropriately scaled x-axis:
        results.(temp.curr_group_name).per_period.epoch_center_seconds = (results.(temp.curr_group_name).per_period.epoch_start_seconds + floor(results.(temp.curr_group_name).per_period.durations ./ 2.0));

         % Leave in terms of the spike rates per unit (14x92 double):
        results.(temp.curr_group_name).spike_rate_per_unit = general_results.per_behavioral_state_period.spike_rate_per_unit(temp.curr_group_indicies, filter_config.filter_active_units);

        % Average Across all of the units
        results.(temp.curr_group_name).spike_rate_all_units.mean = mean(results.(temp.curr_group_name).spike_rate_per_unit, 2); % 14x1 double
        results.(temp.curr_group_name).spike_rate_all_units.stdDev = std(results.(temp.curr_group_name).spike_rate_per_unit, 0, 2); % 14x1 double

        % Compute the average across the REM sessions in each epoch
        results.(temp.curr_group_name).baseline_spike_rate_across_all.mean = mean(results.(temp.curr_group_name).spike_rate_all_units.mean);
        results.(temp.curr_group_name).baseline_spike_rate_across_all.stdDev = std(results.(temp.curr_group_name).spike_rate_all_units.mean);
    
        %% xcorr_all_lags: averaged over lags, preserving all pairs.
        results.(temp.curr_group_name).per_period.xcorr_all_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_lags(temp.curr_group_indicies, filter_config.filter_active_pairs); % 668x7875
        
        %% xcorr_all_pairs:
        results.(temp.curr_group_name).per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(temp.curr_group_indicies, :); % 668x81
        %% xcorr_all_pairs_AND_lags:
        results.(temp.curr_group_name).per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(temp.curr_group_indicies, :); % 668x1

        fprintf('%s: %d periods\n', temp.curr_group_name, results.(temp.curr_group_name).num_behavioral_periods);
        
    end
    
%     fprintf('any_REM: %d periods\n pre_sleep_REM: %d periods\n post_sleep_REM: %d periods\n', results.any_REM.num_behavioral_periods, results.pre_sleep_REM.num_behavioral_periods, results.post_sleep_REM.num_behavioral_periods);

    %%% Plotting Results:
    temp.curr_expt_string = sprintf('experiment[%d]: %s', expt_info.index, expt_info.name);
    plottingOptions.outputs.curr_expt_filename_string = sprintf('%s_', expt_info.name);
    
    
    %% TODO: Plot all others (states that aren't REM at all) as a separate series

    %% Plot Mean Firing Rates across units FOR ALL PERIODS:

    plotResults.figures.meanSpikeRateAllPeriodsFig = figure(119+expt_info.index);
    clf
    hold off;
    
    for i = 1:plottingOptions.num_groups
        
        if plottingOptions.group_plot_types.meanfiringrate{i} 
            temp.curr_group_name = plottingOptions.group_name_list{i};
            temp.curr_group_indicies = filtered.(plottingOptions.group_indexArray_variableName_list{i});

            if strcmpi(plottingOptions.plottingXAxis, 'index')
                plottingOptions.x = [1:results.(temp.curr_group_name).num_behavioral_periods];
            else
                plottingOptions.x = results.(temp.curr_group_name).per_period.epoch_center_seconds;
            end


            [h0] = fnPlotAcrossREMTesting(plottingOptions.plottingMode, plottingOptions.x, ...
                results.(temp.curr_group_name).spike_rate_all_units.mean, ...
                results.(temp.curr_group_name).spike_rate_all_units.stdDev, ...
                results.(temp.curr_group_name).spike_rate_per_unit);

            if ~strcmpi(plottingOptions.plottingMode, 'distributionPlot')
                %     h0(1).Marker = '*';
                h0(1).DisplayName = temp.curr_group_name;
            end

            hold on;
        
        end % end if plottingOptions.group_plot_types.meanfiringrate

    end
    
    
    % Set up axis properties:
    title(sprintf('Firing Rate All periods: %d', results.(temp.curr_group_name).num_behavioral_periods));
    if strcmpi(plottingOptions.plottingXAxis, 'index')
        xlabel('Filtered Period Index')
    else
        xlabel('Period Timestamp Offset (Seconds)')
    end
    ylabel('mean spike rate')
    if ~isempty(plottingOptions.plottingYlim)
        ylim(plottingOptions.plottingYlim)
    end

    legend();
    sgtitle([temp.curr_expt_string ' : Spike Rates - All Periods - Period Index - All Cells'])
    
    
    
    
%     
%     %% Plot Mean Firing Rates across units FOR REM PERIODS ONLY:
%     % Error bars are across units:
%     plotResults.figures.meanSpikeRateFig = figure(9+expt_info.index);
%     clf
%     subplot(2,1,1);
% 
%     if strcmpi(plottingOptions.plottingXAxis, 'index')
%         plottingOptions.x = [1:results.pre_sleep_REM.num_behavioral_periods];
%     else
%         plottingOptions.x = results.pre_sleep_REM.per_period.epoch_center_seconds;
%     end
% 
%     [h1] = fnPlotAcrossREMTesting(plottingOptions.plottingMode, plottingOptions.x, ...
%         results.pre_sleep_REM.spike_rate_all_units.mean, ...
%         results.pre_sleep_REM.spike_rate_all_units.stdDev, ...
%         results.pre_sleep_REM.spike_rate_per_unit);
% 
%     title(sprintf('PRE sleep REM periods: %d', results.pre_sleep_REM.num_behavioral_periods));
%     if strcmpi(plottingOptions.plottingXAxis, 'index')
%         xlabel('Filtered Period Index')
%     else
%         xlabel('Period Timestamp Offset (Seconds)')
%     end
%     ylabel('mean spike rate')
%     if ~isempty(plottingOptions.plottingYlim)
%         ylim(plottingOptions.plottingYlim)
%     end
% 
%     subplot(2,1,2);
% 
%     if strcmpi(plottingOptions.plottingXAxis, 'index')
%         plottingOptions.x = [1:results.post_sleep_REM.num_behavioral_periods];
%     else
%         plottingOptions.x = results.post_sleep_REM.per_period.epoch_center_seconds;
%     end
%     [h2] = fnPlotAcrossREMTesting(plottingOptions.plottingMode, plottingOptions.x, ...
%         results.post_sleep_REM.spike_rate_all_units.mean, ...
%         results.post_sleep_REM.spike_rate_all_units.stdDev, ...
%         results.post_sleep_REM.spike_rate_per_unit);
% 
%     title(sprintf('POST sleep REM periods: %d', results.post_sleep_REM.num_behavioral_periods));
%     if strcmpi(plottingOptions.plottingXAxis, 'index')
%         currPlotConfig.xlabel = 'Filtered Period Index';
%         
%     else
%         currPlotConfig.xlabel = 'Period Timestamp Offset (Seconds)';
% 
%     end
%     currPlotConfig.ylabel = 'mean spike rate';
%     
%     xlabel(currPlotConfig.xlabel)
%     ylabel(currPlotConfig.ylabel)
%     if ~isempty(plottingOptions.plottingYlim)
%         ylim(plottingOptions.plottingYlim)
%     end
%     sgtitle([temp.curr_expt_string ' : Spike Rates - PRE vs Post Sleep REM Periods - Period Index - Pyramidal Only'])
%     % Figure Name:
%     %'Spike Rates - PRE vs Post Sleep REM Periods - Period Index';
%     %'Spike Rates - PRE vs Post Sleep REM Periods - Timestamp Offset';
% 
% %     'Spike Rates - PRE vs Post Sleep REM Periods - Period Index - Pyramidal Only'
%     % MainPlotVariableBeingCompared - PurposeOfComparison - IndependentVariable - FiltersAndConstraints
%     
%     % Build Figure Export File path:
%     currPlotConfig.curr_expt_filename_string = sprintf('%s - %s - %s - %s - %s', ...
%         plottingOptions.outputs.curr_expt_filename_string, ...
%         'Spike Rates', ...
%         'PRE vs Post Sleep REM Periods', ...
%         currPlotConfig.xlabel, ...
%         'Spike Rates'...
%         );
% 
%     currPlotConfig.curr_expt_parentPath = plottingOptions.outputs.rootPath;
%     currPlotConfig.curr_expt_path = fullfile(currPlotConfig.curr_expt_parentPath, currPlotConfig.curr_expt_filename_string);
% 
%     % Perform the export:
%     [plotResults.exports{end+1}.export_result] = fnSaveFigureForExport(plotResults.figures.meanSpikeRateFig, currPlotConfig.curr_expt_path, true, false, false, true);
%     plotResults.configs{end+1} = currPlotConfig;
    
    %% Display the Correlational Results:
    plottingOptions.curr_expt_string = temp.curr_expt_string;
    
    %% Do for all indexes:
    temp.curr_units_to_test = filter_config.filter_active_unit_original_indicies;

%     temp.curr_units_to_test = filter_config.filter_active_unit_original_indicies(1:4);
    [temp.is_pair_included, temp.original_pair_index] = fnFilterPairsWithCriteria(general_results, temp.curr_units_to_test);
    % original_pair_index: 126x86 double
    
    filtered.is_pair_included = temp.is_pair_included(filter_config.filter_active_pairs, :); % 3655x86 logical
    filtered.original_pair_index = temp.original_pair_index(filter_config.filter_active_units, :); % 86x86 double
    
        % Sanity Check:
%     sum(temp.is_pair_included, 1);
%     sum(filtered.is_pair_included, 1);
    
   
    for i = 1:length(temp.curr_units_to_test)
        temp.curr_reference_unit_index = temp.curr_units_to_test(i);
        plotResults.figures.xcorr = fnPlotCellPairsByPeriodHeatmap(active_processing, general_results, results, filter_config, plottingOptions, temp.curr_reference_unit_index);
        
        % Build Figure Export File path:
        currPlotConfig.curr_expt_filename_string = sprintf('%s - %s - %s - %s - %s', ...
            plottingOptions.outputs.curr_expt_filename_string, ...
            'XCorr for all lags', ...
            'All Periods', ...
            'cellPairs', ...
            sprintf('unit[%d]', temp.curr_reference_unit_index));
    
        currPlotConfig.curr_expt_parentPath = fullfile(plottingOptions.outputs.rootPath, 'png');
        if ~exist(currPlotConfig.curr_expt_parentPath, 'dir')
           mkdir(currPlotConfig.curr_expt_parentPath); 
        end
        currPlotConfig.curr_expt_path = fullfile(currPlotConfig.curr_expt_parentPath, currPlotConfig.curr_expt_filename_string);

        % Perform the export:
        [plotResults.exports{end+1}.export_result] = fnSaveFigureForExport(plotResults.figures.xcorr, currPlotConfig.curr_expt_path, false, false, false, true);
        plotResults.configs{end+1} = currPlotConfig;

    end
    
end


function [xcorr_fig] = fnPlotCellPairsByPeriodHeatmap(active_processing, general_results, results, filter_config, plottingOptions, reference_unit_index)
    %% fnPlotCellPairsByPeriodHeatmap: plots a heatmap for a given unit with cell pair on its x-axis and period on the y-axis
        % reference_unit_index: a valid unit index in the set of pairs
    
%     plotResults.figures.xcorr = figure(expt_info.index);
   
    temp.curr_ref_unit_index_string = num2str(reference_unit_index);
    
    temp.curr_active_unit_index = reference_unit_index;
%     temp.curr_active_unit_index = filter_config.filter_active_pair_values(reference_unit_index, 1);
%     temp.curr_active_pair_indicies = (temp.curr_active_unit_index == filter_config.filter_active_pair_values(:,1));
%     temp.curr_active_pair_values = filter_config.filter_active_pair_values(temp.curr_active_pair_indicies,:);
%     temp.curr_active_pair_labels = num2str(temp.curr_active_pair_values(:,2));
    
    %% New
    [temp.is_pair_included, temp.original_pair_index] = fnFilterPairsWithCriteria(general_results, temp.curr_active_unit_index);
    % original_pair_index: 126x1 double
    filtered.is_pair_included = temp.is_pair_included(filter_config.filter_active_pairs, 1); % 3655x1 logical
    filtered.original_pair_index = temp.original_pair_index(filter_config.filter_active_units, 1); % 86x1 double
    
    temp.curr_active_pair_indicies = find(filtered.is_pair_included);
    temp.curr_active_pair_values = filter_config.filter_active_pair_values(temp.curr_active_pair_indicies, :);
    temp.curr_active_pair_other_unit_index = ones([size(temp.curr_active_pair_values, 1), 1]);
    temp.curr_active_pair_other_unit_values = zeros([size(temp.curr_active_pair_values, 1), 1]);
    
    % Find the value in the pair that DOESN'T correspond to the active index (to get the other relevant unit)
    for i = 1:size(temp.curr_active_pair_values, 1)
        if (temp.curr_active_pair_values(i,1) == temp.curr_active_unit_index)
            temp.curr_active_pair_other_unit_index(i) = 2; % Use the other index for the value;
            temp.curr_active_pair_other_unit_values(i) = temp.curr_active_pair_values(i,2); 
        else
            temp.curr_active_pair_other_unit_index(i) = 1; % Use this index for the value;
            temp.curr_active_pair_other_unit_values(i) = temp.curr_active_pair_values(i,1);
        end
    end
    
    temp.curr_active_pair_labels = num2str(temp.curr_active_pair_other_unit_values);
    
    [xcorr_fig, heatmap_handle] = fnPlotAcrossREMXcorrHeatmap(results.all.per_period.xcorr_all_lags(:, temp.curr_active_pair_indicies)); % 668x3655 double

    heatmap_axis = gca;    
    xlabel('Cell Pairs')
    
%     xticklabels('manual');
    xticks(1:length(temp.curr_active_pair_labels));
    xticklabels(temp.curr_active_pair_labels);
    xticklabels('manual')
    xtickangle(90)

    [state_ax, epoch_ax] = fnPlotHelper_AddStateMapSubplot(active_processing, heatmap_axis);

    temp.curr_ref_unit_index_string = sprintf('unit[%d]', reference_unit_index);
    sgtitle([plottingOptions.curr_expt_string ' : XCorr for ' temp.curr_ref_unit_index_string ' and all units - All Periods - Good Units'])
    

end



function [xcorr_fig, handles] = fnPlotAcrossREMXcorrHeatmap(varargin)
    % Plot a heatmap of the xcorr    
%     xcorr_fig = figure(15);
    xcorr_fig = gcf;
    clf
 
    num_heatmaps = nargin;   
    
    for i = 1:num_heatmaps
        subplot(num_heatmaps,1,i);
        
%         handles(i) = heatmap(varargin{i},'GridVisible', false);
        num_timesteps = size(varargin{i},2);
        zero_timestep = num_timesteps / 2;
        
        handles(i) = imagesc(varargin{i});
        
        ylabel('Filtered Trial Index')
        
%         title(sprintf('Periods: %d', size(varargin{i},1)));
        
%         xlabel('Time Lag')
%         xline(zero_timestep, 'red');



        
        
        
    end
    
%     h1 = heatmap(v1);
%     ylabel('Filtered Trial Index')
%     xlabel('Time Lag')
%     title(sprintf('PRE sleep REM periods: %d', size(v1,1)));
% 
% 
%     subplot(num_heatmaps,1,2);
%     h2 = heatmap(v2);
%     ylabel('Filtered Trial Index')
%     xlabel('Time Lag')    
%     title(sprintf('POST sleep REM periods: %d', size(v2,1)));

%     sgtitle('XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units')

end

function [h] = fnPlotAcrossREMTesting(mode, v1, v2, v3, v4)
    % v1: 36x1
    % v4: 36x86 double
    if strcmpi(mode, 'errorbar')
        h = errorbar(v1, ...
            v2, ...
            v3);

    elseif strcmpi(mode, 'scatter')

        num_repeats = size(v4, 2);
%         x = repmat(v1,[num_repeats 1]);
        x = repelem(v1, num_repeats);
        
        y = reshape(v4,[],1);
        h(2) = scatter(x, y, '.','k','MarkerFaceColor','k');
        h(2).AlphaData = repelem(0.3, length(x));
        h(2).MarkerFaceAlpha = 'flat';
        
%         h(2) = errorbar(v1, v2, v3, 'LineStyle','none');

        hold on;
        
        h(1) = scatter(v1, ...
            v2, 'filled');
        
        h(1).AlphaData = repelem(0.8, length(v1));
        h(1).MarkerFaceAlpha = 'flat';
        
        
    elseif strcmpi(mode, 'distributionPlot')
        % v4: 36x86 double
        h = distributionPlot(v4', 'xValues', v1); % defaults 

    elseif strcmpi(mode, 'bar')
        h = bar(v1, v2);
        
    elseif strcmpi(mode, 'stem')
        h = stem(v1, v2);
    else
       error('Invalid mode input!') 
    end
end

