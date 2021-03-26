% Meant to be ran after March24_BuzcodeCCGTest.m produces its results
%% Coming in with ccg_results.by_behavioral_period.ccg.repaired: [668   201   126   126]





%% Get filter info for active units
% [filter_config.filter_active_units, filter_config.filter_active_unit_original_indicies] = fnFilterUnitsWithCriteria(across_experiment_results{expt_index}.active_processing,...
%     across_experiment_results{expt_index}.processing_config.showOnlyAlwaysStableCells,...
%     filter_config.filter_included_cell_types, ...
%     filter_config.filter_maximum_included_contamination_level);
% fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));


%%

% Use all units:
%     temp.activeData = transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts;
% Use filtered units only:
% temp.activeData = transitionAnalysis.results.(temp.currTransitionResultName).byUnit.alignedBinnedMeanSpikeCounts(:, filter_config.filter_active_units);

% A single analysis group represents all the periods whose indicies will be combined (averaged over).
i = 1;
analysisGroups(i).included_epochs = {'pre_sleep'};
analysisGroups(i).included_states = {'rem'};
analysisGroups(i).name = 'pre_sleep_REM';
i = i + 1;
analysisGroups(i).included_epochs = {'pre_sleep'};
analysisGroups(i).included_states = {'rem','nrem'};
analysisGroups(i).name = 'pre_sleep_ALL';
% i = i + 1;
% analysisGroups(i).included_epochs = {'track'};
% analysisGroups(i).included_states = {'rem','nrem'};
i = i + 1;
analysisGroups(i).included_epochs = {'track'};
analysisGroups(i).included_states = {};
analysisGroups(i).name = 'track_ALL';
i = i + 1;
analysisGroups(i).included_epochs = {'post_sleep'};
analysisGroups(i).included_states = {'rem'};
analysisGroups(i).name = 'post_sleep_REM';
i = i + 1;
analysisGroups(i).included_epochs = {'post_sleep'};
analysisGroups(i).included_states = {'rem','nrem'};
analysisGroups(i).name = 'post_sleep_ALL';
num_groups = length(analysisGroups);


[is_period_included, numIncludedPeriods] = fnComputeAnalysisGroupIndexMask(active_processing, analysisGroups);
for i = 1:num_groups
   analysisGroups(i).numIncludedPeriods = numIncludedPeriods(i);
   fprintf('group[%d] %s: periods: %d\n', i, analysisGroups(i).name, numIncludedPeriods(i));
end

% plottingOptions.group_name_list = {'pre_sleep_REM', 'post_sleep_REM', 'any_REM', 'all_except_REM'};
% plottingOptions.group_indexArray_variableName_list = strcat(plottingOptions.group_name_list, '_indicies');
% plottingOptions.num_groups = length(plottingOptions.group_name_list);
% plottingOptions.group_included_epochs = {{'pre_sleep'}, {'post_sleep'}, {}, {}};
% plottingOptions.group_included_states = {{'rem'}, {'rem'}, {'rem'}, {'nrem', 'quiet', 'active'}};
% [plotResults, filtered, results] = fnPerformAcrossPeriodTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions);
    


%% Compute fields for all groups:
for i = 1:num_groups
    temp.curr_group_name = analysisGroups(i).name;
    %% Group Indexing:
    analysisGroups(i).is_period_included = is_period_included(:,i);
    
    % Computations:
    analysisGroupResults.(temp.curr_group_name).num_behavioral_periods = sum(is_period_included(:,i),'all');
    analysisGroupResults.(temp.curr_group_name).per_period.durations = active_processing.behavioral_periods_table.duration(is_period_included(:,i));
    analysisGroupResults.(temp.curr_group_name).per_period.epoch_start_seconds = active_processing.behavioral_periods_table.epoch_start_seconds(is_period_included(:,i));
    analysisGroupResults.(temp.curr_group_name).per_period.epoch_end_seconds = active_processing.behavioral_periods_table.epoch_end_seconds(is_period_included(:,i));

    % Get the CCG results from this group:
    % Average Across all of behavioral periods in this group:
%     analysisGroupResults.(temp.curr_group_name).ccg.mean = squeeze(mean(ccg_results.by_behavioral_period.ccg.repaired(is_period_included(:,i), :, :, :), 1,'omitnan')); % 201x126x126 double
%     analysisGroupResults.(temp.curr_group_name).ccg.stdDev = squeeze(std(ccg_results.by_behavioral_period.ccg.repaired(is_period_included(:,i), :, :, :), 0, 1,'omitnan')); % 14x1 double
    
    
%     %% xcorr_all_lags: averaged over lags, preserving all pairs.
%     analysisGroupResults.(temp.curr_group_name).per_period.xcorr_all_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_lags(is_period_included(:,i), filter_config.filter_active_pairs); % 668x7875
% 
%     %% xcorr_all_pairs:
%     analysisGroupResults.(temp.curr_group_name).per_period.xcorr_all_pairs = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs(is_period_included(:,i), :); % 668x81
%     %% xcorr_all_pairs_AND_lags:
%     analysisGroupResults.(temp.curr_group_name).per_period.xcorr_all_pairs_AND_lags = active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_pairs_AND_lags(is_period_included(:,i), :); % 668x1

    [analysisGroupResults.(temp.curr_group_name).temporalBias.B, analysisGroupResults.(temp.curr_group_name).temporalBias.integration_info] = fnTemporalBias_SkaggsMcNaughton(ccg_results.lag_offsets,...
            analysisGroupResults.(temp.curr_group_name).ccg.mean, ...
            [-0.2, 0.2]); %201x126 double
    
    fprintf('%s: %d periods\n', temp.curr_group_name, analysisGroupResults.(temp.curr_group_name).num_behavioral_periods);
end


% analysisGroupResults.(temp.curr_group_name).ccg.mean: 201x126x126


    
function [is_period_included, numIncludedPeriods] = fnComputeAnalysisGroupIndexMask(active_processing, analysisGroups)
    %% fnComputeAnalysisGroupIndexMask: Takes a struct array of groups, which specify which epochs and states that should included, and returns the index mask for whether each of the behavioral periods should be included in that group.
    % analysisGroups: a struct array with fields 'included_epochs' and 'included_states', one entry for each group.

    %% Output:
    % is_period_include: each column contains whether that behavioral period was included in this group.
    num_groups = length(analysisGroups);
    num_of_behavioral_state_periods = height(active_processing.behavioral_periods_table);
    is_period_included = true(num_of_behavioral_state_periods, num_groups);    
    for i = 1:num_groups
        [is_period_included(:,i)] = fnFilterPeriodsWithCriteria(active_processing, analysisGroups(i).included_epochs, analysisGroups(i).included_states);
        
    end
    numIncludedPeriods = sum(is_period_included, 1);
end
    

function [plotResults, filtered, results] = fnPerformAcrossREMTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions)
    % A REM specific setup for fnPerformAcrossPeriodTesting 
    plottingOptions.group_name_list = {'pre_sleep_REM', 'post_sleep_REM', 'any_REM', 'all_except_REM'};
    plottingOptions.group_indexArray_variableName_list = strcat(plottingOptions.group_name_list, '_indicies');
    plottingOptions.num_groups = length(plottingOptions.group_name_list);
    plottingOptions.group_included_epochs = {{'pre_sleep'}, {'post_sleep'}, {}, {}};
    plottingOptions.group_included_states = {{'rem'}, {'rem'}, {'rem'}, {'nrem', 'quiet', 'active'}};
    
    [plotResults, filtered, results] = fnPerformAcrossPeriodTesting(active_processing, general_results, active_results, processing_config, filter_config, expt_info, plottingOptions);
end
