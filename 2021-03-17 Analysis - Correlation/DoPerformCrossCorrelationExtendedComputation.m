%% All content of Mar17Analysis.mlx
% Pho Hale - 2021-03-22
% Performs extended cross-correlation computations. Shows nothing. Requisite to 'temp_PlotAllXCorrsAllPairs.m'

active_expt_index = 1;
current_binning_index = 1;

active_binning_resolution = across_experiment_results{active_expt_index}.processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = across_experiment_results{active_expt_index}.timesteps_array{current_binning_index};
temp.curr_processed = across_experiment_results{active_expt_index}.active_processing.processed_array{current_binning_index};
active_results = across_experiment_results{active_expt_index}.results_array{current_binning_index};


% Unit Filtering:
filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};
[filter_config.filter_active_units, filter_config.original_unit_index] = fnFilterUnitsWithCriteria(across_experiment_results{active_expt_index}.active_processing, across_experiment_results{active_expt_index}.processing_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);
fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));

% Xcorr lag offsets are given by -9
% active_results.all.pairwise_xcorrelations.lag_offsets
% Find the range 200ms prior and after the 0 point.
temp.curr_xcorr_integration_range.lowerIndicies = find((0 >= active_results.all.pairwise_xcorrelations.lag_offsets) & (active_results.all.pairwise_xcorrelations.lag_offsets >= -0.2));
temp.curr_xcorr_integration_range.upperIndicies = find((0 <= active_results.all.pairwise_xcorrelations.lag_offsets) & (active_results.all.pairwise_xcorrelations.lag_offsets <= 0.2));

temp.curr_xcorr_lag_zero_offset = find(active_results.all.pairwise_xcorrelations.lag_offsets == 0);

[temp.filtered.pre_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(across_experiment_results{active_expt_index}.active_processing, {'pre_sleep'}, {'rem'}); % 668x1
[temp.filtered.post_sleep_REM_indicies] = fnFilterPeriodsWithCriteria(across_experiment_results{active_expt_index}.active_processing, {'post_sleep'}, {'rem'}); % 668x1

[temp.filtered.any_REM_indicies] = fnFilterPeriodsWithCriteria(across_experiment_results{active_expt_index}.active_processing, [], {'rem'}); % 668x1
temp.filtered.all_except_REM_indicies = ~temp.filtered.any_REM_indicies; % 668x1, the complement of the REM indicies


temp.results.any_REM.num_behavioral_periods = sum(temp.filtered.any_REM_indicies,'all');
temp.results.pre_sleep_REM.num_behavioral_periods = sum(temp.filtered.pre_sleep_REM_indicies,'all');
temp.results.post_sleep_REM.num_behavioral_periods = sum(temp.filtered.post_sleep_REM_indicies,'all');
fprintf('any_REM: %d periods\n pre_sleep_REM: %d periods\n post_sleep_REM: %d periods\n', temp.results.any_REM.num_behavioral_periods, temp.results.pre_sleep_REM.num_behavioral_periods, temp.results.post_sleep_REM.num_behavioral_periods);


temp.valid_pairs_reverse_lookup = across_experiment_results{active_expt_index}.general_results.indicies.reverse_lookup_unique_electrode_pairs(filter_config.filter_active_units, filter_config.filter_active_units);
[temp.valid_linear_pair_indicies, ia, ic] = unique(temp.valid_pairs_reverse_lookup);
temp.num_valid_linear_pair_indicies = length(temp.valid_linear_pair_indicies);

temp.valid_units_list = filter_config.original_unit_index;
temp.num_valid_units = length(temp.valid_units_list);

% temp.xcorr_agg_fcn = @(xcorrs_full) sum(xcorrs_full, 1);
% temp.xcorr_agg_fcn = @(xcorrs_full) mean(xcorrs_full, 1);
temp.xcorr_agg_fcn = @(xcorrs_full) max(xcorrs_full, [], 1);


%% Stage 2 - Perform Computations:
% Sum over all num_of_behavioral_state_periods:
temp.by_behavioral_period.curr_xcorr_allPairs.raw = squeeze(temp.xcorr_agg_fcn(active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_full)); % [num_unique_pairs x num_lag_steps] array
% temp.by_behavioral_period.curr_xcorr_allPairs.raw
temp.global_xcorr_pair_max = max(abs(temp.by_behavioral_period.curr_xcorr_allPairs.raw), [], 2); % Max for each pair, should be [num_unique_pairs 1] array
temp.global_xcorr_max = max(temp.global_xcorr_pair_max, [], 'all'); % Absolute global max across all pairs
temp.global_xcorr_zero_offset_values = temp.by_behavioral_period.curr_xcorr_allPairs.raw(:, temp.curr_xcorr_lag_zero_offset); % The zero offset value for each pair, [num_unique_pairs 1] array
temp.global_xcorr_zero_offset_max = max(abs(temp.global_xcorr_zero_offset_values), [], 'all');

% %normalize so that zero-lag has a height of 1
% 

% Normalize by global maximum
temp.by_behavioral_period.curr_xcorr_allPairs.globally_normalized = temp.by_behavioral_period.curr_xcorr_allPairs.raw ./ temp.global_xcorr_max; 
% % Normalize for each pair:
% temp.by_behavioral_period.curr_xcorr_allPairs.normalized = temp.by_behavioral_period.curr_xcorr_allPairs.raw ./ temp.global_xcorr_max; % Normalize by global maximum
% temp.plotVals = 
% temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1



%% Temporal Bias and ByEpoch analysis:

%%% Compute for all **experiment epochs**:
for i = 1:length(active_processing.definitions.behavioral_epoch.classNames)
    temp.curr_epoch_name = active_processing.definitions.behavioral_epoch.classNames{i};
   
    % Perform required normalizations:
    % Get active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr and normalize it if it hasn't been done before.
    if ~isstruct(active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr)
        [temp.output] = fnNormalizeAllPairsXCorr(active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr);
         % Replaces xcorr matrix with structure containing both raw and normalized versions. Original can be accessed by using .raw version.
        active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr = temp.output.xcorr;
    end
    
    temp.xcorr_input =  active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr.raw;
%     temp.xcorr_input =  active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.xcorr.globally_normalized;
    temp.xcorr_input_transformed = fnConvertLinearPairQuantityToUnitByUnit(temp.xcorr_input, general_results, filter_config);
    % Compute the Temporal Bias "B" as defined in SkaggsMcNaughtonScience1996.pdf for this epoch:
    [active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.B, active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.integration_info] = fnTemporalBias_SkaggsMcNaughton(active_results.all.pairwise_xcorrelations.lag_offsets,...
        temp.xcorr_input_transformed, [-0.2, 0.2]);
    
end


% figure
% plot(B)

% Pre-sleep:
% B(temp.filtered.pre_sleep_REM_indicies)
% 
% B(temp.filtered.post_sleep_REM_indicies)

%active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_all_periods = squeeze(mean(active_results.by_behavioral_period.pairwise_xcorrelations.xcorr_full, 1)); % [num_unique_pairs x num_lag_steps] array

% active_results.all.pairwise_xcorrelations.xcorr


%% Testing: Plot the Temporal Bias results:
% active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.B

% temp.B = triu(active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.B);
% temp.B(temp.B == 0) = nan;

temp.plot_pre_track_data.x = [];
temp.plot_pre_track_data.y = [];

temp.plot_post_track_data.x = [];
temp.plot_post_track_data.y = [];

for active_unit_A_index = 1:dim_1.num_valid_units
        
        for active_unit_B_index = 1:dim_2.num_valid_units
            
            % Get current unit bias for each epoch:
             
            temp.curr_pre = active_results.by_epoch.pre_sleep.pairwise_xcorrelations.temporalBias.B(active_unit_A_index, active_unit_B_index);
            temp.curr_track = active_results.by_epoch.track.pairwise_xcorrelations.temporalBias.B(active_unit_A_index, active_unit_B_index);
            temp.curr_post = active_results.by_epoch.post_sleep.pairwise_xcorrelations.temporalBias.B(active_unit_A_index, active_unit_B_index);
            
            if ~isnan(temp.curr_track) & ~isnan(temp.curr_pre)
                temp.plot_pre_track_data.x(end+1) = temp.curr_track;
                temp.plot_pre_track_data.y(end+1) = temp.curr_pre;
            end
            
            if ~isnan(temp.curr_track) & ~isnan(temp.curr_post)
                temp.plot_post_track_data.x(end+1) = temp.curr_track;
                temp.plot_post_track_data.y(end+1) = temp.curr_post;
            end
           
        end % end for dim_2
end % end for dim_1
    
temp.plot_pre_track_data.M = [temp.plot_pre_track_data.x; temp.plot_pre_track_data.y];
temp.plot_post_track_data.M = [temp.plot_post_track_data.x; temp.plot_post_track_data.y];


%% Customize the plotting command to use (stem, plot, area, etc):
% active_plot_cmd = @(ax,x,y) stem(ax, x, y);
% active_plot_cmd = @(ax,x,y) plot(ax, x, y);
active_plot_cmd = @(ax,x,y) scatter(ax, x, y, 'Marker','.','MarkerFaceColor','red');
figure(9);
% clf;
% scatter(temp.plot_pre_track_data.y, ...
%      temp.plot_post_track_data.y, ...
%     'filled')
% 
temp.plot_pre_track_data.ax = subplot(2,1,1);
% plot(reshape(active_results.by_epoch.track.pairwise_xcorrelations.temporalBias.B, [], 1), ...
%     reshape(active_results.by_epoch.pre_sleep.pairwise_xcorrelations.temporalBias.B, [], 1), ...
%     'filled')

% scatter(temp.plot_pre_track_data.ax, temp.plot_pre_track_data.x', ...
%      temp.plot_pre_track_data.y', ...
%     'filled')

active_plot_cmd(temp.plot_pre_track_data.ax, temp.plot_pre_track_data.x', ...
     temp.plot_pre_track_data.y')


temp.plot_pre_track_data.ax.XAxisLocation = 'origin'; 
temp.plot_pre_track_data.ax.YAxisLocation = 'origin';
% axis equal
% daspect([1 1 1])
% set(temp.plot_pre_track_data.ax,'DataAspectRatio',[1 1 1])

xlabel('Bias on track','Interpreter','none')
ylabel('Bias in pre_sleep','Interpreter','none')


temp.plot_post_track_data.ax = subplot(2,1,2);
% scatter(reshape(active_results.by_epoch.track.pairwise_xcorrelations.temporalBias.B, 1,[]), ...
%     reshape(active_results.by_epoch.post_sleep.pairwise_xcorrelations.temporalBias.B, 1,[]), ...
%     'filled')

active_plot_cmd(temp.plot_post_track_data.ax, temp.plot_post_track_data.x, ...
    temp.plot_post_track_data.y)

temp.plot_post_track_data.ax.XAxisLocation = 'origin'; 
temp.plot_post_track_data.ax.YAxisLocation = 'origin';
% axis equal
% daspect([1 1 1])
% set(temp.plot_post_track_data.ax,'DataAspectRatio',[1 1 1])

xlabel('Bias on track','Interpreter','none')
ylabel('Bias in post_sleep','Interpreter','none')

sgtitle('Temporal Bias B')

% %%% Perform for all **experiment epochs**:
% for i = 1:length(active_processing.definitions.behavioral_epoch.classNames)
%     temp.curr_epoch_name = active_processing.definitions.behavioral_epoch.classNames{i};
%    
%     subplot(3,1,i);
%     active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.B
%     
%     % Flatten to a row vector:
%     reshape(active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations.temporalBias.B, 1,[])
% 
%     % Plot each point:
%     plot(
% end