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
analysisGroups(i).included_states = {'nrem'};
analysisGroups(i).name = 'pre_sleep_NREM';
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
analysisGroups(i).included_states = {'nrem'};
analysisGroups(i).name = 'post_sleep_NREM';
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
    if ~isfield(analysisGroupResults.(temp.curr_group_name), 'ccg')
        % Don't overwrite extant results:
        analysisGroupResults.(temp.curr_group_name).ccg.mean = squeeze(mean(ccg_results.by_behavioral_period.ccg.repaired(is_period_included(:,i), :, :, :), 1,'omitnan')); % 201x126x126 double
        analysisGroupResults.(temp.curr_group_name).ccg.stdDev = squeeze(std(ccg_results.by_behavioral_period.ccg.repaired(is_period_included(:,i), :, :, :), 0, 1,'omitnan')); % 14x1 double
    else
        warning('skipping computation of ccg.mean because the value already exists');
    end
    
    [analysisGroupResults.(temp.curr_group_name).temporalBias.B, analysisGroupResults.(temp.curr_group_name).temporalBias.integration_info] = fnTemporalBias_SkaggsMcNaughton(ccg_results.lag_offsets,...
            analysisGroupResults.(temp.curr_group_name).ccg.mean, ...
            [-0.2, 0.2],...
            1); %201x126 double
    
    fprintf('%s: %d periods\n', temp.curr_group_name, analysisGroupResults.(temp.curr_group_name).num_behavioral_periods);
end

% analysisGroupResults.(temp.curr_group_name).ccg.mean: 201x126x126

%% Plot the Temporal Bias information:

% Unit Filtering:
% filter_config.filter_included_cell_types = {};
% filter_config.filter_included_cell_types = {'pyramidal'};
filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};
[filter_config.filter_active_units, filter_config.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, processing_config.showOnlyAlwaysStableCells, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);
fprintf('Filter: Including %d of %d total units\n', sum(filter_config.filter_active_units, 'all'), length(filter_config.filter_active_units));

temp.dest_dir = '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/ResultsTemp/March 26/';



% '/Users/pho/Dropbox/Classes/Spring 2021/PIBS 600 - Rotations/Rotation_3_Kamran Diba Lab/DataProcessingProject/ResultsTemp/March 26/Temporal Bias B - All Cells.fig'


%% Build plotting options:
plotting_options.filter_config = filter_config;
% plotting_options.active_plot_cmd = @(ax,x,y) stem(ax, x, y);
% plotting_options.active_plot_cmd = @(ax,x,y) plot(ax, x, y);
plotting_options.active_plot_cmd = @(ax,x,y) scatter(ax, x, y, 'Marker','.','MarkerFaceColor','red');

[filter_descriptions] = fnGenerateFilterDescriptionString(plotting_options.filter_config);
% Add the filter_descriptions strings to the plotting_options structure:
plotting_options = fnSetStructureFields({'figure_name', 'additional_title'}, ...
        {sprintf('March 26 - %s', filter_descriptions.cells), filter_descriptions.cells}, ...
        plotting_options);
    
    
% REM Only Sleep:
[temp.fig, temp.h, temp.info] = fnBuildTemporalBiasPlot(analysisGroupResults.pre_sleep_REM.temporalBias, ...
        analysisGroupResults.track_ALL.temporalBias, ...
        analysisGroupResults.post_sleep_REM.temporalBias, ...
        plotting_options);

% %% NREM Only Sleep:
% [temp.fig, temp.h, temp.info] = fnBuildTemporalBiasPlot(analysisGroupResults.pre_sleep_NREM.temporalBias, ...
%         analysisGroupResults.track_ALL.temporalBias, ...
%         analysisGroupResults.post_sleep_NREM.temporalBias, ...
%         plotting_options);

    
% %% All Sleep (REM and NREM):    
% [temp.fig, temp.h, temp.info] = fnBuildTemporalBiasPlot(analysisGroupResults.pre_sleep_ALL.temporalBias, ...
%         analysisGroupResults.track_ALL.temporalBias, ...
%         analysisGroupResults.post_sleep_ALL.temporalBias, ...
%         plotting_options);
    

output_filename = strrep(temp.info.fig_title,':','-');
output_filepath = fullfile(temp.dest_dir, output_filename);



savefig(temp.fig, output_filepath, 'compact');
% Requires R2020a or later
exportgraphics(temp.fig,[output_filepath '.pdf'],'ContentType','vector')    

    
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
    


function [fig, h, info] = fnBuildTemporalBiasPlot(preTemporalBias, trackTemporalBias, postTemporalBias, plotting_options, extantFigH)
    % preTemporalBias, trackTemporalBias, postTemporalBias: structs containing at least the field 'B' which contains the temporal bias matrix computed via fnTemporalBias_SkaggsMcNaughton
    % These are each of size numUnits x numUnits
    
    if ~exist('extantFigH','var')
        extantFigH = [];
    end
    if isfield(plotting_options, 'filter_config')
        % Include all units if no filter
        preTemporalBias.B = preTemporalBias.B(plotting_options.filter_config.filter_active_units, plotting_options.filter_config.filter_active_units);
        trackTemporalBias.B = trackTemporalBias.B(plotting_options.filter_config.filter_active_units, plotting_options.filter_config.filter_active_units);
        postTemporalBias.B = postTemporalBias.B(plotting_options.filter_config.filter_active_units, plotting_options.filter_config.filter_active_units);
    end
    
    % size(analysisGroupResults.track_ALL.temporalBias.B): 201x126
    
    
    plot_pre_track_data.x = [];
    plot_pre_track_data.y = [];

    plot_post_track_data.x = [];
    plot_post_track_data.y = [];

    dim_1.num_valid_units = size(preTemporalBias.B, 1);
    dim_2.num_valid_units = size(preTemporalBias.B, 2);

%     if ~isfield(plotting_options, 'filter_config')
%         % Include all units if no filter
%         plotting_options.filter_config.filter_active_units = ones([dim_2.num_valid_units 1]);
%     else
%         
%         dim_2.num_valid_units = sum(plotting_options.filter_config.filter_active_units, 'all');
%     end
    

%     [info.rho.plot_pre_track_data, info.p.plot_pre_track_data] = corr(trackTemporalBias.B, preTemporalBias.B, 'type', 'Spearman','Rows','complete')
%     [info.rho.plot_post_track_data, info.p.plot_post_track_data] = corr(trackTemporalBias.B, postTemporalBias.B, 'type', 'Spearman','Rows','complete')
%     
%     plot_bias_delta_data.y = (plot_post_track_data.y - plot_pre_track_data.y);
%     [info.rho.plot_bias_delta_data, info.p.plot_bias_delta_data] = corr(trackTemporalBias.B, plot_bias_delta_data.y, 'type', 'Spearman','Rows','complete')
    
    
    for active_unit_A_index = 1:dim_1.num_valid_units

            for active_unit_B_index = 1:dim_2.num_valid_units

                % Get current unit bias for each epoch:
                curr_pre = preTemporalBias.B(active_unit_A_index, active_unit_B_index);
                curr_track = trackTemporalBias.B(active_unit_A_index, active_unit_B_index);
                curr_post = postTemporalBias.B(active_unit_A_index, active_unit_B_index);

                if ~isnan(curr_track) & ~isnan(curr_pre)
                    plot_pre_track_data.x(end+1) = curr_track;
                    plot_pre_track_data.y(end+1) = curr_pre;
                end

                if ~isnan(curr_track) & ~isnan(curr_post)
                    plot_post_track_data.x(end+1) = curr_track;
                    plot_post_track_data.y(end+1) = curr_post;
                end

            end % end for dim_2
    end % end for dim_1
    
    %% Determine ranges:
    plotting_options.lims.xrange = [min(plot_pre_track_data.x), max(plot_pre_track_data.x)];
    plotting_options.lims.yrange = [min([plot_post_track_data.y, plot_pre_track_data.y]), max([plot_post_track_data.y, plot_pre_track_data.y])];
    
%     [info.rho.plot_pre_track_data, info.p.plot_pre_track_data] = corr(plot_pre_track_data.x, plot_pre_track_data.y, 'type', 'Spearman','Rows','complete')
%     [info.rho.plot_post_track_data, info.p.plot_post_track_data] = corr(plot_post_track_data.x, plot_post_track_data.y, 'type', 'Spearman','Rows','complete')
% %     
%     plot_bias_delta_data.y = (plot_post_track_data.y - plot_pre_track_data.y);
%     [info.rho.plot_bias_delta_data, info.p.plot_bias_delta_data] = corr(plot_post_track_data.x, plot_bias_delta_data.y, 'type', 'Spearman','Rows','complete')
    
%     disp(info.rho);
%     title('spearman \rho = %d', info.rho.plot_bias_delta_data)

    [fig, h, info] = fnPerformPlotTemporalBias(plot_pre_track_data, plot_post_track_data, plotting_options, extantFigH);
end


function [figH, h, info] = fnPerformPlotTemporalBias(plot_pre_track_data, plot_post_track_data, plotting_options, extantFigH)
% fnPlotTemporalBias: Plot the Temporal Bias results
    % Called only by fnBuildTemporalBiasPlot(...) above
    % Customize the plotting command to use (stem, plot, area, etc):
    num_subplot_rows = 3;
    
    if ~exist('plotting_options','var')
        plotting_options = struct();
    end
    
    plotting_options = fnAddDefaultOptionalArgs({'figure_name', 'additional_title'}, ...
        {'Untitled', ''}, ...
        plotting_options);
    
    
    % Either reuse ths specified figure or create a new figure
    if ~exist('extantFigH','var') | isempty(extantFigH)
        figH = createFigureWithNameIfNeeded(['Temporal Bias B Figure - ' plotting_options.figure_name]); % generate a new figure to plot the sessions.
    else
        figH = extantFigH; % use the existing provided figure    
        figure(figH);
    end
    
    clf(figH);
    
    [info.rho.plot_pre_track_data, info.p.plot_pre_track_data] = corr(plot_pre_track_data.x', plot_pre_track_data.y', 'type', 'Spearman','Rows','complete');
    [info.rho.plot_post_track_data, info.p.plot_post_track_data] = corr(plot_post_track_data.x', plot_post_track_data.y', 'type', 'Spearman','Rows','complete');
    
%     %% Direct Pre vs. Post Plot:
%     [info.rho.plot_pre_post_data, info.p.plot_pre_post_data] = corr(plot_pre_track_data.y', plot_post_track_data.y', 'type', 'Spearman','Rows','complete');
%     
%     h.plot_pre_post_data.ax = subplot(1,1,1);
%     plotting_options.active_plot_cmd(h.plot_pre_post_data.ax, plot_pre_track_data.y', ...
%          plot_post_track_data.y');
%      
%     h.plot_pre_post_data.ax.XAxisLocation = 'origin'; 
%     h.plot_pre_post_data.ax.YAxisLocation = 'origin';
%     xlabel('Bias on pre_sleep','Interpreter','none')
%     ylabel('Bias in post_sleep','Interpreter','none')
%     title('spearman \rho = %d', info.rho.plot_pre_post_data)
    
    

    %% Conventional vs. Track plots:
    h.plot_pre_track_data.ax = subplot(num_subplot_rows,1,1);

    plotting_options.active_plot_cmd(h.plot_pre_track_data.ax, plot_pre_track_data.x', ...
         plot_pre_track_data.y')


    h.plot_pre_track_data.ax.XAxisLocation = 'origin'; 
    h.plot_pre_track_data.ax.YAxisLocation = 'origin';
    % axis equal
    % daspect([1 1 1])
    % set(temp.plot_pre_track_data.ax,'DataAspectRatio',[1 1 1])

    % Set the axis limits if they're specified:
    if isfield(plotting_options, 'lims')
       xlim(h.plot_pre_track_data.ax, plotting_options.lims.xrange);
       ylim(h.plot_pre_track_data.ax, plotting_options.lims.yrange);
    end
    xlabel('Bias on track','Interpreter','none')
    ylabel('Bias in pre_sleep','Interpreter','none')
    title('spearman \rho = %d', info.rho.plot_pre_track_data)
    
    
    h.plot_post_track_data.ax = subplot(num_subplot_rows,1,2);

    plotting_options.active_plot_cmd(h.plot_post_track_data.ax, plot_post_track_data.x, ...
        plot_post_track_data.y)

    h.plot_post_track_data.ax.XAxisLocation = 'origin'; 
    h.plot_post_track_data.ax.YAxisLocation = 'origin';
    % axis equal
    % daspect([1 1 1])
    % set(temp.plot_post_track_data.ax,'DataAspectRatio',[1 1 1])
    % Set the axis limits if they're specified:
    if isfield(plotting_options, 'lims')
       xlim(h.plot_pre_track_data.ax, plotting_options.lims.xrange);
       ylim(h.plot_pre_track_data.ax, plotting_options.lims.yrange);
    end
    
    xlabel('Bias on track','Interpreter','none')
    ylabel('Bias in post_sleep','Interpreter','none')
    title('spearman \rho = %d', info.rho.plot_post_track_data)
    
    if num_subplot_rows > 2
        %% Change in Bias Subplot:
        h.plot_bias_delta_data.ax = subplot(num_subplot_rows,1,3);

        plot_bias_delta_data.y = (plot_post_track_data.y - plot_pre_track_data.y);
        [info.rho.plot_bias_delta_data, info.p.plot_bias_delta_data] = corr(plot_post_track_data.x', plot_bias_delta_data.y', 'type', 'Spearman');
        plotting_options.active_plot_cmd(h.plot_bias_delta_data.ax, plot_post_track_data.x, ...
            plot_bias_delta_data.y)

        h.plot_bias_delta_data.ax.XAxisLocation = 'origin'; 
        h.plot_bias_delta_data.ax.YAxisLocation = 'origin';
        % axis equal
        % daspect([1 1 1])
        % set(temp.plot_post_track_data.ax,'DataAspectRatio',[1 1 1])
        % Set the axis limits if they're specified:
        if isfield(plotting_options, 'lims')
           xlim(h.plot_bias_delta_data.ax, plotting_options.lims.xrange);
           ylim(h.plot_bias_delta_data.ax, plotting_options.lims.yrange);
        end

        xlabel('Bias on track','Interpreter','none')
        ylabel('\delta(Bias): POST - PRE')
        title('spearman \rho = %d', info.rho.plot_bias_delta_data)
    end
    
    
    %% Build Figure Title:
    if ~isempty(plotting_options.additional_title)
        info.fig_title = sprintf('Temporal Bias B: %s', plotting_options.additional_title);
    else
        info.fig_title = 'Temporal Bias B';
    end
    sgtitle(info.fig_title,'Interpreter','none')
end
