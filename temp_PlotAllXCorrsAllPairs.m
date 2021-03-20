% Currently plots a column of XCorr values that were loaded by running Mar17Analysis.mlx (which itself requires loading data from LoadSingleExperimentCombinedResults.m)


%% Set Plotting Options:
plotting_options.should_use_custom_subplots = true;

% Options for tightening up the subplots:
if plotting_options.should_use_custom_subplots
%     plotting_options.subtightplot.gap = [0.01 0.01]; % [intra_graph_vertical_spacing, intra_graph_horizontal_spacing]
    plotting_options.subtightplot.gap = [0.001 0.01]; % [intra_graph_vertical_spacing, intra_graph_horizontal_spacing]
    plotting_options.subtightplot.width_h = [0.01 0.01]; % Looks like [padding_bottom, padding_top]
    plotting_options.subtightplot.width_w = [0.001 0.001];

    plotting_options.opt = {plotting_options.subtightplot.gap, plotting_options.subtightplot.width_h, plotting_options.subtightplot.width_w}; % {gap, width_h, width_w}
    subplot_cmd = @(m,n,p) subtightplot(m, n, p, plotting_options.opt{:});
else
    subplot_cmd = @(m,n,p) subplot(m, n, p);
end

%% Customize the plotting command to use (stem, plot, area, etc):
% active_plot_cmd = @(ax,x,y) stem(ax, x, y);
% active_plot_cmd = @(ax,x,y) plot(ax, x, y);
active_plot_cmd = @(ax,x,y) reduce_plot(x, y);


%% Plot all at once:
figure(1);
clf;
sgtitle('All Valid Pairs');

% xcorr_input = temp.by_behavioral_period.curr_xcorr_allPairs.raw; % un-normalized
xcorr_input = temp.by_behavioral_period.curr_xcorr_allPairs.globally_normalized;

% plotting_options.new_xcorr_plot.plot_single_row = true;
plotting_options.new_xcorr_plot.plot_single_row = false;

% plotting_options.new_xcorr_plot.included_unit_A_indicies = 1:temp.num_valid_units; % Include all
plotting_options.new_xcorr_plot.included_unit_A_indicies = [4];
% plotting_options.new_xcorr_plot.included_unit_A_indicies = [4 8 9];


if ~plotting_options.new_xcorr_plot.plot_single_row
    xcorr_all_plots.num_subplot_rows = temp.num_valid_units;
    xcorr_all_plots.num_subplot_columns = length(plotting_options.new_xcorr_plot.included_unit_A_indicies);
    
else
    xcorr_all_plots.num_subplot_rows = temp.num_valid_units;
    xcorr_all_plots.num_subplot_columns = 1;
end
plotting_options.new_xcorr_plot.num_subplots = xcorr_all_plots.num_subplot_rows * xcorr_all_plots.num_subplot_columns;
xcorr_all_plots.subplots = gobjects(plotting_options.new_xcorr_plot.num_subplots, 1);

%% Switch rows and columns:
% temp.old_num_subplot_rows = xcorr_all_plots.num_subplot_rows;
% xcorr_all_plots.num_subplot_rows = xcorr_all_plots.num_subplot_columns;
% xcorr_all_plots.num_subplot_columns = temp.old_num_subplot_rows;
    

temp.linear_accumulator = 1;
temp.linear_subplot_accumulator = 1; % combine with temp.linear_accumulator?

for active_unit_A_index = 1:temp.num_valid_units
    
    for active_unit_B_index = 1:temp.num_valid_units
        if plotting_options.new_xcorr_plot.plot_single_row & (active_unit_A_index > 1)
            % do nothing, just continue
        elseif ~ismember(active_unit_A_index, plotting_options.new_xcorr_plot.included_unit_A_indicies)
            % do nothing also
        else
            active_pair_index = across_experiment_results{active_expt_index}.general_results.indicies.reverse_lookup_unique_electrode_pairs(active_unit_A_index, active_unit_B_index);
            % Convert subplot index to incement along a column first (each row for the column) and then move to the next column instead of its default (row -> column) mode
            
            [curr_row, curr_col] = ind2subplot(xcorr_all_plots.num_subplot_rows, xcorr_all_plots.num_subplot_columns, temp.linear_subplot_accumulator);
            corrected_active_subplot_index = temp.linear_subplot_accumulator;
            % 3->2, 5->3
            % 2-> (xcorr_all_plots.num_subplot_rows + 1), 4-> (xcorr_all_plots.num_subplot_rows + 2)
%             corrected_active_subplot_index = temp.linear_subplot_accumulator - curr_row+1;
%             corrected_active_subplot_index = ((curr_col - 1) * xcorr_all_plots.num_subplot_rows)+ curr_row;
%               corrected_active_subplot_index = ((curr_row - 1) * xcorr_all_plots.num_subplot_rows)+ curr_col; % Way wrong
%               corrected_active_subplot_index = ((curr_col - 1) * xcorr_all_plots.num_subplot_columns)+ curr_row; % Way wrong
              corrected_active_subplot_index = ((curr_row - 1) * xcorr_all_plots.num_subplot_columns)+ curr_col;

            xcorr_all_plots.subplots(temp.linear_subplot_accumulator) = subplot_cmd(xcorr_all_plots.num_subplot_rows, xcorr_all_plots.num_subplot_columns, corrected_active_subplot_index);
            if active_pair_index > 0 % don't plot the (xcorr of the unit with itself)
                temp.curr_xcorr_forPair = squeeze(xcorr_input(active_pair_index, :)); % [num_of_behavioral_state_periods x num_lag_steps]
                %% Single Plot for All Time:
    %             fnPlotXCorrStem(active_results.all.pairwise_xcorrelations.lag_offsets, temp.curr_xcorr_forPair, 'all');
%                 axes(xcorr_all_plots.subplots(temp.linear_subplot_accumulator);
                active_plot_cmd(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), active_results.all.pairwise_xcorrelations.lag_offsets, temp.curr_xcorr_forPair);
                xline(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), 0, '-r');
                set(xcorr_all_plots.subplots(temp.linear_subplot_accumulator),'xtick',[],'ytick',[])
                title(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), sprintf('Xcorr(u%d, u%d)', active_unit_A_index, active_unit_B_index),'Interpreter','none') % (first is same for entire column)
%                 title(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), sprintf('(linear: %d, corrected_linear: %d, row: %d, col: %d)', temp.linear_subplot_accumulator, corrected_active_subplot_index, curr_row, curr_col),'Interpreter','none') % second is the same for entire column
            end % end if active_pair_index > 0
            temp.linear_subplot_accumulator = temp.linear_subplot_accumulator + 1;
        end
        temp.linear_accumulator = temp.linear_accumulator + 1;
    end
    drawnow
end


% for linear_accumulator = 1:length(xcorr_all_plots.subplots)
% %    set(xcorr_all_plots.subplots(linear_accumulator),'xtick',[],'ytick',[])
% %    ylim(xcorr_all_plots.subplots(linear_accumulator),[0,1]);
%     axes(xcorr_all_plots.subplots(linear_accumulator));
%     ylim([0,0.2]);
% end

% 
% function fnPlotXCorrStem(lag_offsets, plotVals, curr_fig_title)
%     stem(lag_offsets, plotVals);
%     xline(0);
%     title(curr_fig_title,'Interpreter','none')
% end

