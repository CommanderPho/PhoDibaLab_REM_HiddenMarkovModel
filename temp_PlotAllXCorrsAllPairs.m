%% Plot all at once:
figure(1);
clf;
sgtitle('All Valid Pairs');

% xcorr_input = temp.by_behavioral_period.curr_xcorr_allPairs.raw; % un-normalized
xcorr_input = temp.by_behavioral_period.curr_xcorr_allPairs.globally_normalized;

% plotting_options.new_xcorr_plot.plot_single_row = true;
plotting_options.new_xcorr_plot.plot_single_row = false;

% plotting_options.new_xcorr_plot.included_unit_A_indicies = 1:temp.num_valid_units; % Include all
plotting_options.new_xcorr_plot.included_unit_A_indicies = [1 2];


if ~plotting_options.new_xcorr_plot.plot_single_row
    xcorr_all_plots.num_subplot_rows = temp.num_valid_units;
    xcorr_all_plots.num_subplot_columns = length(plotting_options.new_xcorr_plot.included_unit_A_indicies);
    
else
    xcorr_all_plots.num_subplot_rows = temp.num_valid_units;
    xcorr_all_plots.num_subplot_columns = 1;
end
plotting_options.new_xcorr_plot.num_subplots = xcorr_all_plots.num_subplot_rows * xcorr_all_plots.num_subplot_columns;
xcorr_all_plots.subplots = gobjects(plotting_options.new_xcorr_plot.num_subplots, 1);


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
            xcorr_all_plots.subplots(temp.linear_subplot_accumulator) = subplot(xcorr_all_plots.num_subplot_rows, xcorr_all_plots.num_subplot_columns, temp.linear_subplot_accumulator);
            if active_pair_index > 0 % don't plot the (xcorr of the unit with itself)
                temp.curr_xcorr_forPair = squeeze(xcorr_input(active_pair_index, :)); % [num_of_behavioral_state_periods x num_lag_steps]
                %% Single Plot for All Time:
    %             fnPlotXCorrStem(active_results.all.pairwise_xcorrelations.lag_offsets, temp.curr_xcorr_forPair, 'all');
                stem(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), active_results.all.pairwise_xcorrelations.lag_offsets, temp.curr_xcorr_forPair);
                xline(xcorr_all_plots.subplots(temp.linear_subplot_accumulator), 0);
%                 set(xcorr_all_plots.subplots(temp.linear_accumulator),'xtick',[],'ytick',[])
    %             title(xcorr_all_plots.subplots(temp.linear_accumulator), curr_fig_title,'Interpreter','none')
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

