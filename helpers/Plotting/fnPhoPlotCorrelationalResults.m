function [fig, h] = fnPhoPlotCorrelationalResults(active_results, plottingOptions)
%FNPHOPLOTCORRELATIONALRESULTS Summary of this function goes here
%   Detailed explanation goes here

num_of_electrodes = size(active_results.indicies.reverse_lookup_unique_electrode_pairs, 1);

%% Display the Correlational Results:
%%%%%%%%%%%%%%%%%%%%%
fig = figure(3);
clf
% temp.active_idx = 1:5;
% temp.active_idx = num_of_electrodes-5:num_of_electrodes;
temp.active_idx = 10;

% temp.found_lin_idx = find((active_results.indicies.unique_electrode_pairs(:,1) == temp.active_idx) | (active_results.indicies.unique_electrode_pairs(:,2) == temp.active_idx));
temp.found_lin_idx = active_results.indicies.reverse_lookup_unique_electrode_pairs(temp.active_idx, :); % 1x126 double

% temp.preview_subset = temp.active_idx+1:temp.active_idx+80;
temp.preview_subset = 1:num_of_electrodes;
temp.preview_subset_size = length(temp.preview_subset);
temp.active_found_lin_idx = temp.found_lin_idx(temp.preview_subset);

% Overlay version:
% stem(active_results.pairwise_xcorrelations.lags, active_results.pairwise_xcorrelations.processed.sorted.xcorr(temp.active_found_lin_idx, :)','filled');

% % Subplot version:
% % tiledlayout(temp.preview_subset_size, 1)
% for i = 1:temp.preview_subset_size
%     % Subplot
%     temp.isLastSession = (i == temp.preview_subset_size);
%     ax(i) = subplot_cmd(temp.preview_subset_size, 1, i);
%     stem(ax(i), active_results.pairwise_xcorrelations.lags, active_results.pairwise_xcorrelations.processed.sorted.xcorr(temp.active_found_lin_idx(i), :)')
%     
%     set(gca,'ytick',[],'yticklabel',[])
%     
%     if temp.isLastSession
%        xlabel('t');
%     else
% %         title('');
% %         xlabel('');
% %         ylabel('');
%         set(gca,'xtick',[],'xticklabel',[])
% %         set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
%     end
%     
% end


%% Image Matrix Version:
% fnPhoMatrixPlot(active_results.pairwise_xcorrelations.xcorr');

temp.excluded_indicies = find(temp.active_found_lin_idx < 1);
temp.active_found_lin_idx(temp.excluded_indicies) = 1; % temporarily fill it in with a valid index (referring to row 1), but then replace it afterwards with NaN
temp.plot_matrix = active_results.pairwise_xcorrelations.xcorr(temp.active_found_lin_idx, :);
temp.plot_matrix(temp.excluded_indicies, :) = NaN;

% % Cell formatter: "[unitID - timeOffset]: xcorr_value"
% plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('[%d - %d]: %d', row_index, active_results.pairwise_xcorrelations.lag_offsets(column_index), temp.plot_matrix(row_index, column_index)); 

% Cell formatter: "unitID[@t=timeOffset]"
plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('%d[@t=%d]', row_index, active_results.pairwise_xcorrelations.lag_offsets(column_index)); 
plottingOptions.xlabel = '[seconds]';
plottingOptions.ylabel = 'xcorr';
plottingOptions.title = sprintf('Pairwise XCorr for Unit %d',  temp.active_idx);

[h, temp.plot_info] = fnPhoMatrixPlotDetailed(temp.plot_matrix, plottingOptions);
xlabel(plottingOptions.xlabel);
ylabel(plottingOptions.ylabel);
title(plottingOptions.title);


end

