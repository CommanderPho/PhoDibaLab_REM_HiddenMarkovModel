function [fig, h] = fnPhoPlotCorrelationalResults(active_processing, active_results, plottingOptions)
%FNPHOPLOTCORRELATIONALRESULTS Summary of this function goes here
%   Detailed explanation goes here

% temp.active_idx = 1:5;
% temp.active_idx = num_of_electrodes-5:num_of_electrodes;
temp.active_idx = 2; % used only for xcorr mode

num_of_electrodes = size(active_results.indicies.reverse_lookup_unique_electrode_pairs, 1);

%% Display the Correlational Results:
%%%%%%%%%%%%%%%%%%%%%
fig = figure(3);
clf



% % Cell formatter: "[unitID - timeOffset]: xcorr_value"
% plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('[%d - %d]: %d', row_index, active_results.pairwise_xcorrelations.lag_offsets(column_index), temp.plot_matrix(row_index, column_index)); 

% Cell formatter: "unitID[@t=timeOffset]"
if ~exist('plottingOptions','var')
    plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('%d[@t=%d]', row_index, active_results.all.pairwise_xcorrelations.lag_offsets(column_index)); 
    plottingOptions.xlabel = '[seconds]';
    plottingOptions.ylabel = 'xcorr';
    plottingOptions.title = sprintf('Pairwise XCorr for Unit %d',  temp.active_idx);
end

if ~isfield(plottingOptions, 'custom_cell_text_formatter')
    plottingOptions.custom_cell_text_formatter = @(row_index, column_index) sprintf('%d[@t=%d]', row_index, active_results.all.pairwise_xcorrelations.lag_offsets(column_index)); 
end

if ~isfield(plottingOptions, 'xlabel')
    plottingOptions.xlabel = '[seconds]';
end

if ~isfield(plottingOptions, 'showOnlyAlwaysStableCells')
    plottingOptions.showOnlyAlwaysStableCells = false;
end

if ~isfield(plottingOptions, 'plotMode')
    plottingOptions.plotMode = 'xcorr';
end


if plottingOptions.showOnlyAlwaysStableCells
    isAlwaysStable = (active_processing.spikes.stability_count == 3);
    numAlwaysStableCells = sum(isAlwaysStable, 'all');
    alwaysStableIndicies = find(isAlwaysStable);    
end



if strcmpi(plottingOptions.plotMode,'autocorr')
    
    if ~isfield(plottingOptions, 'title')
        if plottingOptions.showOnlyAlwaysStableCells
            plottingOptions.title = sprintf('autocorrelation for Unit %d (for Always Stable Units (%d of %d total)', temp.active_idx, numAlwaysStableCells, length(active_processing.spikes.time));
        else
            plottingOptions.title = sprintf('autocorrelation for Unit %d',  temp.active_idx);
        end
    end
    
    if ~isfield(plottingOptions, 'ylabel')
        plottingOptions.ylabel = 'autocorrelation'; 
    end
    
    % active_results.all.autocorrelations % 70701x1 double
    temp.plot_matrix = cell2mat(active_results.all.autocorrelations); % 70701x126 double
    % Remove unstable items:
    if plottingOptions.showOnlyAlwaysStableCells
        temp.plot_matrix = temp.plot_matrix(:, isAlwaysStable); % Only get the cells that are always stable
    end
    
    if isfield(plottingOptions, 'timestamps')
        h = stackedplot(plottingOptions.timestamps, temp.plot_matrix);
    else
        h = stackedplot(temp.plot_matrix);
    end
    
elseif strcmpi(plottingOptions.plotMode,'xcorr')
    
    if plottingOptions.showOnlyAlwaysStableCells
        % If the selected index isn't stable, get the stable index
        if ~isAlwaysStable(temp.active_idx)
           temp.original_index = temp.active_idx;
           warning(['cell with index ' num2str(temp.original_index) ' is not always stable. Selecting stable index ' num2str(temp.active_idx) ' instead.']);
           temp.active_idx = alwaysStableIndicies(temp.active_idx);
           fprintf('Stable index %d corresponding to original unit index %d.\n', temp.active_idx, temp.original_index);
        end
    end

    if ~isfield(plottingOptions, 'title')
        if plottingOptions.showOnlyAlwaysStableCells
            plottingOptions.title = sprintf('Pairwise XCorr for Unit %d (for Always Stable Units (%d of %d total)',  temp.active_idx, numAlwaysStableCells, length(active_processing.spikes.time));
        else
            plottingOptions.title = sprintf('Pairwise XCorr for Unit %d',  temp.active_idx);
        end
    end
    
    if ~isfield(plottingOptions, 'ylabel')
        plottingOptions.ylabel = 'xcorr'; 
    end

    temp.found_lin_idx = active_results.indicies.reverse_lookup_unique_electrode_pairs(temp.active_idx, :); % 1x126 double
    temp.preview_subset = 1:num_of_electrodes;
    temp.preview_subset_size = length(temp.preview_subset);
    temp.active_found_lin_idx = temp.found_lin_idx(temp.preview_subset);


    temp.excluded_indicies = find(temp.active_found_lin_idx < 1);
    temp.active_found_lin_idx(temp.excluded_indicies) = 1; % temporarily fill it in with a valid index (referring to row 1), but then replace it afterwards with NaN

    % active_results.all.pairwise_xcorrelations.xcorr: 7875x19 double
    temp.plot_matrix = active_results.all.pairwise_xcorrelations.xcorr(temp.active_found_lin_idx, :);
    temp.plot_matrix(temp.excluded_indicies, :) = NaN;
    
    % Remove unstable items:
    if plottingOptions.showOnlyAlwaysStableCells
        temp.plot_matrix = temp.plot_matrix(isAlwaysStable, :); % Only get the cells that are always stable
    end
    
    [h, temp.plot_info] = fnPhoMatrixPlotDetailed(temp.plot_matrix, plottingOptions);
    ylabel(plottingOptions.ylabel);
    
else    
    error(['invalid plot mode: ' plottingOptions.plotMode]);
end


xlabel(plottingOptions.xlabel);
title(plottingOptions.title);


end

