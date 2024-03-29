% Requires the datastructures from "PhoDibaProcess_Stage2.m" to be loaded
% Stage 3 of the processing pipeline

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

clear temp

if ~exist('data_config','var')
    Config;
end

if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    if exist('across_experiment_results','var')
        % Call the LoadFromMultiExperimentResults to load the single-experiment-style variables (active_processing, results_array, ...)
        fprintf('loading multi-experiment results data into single-experiment-style variables...\n'); 
        LoadFromMultiExperimentResults;
    else
        % No multi-experiment results loaded either, just load from the intermediate results file:
        fprintf('loading data from %s...\n', data_config.output.intermediate_file_paths{2}); 
        load(data_config.output.intermediate_file_paths{2}, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');  
    end
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end

if ~exist('results_array','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading results from %s...\n', data_config.output.results_file_path);
    load(data_config.output.results_file_path, 'results_array', 'general_results');
    fprintf('done. Contains results for %d different bin sizes.\n', length(results_array));
else
    fprintf('results_array already exists in workspace. Contains results for %d different bin sizes. Using extant data.\n', length(results_array));
end

%% Begin:
fprintf('PhoDibaProcess_ShowResults ready to process!\n');

plotting_options.showOnlyAlwaysStableCells = processing_config.showOnlyAlwaysStableCells;
plotting_options.plotMode = 'autocorr';

current_binning_index = 2;
active_binning_resolution = processing_config.step_sizes{current_binning_index};
temp.curr_timestamps = timesteps_array{current_binning_index};
temp.curr_processed = active_processing.processed_array{current_binning_index};
active_results = results_array{current_binning_index};


fprintf('Plotting results with bin resolution set to %d.\n', active_binning_resolution);


% Binned Spike rates per time:
[active_binned_spike_data_matrix] = sparse(fnUnitDataCells2mat(active_processing.processed_array{current_binning_index}.all.binned_spike_counts));  % 35351x126 double


%% Display the Correlational Results:
xcorr_fig = figure(14);
plotting_options.plotMode = 'xcorr';
[temp.fig, temp.h] = fnPhoPlotCorrelationalResults(active_processing, general_results, active_results, plotting_options, xcorr_fig);
    
autocorr_fig = figure(15);
% plotting_options.timestamps = temp.curr_timestamps;
plotting_options.plotMode = 'autocorr';
[temp.fig, temp.h] = fnPhoPlotCorrelationalResults(active_processing, general_results, active_results, plotting_options, autocorr_fig);



%% Display Unit Spike Rate Heatmap:
phoPlotSpikeRateHeatmap;


%% Display like interactive spike raster:
phoPlotInteractiveSpikesScrollableRaster;


%% Display the Correlational Results:
% [temp.fig, temp.h] = fnPhoPlotCorrelationalResults(active_results);

    

% phoPlotSpikeRateHeatmap;