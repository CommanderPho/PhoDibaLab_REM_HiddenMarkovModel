% Requires the datastructures from "PhoDibaPrepare_Stage0.m" to be loaded
% Stage 1 of the processing pipeline

addpath(genpath('helpers'));
addpath(genpath('libraries/buzcode/'));

clear temp

if ~exist('data_config','var')
    Config;
end

if ~exist('active_processing','var') %TEMP: cache the loaded data to rapidly prototype the script
    fprintf('loading from %s...\n', data_config.output.intermediate_file_paths{2});
    load(data_config.output.intermediate_file_paths{2}, 'active_processing', 'processing_config', 'num_of_electrodes', 'source_data', 'timesteps_array');
    fprintf('done.\n');
else
    fprintf('active_processing already exists in workspace. Using extant data.\n');
end


%% Pairwise Indexing:
% Generates all unique pairs of indicies for pairwise comparisons (without replacement or repetition)
general_results.indicies.unique_electrode_pairs = nchoose2([1:num_of_electrodes]);
general_results.indicies.num_unique_pairs = length(general_results.indicies.unique_electrode_pairs);

% Build a reverse lookup matrix
general_results.indicies.reverse_lookup_unique_electrode_pairs = zeros(num_of_electrodes);
for linear_pair_idx = 1:general_results.indicies.num_unique_pairs
	curr_pair = general_results.indicies.unique_electrode_pairs(linear_pair_idx,:);
	general_results.indicies.reverse_lookup_unique_electrode_pairs(curr_pair(1), curr_pair(2)) = linear_pair_idx;
	general_results.indicies.reverse_lookup_unique_electrode_pairs(curr_pair(2), curr_pair(1)) = linear_pair_idx;
end


%% Get the duration of the states 
general_results.GroupedByState.groups = findgroups(active_processing.behavioral_periods_table.type);
general_results.GroupedByState.durations = splitapply(@sum, active_processing.behavioral_periods_table.duration, general_results.GroupedByState.groups);


%% Build the Correlational Results:
processing_config.max_xcorr_lag = 9; % Specified the maximum pairwise cross-correlation lag in seconds, the output ranges from -maxlag to maxlag

results_array = cell([length(processing_config.step_sizes) 1]);


%% Loop through each listed binning/step size and built an "results_array"
for current_binning_index = 1:length(processing_config.step_sizes)
	temp.curr_timestamps = timesteps_array{current_binning_index};
	temp.curr_processed = active_processing.processed_array{current_binning_index};
    temp.curr_step_size = processing_config.step_sizes{current_binning_index};
    
	[~, active_results.all.autocorrelations, active_results.all.partial_autocorrelations, active_results.all.pairwise_xcorrelations] = fnProcessCorrelationalMeasures(temp.curr_processed.all.binned_spike_counts, general_results.indicies, processing_config, temp.curr_step_size);

	%% Aggregate Measures:

	%% Split based on experiment epoch:
	for i = 1:length(data_config.behavioral_epoch_names)
		temp.curr_epoch_name = data_config.behavioral_epoch_names{i};
		% size(temp.curr_processed.by_epoch.(temp.curr_epoch_name).binned_spike_counts): 1 x 126
		[~, active_results.by_epoch.(temp.curr_epoch_name).autocorrelations, active_results.by_epoch.(temp.curr_epoch_name).partial_autocorrelations, active_results.by_epoch.(temp.curr_epoch_name).pairwise_xcorrelations] = fnProcessCorrelationalMeasures(temp.curr_processed.by_epoch.(temp.curr_epoch_name).binned_spike_counts, ...
			general_results.indicies, processing_config, temp.curr_step_size);

		active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes = cell2mat(temp.curr_processed.by_epoch.(temp.curr_epoch_name).binned_spike_counts);
		active_results.aggregates.by_epoch.(temp.curr_epoch_name).across_all_cells.count = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 2);
		active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts = sum(active_results.aggregates.by_epoch.(temp.curr_epoch_name).spikes, 'all');
		
		fprintf('epoch: %s\n total_counts: %d\n', temp.curr_epoch_name, active_results.aggregates.by_epoch.(temp.curr_epoch_name).total_counts);
	end


	%timesteps

	if processing_config.show_graphics
		figure(1);
	end

	%% Split based on behavioral state:
	for i = 1:length(active_processing.behavioral_state_names)
		temp.curr_state_name =  active_processing.behavioral_state_names{i};
		
		[~, active_results.by_state.(temp.curr_state_name).autocorrelations, active_results.by_state.(temp.curr_state_name).partial_autocorrelations, active_results.by_state.(temp.curr_state_name).pairwise_xcorrelations] = fnProcessCorrelationalMeasures(temp.curr_processed.by_state.(temp.curr_state_name).binned_spike_counts, ...
			general_results.indicies, processing_config, temp.curr_step_size);

		active_results.aggregates.by_state.(temp.curr_state_name).spikes = cell2mat(temp.curr_processed.by_state.(temp.curr_state_name).binned_spike_counts);
		active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 2);
		active_results.aggregates.by_state.(temp.curr_state_name).total_counts = sum(active_results.aggregates.by_state.(temp.curr_state_name).spikes, 'all');
		
		fprintf('state: %s\n total_counts: %d\n', temp.curr_state_name, active_results.aggregates.by_state.(temp.curr_state_name).total_counts);
		
		if processing_config.show_graphics
			subplot(4,1,i);
			plot(active_results.aggregates.by_state.(temp.curr_state_name).across_all_cells.count);
			ylabel(temp.curr_state_name);
			xlabel('');
		end
	end

	if processing_config.show_graphics
		xlim([timesteps(1), timesteps(end)]);
		title('behavioral state spike counts')
	end

	results_array{current_binning_index} = active_results;
	
end % end for processing_config.step_sizes loop

fprintf('writing out results to %s...\n', data_config.output.results_file_path);
save(data_config.output.results_file_path, 'general_results', 'results_array');
fprintf('done.\n');

fprintf('PhoDibaProcess_Stage1 complete!\n');