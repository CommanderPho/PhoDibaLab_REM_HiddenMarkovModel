%% Uses the Buzcode CCG function to compute the auto and cross-correlations for each provided spike train:
% Pho Hale, 2021-03-24

% All CCG:
% [ccg, t] = CCG(active_processing.spikes.time, []);

% Returns a 201x126x126 CCG


ccg_options.binSize = 0.01; % Seconds
ccg_options.duration = 2; % Seconds





%% Want the ccg for REM, NREM, etc to be the average of all ccg's for each behavioral substate. Not the concatenated values or anything.

% Each spike has a behavioral_duration_indicies entry to indicate which behavioral_period it belongs to:
% active_processing.spikes.behavioral_duration_indicies

% mystats = @(x)[min(x) median(x) max(x)];
% mystats = @(x)[min(x) median(x) max(x)];

% i = 1;
% split_CCG = splitapply((@(X) CCG(X, [])), active_processing.spikes.time{i}, active_processing.spikes.behavioral_duration_indicies{i});


% test = active_processing.spikes.time{:}(active_processing.spikes.behavioral_duration_indicies{:});
% 
% test = cellfun(@(X,i) X(i), active_processing.spikes.time, active_processing.spikes.behavioral_duration_indicies{,'un',0)

% G = findgroups(active_processing.spikes.behavioral_duration_indicies{1});


temp.curr_num_of_units = size(active_processing.spikes.time, 1);
temp.curr_num_of_behavioral_states = height(active_processing.behavioral_periods_table);


% num2cell(1:temp.curr_num_of_behavioral_states)
% test = cellfun(@(I,s) (I == s), active_processing.spikes.behavioral_duration_indicies, {1:temp.curr_num_of_behavioral_states},'UniformOutput',false);

% Allocate a new array to hold the output: which is an array of spike times for each unit for each behavioral state (columns)
behavioral_epoch_spikes = cell(temp.curr_num_of_units, temp.curr_num_of_behavioral_states);


%% Compute CCG Storage:
ccg_results.num_timesteps = (ccg_options.duration ./ ccg_options.binSize) + 1;
ccg_results.num_groups = temp.curr_num_of_units;

if ~isfield(ccg_results,'by_behavioral_period')
    fprintf('Pre-allocating ccg_results output.... this may take several minutes....\n');
    ccg_results.by_behavioral_period.ccg.raw = zeros([temp.curr_num_of_behavioral_states, ccg_results.num_timesteps, ccg_results.num_groups, ccg_results.num_groups]);
else
    warning('ccg_results.by_behavioral_period.ccg.raw already exists. Resuming');
end
% active_results.by_behavioral_period.

% Loop over behavioral activities
for state_index = 1:temp.curr_num_of_behavioral_states
    
%     cellfun(@(X,I) X(i), active_processing.spikes.time, active_processing.spikes.behavioral_duration_indicies(state_index) ==  ,'un',0)
    
    % Gets the matching indicies for this state:
%     test = cellfun(@(I) (I == state_index), active_processing.spikes.behavioral_duration_indicies,'UniformOutput',false);

    % Returns the times only for the spikes that occur within this region:
    behavioral_epoch_spikes(:, state_index) = cellfun(@(X,I) X(I == state_index), active_processing.spikes.time, active_processing.spikes.behavioral_duration_indicies,'UniformOutput',false);
    
    % Do CCG Here too:
    [ccg_results.by_behavioral_period.ccg.raw(state_index, :, :, :), ~] = CCG(active_processing.spikes.time, [], 'binSize', ccg_options.binSize, 'duration', ccg_options.duration);
    
    % [t x ngroups x ngroups] matrix where ccg(t,i,j) is the
%           number (or rate) of events of group j at time lag t with  
%           respect to reference events from group i


    fprintf('behavioral state progress: %d/%d\n', state_index, temp.curr_num_of_behavioral_states);
        
%     temp.curr_state_start = active_processing.behavioral_periods_table.epoch_start_seconds(state_index);
%     temp.curr_state_end = active_processing.behavioral_periods_table.epoch_end_seconds(state_index);
%     temp.curr_state_type = active_processing.behavioral_periods_table.type(state_index);
%     temp.curr_epoch_type = active_processing.behavioral_periods_table.behavioral_epoch(state_index);
    
%     temp.curr_state_spikes = cell(num_of_electrodes, 1);
%     % Extract the spike train for each electrode
%     for electrode_index = 1:num_of_electrodes
%         % Convert spike times to relative to expt start and scale to seconds.
%         temp.curr_electrode_spikes = active_processing.spikes.time{electrode_index};
%         % Get the spike times that belong to this particular state.
%         temp.curr_state_spikes_idx{electrode_index} = find((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
%         
%         temp.curr_state_spikes{electrode_index} = temp.curr_electrode_spikes((temp.curr_state_start < temp.curr_electrode_spikes) & (temp.curr_electrode_spikes < temp.curr_state_end));
%         
%         temp.spikes_behavioral_states{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = temp.curr_state_type;
%         temp.spikes_behavioral_epoch{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = temp.curr_epoch_type;
%         
%         active_processing.spikes.behavioral_duration_indicies{electrode_index}(temp.curr_state_spikes_idx{electrode_index}) = state_index;
%     end

end