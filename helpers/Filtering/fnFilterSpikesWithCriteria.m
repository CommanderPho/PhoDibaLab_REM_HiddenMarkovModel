function [is_target_entry_included] = fnFilterSpikesWithCriteria(filter_target, included_epochs, included_states, target_options)
%FNFILTERSPIKESWITHCRITERIA Summary of this function goes here
%   Detailed explanation goes here
% included_ripple_index

% included_cell_types
% active_processing.spikes.time(plot_outputs.filter_active_units)
% filtered_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), spikesTable.time, spikesTable.isRippleSpike, 'UniformOutput', false); %% Filtered to only show the ripple spikes

    if ~exist('filter_target','var') || ~istabular(filter_target)
        error('fnFilterSpikesWithCriteria: The filter_target must currently be a table! Aborting.')
    end

    if ~exist('target_options','var')
        target_options = struct();
    end
    
    % For compatibility with the old naming conventions, you might need to use:
%     target_options.behavioral_states_variable_name = 'type';

    target_options = fnAddDefaultOptionalArgs({'behavioral_states_variable_name', 'behavioral_epochs_variable_name'}, ...
        {'type', 'behavioral_epoch'}, ...
        target_options);

    

    num_of_target_entries = height(filter_target);
%     is_period_included = logical(ones([num_of_behavioral_state_periods 1]));
    
    %% Check for exclusion by epoch:
    if exist('included_epochs','var') & ~isempty(included_epochs)
%         is_epoch_included = false([num_of_behavioral_state_periods 1]);
%         curr_conditions.is_behavioral_epoch_type = (@(compare_type) (filter_target.behavioral_epoch == compare_type));
%         temp.curr_is_included = cellfun(curr_conditions.is_behavioral_epoch_type, included_epochs, 'UniformOutput', false);
%         % Iterate through all the conditions:
%         for i = 1:length(temp.curr_is_included)
% %             is_period_included = is_period_included & temp.curr_is_included{i};
%             is_epoch_included = is_epoch_included | temp.curr_is_included{i}; % Should be or (|) as it must meet one of the requirements for epoch to be included.
%         end

       [is_epoch_included] = BehavioralEpoch.filter(filter_target.(target_options.behavioral_epochs_variable_name), included_epochs);

    else
        % When empty, return true for all of this type
        is_epoch_included = true([num_of_target_entries 1]);
    end


    %% Check for exclusion by sleep state:
    
    if exist('included_states','var') & ~isempty(included_states)
%         is_state_included = false([num_of_behavioral_state_periods 1]);
%         curr_conditions.is_behavioral_state_type = (@(compare_type) (filter_target.type == compare_type));
%         temp.curr_is_included = cellfun(curr_conditions.is_behavioral_state_type, included_states, 'UniformOutput', false);
%         for i = 1:length(temp.curr_is_included)
% %             is_period_included = is_period_included & temp.curr_is_included{i};
%             is_state_included = is_state_included | temp.curr_is_included{i};
%         end

        [is_state_included] = BehavioralState.filter(filter_target.(target_options.behavioral_states_variable_name), included_states);
    else
        % When empty, return true for all of this type
        is_state_included = true([num_of_target_entries 1]);
    end

    % Finally, for the period to be included in general it must meet both the epoch and state requirements, so it should be logical AND
    is_target_entry_included = is_epoch_included & is_state_included;

end

