function [is_period_included] = fnFilterPeriodsWithCriteria(active_processing, included_epochs, included_states)
%fnFilterPeriodsWithCriteria Period will be included if it belongs to one of the specified epoch types AND one of the specified included states
%   Detailed explanation goes here

% included_epochs: a cell array of epoch types to include. Period will be included if it belongs to one of the specified epoch types AND one of the specified included states
    % e.g. {'pre_sleep', 'post_sleep'}
% included_states: a cell array of behavioral state types to include.
    % e.g. {'rem'}
    
    num_of_behavioral_state_periods = height(active_processing.behavioral_periods_table);
    is_period_included = logical(ones([num_of_behavioral_state_periods 1]));

    if exist('included_epochs','var') & ~isempty(included_epochs)
        curr_conditions.is_behavioral_epoch_type = (@(compare_type) (active_processing.behavioral_periods_table.behavioral_epoch == compare_type));
        temp.curr_is_included = cellfun(curr_conditions.is_behavioral_epoch_type, included_epochs, 'UniformOutput', false);
        for i = 1:length(temp.curr_is_included)
            is_period_included = is_period_included & temp.curr_is_included{i};
        end

    end

    if exist('included_states','var') & ~isempty(included_states)
        curr_conditions.is_behavioral_state_type = (@(compare_type) (active_processing.behavioral_periods_table.type == compare_type));
        temp.curr_is_included = cellfun(curr_conditions.is_behavioral_state_type, included_states, 'UniformOutput', false);
        for i = 1:length(temp.curr_is_included)
            is_period_included = is_period_included & temp.curr_is_included{i};
        end

    end



end

