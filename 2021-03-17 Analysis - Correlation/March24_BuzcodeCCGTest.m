%% Uses the Buzcode CCG function to compute the auto and cross-correlations for each provided spike train:
% Pho Hale, 2021-03-24

% All CCG:
% [ccg, t] = CCG(active_processing.spikes.time, []);

% Returns a 201x126x126 CCG
addpath(genpath('helpers'));

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
if ~exist('behavioral_epoch_spikes','var')
    behavioral_epoch_spikes = cell(temp.curr_num_of_units, temp.curr_num_of_behavioral_states); % 126x668
end


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

%% TEMP:
% Reshape the 126x668 cell matrix into a (126 * 668)x1 column vector
behavioral_epoch_spikes_flat = reshape(behavioral_epoch_spikes,[],1);
% size(behavioral_epoch_spikes_flat): 84168x1 cell
% Convert Back:
% behavioral_epoch_spikes = reshape(behavioral_epoch_spikes,temp.curr_num_of_units, temp.curr_num_of_behavioral_states);


% [ccg_results.by_behavioral_period.ccg.raw(state_index, :, :, :), ~] = CCG(behavioral_epoch_spikes_flat, [], 'binSize', ccg_options.binSize, 'duration', ccg_options.duration);

% [ccg_results.by_behavioral_period.ccg.raw_flat(:, :, :, :), ~] = CCG(behavioral_epoch_spikes_flat, [], 'binSize', ccg_options.binSize, 'duration', ccg_options.duration);
% [201, 3780, 3780]

% Split back into the appropriate dimensions:
unit_has_no_spikes = false(temp.curr_num_of_units, temp.curr_num_of_behavioral_states);

%% /TEMP

% Loop over behavioral activities
for state_index = 1:temp.curr_num_of_behavioral_states
    
    % Returns the times only for the spikes that occur within this region:
    behavioral_epoch_spikes(:, state_index) = cellfun(@(X,I) X(I == state_index), active_processing.spikes.time, active_processing.spikes.behavioral_duration_indicies,'UniformOutput',false);
    
    % Find any units that have no spikes for this behavioral_epoch. This epoch will be excluded from analysis for those epochs then.
    unit_has_no_spikes(:, state_index) = cellfun(@isempty, behavioral_epoch_spikes(:, state_index)); 
    
    % Do CCG Here too:
    [ccg_results.by_behavioral_period.ccg.raw(state_index, :, :, :), ~, cell_conversion] = fnPhoCCG(behavioral_epoch_spikes(:, state_index), [], 'binSize', ccg_options.binSize, 'duration', ccg_options.duration);
    
    % [t x ngroups x ngroups] matrix where ccg(t,i,j) is the
%           number (or rate) of events of group j at time lag t with  
%           respect to reference events from group i


    fprintf('behavioral state progress: %d/%d\n', state_index, temp.curr_num_of_behavioral_states);
        
end


%% 
% ccg_results.by_behavioral_period.ccg.raw: [668   201   126   126]
%   ccg_results.by_behavioral_period.ccg.raw
temp.isValidCCGResults = ~cellfun(@isempty, behavioral_epoch_spikes); %126x668
% isempty(behavioral_epoch_spikes);


ccg_results.by_behavioral_period.ccg.repaired = ccg_results.by_behavioral_period.ccg.raw;


for i = 1:size(temp.isValidCCGResults, 1)
    % Remove the bad entries:
    ccg_results.by_behavioral_period.ccg.repaired(~temp.isValidCCGResults(i,:),:,i,i) = nan;
end


ccg_results.all.ccg.aggregates.mean = squeeze(mean(ccg_results.by_behavioral_period.ccg.repaired, 1, 'omitnan'));

[ccg_results.lag_offsets, ccg_results.nBins] = fnComputeCCGTimes(ccg_options.binSize, ccg_options.duration);




