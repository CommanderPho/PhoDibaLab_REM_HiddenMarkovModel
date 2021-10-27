%% Working on 2021-10-26

% Load Files:
% load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1.mat');
% load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/toAddVariables.mat');


% 
% % Get only the ripple spikes
% active_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), active_processing.spikes.time(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
% active_spike_ripple_indices = cellfun(@(spike_ripple_indicies, is_spike_ripple) spike_ripple_indicies(is_spike_ripple), active_processing.spikes.RippleIndex(plot_outputs.filter_active_units), active_processing.spikes.isRippleSpike(plot_outputs.filter_active_units), 'UniformOutput', false); %% Filtered to only show the ripple spikes
% 
% [active_spike_ripple_unique_count, active_spike_ripple_unique_indicies] = cellfun(@(spike_ripple_indicies) hist(spike_ripple_indicies, unique(spike_ripple_indicies)), active_spike_ripple_indices, 'UniformOutput', false); %% Filtered to only show the ripple spikes
% [active_spike_flatSpikeTimes, active_spike_flatSpikeUnitIDs] = fnUnitSpikeCells2FlatSpikes(active_spike_times);


%% Filtering Options:
% filter_config.filter_included_cell_types = {};
filter_config.filter_included_cell_types = {'pyramidal'};
% filter_config.filter_included_cell_types = {'interneurons'};
filter_config.filter_maximum_included_contamination_level = {2};

%% Get filter info for active units
[plot_outputs.filter_active_units, plot_outputs.original_unit_index] = fnFilterUnitsWithCriteria(active_processing, true, filter_config.filter_included_cell_types, ...
    filter_config.filter_maximum_included_contamination_level);

temp.num_active_units = sum(plot_outputs.filter_active_units, 'all');
fprintf('Filter: Including %d of %d total units\n', temp.num_active_units, length(plot_outputs.filter_active_units));


%% Flatten subset of spikes table for efficient ripple-related processing:
curr_cell_table = table(active_processing.spikes.time, ...
        active_processing.spikes.isRippleSpike, ...
        active_processing.spikes.RippleIndex, ...
        'VariableNames', {'time','isRippleSpike','RippleIndex'});
% Add the optional variables
curr_cell_table.behavioral_duration_indicies = fnCellContentsTranpose(active_processing.spikes.behavioral_duration_indicies);
curr_cell_table.behavioral_states = fnCellContentsTranpose(active_processing.spikes.behavioral_states);
curr_cell_table.behavioral_epoch = fnCellContentsTranpose(active_processing.spikes.behavioral_epoch);

% Filter the inactive units to save on processing overhead before flattening the cells
curr_cell_table = curr_cell_table(plot_outputs.filter_active_units, :);

% Flatten the cells:
[curr_flattened_table] = fnFlattenCellsToContents(curr_cell_table);

% curr_flattened_table(curr_flattened_table.isRippleSpike)

%% Removes rows with missing values, meaning the rows with RippleIndex == NaN (meaning they aren't ripple spikes)
curr_filtered_table = rmmissing(curr_flattened_table);

%% Exclude sleep states:
% curr_filtered_table = curr_filtered_table(('rem' == curr_filtered_table.behavioral_states), :); % Only those occuring during rem
% curr_filtered_table = curr_filtered_table(('nrem' == curr_filtered_table.behavioral_states), :);
% curr_filtered_table = groupfilter(curr_filtered_table, 'groupID', @(x) all(x == 'rem'), 'behavioral_states');
% 
% curr_filtered_table = curr_filtered_table((('rem' == curr_filtered_table.behavioral_states) | ('nrem' == curr_filtered_table.behavioral_states)), :);

num_active_units = height(curr_cell_table);
target_options.behavioral_states_variable_name = 'behavioral_states';
% [is_target_entry_included] = fnFilterSpikesWithCriteria(curr_filtered_table, [], {'rem','nrem'}, target_options);

% [is_target_entry_included] = fnFilterSpikesWithCriteria(curr_filtered_table, {'pre_sleep'}, {'active'}, target_options);
[track_active_included] = fnFilterSpikesWithCriteria(curr_filtered_table, {'track'}, {'active'}, target_options);
[track_active_outputs.filtered_table, track_active_outputs.eachRipple_filtered_flattened_table, track_active_outputs.eachRipple_Matricies] = fnSplitIntoSpecificRipples(curr_filtered_table, track_active_included, num_active_units);



[track_quiet_included] = fnFilterSpikesWithCriteria(curr_filtered_table, {'track'}, {'quiet'}, target_options);
[track_quiet_outputs.filtered_table, track_quiet_outputs.eachRipple_filtered_flattened_table, track_quiet_outputs.eachRipple_Matricies] = fnSplitIntoSpecificRipples(curr_filtered_table, track_quiet_included, num_active_units);


%% fnReconstructCellsFromFlattenedContents
% Test reconstruction of original test table from the flattened table:
% 
% % We surprisingly need to flatten the output, not reconstruct it, because it's been bound up into a cell array with one cell for each ripple containing a table
% 
% % Flatten the cells:
% 
% % Flatten the changes all back so we can get the extra columns (rippleRelativeSequenceIndex and rippleRelativeTimeOffsets) in the flat filtered table
% track_quiet_outputs.filtered_table = vertcat(track_quiet_outputs.eachRipple_filtered_flattened_table{:});
% 
% 
% [track_quiet_outputs.reconstructedCellTable] = fnFlattenCellsToContents(track_quiet_outputs.eachRipple_filtered_flattened_table);

[track_quiet_outputs.reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(track_quiet_outputs.filtered_table);
% [track_quiet_outputs.reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(track_quiet_outputs.filtered_table);

%%% This is going to be very inefficient
% 
% for unit_idx = 1:height(track_quiet_outputs.reconstructedCellTable)
%     
%     % Get the ripple indexes this belongs in:
%     currUnitRippleIndicies = track_quiet_outputs.reconstructedCellTable.rippleRelativeSequenceIndex{unit_idx};
%     
%     for unit_ripple_idx = 1:length(currUnitRippleIndicies)
%         % All the unitIDs within this ripple
%         currUnitRippleUnitIDs = track_quiet_outputs.reconstructedCellTable.eachRipple_filtered_flattened_table{unit_ripple_idx}.flattened_UnitIDs:
%         currUnitRippleUnitRelativeTimeOffsets = track_quiet_outputs.reconstructedCellTable.eachRipple_filtered_flattened_table{unit_ripple_idx}.rippleRelativeTimeOffsets:
% 
%         %%% NOT YET FINISHED 2021-10-26
%         error('Not yet implemented')
%         %%% NOT YET FINISHED
%         
% 
%         % Find the location of the current unit
%         currUnitFoundIndicies = find(currUnitRippleUnitIDs == unit_idx);
% 
%   
%         %% TODO: The computational complexity here is going to take off, I need to find a smarter implementation.
%         % The concept was to find all the units that preceed our current unit_idx (curr_ripple_leading_units)
%         
%         currUnitRippleUnitRelativeTimeOffsets(curr_ripple_leading_units) % get the offsets for the units that lead our unit
% 
%         master_weights{unit_idx, curr_ripple_leading_units} = -1 .* currUnitRippleUnitRelativeTimeOffsets(curr_ripple_leading_units);
% 
%         %% ?? TODO: Need to set the reciprocal weights too?
%         master_weights{curr_ripple_leading_units, unit_idx} = currUnitRippleUnitRelativeTimeOffsets(curr_ripple_leading_units);
%         
%         currUnitRippleSequenceIndex = track_quiet_outputs.reconstructedCellTable.rippleRelativeSequenceIndex{unit_idx};
%        
%     
% 
%     end
% 
% end



%% Plot the outputs:

% currMatrixName = 'activeSet_Matrix';
currMatrixName = 'relativeSequenceIndex_Matrix';
% currMatrixName = 'relativeProportionalTimeOffset_Matrix';

figure(1);
clf;
fnPhoMatrixPlot(track_quiet_outputs.eachRipple_Matricies.(currMatrixName))
title('quiet wake')
xlabel('SWR Index')
ylabel('Unit ID')

figure(2);
clf;
fnPhoMatrixPlot(track_active_outputs.eachRipple_Matricies.(currMatrixName))
title('active wake')
xlabel('SWR Index')
ylabel('Unit ID')


% %% Clustering exploration
% 
% eachRipple_activeSet_Matrix_active = track_active_outputs.eachRipple_activeSet_Matrix';
% eachRipple_activeSet_Matrix_quiet = track_quiet_outputs.eachRipple_activeSet_Matrix';
% 
% 
% M = track_active_outputs.eachRipple_activeSet_Matrix'; % Transpose MM% To Compare Columns
% D1 = pdist(M,'cityblock'); % Pairwise distance between pairs of observations                                
% D2 = pdist(X,'minkowski');
% 
% Result1 = squareform(D1);
% Result2 = squareform(D2);
% 
% 
% % [coeff,score,pcvar,mu] = pca(y,3);

% track_active_outputs.eachRipple_Matricies.relativeProportionalTimeOffset_Matrix(
% 
% 
% dist = zeros(size(M,2), size(M,2), size(M,1));
% for k1 = 1:size(M,2)
%     for k2 = 1:size(M,2)
%         dist(k1,k2,:) = M(:,k1)-M(:,k2);
%     end
% end
% 
% 
% 
% X = track_quiet_outputs.eachRipple_activeSet_Matrix(plot_outputs.filter_active_units, :);
% [idx,V,D] = spectralcluster(X, 3);
% % figure(3);
% % clf;
% % gscatter(X(:,1),X(:,2),idx);
% 
% % plottingOptions = struct();
% % [h, temp.plot_info] = fnPhoMatrixPlotDetailed(eachRipple_activeSet_Matrix(plot_outputs.filter_active_units, :), plottingOptions);
% 


function [curr_filtered_table, eachRipple_filtered_flattened_table, eachRipple_Matricies] = fnSplitIntoSpecificRipples(target_table, is_target_entry_included, num_active_units)
    
    curr_filtered_table = target_table(is_target_entry_included, :);
    
    %% Filter again to a specific ripple index
    unique_ripple_indices = unique(curr_filtered_table.RippleIndex);
    eachRipple_filtered_flattened_table = cell([length(unique_ripple_indices), 1]);
    for ripple_idx = 1:length(unique_ripple_indices)
        eachRipple_filtered_flattened_table{ripple_idx} = curr_filtered_table(curr_filtered_table.RippleIndex == unique_ripple_indices(ripple_idx), :);
        curr_ripple_num_spikes = height(eachRipple_filtered_flattened_table{ripple_idx});
        eachRipple_filtered_flattened_table{ripple_idx}.rippleRelativeSequenceIndex = [1:curr_ripple_num_spikes]';

        curr_ripple_time_offsets = eachRipple_filtered_flattened_table{ripple_idx}.time;
        curr_ripple_time_offsets = (eachRipple_filtered_flattened_table{ripple_idx}.time - curr_ripple_time_offsets(1)); % convert to relative time by subtracting the start timestamp. This means each spike's time corresponds to the time ellapsed since the start of the ripple
        curr_ripple_normalized_duration_time_offsets = curr_ripple_time_offsets ./ curr_ripple_time_offsets(end); % this should give each relative offset as a value between 0.0 and 1.0

        eachRipple_filtered_flattened_table{ripple_idx}.rippleRelativeTimeOffsets = curr_ripple_time_offsets; 
    end
    
    %% Compute the adirectional active set for each ripple:
    eachRipple_Matricies.activeSet_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)]));
    eachRipple_Matricies.relativeSequenceIndex_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)]));
    eachRipple_Matricies.relativeProportionalTimeOffset_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)])); % eachRipple_relativeProportionalTimeOffset_Matrix: the time
    for ripple_idx = 1:length(unique_ripple_indices)
        curr_ripple_num_spikes = height(eachRipple_filtered_flattened_table{ripple_idx});
        % build a binary mask where an entry is 1 only if that cell is active at some point during this ripple:
        eachRipple_Matricies.activeSet_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = 1;
        % build an increasing series of indicies that gives the spike's position within each given ripple:
        eachRipple_Matricies.relativeSequenceIndex_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = [1:curr_ripple_num_spikes];
        
%         curr_ripple_time_offsets = eachRipple_filtered_flattened_table{ripple_idx}.time;
        curr_ripple_time_offsets = (eachRipple_filtered_flattened_table{ripple_idx}.time - curr_ripple_time_offsets(1)); % convert to relative time by subtracting the start timestamp. This means each spike's time corresponds to the time ellapsed since the start of the ripple
        curr_ripple_normalized_duration_time_offsets = curr_ripple_time_offsets ./ curr_ripple_time_offsets(end); % this should give each relative offset as a value between 0.0 and 1.0

        eachRipple_Matricies.relativeProportionalTimeOffset_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = curr_ripple_time_offsets;        
%         [0, diff(eachRipple_filtered_flattened_table{ripple_idx}.time)]

    end

    % Flatten the changes all back so we can get the extra columns (rippleRelativeSequenceIndex and rippleRelativeTimeOffsets) in the flat filtered table
    curr_filtered_table = vertcat(eachRipple_filtered_flattened_table{:});
    
end




