%% Working on 2021-10-26

% Load Files:
% load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/Results/PhoResults_Expt1_RoyMaze1.mat');
% load('/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/Data/Rotation_3_Kamran Diba Lab/DataProcessingProject/Hiro_Datasets/analysesResults_13-Oct-2021/Roy-maze1/toAddVariables.mat');

%% Windows Loading paths:
% load('C:\Share\data\RoyMaze1\PhoResults_Expt1_RoyMaze1.mat');
% load('C:\Share\data\RoyMaze1\analysesResults_13-Oct-2021\toAddVariables.mat');


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


%% We come in with a table where each row is a unit (active_processing.spikes)
%% From there we build 

%% Temporary copy of spikes table with subset of columns that will be filtered pre-processed to prepare for flattening:
% all table variables are cell arrays of different sizes
curr_unit_rows_table = table(active_processing.spikes.time, ...
        active_processing.spikes.isRippleSpike, ...
        active_processing.spikes.RippleIndex, ...
        fnCellContentsTranpose(active_processing.spikes.behavioral_duration_indicies), fnCellContentsTranpose(active_processing.spikes.behavioral_states), fnCellContentsTranpose(active_processing.spikes.behavioral_epoch), ...
        'VariableNames', {'time','isRippleSpike','RippleIndex', 'behavioral_duration_indicies', 'behavioral_states', 'behavioral_epoch'});
% Filter the inactive units from curr_unit_rows_table to save on processing overhead before flattening the cells
curr_unit_rows_table = curr_unit_rows_table(plot_outputs.filter_active_units, :);
num_active_units = height(curr_unit_rows_table);

%% Flatten over the rows (which are units) subset of spikes table for efficient ripple-related processing:
[curr_flattenedOverUnits_table] = fnFlattenCellsToContents(curr_unit_rows_table); % Flatten the cells:

%% Removes rows with missing values, meaning the rows with RippleIndex == NaN (meaning they aren't ripple spikes)
curr_filtered_spikeFlattened_table = rmmissing(curr_flattenedOverUnits_table);


%% Exclude sleep states:
% curr_filtered_spikeFlattened_table = curr_filtered_spikeFlattened_table(('rem' == curr_filtered_spikeFlattened_table.behavioral_states), :); % Only those occuring during rem
% curr_filtered_spikeFlattened_table = curr_filtered_spikeFlattened_table(('nrem' == curr_filtered_spikeFlattened_table.behavioral_states), :);
% curr_filtered_spikeFlattened_table = groupfilter(curr_filtered_spikeFlattened_table, 'groupID', @(x) all(x == 'rem'), 'behavioral_states');
% 
% curr_filtered_spikeFlattened_table = curr_filtered_spikeFlattened_table((('rem' == curr_filtered_spikeFlattened_table.behavioral_states) | ('nrem' == curr_filtered_spikeFlattened_table.behavioral_states)), :);

target_options.behavioral_states_variable_name = 'behavioral_states';

%% "Quiet" - 1

[track_quiet_included] = fnFilterSpikesWithCriteria(curr_filtered_spikeFlattened_table, {'pre_sleep'}, {'nrem','rem'}, target_options);
% [track_quiet_included] = fnFilterSpikesWithCriteria(curr_filtered_spikeFlattened_table, {'track'}, {'quiet'}, target_options);
[track_quiet_outputs.filtered_table, track_quiet_outputs.eachRipple_filtered_flattened_table, track_quiet_outputs.eachRipple_Matricies] = fnSplitIntoSpecificRipples(curr_filtered_spikeFlattened_table, track_quiet_included, num_active_units);


%% "Active" - 2
[track_active_included] = fnFilterSpikesWithCriteria(curr_filtered_spikeFlattened_table, {'post_sleep'}, {'nrem','rem'}, target_options);
% [track_active_included] = fnFilterSpikesWithCriteria(curr_filtered_spikeFlattened_table, {'track'}, {'active'}, target_options);
[track_active_outputs.filtered_table, track_active_outputs.eachRipple_filtered_flattened_table, track_active_outputs.eachRipple_Matricies] = fnSplitIntoSpecificRipples(curr_filtered_spikeFlattened_table, track_active_included, num_active_units);




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

% [track_quiet_outputs.reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(track_quiet_outputs.filtered_table);
% [track_quiet_outputs.reconstructedCellTable] = fnReconstructCellsFromFlattenedContents(track_quiet_outputs.filtered_table);


%% Pairwise Spiking Events Implementation Attempt: 
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

% %% Another attempt at processing pairwise relations between units spiking events:
% % track_quiet_outputs.eachRipple_filtered_flattened_table{111};
% D = pdist(track_active_outputs.eachRipple_Matricies.activeSet_Matrix,'hamming'); % The percentage of each sequence's coordinates that differ
% % [ist,ind,dst] = findsignal(corr,sgn,'TimeAlignment','dtw');
% D_sqr = squareform(D);



%% Plot the outputs:

plotting_options.ordered_by = 'participation_rate';
plotting_options.ordering_matrix_variable_name = 'activeSet_Matrix';

% plotting_options.currMatrixName = 'activeSet_Matrix';
% plotting_options.currMatrixName = 'relativeSequenceIndex_Matrix';
plotting_options.currMatrixName = 'relativeProportionalTimeOffset_Matrix';

%% "Quiet" - 1
plotting_options.figure_title = 'pre sleep';
figure(1);
clf;
[~, ordering_matrics, ordering_indicies] = fnPlotRippleMatrixComparisonResult(track_quiet_outputs.eachRipple_Matricies, plotting_options);

%% "Active" - 2
plotting_options.figure_title = 'post sleep';
% uses the ordering provided from the previous plot
plotting_options.ordered_by = 'provided_indicies';
plotting_options.provided_indicies = ordering_indicies;
figure(2);
clf;
[~] = fnPlotRippleMatrixComparisonResult(track_active_outputs.eachRipple_Matricies, plotting_options);


% print('-clipboard','-dmeta')

function [fig, ordering_matrics, ordering_indicies] = fnPlotRippleMatrixComparisonResult(results_matrix_struct, plotting_options)
    %% fnPlotRippleMatrixComparisonResult: Plots the spikes during each ripple based on options provided in plotting_options
    % 2021-10-26 Pho Hale
    % note: assumes all ordering is done on the first column of the matrix (corresponding to the units)

    fig = gcf;
    active_mat = results_matrix_struct.(plotting_options.currMatrixName);
    
    % Order the output if specified:
    active_ordering_mat = results_matrix_struct.(plotting_options.ordering_matrix_variable_name);
    % calculate particpation_rate for sorting)
    ordering_matrics.participation_rate = sum(active_ordering_mat, 2); % count the participation of each unit for sorting
    
    if strcmpi(plotting_options.ordered_by, 'participation_rate')
        [~, ordering_indicies] = sort(ordering_matrics.participation_rate, "descend");
        active_mat = active_mat(ordering_indicies, :);
    elseif strcmpi(plotting_options.ordered_by, 'provided_indicies')
        % provided_indicies: uses user-specified indicies provided in plotting_options.provided_indicies
        if ~isfield(plotting_options, 'provided_indicies')
            error('plotting_options.ordered_by was set to "provided_indicies" but no plotting_options.provided_indicies was provided!');
        end
        ordering_indicies = plotting_options.provided_indicies;
        active_mat = active_mat(ordering_indicies, :);

    else
        % no change in ordering
        ordering_indicies = [1:size(active_mat,1)];
    end

    % Plot:
    fnPhoMatrixPlot(active_mat)

    title(plotting_options.figure_title)
    xlabel('SWR Index')
    ylabel('Unit ID')
    xticks([]); yticks([]); % The tick marks are wrong due to a bug with imagesc, so just remove them)    
end



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


function [curr_filtered_spikeFlattened_table, eachRipple_filtered_flattened_table, eachRipple_Matricies] = fnSplitIntoSpecificRipples(target_table, is_target_entry_included, num_active_units)
    %% fnSplitIntoSpecificRipples: splits the tables of entries describing * to a cell array (eachRipple_filtered_flattened_table) of small tables, each containing information about the spikes that occured within the given ripple.
    %   These can be used for easily processing in a ripple-by-ripple manner
    % also returns matricies that can be used to compare unit activity across ripples
    curr_filtered_spikeFlattened_table = target_table(is_target_entry_included, :);
    
    %% Filter again to a specific ripple index
    unique_ripple_indices = unique(curr_filtered_spikeFlattened_table.RippleIndex);
    %% Build a cell array (eachRipple_filtered_flattened_table) of small tables, each containing information about the spikes that occured within the given ripple. This will be returned and also used in the next step to produce the matricies
    eachRipple_filtered_flattened_table = cell([length(unique_ripple_indices), 1]);
    for ripple_idx = 1:length(unique_ripple_indices)
        eachRipple_filtered_flattened_table{ripple_idx} = curr_filtered_spikeFlattened_table(curr_filtered_spikeFlattened_table.RippleIndex == unique_ripple_indices(ripple_idx), :);
        curr_ripple_num_spikes = height(eachRipple_filtered_flattened_table{ripple_idx});
        eachRipple_filtered_flattened_table{ripple_idx}.rippleRelativeSequenceIndex = [1:curr_ripple_num_spikes]';

        curr_ripple_time_offsets = eachRipple_filtered_flattened_table{ripple_idx}.time;
        curr_ripple_time_offsets = (eachRipple_filtered_flattened_table{ripple_idx}.time - curr_ripple_time_offsets(1)); % convert to relative time by subtracting the start timestamp. This means each spike's time corresponds to the time ellapsed since the start of the ripple
%         curr_ripple_normalized_duration_time_offsets = curr_ripple_time_offsets ./ curr_ripple_time_offsets(end); % this should give each relative offset as a value between 0.0 and 1.0
        eachRipple_filtered_flattened_table{ripple_idx}.rippleRelativeTimeOffsets = curr_ripple_time_offsets; 
    end
    
    %% Build a struct eachRipple_Matricies containing fields with several computed matricies that will be returned for plotting:
    eachRipple_Matricies.activeSet_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)])); %% Compute the adirectional active set for each ripple:
    eachRipple_Matricies.relativeSequenceIndex_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)]));
    eachRipple_Matricies.relativeProportionalTimeOffset_Matrix = sparse(zeros([num_active_units, length(unique_ripple_indices)])); % eachRipple_relativeProportionalTimeOffset_Matrix: the time
    for ripple_idx = 1:length(unique_ripple_indices)
        curr_ripple_num_spikes = height(eachRipple_filtered_flattened_table{ripple_idx});
        % build a binary mask where an entry is 1 only if that cell is active at some point during this ripple:
        eachRipple_Matricies.activeSet_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = 1;
        % build an increasing series of indicies that gives the spike's position within each given ripple:
        eachRipple_Matricies.relativeSequenceIndex_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = [1:curr_ripple_num_spikes];
        
        curr_ripple_time_offsets = eachRipple_filtered_flattened_table{ripple_idx}.rippleRelativeTimeOffsets;
        eachRipple_Matricies.relativeProportionalTimeOffset_Matrix(eachRipple_filtered_flattened_table{ripple_idx}.flattened_UnitIDs, ripple_idx) = curr_ripple_time_offsets;        
    end

    % Flatten the changes all back so we can get the extra columns (rippleRelativeSequenceIndex and rippleRelativeTimeOffsets) in the flat filtered table
    curr_filtered_spikeFlattened_table = vertcat(eachRipple_filtered_flattened_table{:});
    
end




