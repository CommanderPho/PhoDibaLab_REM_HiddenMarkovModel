% ClusteringTesting2021 - One line description of what the script performs
% Detailed explanation goes here
% 
% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 29-Oct-2021 ; Last revision: 29-Oct-2021 


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
% ------------- END OF CODE --------------


% ------------- END OF CODE --------------
