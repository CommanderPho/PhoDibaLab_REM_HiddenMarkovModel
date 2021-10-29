%% Partition each ripple


%% for each unique ripple:
    % split the table of spikes on a given unitID of interest, resulting in groups of spikes

   [G,powerLosses] = findgroups(T1);
powerLosses.maxLoss = splitapply(@max,T.Loss,G);


%% Find inactive set for the active periods on the track by setting a threshold for including a unit in the active set by firing rate.

