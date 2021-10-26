function [filtered_spike_times] = fnFilterSpikesWithCriteria(spikesTable, included_ripple_index)
%FNFILTERSPIKESWITHCRITERIA Summary of this function goes here
%   Detailed explanation goes here

if 

included_cell_types

active_processing.spikes.time(plot_outputs.filter_active_units)

filtered_spike_times = cellfun(@(spike_times, is_spike_ripple) spike_times(is_spike_ripple), spikesTable.time, spikesTable.isRippleSpike, 'UniformOutput', false); %% Filtered to only show the ripple spikes

end

