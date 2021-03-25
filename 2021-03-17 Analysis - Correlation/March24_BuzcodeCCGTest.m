%% Uses the Buzcode CCG function to compute the auto and cross-correlations for each provided spike train:
% Pho Hale, 2021-03-24

[ccg, t] = CCG(active_processing.spikes.time, []);
% Returns a 201x126x126 CCG



