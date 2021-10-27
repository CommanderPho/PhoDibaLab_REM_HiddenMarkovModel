function firingRates = computeFiringRates (spikeTimes, tstart, tend, tau)
nCells = length(spikeTimes);
nTimeBins = ceil((tend-tstart)/tau);
win = linspace(tstart, tend, nTimeBins);
firingRates = zeros(nCells, nTimeBins-1);
for i = 1:nCells
    firingRates(i,:) = histcounts(spikeTimes{i}, win)/tau;
end