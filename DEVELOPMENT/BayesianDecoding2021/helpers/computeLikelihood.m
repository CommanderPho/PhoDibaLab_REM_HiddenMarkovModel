function likelihood = computeLikelihood(spikeCounts, placeFields, tau)

    pFields = (placeFields * tau)';
    xyBins = size(placeFields,2);
    nTimeBins = size(spikeCounts,1);
    likelihood = zeros(xyBins, nTimeBins);
    parfor i = 1:nTimeBins
        nSpikes = repmat(spikeCounts(i,:), xyBins, 1);
        maxL = poisspdf(nSpikes,pFields);
        maxL = prod(maxL,2);
        likelihood(:,i) = maxL;
    end 

end