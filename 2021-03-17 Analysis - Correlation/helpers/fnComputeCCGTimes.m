function [t, nBins] = fnComputeCCGTimes(binSize, duration)
    %% fnComputeCCGTimes
    % [t, nBins] = fnComputeCCGTimes(ccg_options.binSize, ccg_options.duration)
    halfBins = round(duration/binSize/2);
    nBins = 2*halfBins+1;
    t = (-halfBins:halfBins)'*binSize;
end