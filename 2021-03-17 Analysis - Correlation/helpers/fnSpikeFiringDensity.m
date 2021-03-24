function [gau_sdf, spk_count, tbin_edges, tbin_centers] = fnSpikeFiringDensity(spike_times, binsize, startEndTimeRange)
% fnSpikeFiringDensity - Computes the spike firing density function for the given spikes
% Detailed explanation goes here
% 
% Syntax:  
%     [gau_sdf, spk_count, tbin_edges, tbin_centers] = fnSpikeFiringDensity(spike_times, binsize, startEndTimeRange)
% 
% Inputs:
%    spike_times - Description
%    binsize - binning size in seconds, so everything else should be seconds too
%    startEndTimeRange - the minimimum and maximum start/stop time range
% 
% Outputs:
%    gau_sdf - gau_sdf
%    spk_count - spk_count
%    tbin_edges -
%    tbin_centers - 
% 

% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 24-Mar-2021 ; Last revision: 24-Mar-2021 

% ------------- BEGIN CODE --------------

    % spike_times
    tbin_edges = startEndTimeRange(1):binsize:startEndTimeRange(2);
    tbin_centers = tbin_edges(1:end-1)+binsize/2;

    if ~iscell(spike_times)
       spike_times = {spike_times}; % Add to cell array for compatibility. 
    end
    
    numSeries = size(spike_times, 1);
    numBins = length(tbin_centers);
    
    spk_count = zeros([numSeries, (numBins+1)]);
    
    for i = 1:numSeries
        spk_count(i,:) = histc(spike_times{i}, tbin_edges);        
    end
    
    % Drop the last column for all series.
    spk_count = spk_count(:, 1:end-1);

    fprintf('Done getting spike count!\n');
    % rw_sdf = conv2(spk_count,rectwin(50),'same'); % convolve with rectangular window
    % plot(tbin_centers,rw_sdf,'b');
    %  
    % gau_sdf = conv2(spk_count,gausswin(50),'same'); % convolve with gaussian window
    % plot(tbin_centers,gau_sdf,'g');


    %% Version 2:
    % binsize = 0.001; % 
    gauss_window = 1./binsize; % 1 second window
    gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
    % gaussKernel = Gauss(gauss_window, gauss_SD);
    gaussKernel = fspecial('gaussian', gauss_window, gauss_SD);
    gaussKernel = gaussKernel./binsize; % normalize by binsize
    
    
    % Begin:
    fprintf('Performing Convolution...\n');
    gau_sdf = conv2(spk_count, gaussKernel,'same'); % convolve with gaussian window
    fprintf('done!\n');

end