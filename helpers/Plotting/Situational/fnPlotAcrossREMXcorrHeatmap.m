function [xcorr_fig, handles] = fnPlotAcrossREMXcorrHeatmap(varargin)
% fnPlotAcrossREMXcorrHeatmap - Plot a heatmap of the xcorr
% Detailed explanation goes here
% 
% Syntax:  
%     [xcorr_fig, handles] = fnPlotAcrossREMXcorrHeatmap(varargin)
% 
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
% 
% Outputs:
%    output1 - Description
%    output2 - Description
% 

% Author: Pho Hale
% PhoHale.com 
% email: halechr@umich.edu
% Created: 26-Mar-2021 ; Last revision: 26-Mar-2021 

% ------------- BEGIN CODE --------------

%     xcorr_fig = figure(15);
    xcorr_fig = gcf;
    clf
 
    num_heatmaps = nargin;   
    
    for i = 1:num_heatmaps
        subplot(num_heatmaps,1,i);
        
%         handles(i) = heatmap(varargin{i},'GridVisible', false);
        num_timesteps = size(varargin{i},2);
        zero_timestep = num_timesteps / 2;
        
        handles(i) = imagesc(varargin{i});
        
        ylabel('Filtered Trial Index')
        
%         title(sprintf('Periods: %d', size(varargin{i},1)));
        
%         xlabel('Time Lag')
%         xline(zero_timestep, 'red');
        
        
        
    end
    
%     h1 = heatmap(v1);
%     ylabel('Filtered Trial Index')
%     xlabel('Time Lag')
%     title(sprintf('PRE sleep REM periods: %d', size(v1,1)));
% 
% 
%     subplot(num_heatmaps,1,2);
%     h2 = heatmap(v2);
%     ylabel('Filtered Trial Index')
%     xlabel('Time Lag')    
%     title(sprintf('POST sleep REM periods: %d', size(v2,1)));
%     sgtitle('XCorr for all pairs - PRE vs Post Sleep REM Periods - Period Index - All Units')
end


% ------------- END OF CODE --------------
