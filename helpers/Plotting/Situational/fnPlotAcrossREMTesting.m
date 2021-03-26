function [h] = fnPlotAcrossREMTesting(mode, v1, v2, v3, v4)
% fnPlotAcrossREMTesting - One line description of what the function or script performs
% Detailed explanation goes here
% 
% Syntax:  
%     [h] = fnPlotAcrossREMTesting(mode, v1, v2, v3, v4)
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


    % v1: 36x1
    % v4: 36x86 double
    if strcmpi(mode, 'errorbar')
        h = errorbar(v1, ...
            v2, ...
            v3);
    elseif strcmpi(mode, 'scatter')
        num_repeats = size(v4, 2);
%         x = repmat(v1,[num_repeats 1]);
        x = repelem(v1, num_repeats);
        
        y = reshape(v4,[],1);
        h(2) = scatter(x, y, '.','k','MarkerFaceColor','k');
        h(2).AlphaData = repelem(0.3, length(x));
        h(2).MarkerFaceAlpha = 'flat';
        
%         h(2) = errorbar(v1, v2, v3, 'LineStyle','none');
        hold on;
        
        h(1) = scatter(v1, ...
            v2, 'filled');
        
        h(1).AlphaData = repelem(0.8, length(v1));
        h(1).MarkerFaceAlpha = 'flat';
        
        
    elseif strcmpi(mode, 'distributionPlot')
        % v4: 36x86 double
        h = distributionPlot(v4', 'xValues', v1); % defaults 
    elseif strcmpi(mode, 'bar')
        h = bar(v1, v2);
        
    elseif strcmpi(mode, 'stem')
        h = stem(v1, v2);
    else
       error('Invalid mode input!') 
    end
end


% ------------- END OF CODE --------------
