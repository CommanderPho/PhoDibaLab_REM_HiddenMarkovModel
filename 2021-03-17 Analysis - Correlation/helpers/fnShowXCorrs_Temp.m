function fnShowXCorrs_Temp(temp, lag_offsets, xcorr_input, active_pair_index)
% fnShowXCorrs_Temp - One line description of what the function or script performs (H1 line)
% Optional file header info (to give more details about the function than in the H1 line)
% 
% Syntax:  
%     fnShowXCorrs_Temp(temp, lag_offsets, xcorr_input, active_pair_index)
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
% Created: 24-Mar-2021 ; Last revision: 24-Mar-2021 

% ------------- BEGIN CODE --------------


    temp.by_behavioral_period.curr_xcorr_forPair = squeeze(xcorr_input(:, active_pair_index, :)); % [num_of_behavioral_state_periods x num_lag_steps]
    
    %% Single Plot for All Time:
    figure
    subplot(1,1,1)
    temp.plotVals = sum(temp.by_behavioral_period.curr_xcorr_forPair(:, :),1);
    temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
    fnPlotXCorrStem(lag_offsets, temp.plotVals, 'all');
    sgtitle(temp.curr_fig_title);
    
    % Separate plots for each state of interest:
    temp.n_subplots.rows = 3;
    figure
    subplot(temp.n_subplots.rows,1,1)
    temp.plotVals = sum(temp.by_behavioral_period.curr_xcorr_forPair(temp.filtered.pre_sleep_REM_indicies, :),1);
    temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
    stem(lag_offsets, temp.plotVals);
    xline(0);
    title('pre_sleep_REM_indicies','Interpreter','none')
    
    subplot(temp.n_subplots.rows,1,2)
    temp.plotVals = sum(temp.by_behavioral_period.curr_xcorr_forPair(temp.filtered.all_except_REM_indicies, :),1);
    temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
    stem(lag_offsets, temp.plotVals);
    xline(0);
    title('all_except_REM_indicies','Interpreter','none')
    
    subplot(temp.n_subplots.rows,1,3)
    temp.plotVals = sum(temp.by_behavioral_period.curr_xcorr_forPair(temp.filtered.post_sleep_REM_indicies, :),1);
    temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
    stem(lag_offsets, temp.plotVals);
    xline(0);
    title('post_sleep_REM_indicies','Interpreter','none')
    
    % subplot(4,1,4)
    % temp.plotVals = sum(temp.by_behavioral_period.curr_xcorr_forPair(:, :),1);
    % temp.plotVals = temp.plotVals ./ temp.plotVals(temp.curr_xcorr_lag_zero_offset); %normalize so that zero-lag has a height of 1
    % stem(lag_offsets, temp.plotVals);
    % xline(0);
    % title('all','Interpreter','none')
    
    
    sgtitle(temp.curr_fig_title)
    
    
end


function fnPlotXCorrStem(lag_offsets, plotVals, curr_fig_title)
    stem(lag_offsets, plotVals);
    xline(0);
    title(curr_fig_title,'Interpreter','none')
end


% ------------- END OF CODE --------------
