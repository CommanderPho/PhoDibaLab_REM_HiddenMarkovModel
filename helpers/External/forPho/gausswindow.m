%Generate Gaussian window of length M in vector of length n, start at j=m go to
%j=m+M-1.  Centers window at j = m+(M-1)/2, truncates to zero for j<m or
%j>m+M-1.  Uses "width" parameter "alpha", scaled according to value of M.
%Generally alpha = 2 to 3 is reasonable. Assumes vector is indexed starting
%at j = 0.

% function y = gausswindow(n, m, M, alpha)
% 
% y = zeros(1,n);
% cen = m + 0.5*(M-1);
% y = [1:n]-cen;
% y = 4*y.^2/(M-1)^2;
% y = exp(-alpha*y);
% y(1:(m-1)) = 0.0; y(m+M:n) = 0.0;


function y = gausswindow(sigma, halfwidth)
    y = gauss((halfwidth*2), sigma);
    % y = zeros(1,n);
    % cen = m + 0.5*(M-1);
    % y = [1:n]-cen;
    % y = 4*y.^2/(M-1)^2;
    % y = exp(-alpha*y);
    % y(1:(m-1)) = 0.0; y(m+M:n) = 0.0;

    