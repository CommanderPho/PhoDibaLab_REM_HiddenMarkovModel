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
%     y = gauss((halfwidth*2), sigma);
%     y = gausswin((halfwidth*2), sigma);

%     incr = 2*sds/(frames-1);
%     y = exp(-(-sds:incr:sds).^2);
     gauss_window = (halfwidth*2); % 1 second window
     gauss_SD = sigma; % 0.02 seconds (20ms) SD

%     gauss_window = 1./binsize; % 1 second window
%     gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
%     % gaussKernel = Gauss(gauss_window, gauss_SD);
%     y = fspecial('gaussian', gauss_window, gauss_SD);
%     y = y./binsize; % normalize by binsize
    
    % y = zeros(1,n);
    % cen = m + 0.5*(M-1);
    % y = [1:n]-cen;
    % y = 4*y.^2/(M-1)^2;
    % y = exp(-alpha*y);
    % y(1:(m-1)) = 0.0; y(m+M:n) = 0.0;

    y = subfnGauss(gauss_window, gauss_SD);


end

function outvec = subfnGauss(window_length, std_devs)
    % subfnGausss() - return a smooth Gaussian window
    %
    % Usage:
    %   >> outvector = gauss(frames,sds);
    %
    % Inputs:
    %   frames = window length
    %   sds    = number of +/-std. deviations = steepness 
    %            (~0+ -> flat; >>10 -> spike)
    
    outvec = [];
    if nargin < 2
      help gauss
      return
    end
    if std_devs <=0 | window_length < 1
      help gauss
      return
    end
    
    incr = 2*std_devs/(window_length-1);
    outvec = exp(-(-std_devs:incr:std_devs).^2);

end
