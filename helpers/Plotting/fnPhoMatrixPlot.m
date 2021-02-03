function [h, info] = fnPhoMatrixPlot(data)
%FNPHOMATRIXPLOT Plots a 2D matrix of unnormalized data.
%   data should be a 2D matrix.
info.dim.x = size(data, 1);
info.dim.y = size(data, 2);

% [xx, yy] = meshgrid(1:dim.x, 1:dim.y);
% h = plot3(xx, yy, data);

xx = [1:info.dim.x];
yy = [1:info.dim.y];

info.pixel_width = (xx(2)-xx(1))/(size(data,2)-1);
info.pixel_height = (yy(2)-yy(1))/(size(data,1)-1);

h = imagesc(xx, yy, data);
colorbar

end

