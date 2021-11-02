% function linearPos = linearizePosition2(fileinfo, behavior, mazeShape)
function [linearPos, userSelectedCenters] = linearizePosition2(xyt, maze_start_time, maze_end_time, mazeShape, userSelectedCenters)
% Interactive: requires selection of the centers
%% Conversion from old call style: linearizePosition2(fileinfo.xyt, behavior.time(2,1), behavior.time(2,2), mazeShape);

time = xyt(:, 3);
xpos = xyt(:, 1);
ypos = xyt(:, 2);

linearPos = nan(length(xyt), 1);

runIdx = time > maze_start_time & time < maze_end_time; % only the positions on the maze

xpos = xpos(runIdx);
ypos = ypos(runIdx);

nTimeBins = length(xpos);

runLinearPos = nan(nTimeBins, 1);

nonNanBins = find(~isnan(xpos) & ~isnan(ypos));


figure(1);
plot(xpos, ypos, '.', 'markersize', 2, 'color', [.7 .7 .7])


switch mazeShape
	case 'linear'
		if ~exist('userSelectedCenters','var') | (length(userSelectedCenters) ~= 2)
			userSelectedCenters = ginput(2);
		end
		centers = userSelectedCenters;
		
		projectedPnts = zeros(nTimeBins, 2);
		for ii = 1: numel(nonNanBins)
			bin = nonNanBins(ii);            
			projectedPnts(bin, :) = projectPnt([xpos(bin) ypos(bin)], centers(1, :), centers(2, :), 'none');
		end
		
		minX = nanmin(projectedPnts(:,1));
		minY = nanmin(projectedPnts(:,2));
		
		projectedPnts(:,1) = projectedPnts(:,1) - minX;
		projectedPnts(:,2) = projectedPnts(:,2) - minY;
		
		runLinearPos = sqrt(sum(projectedPnts.^2, 2));
	case 'L-shape' % we can add an imaginary arc track where the platform is to avoid omitting the position data on the platform
		if ~exist('userSelectedCenters','var') | (length(userSelectedCenters) ~= 3)
			userSelectedCenters = ginput(3);
		end
		centers = userSelectedCenters;

		projectedPnts1 = zeros(nTimeBins, 2);
		projectedPnts2 = zeros(nTimeBins, 2);
		
		dist1 = zeros(nTimeBins, 1);
		dist2 = zeros(nTimeBins, 1);
		
		for ii = 1: numel(nonNanBins)
			
			 bin = nonNanBins(ii);
				
			% prjojection on line segment 1
			[projectedPnts1(bin, :), dist1(bin)] = projectPnt([xpos(bin) ypos(bin)], centers(1, :), centers(2, :), 'second');


			% projection on line segment 2
			[projectedPnts2(bin, :), dist2(bin)] = projectPnt([xpos(bin) ypos(bin)], centers(2, :), centers(3, :), 'first');

		end
		
		segment1_Projections = projectedPnts1(dist1 < dist2, :);
		[segment1_length, maxIndex] = max(sqrt(sum((segment1_Projections - centers(2, :)).^2, 2)));
		segment1_firstPnt = segment1_Projections(maxIndex, :);
		
		segment2_firstPnt = centers(2, :);
		
		projectedPnts = zeros(nTimeBins, 2);
		for ii = 1:numel(nonNanBins)
			
			bin = nonNanBins(ii);
			
			[~, idx] = min([dist1(bin) dist2(bin)]);
			
			switch idx
				
				case 1
					projectedPnts(bin, :) = projectedPnts1(bin, :);
					runLinearPos(bin)     = norm(projectedPnts(bin, :)- segment1_firstPnt);
					
				case 2
					projectedPnts(bin, :) = projectedPnts2(bin, :);
					runLinearPos(bin)     = norm(projectedPnts(bin, :) - segment2_firstPnt)+ segment1_length;
			end
			
		end
					
		
			  
		
	case 'U-shape' % could be z-shape too, whatever maze with 3 arms/segments 
		if ~exist('userSelectedCenters','var') | (length(userSelectedCenters) ~= 4)
			userSelectedCenters = ginput(4);
		end
		centers = userSelectedCenters;
		
		projectedPnts1 = zeros(nTimeBins, 2);
		projectedPnts2 = zeros(nTimeBins, 2);
		projectedPnts3 = zeros(nTimeBins, 2);
		
		dist1 = zeros(nTimeBins, 1);
		dist2 = zeros(nTimeBins, 1);
		dist3 = zeros(nTimeBins, 1);
		
		for ii = 1: numel(nonNanBins)
			
			 bin = nonNanBins(ii);
				
			% prjojection on line segment 1
			[projectedPnts1(bin, :), dist1(bin)] = projectPnt([xpos(bin) ypos(bin)], centers(1, :), centers(2, :), 'second');


			% projection on line segment 2
			[projectedPnts2(bin, :), dist2(bin)] = projectPnt([xpos(bin) ypos(bin)], centers(2, :), centers(3, :), 'both');
			
			% projection on line segment 3
			[projectedPnts3(bin, :), dist3(bin)] = projectPnt([xpos(bin) ypos(bin)], centers(3, :), centers(4, :), 'first');

		end
		
		segment1_Projections = projectedPnts1(dist1 < dist2 & dist1 < dist3, :);
		[segment1_length, maxIndex] = max(sqrt(sum((segment1_Projections - centers(2, :)).^2, 2)));
		segment1_firstPnt = segment1_Projections(maxIndex, :);
		
		segment2_length = norm(centers(3,:)-centers(2,:));
		segment2_firstPnt = centers(2,:);
		
		segment3_firstPnt = centers(3,:);
		
		projectedPnts = zeros(nTimeBins, 2);
		for ii = 1:numel(nonNanBins)
			bin = nonNanBins(ii);
			[~, idx] = min([dist1(bin) dist2(bin) dist3(bin)]);
			
			switch idx
				case 1
					projectedPnts(bin, :) = projectedPnts1(bin, :);
					runLinearPos(bin)     = norm(projectedPnts(bin, :)- segment1_firstPnt);
				case 2
					projectedPnts(bin, :) = projectedPnts2(bin, :);
					runLinearPos(bin)     = norm(projectedPnts(bin, :)- segment2_firstPnt) + segment1_length;
				case 3
					projectedPnts(bin, :) = projectedPnts3(bin, :);
					runLinearPos(bin)     = norm(projectedPnts(bin, :)- segment3_firstPnt) + segment1_length + segment2_length;
			end
			
		end
		
		
	case 'circular'
		
		cx = min(xpos) + range(xpos)/2;
		cy = min(ypos) + range(ypos)/2;
		
		radius = nanmean(sqrt((xpos - cx).^2 + (ypos - cy).^2));
		
		xgrid = transpose(linspace(cx-radius, cx+radius, 100));
		
		trackX = [xgrid(end:-1:1); xgrid];
		trackY = [cy+sqrt(radius^2 - (xgrid - cx).^2); cy-sqrt(radius^2 - (xgrid - cx).^2)];
		
		% calcualte the projections
		
		projectedPnts = zeros(nTimeBins, 2);
		for ii = 1: numel(nonNanBins)
			
			bin = nonNanBins(ii);
		
			distFromCntr = sqrt((xpos(bin) - cx)^2 + (ypos(bin) - cy)^2);
			
			projectedPnts(bin, 1) = radius * (xpos(bin) - cx)/distFromCntr + cx;
			projectedPnts(bin, 2) = radius * (ypos(bin) - cy)/distFromCntr + cy;
			
			px = projectedPnts(bin, 1);
			py = projectedPnts(bin, 2);
			
			if px >= cx && py >= cy
			   theta = atan(abs(py-cy)/abs(px-cx));
			elseif px < cx && py >= cy
			   theta = (pi/2-atan(abs(py-cy)/abs(px-cx)))+pi/2;
			elseif px < cx && py < cy
			   theta = atan(abs(py-cy)/abs(px-cx))+pi;
			elseif px >= cx && py < cy
			   theta = (pi/2-atan(abs(py-cy)/abs(px-cx))) + 3*pi/2;
			end
   
			runLinearPos(bin) = radius* (2*pi-theta); % if the animal is moving clockwise, use (2*pi-theta) and if ccw use theta

		end % end for ii = ...
		

end % end switch mazeShape


hold on
if exist('centers','var')
    scatter(centers(:,1), centers(:,2), 42.0, 'b', 'filled')
%     plot(projectedPnts(:,1)+centers(1,1), projectedPnts(:,2)+centers(1,2), '.r', 'linewidth', 2)
end
plot(projectedPnts(:,1), projectedPnts(:,2), '.r', 'linewidth', 2)

linearPos(runIdx) = runLinearPos;

figure(2); 
plot(runLinearPos, '.r')


end % end function 


function [Q, dist] = projectPnt(pnt, V1, V2, removedPlatforms)

	v = (V2- V1)/norm(V2-V1);
	Q = dot(pnt-V1, v)*v + V1;

	dist = norm(pnt-Q);

	flag = 0;
	segmentLen = norm(V2 - V1);
	switch removedPlatforms
		
		case 'first'
			flag = norm(Q - V2) > segmentLen;
			
		case 'second'
			flag = norm(Q - V1) > segmentLen;
			
		case 'both'
			flag = (norm(Q - V1) > segmentLen) || (norm(Q - V2) > segmentLen);
			
		case 'none'
			
			flag = 0;
	end

	if flag
		Q = [nan nan];
	end

end % end function 


