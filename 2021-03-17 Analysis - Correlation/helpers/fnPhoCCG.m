%fnPhoCCG - Compute multiple cross- and auto-correlograms
%	- Refactored from Buzcode/CCG.m

%  USAGE
%
%    [ccg,t] = fnPhoCCG(times,groups,<options>)
%
%    times          times of all events
%                   NOTE: spiketimes in SECONDS.
%    groups         group IDs for each event in time list (should be
%                   integers 1:nGroups)
%
%   alternate buzcode usage: 
%   times           {Ncells} array of [Nspikes]  spiketimes for each cell.
%   groups          []
%
%                   spikes = bz_GetSpikes
%                   [ccg,t] = CCG(spikes.times,[])  
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'norm'        normalization of the CCG, 'counts' or 'rate' (DL added 8/1/17)
%                   'counts' gives raw event/spike count,
%                   'rate' returns CCG in units of spks/second (default: counts)
%    =========================================================================
%
%
%  OUTPUT
%   ccg     [t x ngroups x ngroups] matrix where ccg(t,i,j) is the
%           number (or rate) of events of group j at time lag t with  
%           respect to reference events from group i
%   t       time lag vector (units: seonds)
%
%  SEE
%
%    See also ShortTimeCCG.

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [ccg,t, cell_conversion] = fnPhoCCG(times,groups,varargin)

% Default values
duration = 2;
binSize = 0.01;
Fs = 1/20000;
normtype = 'counts';

% Option for spike times to be in {Ncells} array of spiketimes DL2017
if iscell(times) && isempty(groups)
    numcells = length(times);
    cell_conversion.has_empty_groups = false;
    cell_conversion.group_times_size = zeros([numcells 1]);
    cell_conversion.group_is_empty = false([numcells 1]);
    
    for cc = 1:numcells
        cell_conversion.group_times_size(cc) = length(times{cc});
        cell_conversion.group_is_empty(cc) = isempty(times{cc});
        
        if cell_conversion.group_is_empty(cc)
          groups{cc}=cc;   % add a single placeholder entry to be zeroed out later.
        else
        groups{cc}=cc.*ones(size(times{cc}));
        end
    end
    times = cat(2,times{:}); 
    groups = cat(2,groups{:}); % Unwraps the cell array of doubles to a single array of doubles. The membership (cell) they belong to is tracked instead by 'groups'
end

%Sort
[times,sortidx] = sort(times);
groups = groups(sortidx);

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
%if ~isdvector(times),
%	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
%end
if ~isdscalar(groups) && ~isdvector(groups),
	error('Parameter ''groups'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(groups) && length(times) ~= length(groups),
	error('Parameters ''times'' and ''groups'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
groups = groups(:);
times = times(:);

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			%if ~isdscalar(binSize,'>0'),
		%		error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
	%		end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
            
            case 'Fs',
			Fs = varargin{i+1};
			if ~isdscalar(Fs,'>0'),
				error('Incorrect value for property ''Fs'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
        case 'norm'
            normtype = varargin{i+1};
     
  end
end



% Number of groups, number of bins, etc.
if length(groups) == 1,
	groups = ones(length(times),1);
	nGroups = 1;
else
    % Error: The input cell array's groups aren't preserved if it contains empty elements.
	nGroups = max(unique(groups));
    
    if exist('numcells','var')
        if nGroups < numcells
            fprintf('nGroups < numCells: fixing.\n');
            cell_conversion.has_empty_groups = true;
           nGroups = max(nGroups, numcells); 
        end
    end
end


halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;
times = round(times/Fs);
binSize_Fs = round(binSize/Fs);



if length(times) <= 1,
	% ---- MODIFIED BY EWS, 1/2/2014 ----
    % *** Use unsigned integer format to save memory ***
    ccg = uint16(zeros(nBins,nGroups,nGroups));
    % -----------------------------------
    return
end

% Compute CCGs
nEvents = length(times);
% 

% gotta see how this function determines what to return.
counts = double(CCGHeart(times,uint32(groups),binSize_Fs,uint32(halfBins)));
% -----------------------------------
% 
% Reshape the results



n = max(groups);
counts = reshape(counts,[nBins n n]);


% if exist('cell_conversion')
%     
%     if cell_conversion.has_empty_groups
%          % Here is where we must fix cell array stuff:
%         fprintf('Fixing groups:\n');
%         
%         % The groups that these belong to are those that aren't empty:
%         non_empty_indicies = find(~cell_conversion.group_is_empty);
%         
%         full_counts = uint16(zeros(nBins,nGroups,nGroups));
%         full_counts(:, non_empty_indicies, non_empty_indicies) = counts;
%         counts = full_counts; % copy back to counts
%         fprintf('done.\n');
%     
%     end
% end

if n < nGroups,
	counts(nBins,nGroups,nGroups) = 0;
end

%Rate normalization: counts/numREFspikes/dt to put in units of spikes/s. DL
switch normtype
    case 'rate'
        for gg = 1:nGroups
            numREFspikes = sum(groups==gg);%number of reference events for group
            counts(:,gg,:) = counts(:,gg,:)./numREFspikes./binSize;
        end
end
        
        
ccg = flipud(counts);