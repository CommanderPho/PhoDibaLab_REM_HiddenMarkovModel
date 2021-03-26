function [formattedOutputString] = fnGenerateBackupTimestampString(overrideDatetime, overrideOutputFormat)
%fnGenerateBackupTimestampString Outputs the current date/time in a format suitable for building backup file names
%   overrideDatetime: Accepts a datetime value to be returned as a formatted string. Defaults to datetime('now') if not specified.
%   overrideOutputFormat: Unless overriden by specifying overrideOutputFormat:
%       the time component is expressed to hour precision, 12-hour (AM/PM) fixed-width format, in Ann Arbor, MI local time. - 'HH PM'
%       The date component is expressed as fixed-width 'yyyy-mm-dd' (Hyphenated ISO 8601) format
%
%   formattedOutputString: the datetime expressed in 'yyyy-mm-dd_T_HH PM' format.

if ~exist('overrideDate','var')
   overrideDatetime = datetime('now'); 
end

if ~exist('overrideOutputFormat','var')
   overrideOutputFormat = 'yyyy-mm-dd_T_HH PM'; 
end

formattedOutputString = datestr(overrideDatetime, overrideOutputFormat);

end

