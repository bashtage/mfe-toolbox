function seconds = wall2seconds(wall)
% Convert wall times to seconds past midnight
%
% USAGE:
%   [SECONDS] = wall2seconds(WALL)
%
% INPUTS:
%   WALL      - m by 1 column vector of wall times (e.g. 93047, 134529, etc.)
%
% OUTPUTS:
%   SECONDS   - m by 1 column vector of times measured as seconds past midnight
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also SECONDS2WALL, REALIZED_KERNEL, REALIZED_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=1
    error('One inputs required.')
end
if any(wall>=240000) || any(wall<0)
    error('WALL does not contain valid wall times.  Wall times should be of the form HHMMSS (e.g. 101534).');
end
% Inserted to protect against inputing integer times
wall = double(wall);

hr=floor(wall/10000);
mm=floor(wall/100)-hr*100;
ss=rem(wall,100);

if any(mm>60) || any(ss>60)
    error('WALL does not contain valid wall times.  Wall times should be of the form HHMMSS (e.g. 101534).');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


seconds = 3600 * hr + 60 * mm + ss;

