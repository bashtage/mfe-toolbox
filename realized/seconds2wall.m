function wall = seconds2wall(seconds)
% Convert seconds past midnight to wall times
%
% USAGE:
%   [WALL] = wall2seconds(SECONDS)
%
% INPUTS:
%   SECONDS   - m by 1 column vector of times measured as seconds past midnight
%
% OUTPUTS:
%   WALL      - m by 1 column vector of wall times (e.g. 93047, 134529, etc.)
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also WALL2SECONDS, REALIZED_KERNEL, REALIZED_VARIANCE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=1
    error('One inputs required.')
end
if any(seconds>=(24*3600)) || any(seconds<0)
    error('SECONDS does not contain valid numerical times.  Numerical times must satisft 0<=SECONDS<86400');
end
seconds = double(seconds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Parse and recombine
hr = floor(seconds/3600);
mm = floor((seconds - hr*3600)/60);
ss = rem(seconds,60);
wall = hr*10000 + mm*100 + ss;

