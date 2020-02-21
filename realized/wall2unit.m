function unit = wall2unit(wall,wall0,wall1)
% Convert wall times to unit times
%
% USAGE:
%   [UNIT] = wall2unit(WALL,WALL0,WALL1)
%
% INPUTS:
%   WALL    - m by 1 column vector of wall times (e.g. 93047, 134529, etc.)
%               where m>=2
%   WALL0   - Base time, maps to 0 in the unit interval
%   WALL1   - End time, maps to 1 in the unit interval
%
% OUTPUTS:
%   UNIT    - m by 1 column vector of times measured as fraction of time
%               between min(WALL) and max(WALL)
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also WALL2SECONDS, SECONDS2WALL, REALIZED_KERNEL, REALIZED_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('Three inputs required.')
end
if any(wall>=240000) || any(wall<0)
    error('WALL does not contain valid numerical times.  Numerical times should be of the form HHMMSS (e.g. 101534).');
end

if size(wall,2)>size(wall,1)
    wall=wall';
end

if size(wall,2)>1
    error('WALL must be an m by 1 column vector');
end

if length(wall)<2
    error('WALL must contain at least 2 elements');
end

% Inserted to protect against inputing integer times
wall = double(wall);

hr=floor(wall/10000);
mm=floor(wall/100)-hr*100;
ss=rem(wall,100);

if any(mm>60) || any(ss>60)
    error('WALL does not contain valid numerical times.  Numerical times should be 24-hour and of the form HHMMSS (e.g. 101534, 221313).');
end

if wall0<0 || wall0>=240000 || rem(wall0,100)>=60 || ((rem(wall0,10000)-rem(wall0,100))/100)>=60
    error('WALL0 must be a valid numerical time');
end

if wall1<0 || wall1>=240000 || rem(wall1,100)>=60 || ((rem(wall1,10000)-rem(wall1,100))/100)>=60
    error('WALL1 must be a valid numerical time');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wall0=wall2seconds(wall0);
wall1=wall2seconds(wall1);
seconds = wall2seconds(wall);
unit=(seconds-wall0)/(wall1 - wall0);

