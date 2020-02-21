function wall = unit2wall(unit, wall0, wall1)
% Convert unit times to wall times 
%
% USAGE:
%   [WALL] = wall2unit(UNIT,WALL0,WALL1)
%
% INPUTS:
%   UNIT    - m by 1 column vector of times measured as fraction of time
%               between min(WALL0) and max(WALL1)
%   WALL0   - Base time, maps to 0 in the unit interval
%   WALL1   - End time, maps to 1 in the unit interval
%
% OUTPUTS:
%   WALL    - m by 1 vector of times measures in 24-hour wall time, such as
%               103421 for 10:34:21 AM or 145311 for 2:53:11 PM
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also WALL2UNIT, WALL2SECONDS, SECONDS2WALL, REALIZED_KERNEL, REALIZED_VARIANCE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('Three inputs required.')
end
 
if size(unit,2)>size(unit,1)
    unit=unit';
end
 
if size(unit,2)>1
    error('UNIT must be an m by 1 column vector');
end
if length(unit)<2
    error('UNIT must contain at least 2 elements');
end

unit = double(unit);
wall0 = double(wall0);
wall1 = double(wall1);

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

seconds=wall0+(wall1-wall0)*unit;
wall = round(100000*seconds2wall(seconds))/100000;