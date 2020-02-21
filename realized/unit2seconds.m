function seconds = unit2seconds(unit, seconds0, seconds1)
% Convert unit times to seconds past midnight 
%
% USAGE:
%   [SECONDS] = seconds2unit(UNIT,SECONDS0,SECONDS1)
%
% INPUTS:
%   UNIT       - m by 1 column vector of times measured as fraction of time
%                 between min(SECONDS) and max(SECONDS)
%   SECONDS0   - Base time, maps to 0 un the unit interval
%   SECONDS1   - End time, maps to 1 un the unit interval
%
% OUTPUTS:
%   SECONDS - m by 1 column vector of times expressed as seconds past midnight
%               (e.g. 1:00:00 is 3600, 12:00:15 is 43215) where m>=2
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also SECONDS2UNIT, WALL2UNIT, WALL2SECONDS, SECONDS2WALL, REALIZED_KERNEL, REALIZED_VARIANCE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('Three inputs required.')
end
if seconds0<0
    error('SECONDS0 must be non-negative');
end
if seconds1<0
    error('SECONDS1 must be non-negative');
end

if (seconds1-seconds0)<eps 
    error('SECONDS1 must be larger than SECONDS0.');
end


if size(unit,2)>size(unit,1)
    unit=unit';
end
 
unit = double(unit);
seconds0 = double(seconds0);
seconds1 = double(seconds1);

if size(unit,2)>1
    error('SECONDS must be an m by 1 column vector');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



seconds=seconds0+(seconds1-seconds0)*unit;
