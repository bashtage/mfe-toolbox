function unit = seconds2unit(seconds,seconds0,seconds1)
% Convert seconds past midnight to unit times
%
% USAGE:
%   [UNIT] = seconds2unit(SECONDS,SECONDS0,SECONDS1)
%
% INPUTS:
%   SECONDS - m by 1 column vector of times expressed as seconds past midnight
%               (e.g. 1:00:00 is 3600, 12:00:15 is 43215) where m>=2
%   SECONDS0   - Base time, maps to 0 un the unit interval
%   SECONDS1   - End time, maps to 1 un the unit interval
%
% OUTPUTS:
%   UNIT    - m by 1 column vector of times measured as fraction of time
%               between min(SECONDS) and max(SECONDS)
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
if any(seconds<0) 
    error('SECONDS must be non-negative');
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

if size(seconds,2)>size(seconds,1)
    seconds=seconds';
end

seconds = double(seconds);
seconds0 = double(seconds0);
seconds1 = double(seconds1);

if size(seconds,2)>1
    error('SECONDS must be an m by 1 column vector');
end
if length(seconds)<2
    error('SECONDS must contain at least 2 elements');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit=(seconds-seconds0)/(seconds1 - seconds0);