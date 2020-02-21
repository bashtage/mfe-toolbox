function mldate = c2mdate(crspdate)
% C2MDATE provides a simple method to convert between CRSP dates provided by WRDS and MATLAB dates.
%
% USAGE:
%   [MLDATE] = c2mdate(CRSPDATE)
%
% INPUTS:
%   CRSPDATE  - A scalar or vector of CRSP dates.
%
% OUTPUTS:
%   MLDATE    - A vector with the same size as CRSPDATE consisting of MATLAB dates.
%
% EXAMPLE:
%   crspdate = [19951028  20090706 20120401]'
%   mldate  = c2mdate(crspdate);
%   datestr(mldate)
%       28-Oct-1995
%       06-Jul-2009
%       01-Apr-2012
%
% COMMENTS:
%   This is provided to make it easy to move between CRSP and MATLAB dates.
%
% See also X2MDATE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 10/17/2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(ischar(crspdate))
    error('CRSPDATE must be numeric')
end

if ndims(crspdate)~=2 &&  min(size(crspdate))~=1
    error('CRSPDATE must be a T by 1 or 1 by T vector');
end

if nargin~=1
    error('1 inputs only');
end

if any(crspdate<18000000)
    error('CRSPDATE appears to be incorrectly formatted.  Only supported format is YYYYMMDD');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yr = floor(crspdate/10000);
mo = floor(mod(crspdate,10000)/100);
dd = mod(crspdate,100);

mldate = datenum(yr,mo,dd);



