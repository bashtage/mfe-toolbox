function  crspdate = m2cdate(mldate)
% M2CDATE provides a simple method to convert between MATLAB dates and CRSP
% dates provided by WRDS.
%
% USAGE:
%   [CRSPDATE] = c2mdate(MLDATE)
%
% INPUTS:
%   MLDATE    - A vector with the same size as CRSPDATE consisting of MATLAB dates.
%
% OUTPUTS:
%   CRSPDATE  - A scalar or vector of CRSP dates.
%
% EXAMPLE:
%   mldate   = [728960 733960 734960]';
%   datestr(mldate)
%       28-Oct-1995
%       06-Jul-2009
%       01-Apr-2012
%   crspdate = c2mdate(crspdate);
%
% COMMENTS:
%   This is provided to make it easy to move between CRSP and MATLAB dates.
%
% See also C2MDATE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 10/17/2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(ischar(mldate))
    error('CRSPDATE must be numeric')
end

if ndims(mldate)~=2 &&  min(size(mldate))~=1
    error('CRSPDATE must be a T by 1 or 1 by T vector');
end

if nargin~=1
    error('1 inputs only');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=datevec(mldate);
crspdate = 10000*v(:,1) + 100*v(:,2) + v(:,3);



