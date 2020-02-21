function mldate = x2mdate(xlsdate, type)
% X2MDATE provides a simple method to convert between excel dates and MATLAB dates.  
%
% USAGE:
%   [MLDATE] = x2mdate(XLSDATE)
%   [MLDATE] = x2mdate(XLSDATE, TYPE)
%
% INPUTS:
%   XLSDATE   - A scalar or vector of Excel dates. 
%   TYPE      - [OPTIONAL] A scalar or vector of the same size as XLSDATE that describes the Excel
%                 basedate.  Can be either 0 or 1.  If 0 (default), the base date of Dec-31-1899 is
%                 used.  If 1, the base date is Jan 1, 1904.
%
% OUTPUTS:
%   MLDATE    - A vector with the same size as XLSDATE consisting of MATLAB dates.
%
% EXAMPLE:
%   XLSDATE = [35000 40000 41000];
%   MLDATE  = x2mdate(XLSDATE);
%   datestr(MLDATE)
%       28-Oct-1995
%       06-Jul-2009
%       01-Apr-2012
%
% COMMENTS:
%   This is a reverse engineered clone of the MATLAB function x2mdate and should behave the same.
%   You only need it if you do not have the financial toolbox installed. 
%
% See also C2MDATE


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 10/27/2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==1
    type=0;
end

if any(ischar(xlsdate))
    error('XLSDATE must be numeric')
end

if ndims('xlsdate')~=2 &&  min(size(type))~=1
    error('XLSDATE must be a T by 1 or 1 by T vector');
end

if nargin>1
    if nargin>2
        error('1 or 2 inputs only');
    end
    if min(size(type))~=1 || ndims(type)~=2
        error('TYPE must be either a scalar or a vector conformable to xlsdate');
    end
    if max(size(type))~=1
        if max(size(type))~=max(size(xlsdate))
            error('TYPE must be either a scalar or a vector conformable to xlsdate');
        end
    end
    if any(~ismember(type,[0 1]))
        error('TYPE must be either 0 or 1.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isscalar(type)
    type=type*ones(size(xlsdate));
end
type = logical(type);

mldate=zeros(size(xlsdate));
mldate(~type) = xlsdate(~type) + datenum('30-Dec-1899') ;
mldate(type) = xlsdate(type) + datenum('1-Jan-1904') ;
