function p=chi2cdf(x,v)
% Cumulative Distribution Function (CDF) of the Chi-square distribution
% Work-a-round for the Stat's Toolbox chi2cdf
%
% USAGE:
%   P = chi2cdf(X,V)
%
% INPUTS:
%   X     - Chi-square random variables
%   V     - Degree of freedom parameters, either scalar or size(X)
%
% OUTPUTS:
%   P     - Cumulative distribution evaluated at X
%
% COMMENTS:
%   V>0
%
% REFERENCES:
%

% Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 8/1/2005

if nargin~=2
    error('2 inputs required')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(x));
if err
    error(errtext)
end

%Compute the st
p=repmat(NaN,sizeOut);
pl = v>0;
temp=gammainc(x(pl)/2,v(pl)*0.5);
p(pl)=temp;
