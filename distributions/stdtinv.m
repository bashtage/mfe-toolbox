function x=stdtinv(p,v)
% Inverse Cumulative Distribution Function (CDF) of the Standardized T
% distribution; Maps [0,1] to a standardized Students-t with V degrees of freedom
%
% USAGE:
%   X = stdtinv(P,V)
%
% INPUTS:
%   P     - Values to be inverted, P in [0,1]
%   V     - Degree of freedom parameters, either scalar or size(X)
%
% OUTPUTS:
%   X     - Standardized T distributed random variables corresponding to P
%
% COMMENTS:
%   V>2
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Interence'
%
% See also STDTCDF, STDTINV, STDTRND, STDTLOGLIK, TPDF

% Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


if nargin~=2 
    error('2 inputs required')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(p));
if err
    error(errtext)
end

x=tinv(p,v);

stdev=sqrt(v./(v-2));
stdev(v<=2)=NaN;
x=x./stdev;