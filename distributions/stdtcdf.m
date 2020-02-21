function p=stdtcdf(x,v)
% Cumulative Distribution Function (CDF) of the Standardized T distribution
%
% USAGE:
%   P = stdtcdf(X,V)
%
% INPUTS:
%   X     - Standardized T random variables
%   V     - Degree of freedom parameters, either scalar or size(X)
%
% OUTPUTS:
%   P     - Cumulative distribution evaluated at X
%
% COMMENTS:
%   V>2
%   Uses TCDF
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Inference'
%
% See also STDTPDF, STDTINV, STDTRND, STDTLOGLIK, TCDF

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
stdev=sqrt(v./(v-2));
stdev(v<=2)=NaN;
x=x.*stdev;

p=tcdf(x,v);
