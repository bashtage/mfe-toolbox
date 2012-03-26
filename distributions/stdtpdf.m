function y=stdtpdf(x,v)
% Probability Density Function (PDF) for the Standardized T distribution
%
% USAGE:
%   Y = stdtpdf(X,V)
%
% INPUTS:
%   X     - Standardized T random variables
%   V     - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   Y     - Probability density evaluated at x
%
% COMMENTS:
%   V>2
%   Uses TPDF
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Inference'
%
% See also STDTCDF, STDTINV, STDTRND, STDTLOGLIK, TPDF  

% Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 8/1/2005

% v must either be scalar or the same size as x
if nargin~=2
    error('2 inputs required')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(x));

if err
    error(errtext)
end

%stdev=sqrt(v./(v-2));
%stdev(v<=2)=NaN;
%x=x.*stdev;

constant = exp(gammaln((v+1)/2)-gammaln(v/2));
y = constant./sqrt(pi*(v-2)).*(1+x.^2./(v-2)).^(-(v+1)/2);
%y=tpdf(x,v);