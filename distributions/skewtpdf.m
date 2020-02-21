function y = skewtpdf(x, v, lambda)
% Probability Density Function (PDF) of Hansen's (1994) Skewed T distribution
%
% USAGE:
%   Y = skewtpdf(X,V,LAMBDA)
%
% INPUTS:
%   X      - Standardized T random variables
%   V      - Degree of freedom parameters, either scalar or size(x)
%   LAMBDA - Asymmetry Parameter
%
% OUTPUTS:
%   Y     - Probability density evaluated at X
%
% COMMENTS:
%   V>2
%   -1<LAMBDA<1
%   Uses TCDF
%
% REFERENCES:
%   [1] Hansen (1994), Intl.Econ.Rev. (35)
%
% See also SKEWTCDF, SKEWTINV, SKEWTRND, SKEWTLOGLIK, TCDF

% Copyright: Andrew Patton
% a.patton@lse.ac.uk
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 12/31/2001

if nargin~=3
    error('Requires 3 input arguments.')
end

v(v<2)=NaN;
lambda(lambda<-1 | lambda>1)=NaN;

[err, errtext, sizeOut, v, lambda] = iscompatible(2,v,lambda,size(x));

if err
    error(errtext)
end

c = gamma((v+1)/2)./(sqrt(pi*(v-2)).*gamma(v/2));
a = 4*lambda.*c.*((v-2)./(v-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

y1 = b.*c.*(1 + 1./(v-2).*((b.*x+a)./(1-lambda)).^2).^(-(v+1)/2);
y2 = b.*c.*(1 + 1./(v-2).*((b.*x+a)./(1+lambda)).^2).^(-(v+1)/2);
y  = y1.*(x<(-a./b)) + y2.*(x>=(-a./b));