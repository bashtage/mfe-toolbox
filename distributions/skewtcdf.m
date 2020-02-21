function p = skewtcdf(x,v,lambda)
% Cumulative Distribution Function (CDF) of Hansen's (1994) Skewed T distribution
%
% USAGE:
%   P = skewtcdf(X,V,LAMBDA)
%
% INPUTS:
%   X      - Standardized T random variables
%   V      - Degree of freedom parameters, either scalar or size(X)
%   LAMBDA - Asymmetry Parameter
%
% OUTPUTS:
%   P      - Cumulative distribution evaluated at x
%
% COMMENTS:
%   V>2
%   -1<LAMBDA<1
%   Uses TCDF
%
% REFERENCES:
%   [1] Hansen (1994), Intl.Econ.Rev. (35)
%
% See also SKEWTPDF, SKEWTINV, SKEWTRND, SKEWTLOGLIK, TCDF

% Copyright: Andrew Patton
% a.patton@lse.ac.uk
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 12/31/2001

if nargin~=3
    error('Requires 3 input arguments.')
end

[err, errtext, sizeOut, v, lambda] = iscompatible(2,v,lambda,size(x));

if err
    error(errtext)
end

c = gamma((v+1)/2)./(sqrt(pi*(v-2)).*gamma(v/2));
a = 4*lambda.*c.*((v-2)./(v-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

y1 = (b.*x+a)./(1-lambda).*sqrt(v./(v-2));
y2 = (b.*x+a)./(1+lambda).*sqrt(v./(v-2));

p = (1-lambda).*tcdf(y1,v).*(x<-a./b);			% this method seems to cause error message - work it out later...
p = p + (x>=-a./b).*((1-lambda)/2 + (1+lambda).*(tcdf(y2,v)-0.5));

%Original Nonvector code
%cdf = -999.99*ones(T,1);
%for tt = 1:T;
%   if x(tt)<(-a(tt)/b(tt))
%      cdf(tt) = (1-lambda(tt)).*tdis_cdf(y1(tt),v(tt));
%   else
%      cdf(tt) = ((1-lambda(tt))/2 + (1+lambda(tt)).*(tdis_cdf(y2(tt),v(tt))-0.5));
%   end
%end