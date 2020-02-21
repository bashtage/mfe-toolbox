function x=skewtinv(p,v,lambda)
% Inverse Cumulative Distribution Function (CDF) of Hansen's (1994) 'skewed t'
% distribution; Maps [0,1] to a Skewed-t with V degrees of freedom, shape LAMBDA
%
% USAGE:
%   X=skewtinv(P,V,LAMBDA)
%
% INPUTS:
%   P      - Values to be inverted, P in [0,1]
%   V      - Degree of freedom parameters, either scalar or size(X)
%   LAMBDA - Degree of freedom parameters, either scalar or size(X)
%
% OUTPUTS:
%   X     - Skewed T distributed random variables corresponding to P
%
% COMMENTS:
%   V>2
%   -.99<LAMBDA<.99
%
% REFERENCES:
%   [1] Hansen (1994), Intl.Econ.Rev. (35)
%
% See also SKEWTCDF, SKEWTINV, SKEWTRND, SKEWTLOGLIK

% Copyright: Andrew Patton
% a.patton@lse.ac.uk
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin~=3
    error('3 inputs required')
end

[err, errtext, sizeOut, v, lambda] = iscompatible(2,v,lambda,size(p));

if err
    error(errtext)
end

v(v<=2)=NaN;
lambda(lambda<-.99 | lambda>.99)=NaN;

c = gamma((v+1)/2)./(sqrt(pi*(v-2)).*gamma(v/2));
a = 4*lambda.*c.*((v-2)./(v-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

f1 = find(p<((1-lambda)/2));
f2 = find(p>=((1-lambda)/2));

inv1 = (1-lambda(f1))./b(f1).*sqrt((v(f1)-2)./v(f1)).*tinv(p(f1)./(1-lambda(f1)),v(f1))-a(f1)./b(f1);
inv2 = (1+lambda(f2))./b(f2).*sqrt((v(f2)-2)./v(f2)).*tinv(0.5+1./(1+lambda(f2)).*(p(f2)-(1-lambda(f2))./2),v(f2))-a(f2)./b(f2);
x = repmat(NaN,sizeOut);
x(f1) = inv1;
x(f2) = inv2;