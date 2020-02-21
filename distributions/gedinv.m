function x = gedinv(p,v)
% Inverse Cumulative Distribution Function (CDF) of the Generalized Error
% Distribution (GED); Maps [0,1] to a GED with shape parameter V
%
% USAGE:
%   X = gedinv(P,V)
%
% INPUTS:
%   P     - Values to be inverted, P in [0,1]
%   V     - Shape parameters, either scalar or size(X)
%
% OUTPUTS:
%   X     - GED random variables corresponding to P
%
% COMMENTS:
%   A scalar GED r.v. with variance normalized to 1 has probability
%   density given by:
%       f(x,V) = [V/(lda*2^(1+1/V)*gamma(1/V))]*exp(-0.5*|x/lda|^V)
%       lda = [2^(-2/V)*gamma(1/V)/gamma(3/V)]^0.5
%   GAMINV does the computational work
%
% REFERENCES:
%   [1] Tadikamalla (1980), J.Am.Stat.Assoc. (75)
%   [2] Nelson (1991), Econometrica
%
% See also GEDCDF, GEDINV, GEDRND, GEDLOGLIK

% Copyright: Ivana Komunjer
% komunjer@hss.caltech.edu
% Revision: 1    Date: 23/08/2002
% Modifications Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin~=2
    error('Requires two input arguements.')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(p));

if err
    error(errtext)
end

%   Initialize X to zero.
x = zeros(size(p));
scalex = (gamma(3./v)./gamma(1./v)).^0.5;

k = find(p<0 | p>1 | v < 1);
if any(k),
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.
k0 = find(p == 0 & v >= 1);
if any(k0),
    x(k0) = zeros(size(k0));
end

k1 = find(p == 1 & v >= 1);
if any(k1),
    tmp = Inf;
    x(k1) = tmp(ones(size(k1)));
end

% Tadikamalla's Method: uses the fact that Y=abs(X)^v has
% a Gamma distribution with shape parameter 1/v
kl = find(p > 0  &  p < 0.5 & v >= 1);
if any(kl(:))
    pkl = p(kl); 
    vkl = v(kl);
    xkl = gaminv(1-2*pkl,1./vkl);
    xkl = (-1)*xkl.^(1./vkl);
    x(kl) = xkl;
end

ku = find(p >= 0.5  &  p < 1 & v >= 1);
if any(ku(:))
    pku = p(ku); 
    vku = v(ku);
    xku = gaminv(2*pku-1,1./vku); 
    xku = xku.^(1./vku);
    x(ku) = xku;
end

% need to standardize the quantiles
x = x./scalex;

