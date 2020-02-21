function p = gedcdf(x,v)
% Cumulative Distribution Function (CDF) of the Generalized Error Distribution(GED)
%
% USAGE:
%   P = gedcdf(X,V)
%
% INPUTS:
%   X     - Generalized Error Distribution distributed random variables
%   V     - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   P     - Cumulative distribution evaluated at x
%
% COMMENTS:
%   V>1
%   A scalar GED r.v. with variance normalized to 1 has probability
%   density given by:
%       f(x,V) = [V/(lda*2^(1+1/V)*gamma(1/V))]*exp(-0.5*|x/lda|^V)
%       lda = [2^(-2/V)*gamma(1/V)/gamma(3/V)]^0.5
%   Uses GAMCDF
%
% REFERENCES:
%   [1] Tadikamalla (1980), J.Am.Stat.Assoc. (75)
%   [2] Nelson (1991), Econometrica
%
% See also GEDPDF, GEDINV, GEDRND, GEDLOGLIK GAMCDF

% Copyright:
% Ivana Komunjer
% komunjer@hss.caltech.edu
% Revision: 1    Date: 23/08/2002
% Modifications Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 8/1/2005
if nargin~=2
    error('Requires two input arguments.')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(x));

if err
    error(errtext)
end

% Initialize P to zero.
p = zeros(size(x));
scalex = (gamma(3./v)./gamma(1./v)).^0.5;

%   Return NaN if the argument V is outside its limits.
p(v < 1) = NaN;

ku = find(x >= 0 & ~(v < 1));
vku = v(ku);
xku = x(ku);
scalexku = scalex(ku);
if any(ku),
    p(ku) = 0.5*(1 + gamcdf((abs(xku.*scalexku)).^vku,1./vku));
end

kl = find(x < 0 & ~(v < 1));
vkl = v(kl);
xkl = x(kl);
scalexkl = scalex(kl);
if any(kl),
    p(kl) = 0.5*(1 - gamcdf((abs(xkl.*scalexkl)).^vkl,1./vkl));
end

% Make sure that round-off errors never make P less than 0 or greater than 1.
p(p < 0) = 0;
p(p > 1) = 1;