function y = gedpdf(x,v)
% Probability Density Function (PDF) for the Generalized Error Distribution (GED)
%
% USAGE:
%   Y = gedpdf(X,V)
%
% INPUTS:
%   X     - Generalized Error Distribution distributed random variables
%   V     - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   Y     - Probability density evaluated at x
%
% COMMENTS:
%   V>=1
%   A scalar GED r.v. with variance normalized to 1 has probability
%   density given by:
%       f(x,V) = [V/(lda*2^(1+1/V)*gamma(1/V))]*exp(-0.5*|x/lda|^V)
%       lda = [2^(-2/V)*gamma(1/V)/gamma(3/V)]^0.5
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
% Revision: 3    Date: 8/1/2005

if nargin~=2
    error('Requires two input arguments.')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,size(x));

if err
    error(errtext)
end

% Initialize D to zero.
y = zeros(size(x));

%   Return NaN if the argument V is outside its limit.
y(v < 1) = NaN;

k = find(v >= 1);
vk = v(k); xk = x(k);
if any(k)
    ldak = ((2.^(-2./vk)).*(gamma(1./vk))./(gamma(3./vk))).^(0.5);
    y(k) = vk.*exp(-0.5*((abs(xk./ldak)).^vk))./(ldak.*gamma(1./vk).*(2.^(1+1./vk)));
end