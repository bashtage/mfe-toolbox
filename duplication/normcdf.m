function p = normcdf(x,mu,sigma)
% Cumulative Distribution Function (CDF) of the Normal Distribution (Gaussian)
%
% USAGE:
%   P = normcdf(X,MU,SIGMA)
%
% INPUTS:
%   X     - Normally distributed random variables
%   MU    - Mean parameter
%   SIGMA - Standard deviation parameter
%
% OUTPUTS:
%   P     - Cumulative distribution evaluated at x
%
% COMMENTS:
%   SIGMA>0
%
%   Uses ERFC
%
% See also NORMPDF, NORMINV, NORMRND, NORMLOGLIK

%%%%%%%%%%%%%%%%%
% Error Checking
%%%%%%%%%%%%%%%%%
if nargin==1
    mu=0;
    sigma=1;
elseif nargin==2
    sigma=1;
elseif nargin<1 || nargin>3
    error('Requires one to three input arguments.')
end

[err, errtext, sizeOut, mu, sigma] = iscompatible(2,mu,sigma,size(x));

if err
    error(errtext);
end
%%%%%%%%%%%%%%%%%
% Error Checking
%%%%%%%%%%%%%%%%%


% Initialize P to zero.
e = (x-mu)./sigma;
% Find any nan values
good = sigma>0;
p=repmat(NaN,size(x));
p(good)=1-erfc(e(good)./sqrt(2))/2;