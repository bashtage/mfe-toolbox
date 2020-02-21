function x = norminv(p,mu,sigma)
% Inverse Cumulative Distribution Function (CDF) of the Normal Distribution (Gaussian)
%
% USAGE:
%   X = normcdf(P,MU,SIGMA)
%
% INPUTS:
%   P     - Cumulative distribution values
%   MU    - Mean parameter
%   SIGMA - Standard deviation parameter
%
% OUTPUTS:
%   X     - Normally distributed random variables
%
% COMMENTS:
%   SIGMA>0
%
%   Uses ERFINV
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

[err, errtext, sizeOut, mu, sigma] = iscompatible(2,mu,sigma,size(p));

if err
    error(errtext);
end
%%%%%%%%%%%%%%%%%
% Error Checking
%%%%%%%%%%%%%%%%%

% Find any nan values
good = sigma>0 & p>=0 & p<=1;
x=repmat(NaN,size(p));
x(p==.5)=0;
x(good & p>.5) = erfinv(2*(p(good & p>.5)-0.5))*sqrt(2);
x(good & p<.5) = -erfinv(2*(-p(good & p<.5)+0.5))*sqrt(2);
x=x.*sigma+mu;