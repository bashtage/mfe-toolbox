function y = normpdf(x,mu,sigma)
% Probability Density Function (PDF) of the Normal Distribution (Gaussian)
%
% USAGE:
%   Y = normpdf(X,MU,SIGMA)
%
% INPUTS:
%   X     - Normally distributed random variables
%   MU    - Mean parameter
%   SIGMA - Standard deviation parameter
%
% OUTPUTS:
%   Y     - Probability density evaluated at x
%
% COMMENTS:
%   SIGMA>0
%
%   This is a clean room implementation that mimics the core function of
%   MATLAB's normpdf.  
%
%
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


% Initialize y to nan.
% Find any nan values
good = sigma>0;
y=repmat(NaN,size(x));
y(good)=(1./(sqrt(2*pi)*sigma)).*exp(-(x(good)-mu(good)).^2./(2*sigma(good).^2));
