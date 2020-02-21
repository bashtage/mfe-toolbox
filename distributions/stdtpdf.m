function y=stdtpdf(x,mu,sigma2,nu)
% Probability Density Function (PDF) for the Standardized T distribution
%
% USAGE:
%   Y = stdtpdf(X,MU,SIGMA2,NU)
%
% INPUTS:
%   X      - Standardized T random variables
%   MU     - Mean of X, either scalar or size(x) 
%   SIGMA2 - Variance of X, either scalar or size(x)
%   NU     - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   Y     - Probability density evaluated at X
%
% COMMENTS:
%   NU>2
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Inference'
%
% See also STDTCDF, STDTINV, STDTRND, STDTLOGLIK, TPDF  

% Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 6    Date: 8/21/2014

[T,K]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K~=1
    error('X must be a column vector');
end

if nargin==4
    if length(mu)~=1 && ~all(size(mu)==[T K])
        error('mu must be either a scalar or the same size as X');
    end
    if any(sigma2<=0)
        error('sigma2 must contain only positive elements')
    end
    if length(sigma2)==1
        sigma2=sigma2*ones(T,K);
    elseif size(sigma2,1)~=T || size(sigma2,2)~=1
        error('sigma2 must be a scalar or a vector with the same dimensions as X');
    end
    if length(nu)>1 || nu<=2
        error('nu must be a scalar greater than 2');
    end
    x=x-mu;
else
    error('Only 4 inputs supported');
end


constant = exp(gammaln( 0.5 * (nu + 1)) - gammaln(0.5 * nu));
y = constant ./ sqrt(pi * (nu - 2) * sigma2) .* (1 + (x-mu) .^ 2.0 / (sigma2 * (nu - 2))) .^ (-(nu + 1) / 2);
