function [ll,lls]=stdtloglik(x,mu,sigma2,nu)
% Log likelihood of the Standardized T distribution
%
% USAGE:
%   [LL,LLS]=stdtloglik(X,MU,SIGMA2,NU)
%
% INPUTS:
%   X      - Standardized T random variables, either scalar or column vector
%   MU     - Mean of X, either scalar or size(x) 
%   SIGMA2 - Variance of X, either scalar or size(x)
%   V      - Degree of freedom parameters, either scalar or size(x)
%
% OUTPUTS:
%   LL    - Log-likelihood evaluated at X
%   LLS   - Vector of log-likelihoods corresponding to X
%
% COMMENTS:
%   V>2
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Inference'
%
% See also STDTCDF, STDTINV, STDTRND, STDTPDF  

% Copyright: 
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2004

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the log likelihood
ll = T*gammaln(0.5*(nu+1)) - T*gammaln(nu/2) - T/2*log(pi*(nu-2)) ...
    - 0.5*sum(log(sigma2)) - ((nu+1)/2)*sum(log(1 + (x.^2)./(sigma2*(nu-2))));

%Compute the individual log likelihoods if needed
if nargout>1
    lls = gammaln(0.5*(nu+1)) - gammaln(nu/2) - 1/2*log(pi*(nu-2))...
        - 0.5*(log(sigma2)) - ((nu+1)/2)*(log(1 + (x.^2)./(sigma2*(nu-2))));
end
