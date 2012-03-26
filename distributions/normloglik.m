function [LL,lls]=normloglik(x,mu,sigma2)
% Log-likelihood for the normal distribution
%
% USAGE:
%   [LL,LLS]=normloglik(X,MU,SIGMA2)
%
% INPUTS:
%   X      - Normal random variables, either scalar or column vector
%   MU     - Mean of X, either scalar or size(x)
%   SIGMA2 - Variance of X, either scalar or size(x)
%
% OUTPUTS:
%   LL    - Log-likelihood evaluated at X
%   LLS   - Vector of log-likelihoods corresponding to X
%
% COMMENTS:
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Interence'
%
% See also NORMCDF, NORMINV, NORMRND, NORMPDF

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2005

[T,K]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K~=1
    error('x must be a column vector');
end

if nargin==1
    sigma2=ones(T,K);
elseif nargin==2
    %Check mu's conformability
    if length(mu)~=1 && ~all(size(mu)==[T K])
        error('mu must be either a scalar or the same size as X');
    end
    x=x-mu;
    sigma2=ones(T,K);
elseif nargin==3
    if length(mu)~=1 && ~all(size(mu)==[T K])
        error('mu must be either a scalar or the same size as X');
    end
    if any(sigma2<=0)
        error('sigma2 must contain only positive elements')
    end
    if length(sigma2)==1
        sigma2=sigma2*ones(T,K);
    elseif size(sigma2,1)~=T || size(sigma2,2)~=1
        error('sigma2 must be a scalar or a vector with the same dimensions as x');
    end
    x=x-mu;
else
    error('Only 1 to 3 inputs supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargout>1
%Compute the individual log likelihoods if needed
    lls = -0.5*(log(2*pi) + log(sigma2) + x.^2./sigma2);
    % Use these to comput the LL
    LL = sum(lls);
else
    %Compute the log likelihood
    LL  =  -0.5 * (sum(log(sigma2)) + sum((x.^2)./sigma2)  +  T*log(2*pi));
end
