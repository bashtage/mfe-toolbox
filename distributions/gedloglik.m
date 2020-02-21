function [LL,LLS]=gedloglik(x,mu,sigma2,v)
% Log likelihood for the Generalized Error Distribution (GED) distribution
%
% USAGE:
%   [LL,LLS]=gedloglik(X,MU,SIGMA2,V)
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
%   V>1
%   f(x,V) = [v/(lda*2^(1+1/v)*gamma(1/v))]*exp(-0.5*|x/lda|^v)
%      lda = [2^(-2/v)*gamma(1/v)/gamma(3/v)]^0.5
%
% REFERENCES:
%   [1] Tadikamalla (1980), J.Am.Stat.Assoc. (75)
%   [2] Nelson (1991), Econometrica
%
% See also GEDCDF, GEDINV, GEDRND, GEDPDF

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2005

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
    if length(v)>1 || v<=1
        error('V must be a scalar greater than 1');
    end
    x=x-mu;
else
    error('Only 4 inputs supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the log likelihood
logl=0.5*(-2/v*log(2)+gammaln(1/v)-gammaln(3/v));
l=exp(logl);
LL = (T * log(v)) - (T*logl) - (T*gammaln(1/v)) - T*(1+1/v)*log(2);
LL = LL - 0.5 * sum(log(sigma2)) - 0.5 * sum((abs(x./(sqrt(sigma2)*l))).^v);

%Compute the individual LLS if needed
if nargout>1
    LLS = log(v) - logl - gammaln(1/v) - (1+1/v)*log(2) - 0.5 * log(sigma2) - 0.5 * (abs(x./(sqrt(sigma2)*l))).^v;
    %    LLS=log(v)-logl-0.5*abs(x./l).^v-logl-(1+1/v)*log(2)-gammaln(1/v);
end
