function [LL,lls]=mvnormloglik(x,mu,sigma)
% Log-likelihood for the multivariate normal distribution
%
% USAGE:
%   [LL,LLS]=mvnormloglik(X,MU,SIGMA)
%
% INPUTS:
%   X      - Normal random variables, T by K
%   MU     - Mean of X, either K by 1 or T by K
%   SIGMA  - Covariance of X, either K by K or K by K by T
%
% OUTPUTS:
%   LL     - Log-likelihood evaluated at X
%   LLS    - Vector of log-likelihoods corresponding to X
%
% COMMENTS:
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Interence'
%
% See also NORMLOGLIK

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 2/1/2007

[T,K]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(x)>2
    error('X must be a T by K matrix');
end

if nargin==1
    sigma = eye(K);
    sigma_is_2D=true;
elseif nargin==2
    %Check mu's conformability
    if ndims(mu)>2 || ~(all(size(mu)==[K 1]) || all(size(mu)==[T K]))
        error('MU must be either a K by 1 vector or the same size as X');
    end
    if size(mu,1)==K
        mu=repmat(mu',T,1);
    end
    x=x-mu;
    sigma=eye(K);
    sigma_is_2D=true;
elseif nargin==3
    if ndims(mu)>2 || ~(all(size(mu)==[K 1]) || all(size(mu)==[T K]))
        error('MU must be either a K by 1 vector or the same size as X');
    end
    if size(mu,1)==K
        mu=repmat(mu',T,1);
    end
    x=x-mu;
    switch ndims(sigma)
        case 2
            % Must by K by K and PD
            if size(sigma,1)~=K || size(sigma,2)~=K || min(eig(sigma))<=0 || any(any(sigma~=sigma'))
                error('SIGMA must by either K by K and positive definite, or K by K by T and contain only positive definite matrices.')
            end
            sigma_is_2D=true;
        case 3
            if size(sigma,1)~=K || size(sigma,2)~=K || size(sigma,3)~=T
                error('SIGMA must by either K by K and positive definite, or K by K by T and contain only positive definite matrices.')
            end
            for t=1:T
                s=sigma(:,:,t);
                if min(eig(s)<=0) || any(any(s~=s'))
                    error('SIGMA must by either K by K and positive definite, or K by K by T and contain only positive definite matrices.')
                end
            end
            sigma_is_2D=false;
        otherwise
            error('SIGMA must by either K by K and positive definite, or K by K by T and contain only positive definite matrices.')
    end
else
    error('Only 1 to 3 inputs supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the log likelihood
if sigma_is_2D
    stdx=x*sigma^(-0.5);
    LL  = -0.5 * (T*K*log(2*pi) + T*log(det(sigma)) + sum(sum(stdx.^2)));
    %Compute the individual log likelihoods if needed
    if nargout>1
        lls = -0.5*(K*log(2*pi) + log(det(sigma)) + sum(stdx.^2,2));
    end
else
    lls=zeros(T,1);
    for t=1:T
        s=sigma(:,:,t);
        lls(t) = -0.5*(K*log(2*pi) + log(det(s)) + x(t,:)*inv(s)*x(t,:)');
    end
    LL=sum(lls);
end

