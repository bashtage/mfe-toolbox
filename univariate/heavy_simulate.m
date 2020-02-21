function [data,ht] = heavy_simulate(T,K,parameters,p,q,m,R)
% Simulation of HEAVY volatility model of Shephard and Sheppard.  Also simulates general volatility
% spillover models with more than 2 dimensions.
%
% USAGE:
%  [DATA,HT] = heavy_simulate(T,K,PARAMETERS,P,Q,M,R)
%
% INPUTS:
%   T          - Either a scalar containing the length of the series to simulate, or a T by K matrix
%                  of simulated random variables.  The default is to use standard normal random
%                  variables.  Providing a T by K matrix allows other distributions to be used.
%   K          - Number of series to simulate
%   PARAMETERS - A vector with K+sum(sum(P))+sum(sum(Q)) elements. See COMMENTS.
%   P          - A K by K matrix containing the lag length of model innovations.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   Q          - A K by K matrix containing the lag length of conditional variances.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   M          - [OPTIONAL] K by 1 vector containing the number of samples to use when computing
%                  realized-like data.  Default is ones(K,1).  When M is 1, the returned data is
%                  like a return, so that data(t,i) has mean 0 and variance h(t,i).  When M>1, the
%                  returned data is like a realized measure, so that the mean of data(t,i) is h(t,i).
%   R          - [OPTIONAL] Correlation matrix for Gaussian copula of the underlying data.
%                  Only used if T is scalar.  Default is eye(K).
%
% OUTPUTS:
%   DATA   - A T by K matrix of simulated data
%   HT     - A T by K matrix of conditional variances
%
% COMMENTS:
%   Dynamics are given by:
%     h(t,:)' = O + A(:,:,1)*f(data(t-1,:))' + ... + A(:,:,maxP)*f(data(t-maxP,:))' + ...
%                 + B(:,:,1)*h(t-1,:)' + ... + B(:,:,maxQ)*h(t-maxQ,:)'
%
%   PARAMETERS are ordered:
%   [O' A(1,1,1:p(1,1)) A(1,2,1:p(1,2)) ... A(1,K,1:p(1,K)) A(2,1,1:p(2,1)) ... A(2,K,1:p(2,K)) 
%       ... A(K,1,1:p(K,1)) ... A(K,K,1:p(K,K)) ... B(1,1,1:q(1,1)) ... B(1,K,1:q(1,K)) 
%       ... B(K,1,1:q(K,1)) ... B(K,K,1:q(K,K)) ]
%
% EXAMPLES:
%   % Standard HEAVY simulation
%   parameters = [.15 .05 .2 .4 .7 .55]
%   p = [0 1;0 1]
%   q = eye(2)
%   m = [1 78];
%   [data,ht] = heavy_simulate(1000,2,parameters,p,q,m)
%
% See also HEAVY

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/2/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    error('5 to 7 inputs only')
end

if ~isscalar(T)
    e = T;
    [T,K] = size(e);
    R = eye(K);
else
    e = [];
end

switch nargin
    case 5
        m = ones(K,1);
        R = eye(K);
    case 6
        R = eye(K);
    case 7
        % Nothing
    otherwise
        error('5 to 7 inputs only')
end

if any(size(p)~=K) || any(any(floor(p)~=p)) || any(any(p<0))
    error('P must be a K by K matrix of non-negative integers.')
end
if any(size(q)~=K) || any(any(floor(q)~=q)) || any(any(q<0))
    error('Q must be a K by K matrix of non-negative integers.')
end
% Parameters
count = K+sum(sum(p))+sum(sum(q));
if length(parameters)~=count
    error('PARAMETERS has the wrong number of inputs')
end
% M
if length(m)<K || any(floor(m)~=m)  || any(m<1)
    error('M must be a K by 1 vector of positive integers.')
end
% R
if all(size(R)~=K) || any(any(R~=R')) || min(eig(R))<=0
    error('R must be a K by K positive definite matrix')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if T<=500
    tau = 2*T;
else
    tau = T+500;
end

if isempty(e)
    u = R^(0.5)*randn(K,tau);
    e = u;
    u = normcdf(e);
    signs = zeros(K,tau);
    for i=1:K
        if m(i)~=1
            e(i,:) = chi2inv(u(i,:),m(i))/m(i);
        else
            signs(i,:) = 2*(e(i,:)>0)-1;
            e(i,:) = e(i,:).^2;
        end
    end
end

[O,A,B] = heavy_parameter_transform(parameters,p,q,K);
% Determine if stationary
pq = max([p(:)',q(:)']);
pMax = max(p(:));
qMax = max(q(:));

if pq==1
    companion = A+B;
else
    companion = zeros(K*pq,K*pq);
    companion(K+1:K*pq,1:K*(pq-1)) = eye(K*(pq-1));
    companion(1:K,1:K*pMax) = reshape(A,[K  K*pMax]);
    companion(1:K,1:K*qMax) = companion(1:K,1:K*qMax) + reshape(B,[K  K*qMax]);
end

lambda = eig(companion);
if max(abs(lambda))<1
    uncond = ((eye(K)-sum(A,3)-sum(B,3))\eye(K))*O;
else
    uncond = O/.001;
    warning('MFE:nonstationary','The HEAVY process is non-stationary');
end
data = zeros(K,tau);
ht = ones(K,tau);
backCast = uncond;

for t=1:tau
    ht(:,t)= O;
    for j=1:pMax
        if (t-j)<=0
            ht(:,t) = ht(:,t) + A(:,:,j)*backCast;
        else
            ht(:,t) = ht(:,t) + A(:,:,j)*data(:,t-j);
        end
    end
    for j=1:qMax
        if (t-j)<=0
            ht(:,t) = ht(:,t) + B(:,:,j)*backCast;
        else
            ht(:,t) = ht(:,t) + B(:,:,j)*ht(:,t-j);
        end
    end
    data(:,t) = e(:,t).*ht(:,t);
end

data(m==1,:) = signs(m==1,:).*sqrt(data(m==1,:));
ht = ht';
data= data';

data = data(tau-T+1:tau,:);
ht = ht(tau-T+1:tau,:);