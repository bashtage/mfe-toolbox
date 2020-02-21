function [simulatedata,ht,pseudorc]=scalar_vt_vech_simulate(t,parameters,c,p,o,q,m)
% Symmetric and Asymmetric SCALAR_VT_VECH(P,O,Q) time series simulation
%
% USAGE:
%   [SIMULATEDATA, HT, PSEUDORC] = scalar_vt_vech_simulate(T, PARAMETERS, C, P, O, Q, M)
%
% INPUTS:
%   T            - Length of the time series to be simulated
%   PARAMETERS   - a P+O+Q x 1 parameter vector
%                    [alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q)]'.
%   C            - K x K matrix containing the intercept
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative, scalar integer representing the number of asymmetric lags
%   Q            - Non-negative, scalar integer representing the number of conditional variance lags 
%   M            - [OPTIONAL] Number of ``intradaily'' returns to simulate to pseudo-Realized Covariance
%
% OUTPUTS:
%   SIMULATEDATA - A time series with scalar variance-targeting vech covariances
%   HT           - A [k k t] matrix of simulated conditional covariances
%   PSEUDORC     - A [k k t] matrix of pseudo-Realized Covariances
%
% COMMENTS:
%   The conditional variance, H(t), of a scalar vech is modeled as follows:
%      H(t) = C + 
%             alpha(1)*r_{t-1}'*r_{t-1} + ... + alpha(p)*r_{t-p}'*r_{t-p}+...
%             gamma(1)*n_{t-1}'*n_{t-1} + ... + gamma(o)*n_{t-p}'*n_{t-p}+...
%             beta(1)*H(t-1) +...+ beta(q)*H(t-q)
%
%       where n_{t-1} = r_{t-1} * (r_{t-1}<0).
%
%   Pseudo Realized Covariances are simulated by generating m-intra daily returns from a N(0,1/m)
%   and computing the Realized Covariance of these. These were used in Patton and Sheppard (2009)
%   when evaluating variance and covariance specifications in a Monte Carlo.  If M=1, then PSEUDORC
%   is just the other product of the SIMULATEDATA.
%
%   NOTE: This program generates 2000 more than required to minimize any start-up bias
%
% See also DIAGONAL_BEKK_SIMULATE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 4/1/2004


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check nargin
switch nargin
    case 6
        m = 72;
    case 7
    % Nothing
    otherwise
    error('6 or 7 arguments required.');
end

%Check t
if ~isscalar(t) && t>0
    error('T must be a positive scalar value.');
end

%p, q much be non-negative scalars
if length(p)>1 || any(p<1) || floor(p)~=p
    error('P must be a positive scalar');
end

if length(o)>1 || any(o<0) || floor(o)~=o
    error('O must be a non-negative scalar');
end

if length(q)>1 || any(q<0) || floor(q)~=q
    error('Q must be a non-negative scalar');
end

%Check c
if ndims(c)~=2 || size(c,1)~=size(c,2) || min(eig(c))<0
    error('C must be a positive definite covariance matrix.')
end

%Check parameters

if any(parameters<0) || length(parameters)<(p+o+q)
    error('All PARAMETERS must be nonnegative and sum to less than 1 (sum(PARAMETERS)<1).');
end
if (sum(parameters(1:p)) + 0.5*sum(parameters(p+1:p+o)) + sum(parameters(p+o+1:p+o+q)))>=1
    warning('MFEToolbox:Stationarity','PARAMETERS are not compatible with covariance stationary process.  Please check HT for problems.')
    initialValue = c./(1-.99);
else
    initialValue = c./(1-(sum(parameters(1:p)) + 0.5*sum(parameters(p+1:p+o)) + sum(parameters(p+o+1:p+o+q))));
end
if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set the burnin amount, 2000 is probably reasonable
burnin=2000;
%Up t by the burnin amount
t=t+burnin;

%Determine K from the size of C
k=length(c);

% Draw some normal random numbers
intraRandomNums=randn(m*t,k)*sqrt(1/m);
randomNums = cumsum(intraRandomNums);
randomNums = diff([zeros(1,k);randomNums(m:m:m*t,:)]);

%Determine the back cast length
n=max([p,o,q]);

%Initialize the random numbers in the back cast range to the unconditional covariance
randomNums(1:n,:)=randomNums(1:n,:)*initialValue^(0.5);

%Initialize the covariance
ht=repmat(initialValue,[1 1 t]);

%Initialize the covariance
pseudorc=repmat(initialValue,[1 1 t]);

%Parse the parameters
alpha=parameters(1:p);
gamma=parameters(p+1:p+o);
beta=parameters(p+o+1:p+o+q);

%Set the data equal to the randomnums to make things easiest
data = randomNums;
eta  = data.*sqrt(0.5);
%Perform the recursion
for i=n+1:t
    ht(:,:,i)=c;
    for j=1:p
        ht(:,:,i)=ht(:,:,i)+alpha(j)*(data(i-j,:)'*data(i-j,:));
    end
    for j=1:o
        ht(:,:,i)=ht(:,:,i)+gamma(j)*(eta(i-j,:)'*eta(i-j,:));
    end
    for j=1:q
        ht(:,:,i)=ht(:,:,i)+beta(j)*ht(:,:,i-j);
    end
    ht(:,:,i) = (ht(:,:,i) + ht(:,:,i)')/2;
    ht12 = ht(:,:,i)^(0.5);
    r = intraRandomNums((i-1)*m+1:i*m,:)*ht12; 
    pseudorc(:,:,i) = r'*r;
    data(i,:)=data(i,:)*ht12;
    eta(i,:)  = data(i,:).*(data(i,:)<0);
end
%Truncate the data and the covariance
simulatedata=data(burnin+1:t,:);
pseudorc = pseudorc(:,:,burnin+1:t);
ht=ht(:,:,burnin+1:t);