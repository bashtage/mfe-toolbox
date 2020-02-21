function [simulatedata, ht] = agarch_simulate(t,parameters,p,q,model_type,error_type)
% AGARCH(P,Q) and NAGARCH(P,Q) time series simulation with multiple error distributions
%
% USAGE:
%   [SIMULATEDATA, HT] = agarch_simulate(T,PARAMETERS,P,Q,MODEL_TYPE,ERROR_TYPE)
%
% INPUTS:
%   T             - Length of the time series to be simulated  OR
%                   T by 1 vector of user supplied random numbers (i.e. randn(1000,1))
%   PARAMETERS    - a 2+P+Q (+1 or 2, depending on error distribution) x 1 parameter vector
%                   [omega alpha(1) ... alpha(p) gamma beta(1) ... beta(q) [nu lambda]]'.
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance
%   MODEL_TYPE    - [OPTIONAL] The type of variance process, either
%                     'AGARCH'  - Asymmetric GARCH, Engle (1990) [DEFAULT]
%                     'NAGARCH' - Nonlinear Asymmetric GARCH, Engle & Ng (1993)
%   ERROR_TYPE    - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%
% OUTPUTS:
%   SIMULATEDATA  - A time series with AGARCH or NAGARCH variances
%   HT            - A vector of conditional variances used in making the time series
%
% COMMENTS:
%    The conditional variance, h(t), of a AGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma)^2 + ... + alpha(p)*(r_{t-p}-gamma)^2
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%    The conditional variance, h(t), of a NAGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma*sqrt(h(t-1)))^2 + ... + alpha(p)*(r_{t-p}-gamma*sqrt(h(t-p)))^2 
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%   NOTE: This program generates 2000 more than required to minimize any starting bias
%
%  See also AGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==4
    error_type='NORMAL';
    model_type='AGARCH';
elseif nargin==5
    error_type='NORMAL';
elseif nargin==6
    %nothing
else
    error('4 to 6 inputs only.');
end

if isempty(error_type)
    error_type='NORMAL';
end

if strcmp(error_type,'NORMAL')
    extrap=0;
elseif strcmp(error_type,'STUDENTST')
    extrap=1;
elseif strcmp(error_type,'GED')
    extrap=1;
elseif strcmp(error_type,'SKEWT')
    extrap=2;
else
    error('Unknown error type')
end

if ~ismember(model_type,{'AGARCH','NAGARCH'})
    error('MODEL_TYPE must be either ''AGARCH'' or ''NAGARCH''')
end

if size(parameters,1)~=(2+p+q+extrap) || size(parameters,2)>1
    error('PARAMETERS must be a column vector with the correct number of parameters.');
end

if (sum(parameters(2:p+1)) + sum(parameters(p+3:p+q+2))) >=1
    warning('UCSD_GARCH:nonstationary','PARAMETERS are in the non-stationary space, be sure to check that H is not inf.');
end

if length(p)>1 || length(q)>1 || any(q<0) || any(p<1)
    error('P and Q must be scalars with P positive and Q non-negative')
end

if size(t,2)~=1
    error('T must be either a positive scalar or a vector of random numbers.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Separate the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
gamma=parameters(p+2);
beta=parameters(p+3:p+q+2);

%Initialize the random numbers
if isscalar(t)
    t=t+2000;
    if strcmp(error_type,'NORMAL')
        RandomNums=randn(t,1);
    elseif strcmp(error_type,'STUDENTST')
        nu=parameters(3+p+q);
        RandomNums=stdtrnd(nu,t,1);
    elseif strcmp(error_type,'GED')
        nu=parameters(3+p+q);
        RandomNums=gedrnd(nu,t,1);
    elseif strcmp(error_type,'SKEWT')
        nu=parameters(3+p+q);
        lambda=parameters(4+p+q);
        RandomNums=skewtrnd(nu,lambda,t,1);
    else
        error('Unknown error type')
    end
else
    RandomNums=t;
    t=length(RandomNums);
    seeds=ceil(rand(2000,1)*t);
    RandomNums=[RandomNums(seeds);RandomNums];
    t=length(RandomNums);
end

m  =  max([p,q]);
%Back casts, zeros are fine since we are throwing away over 2000
RandomNums=[zeros(m,1);RandomNums];

if (1-sum(alpha)-sum(beta))>0
    % Approximate, shoudln't matter much
    UncondStd =  sqrt(omega/(1-sum(alpha)-sum(beta)));
else %non stationary
    UncondStd=1;
end

h=UncondStd.^2*ones(t+m,1);
data=UncondStd*ones(t+m,1);
T=size(data,1);
shock=zeros(t+m,1);


parameters=[omega;alpha;beta];
if strcmp(model_type,'AGARCH')
    for t=1:m
        shock(t) = (UncondStd-gamma)^2;
    end
    for t = (m + 1):T
        h(t)     = parameters' * [1 ; shock(t-(1:p)); h(t-(1:q)) ];
        data(t)  = RandomNums(t)*sqrt(h(t));
        shock(t) = (data(t)-gamma)^2;
    end
else
    for t=1:m
        shock(t) = (UncondStd-gamma*UncondStd)^2;
    end    
    for t = (m + 1):T
        h(t)     = parameters' * [1 ; shock(t-(1:p)); h(t-(1:q)) ];
        data(t) = RandomNums(t)*sqrt(h(t));
        shock(t) = (data(t)-gamma*sqrt(h(t)))^2;
    end
end
simulatedata=data((m+1+2000):T);
ht=h(m+1+2000:T);