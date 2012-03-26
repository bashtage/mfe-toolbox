function [simulatedata, ht] = tarch_simulate(t,parameters,p,o,q,error_type,tarch_type)
% TARCH(P,O,Q) time series simulation with multiple error distributions
%
% USAGE:
%   [SIMULATEDATA, HT] = tarch_simulate(T, PARAMETERS, P, O, Q, ERROR_TYPE, TARCH_TYPE)
%
% INPUTS:
%   T             - Length of the time series to be simulated  OR
%                   T by 1 vector of user supplied random numbers (i.e. randn(1000,1))
%   PARAMETERS    - a 1+P+O+Q (+1 or 2, depending on error distribution) x 1 parameter vector
%                   [omega alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q) [nu lambda]]'.
%   P             - Positive, scalar integer representing the number of symmetric innovations
%   O             - Non-negative scalar integer representing the number of asymmetric innovations (0
%                     for symmetric processes) 
%   Q             - Non-negative, scalar integer representing the number of lags of conditional
%                     variance (0 for ARCH) 
%   ERROR_TYPE    - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   TARCH_TYPE    - [OPTIONAL] The type of variance process, either
%                     1 - Model evolves in absolute values
%                     2 - Model evolves in squares [DEFAULT]
%
% OUTPUTS:
%   SIMULATEDATA  - A time series with ARCH/GARCH/GJR/TARCH variances
%   HT            - A vector of conditional variances used in making the time series
%
% COMMENTS:
% The conditional variance, h(t), of a TARCH(P,O,Q) process is modeled as follows:
%     g(h(t)) = omega
%             + alpha(1)*f(r_{t-1}) + ... + alpha(p)*f(r_{t-p})
%             + gamma(1)*I(t-1)*f(r_{t-1}) +...+ gamma(o)*I(t-o)*f(r_{t-o})
%             + beta(1)*g(h(t-1)) +...+ beta(q)*g(h(t-q))
%
%     where f(x) = abs(x)  if tarch_type=1
%           g(x) = sqrt(x) if tarch_type=1
%           f(x) = x^2     if tarch_type=2
%           g(x) = x       if tarch_type=2
%
% NOTE: This program generates 2000 more than required to minimize any starting bias
%
% See also TARCH, EGARCH_SIMULATE, APARCH_SIMULATE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==5
    error_type='NORMAL';
    tarch_type=2;
elseif nargin==6
    tarch_type=2;
elseif nargin==7
    %nothing
else
    error('5 to 7 inputs only.');
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
    error('Unknown ERROR_TYPE')
end

if ~(tarch_type==1 || tarch_type==2)
    error('TARCH_TYPE must be either 1 or 2')
end

if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end

if size(parameters,1)~=(1+p+o+q+extrap) || size(parameters,2)>1
    error('PARAMETERS must be a column vector with the correct number of parameters.');
end

if (sum(parameters(2:p+1)) + 0.5*sum(parameters(p+2:p+o+1)) + sum(parameters(p+o+2:p+o+q+1))) >=1
    warning('UCSD_GARCH:nonstationary','PARAMETERS are in the non-stationary space, be sure to check that H is not inf.');
end

if length(p)>1 || length(o)>1 || length(q)>1 || any(q<0) || any(p<1) || any(o<0)
    error('P, O and Q must be scalars with P positive and O and Q non-negative')
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
gamma=parameters(p+2:p+o+1);
beta=parameters(p+o+2:p+o+q+1);

%Initialize the random numbers
if isscalar(t)
    t=t+2000;
    if strcmp(error_type,'NORMAL')
        RandomNums=randn(t,1);
    elseif strcmp(error_type,'STUDENTST')
        nu=parameters(2+o+p+q);
        RandomNums=stdtrnd(nu,t,1);
    elseif strcmp(error_type,'GED')
        nu=parameters(2+o+p+q);
        RandomNums=gedrnd(nu,t,1);
    elseif strcmp(error_type,'SKEWT')
        nu=parameters(2+o+p+q);
        lambda=parameters(3+o+p+q);
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

m  =  max([p,o,q]);
%Back casts, zeros are fine since we are throwing away over 2000
RandomNums=[zeros(m,1);RandomNums];

if (1-sum(alpha)-sum(beta)-0.5*sum(gamma))>0
    UncondStd =  sqrt(omega/(1-sum(alpha)-sum(beta)-0.5*sum(gamma)));
else %non stationary
    UncondStd=1;
end

h=UncondStd.^2*ones(t+m,1);
data=UncondStd*ones(t+m,1);
T=size(data,1);
Idata=zeros(size(data));

parameters=[omega;alpha;gamma;beta];
if tarch_type==1
    for t = (m + 1):T
        h(t) = parameters' * [1 ; abs(data(t-(1:p)));  Idata(t-(1:o)).*abs(data(t-(1:o))); h(t-(1:q)) ];
        data(t)=RandomNums(t)*h(t);
        Idata(t)=data(t)<0;
    end
    h=h.^2;
else
    for t = (m + 1):T
        h(t) = parameters' * [1 ; data(t-(1:p)).^2;  Idata(t-(1:o)).*data(t-(1:o)).^2; h(t-(1:q)) ];
        data(t)=RandomNums(t)*sqrt(h(t));
        Idata(t)=data(t)<0;
    end
end
simulatedata=data((m+1+2000):T);
ht=h(m+1+2000:T);