function [simulatedata, ht] = egarch_simulate(t,parameters,p,o,q,error_type)
% EGARCH(P,O,Q) time series simulation with multiple error distributions
%
% USAGE:
%   [SIMULATEDATA, HT] = egarch_simulate(T, PARAMETERS, P, O, Q, ERROR_TYPE)
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
%
% OUTPUTS:
%   SIMULATEDATA  - A time series with EGARCH variances
%   HT            - A vector of conditional variances used in making the time series
%
% COMMENTS:
%   The conditional variance, h(t), of a EGARCH(P,O,Q) process is modeled as follows:
%   ln(h(t)) = omega
%            + alpha(1)*(abs(e_{t-1})-C) + ... + alpha(p)*(abs(e_{t-p})-C)+...
%            + gamma(1)*e_{t-1} +...+ e_{t-o} +...
%              beta(1)*ln(h(t-1)) +...+ beta(q)*ln(h(t-q))
%
%       where: ln is natural log
%              e_t = r_t/sqrt(h_t)
%              C   = 1/sqrt(pi/2)
%
%   NOTE: This program generates 2000 more than required to minimize any starting bias
%
% EXAMPLES:
%   Simulate a symmetric EGARCH(1,0,1) process                             
%       simulatedDate = egarch_simulate(1000,[0 .1  .95],1,0,1);                 
%   Simulate a standard EGARCH(1,1,1) process                              
%       simulatedDate = egarch_simulate(1000,[0 .1 -.1 .95],1,1,1);              
%   Simulate a standard EGARCH(1,1,1) process with Student's T innovations 
%       simulatedDate = egarch_simulate(1000,[0 .1 -.1 .95 6],1,1,1,'STUDENTST');
%   Simulate a standard EGARCH(1,1,1) process with GED innovations         
%       simulatedDate = egarch_simulate(1000,[0 .1 -.1 .95 1.5],1,1,1,'GED');    
%   
% See also EGARCH, TARCH_SIMULATE, APARCH_SIMULATE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==5
    error_type='NORMAL';
elseif nargin==6
    %nothing
else
    error('5 or 6 inputs only.');
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

if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end
if size(parameters,1)~=(1+p+o+q+extrap) || size(parameters,2)>1
    error('parameters must be a column vector with the correct number of parameters.');
end


if  ~all(egarch_nlcon(parameters,1,p,o,q,1,1,1,1)<0)
    warning('parameters are in the non-stationary space, be sure to make sure H isn''t inf.');
end

if length(p)>1 || length(o)>1 || length(q)>1 || any(q<0) || any(p<1) || any(o<0)
    error('p and q must be scalars with p positive and q non-negative')
end

if size(t,2)~=1  && (size(t,1)>1 || floor(t)==t);
    error('T must be either a positive scalar or a vector of random numbers.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Seperate the parameters
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
%Back-casts, zeros are fine since we are throwing away over 2000
RandomNums=[zeros(m,1);RandomNums];

if all(egarch_nlcon(parameters,1,p,o,q,1,1,1,1)<0) && sum(beta)<1 %Bad guesstimate
    UncondStd =  exp(omega./(1-sum(beta)));
else %non stationary
    UncondStd=1;
end

h=UncondStd.^2*ones(t+m,1);
absRandomNums=abs(RandomNums);
T=size(h,1);

parameters=[omega;alpha;gamma;beta];
const=1/sqrt(pi/2);
for t = (m + 1):T
    h(t) = parameters' * [1 ; absRandomNums(t-(1:p))-const;  RandomNums(t-(1:o)); h(t-(1:q)) ];
end
h=exp(h);
data=RandomNums.*sqrt(h);

simulatedata=data((m+1+2000):T);
ht=h(m+1+2000:T);
