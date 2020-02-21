function [simulatedata, ht] = aparch_simulate(t,parameters,p,o,q,error_type)
% APARCH(P,O,Q) time series simulation with multiple error distributions
%
% USAGE:
%   [SIMULATEDATA, HT] = aparch_simulate(T, PARAMETERS, P, O, Q, ERROR_TYPE)
%
% INPUTS:
%   T            - Length of the time series to be simulated  OR
%                  T by 1 vector of user supplied random numbers (i.e. randn(1000,1))
%   PARAMETERS   - a 1+P+O+Q (+1 or 2, depending on error distribution) x 1 parameter vector
%                  [omega alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q) delta 
%                  [nu lambda]]'
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative scalar integer representing the number of asymmetric innovations (0
%                    for symmetric processes).  Must be less than or equal to P
%   Q            - Non-negative, scalar integer representing the number of lags of conditional
%                    variance (0 for ARCH) 
%   ERROR_TYPE   - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%
% OUTPUTS:
%   SIMULATEDATA - A time series with APARCH variances
%   HT           - A vector of conditional variances used in making the time series
%
% COMMENTS:
%   The conditional variance, h(t), of a APARCH(P,O,Q) process is modeled as follows:
%
%    h(t)^(delta/2) = omega
%             + alpha(1)*(abs(r(t-1))+gamma(1)*r(t-1))^delta + ...
%               alpha(p)*(abs(r(t-p))+gamma(p)*r(t-p))^delta +
%               beta(1)*h(t-1)^(delta/2) +...+ beta(q)*h(t-q)^(delta/2)
%
%    Required restrictions on parameters:
%    delta > 0
%    -1<gamma<1
%    -1<lambda<1
%    nu>2 for T
%    nu>1 for GED
%    alpha(i) > 0
%
% NOTE: This program generates 2000 more than required to minimize any starting bias
%
% EXAMPLES:
%   Simulate a GARCH(1,1)
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 .85 2], 1, 0, 1)
%   Simulate an AVARCH(1,1)
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 .85 1], 1, 0, 1)
%   Simulate a GJR-GARCH(1,1,1)
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 -.1 .8 2], 1, 1, 1)
%   Simulate a TARCH(1,1,1)
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 -.1 .8 1], 1, 1, 1)
%   Simulate an APARCH(1,1,1)
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 -.1 .8 .8], 1, 1, 1)
%   Simulate an APARCH(1,1,1) with Student's T innovations
%       [SIMULATEDATA, HT] = aparch_simulate(1000, [.1 .1 -.1 .85 2 6], 1, 1, 1, 'STUDENTST')
%
% See also APARCH, TARCH_SIMULATE, EGARCH_SIMULATE

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
    error('5 to 6 inputs only.');
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

if o>p
    error('O must be less than P')
end

if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end
if size(parameters,1)~=(2+p+o+q+extrap) || size(parameters,2)>1
    error('parameters must be a column vector with the correct number of parameters.');
end

%Separate the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
gamma=parameters(p+2:p+o+1);
beta=parameters(p+o+2:p+o+q+1);
delta = parameters(p+o+q+2);

if o>0
    gamma2=[gamma;zeros(p-o,1)];
    sumag = 0.5*sum(alpha.*gamma2);
else
    sumag=0;
end

if (sum(alpha)+ sumag + sum(beta)) >=1
    warning('UCSD_GARCH:nonstationary','PARAMETERS may be in the non-stationary space, be sure to check that H is not inf.');
    nonstationary =true;
else
    nonstationary =false;
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

%Initialize the random numbers
if isscalar(t)
    t=t+2000;
    if strcmp(error_type,'NORMAL')
        RandomNums=randn(t,1);
    elseif strcmp(error_type,'STUDENTST')
        nu=parameters(3+o+p+q);
        RandomNums=stdtrnd(nu,t,1);
    elseif strcmp(error_type,'GED')
        nu=parameters(3+o+p+q);
        RandomNums=gedrnd(nu,t,1);
    elseif strcmp(error_type,'SKEWT')
        nu=parameters(3+o+p+q);
        lambda=parameters(4+o+p+q);
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
%Back casts, zeros are fine since we are throwing away 2000
RandomNums=[zeros(m,1);RandomNums];

if nonstationary
    UncondStd =  sqrt(omega/(1-sum(alpha)-sum(beta)));
else %non stationary
    UncondStd=1;
end

hdelta=UncondStd.^delta*ones(t+m,1);
h=UncondStd.^delta*ones(t+m,1);
data=UncondStd*ones(t+m,1);
T=size(data,1);


for t = (m + 1):T
    hdelta(t) = omega;
    for j=1:p
        if o>=j
            hdelta(t) = hdelta(t) + alpha(j)*(abs(data(t-j))+gamma(j)*data(t-j))^delta;
        else
            hdelta(t) = hdelta(t) + alpha(j)*abs(data(t-j))^delta;
        end
    end
    for j=1:q
        hdelta(t) = hdelta(t) + beta(j)*hdelta(t-j);
    end
    h(t) = hdelta(t)^(2/delta);
    data(t) = sqrt(h(t))*RandomNums(t);
end

simulatedata=data((m+1+2000):T);
ht=h(m+1+2000:T);