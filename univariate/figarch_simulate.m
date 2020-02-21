function [simulatedata, ht, lambda] = figarch_simulate(t,parameters,p,q,errorType,truncLag,bcLength)
% FIGARCH(Q,D,P) time series simulation with multiple error distributions for P={0,1} and Q={0,1}
%
% USAGE:
%   [SIMULATEDATA, HT, LAMBDA] = figarch_simulate(T, PARAMETERS, P, Q, ERRORTYPE, TRUNCLAG, BCLENGTH)
%
% INPUTS:
%   T             - Length of the time series to be simulated  OR
%                     T by 1 vector of user supplied random numbers (i.e. randn(1000,1))
%   PARAMETERS    - a 2+P+Q (+1 or 2, depending on error distribution) x 1 parameter vector
%                     [omega phi d beta [nu lambda]]'.  Parameters should satisfy conditions in
%                     FIGARCH_ITRANSFORM 
%   P             - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q             - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERRORTYPE    - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   TRUNCLAG     - [OPTIONAL] Truncation lag for use in the construction of lambda. Default value is
%                    2500. 
%   BCLENGTH     - [OPTIONAL] Number of extra observations to produce to reduce start up bias.
%                    Default value is 2500. 
%
% OUTPUTS:
%   SIMULATEDATA  - A time series with ARCH/GARCH/GJR/TARCH variances
%   HT            - A vector of conditional variances used in making the time series
%   LAMBDA        - TRUNCLAG by 1 vector of weights used when computing the conditional variances
%
% COMMENTS:
%    The conditional variance, h(t), of a FIGARCH(1,d,1) process is modeled as follows:
%
%    h(t) = omega + [1-beta L - phi L (1-L)^d] epsilon(t)^2 + beta * h(t-1)
%    
%    which is estimated using an ARCH(oo) representation, 
%
%    h(t) = omega + sum(lambda(i) * epsilon(t-1)^2)
%    
%    where lambda(i) is a function of the fractional differencing parameter, phi and beta
%
% EXAMPLES:
%   FIGARCH(0,d,0) simulation
%       simulatedData = figarch_simulate(2500, [.1 .42],0,0)
%   FIGARCH(1,d,1) simulation
%       simulatedData = figarch_simulate(2500, [.1 .1 .42 .4],1,1)
%   FIGARCH(0,d,0) simulation with Student's T errors
%       simulatedData = figarch_simulate(2500, [.1 .42],0,0,'STUDENTST')
%   FIGARCH(0,d,0) simulation with a truncation lag of 5000
%       simulatedData = figarch_simulate(2500, [.1 .42],0,0,[],5000)
%
% See also FIGARCH, FIGARCH_TRANSFORM, FIGARCH_ITRANSFORM

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/18/2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 4
    errorType='NORMAL';
    truncLag = 2500;
    bcLength = 2500;
    case 5
    truncLag = 2500;
    bcLength = 2500;
    case 6
    bcLength = 2500;
    case 7
        % Nothing
    otherwise 
    error('4 to 7 inputs only.');
end

if isempty(errorType)
    errorType='NORMAL';
end

if strcmp(errorType,'NORMAL')
    extrap=0;
elseif strcmp(errorType,'STUDENTST')
    extrap=1;
elseif strcmp(errorType,'GED')
    extrap=1;
elseif strcmp(errorType,'SKEWT')
    extrap=2;
else
    error('Unknown error type')
end

if size(parameters,2)>size(parameters,1);
    parameters = parameters';
end
if size(parameters,1)~=(2+p+q+extrap) || size(parameters,2)>1
    error('PARAMETERS must be a column vector with the correct number of parameters.');
end

if length(p)>1 || length(q)>1 || ~ismember(p,[0 1]) || ~ismember(q,[0 1])
    error('P and Q must be scalars with P positive and O and Q non-negative')
end

if size(t,2)~=1
    error('T must be either a positive scalar or a vector of random numbers.');
end

%Separate the parameters
omega = parameters(1);
if p
    phi = parameters(2);
    d = parameters(3);
else
    phi = 0;
    d = parameters(2);
end
if q
    beta = parameters(3+p);
else
    beta = 0;
end
if omega<0
    error('omega must be positive');
end
if phi<0 || beta<0
    error('phi and beta must be non-negative')
end
if phi + d - beta<0
    error('phi + d - beta must be non-negative')
end
lambda = figarch_weights(d,0,0,2);
if phi>lambda(2)/lambda(1);
    error('phi must be less than lambda(2)/lambda(1) from a call to figarch_weights(d,0,0,2)');
end

if isempty(truncLag)
    truncLag = 2500;
end

if isempty(bcLength)
    bcLength = 2500;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize the random numbers
if isscalar(t)
    t=t+bcLength;
    if strcmp(errorType,'NORMAL')
        RandomNums=randn(t,1);
    elseif strcmp(errorType,'STUDENTST')
        nu=parameters(3+p+q);
        RandomNums=stdtrnd(nu,t,1);
    elseif strcmp(errorType,'GED')
        nu=parameters(3+p+q);
        RandomNums=gedrnd(nu,t,1);
    elseif strcmp(errorType,'SKEWT')
        nu=parameters(3+p+q);
        lambda=parameters(4+p+q);
        RandomNums=skewtrnd(nu,lambda,t,1);
    else
        error('Unknown error type')
    end
else
    RandomNums=t;
    t=length(RandomNums);
    seeds=ceil(rand(bcLength,1)*t);
    RandomNums=[RandomNums(seeds);RandomNums];
    t=length(RandomNums);
end
lambda = figarch_weights(parameters(2:2+p+q),p,q,truncLag);
%Back casts, zeros are fine since we are throwing away over 2000
if sum(lambda)<1
    backCast = omega/(1-sum(lambda));
else
    backCast = omega/.01;
end
% Setup
r2=[backCast*ones(truncLag,1);zeros(size(RandomNums))];
e2=[zeros(truncLag,1);RandomNums.^2];
e=[zeros(truncLag,1);RandomNums];
data=[zeros(truncLag,1);zeros(size(RandomNums))];
h = zeros(size(e2));
% Loop to construct variances
for i=truncLag+1: t + truncLag
    h(i) = omega + lambda'*r2(i-1:-1:i-truncLag);
    r2(i) = e2(i) * h(i);
    data(i) = e(i) * sqrt(h(i));
end

% Discard unused data
tau = (truncLag+1+bcLength):truncLag+t;
simulatedata=data(tau);
ht=h(tau);
