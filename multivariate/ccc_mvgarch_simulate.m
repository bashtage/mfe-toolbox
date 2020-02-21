function [simulatedata,ht,pseudorc]=ccc_mvgarch_simulate(t,k,parameters,p,o,q,m)
% TARCH(p,o,q) - Constant Conditional Correlation MV GARCH simulation
%
% USAGE:
%   [SIMULATEDATA, HT, PSEUDORC] = ccc_mvgarch_simulate(T, K, PARAMETERS, P, O, Q, M)
%
% INPUTS:
%   T            - Length of the time series to be simulated
%   K            - Cross-sectional dimension
%   PARAMETERS   - a parameter vector of the form
%                  [tarch(1)' tarch(2)'  ... tarch(k)' corr_vech(R)]
%                    where each set of TARCH parameters is
%                    tarch(i) =
%                    [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                   ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
%                    and where R is the constant conditional correlation.
%   P            - Positive, scalar integer representing the number of symmetric innovations -OR-
%                     K by 1 vector of individual symmetric innovations order
%   O            - Non-negative, scalar integer representing the number of asymmetric lags -OR-
%                     K by 1 vector of individual asymmetric innovations order
%   Q            - Non-negative, scalar integer representing the number of conditional variance lags -OR-
%                     K by 1 vector of individual conditional variance lags
%   M            - [OPTIONAL] Number of ``intradaily'' returns to simulate to pseudo-Realized
%                    Covariance. If omitted, set to 72.
%
% OUTPUTS:
%   SIMULATEDATA - A time series with constant conditional correlation covariance
%   HT           - A [k k t] matrix of simulated conditional covariances
%   PSEUDORC     - A [k k t] matrix of pseudo-Realized Covariances
%
% COMMENTS:
%   The conditional variance, H(t), of a constant conditional correlation model is
%      H(t) = Sigma(t) * R * Sigma(t)
%
%       where Sigma(t) is a diagonal matrix with TARCH(P,O,Q) volatilities on its diagonal.
%
%   Pseudo Realized Covariances are simulated by generating m-intra daily returns from a N(0,1/m)
%   and computing the Realized Covariance of these. These were used in Patton and Sheppard (2009)
%   when evaluating variance and covariance specifications in a Monte Carlo.  If M=1, then PSEUDORC
%   is just the outer product of the SIMULATEDATA.
%
%   NOTE: This program generates 2000 more than required to minimize any start-up bias
%
% EXAMPLES:
%    % 3 by 3 CCC model
%    garch_parameters = [0.1 0.1 0.8]';
%    R = [1 .2 .5;.2 1 .3;.5 .3 1];
%    p = 1; o = 0; q = 1;
%    parameters = [garch_parameters; garch_parameters; garch_parameters];
%    parameters = [parameters; corr_vech(R)]
%    [data,Ht] = ccc_mvgarch_simulate(1000, 3, parameters, p, o, q) 
%
% See also TARCH_SIMULATE, SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/28/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 6
        m=[];
    case 7
        % nothing
    otherwise
        error('MFEToolbox:Input','6 or 7 inputs required.');
end

% t
if ~isscalar(t) || t<0 || floor(t)~=t
    error('T must be a positive integer.');
end
% k
if ~isscalar(k) || k<2 || floor(k)~=k
    error('K must be a positive integer greater than or equal to 2.');
end
% p
if length(p)==1
    if p<1 || floor(p)~=p
        error('P must be a positive integer if scalar.');
    end
    p = ones(k,1) * p;
else
    if length(p)~=k || min(size(p))~=1 || any(p<1) || any(floor(p)~=p)
        error('P must contain K positive integer elements if a vector.');
    end
end
% o
if length(o)==1
    if o<0 || floor(o)~=o
        error('O must be a non-neagative integer if scalar.');
    end
    o = ones(k,1) * o;
else
    if length(o)~=k || min(size(o))~=1 || any(o<0) || any(floor(o)~=o)
        error('O must contain K non-negative integer elements if a vector.');
    end
end
% q
if length(q)==1
    if q<0 || floor(q)~=q
        error('Q must be a non-neagative integer if scalar.');
    end
    q = ones(k,1) * q;
else
    if length(q)~=k || min(size(q))~=1 || any(q<0) || any(floor(q)~=q)
        error('Q must contain K non-negative integer elements if a vector.');
    end
end
% m
if isempty(m)
    m = 72;
end
if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
% parameters
parameterCount = k + sum(p) + sum(o) + sum(q) + k*(k-1)/2;
if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end
n = length(parameters);
if n~=parameterCount
    error('PARAMETERS must have K + sum(P) + sum(O) + sum(Q) + K*(K-1)/2 elements.');
end
R = corr_ivech(parameters(n-k*(k-1)/2+1:n));
if min(eig(R))<=0 || any(any(R~=R'))
    error('The correlation matrix R must be postive definite.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the correlation matrix


%Set the burnin amount, 2000 is probably reasonable
burnin=2000;
%Up t by the burnin amount
t=t+burnin;



%Draw some normal random numbers
intraRandomNums=randn(m*t,k)*sqrt(1/m);
intraRandomNums = intraRandomNums * R^(0.5);
randomNums = cumsum(intraRandomNums);
randomNums = diff([zeros(1,k);randomNums(m:m:m*t,:)]);

parameterStart = 1;
data = zeros(t,k);
htMat = zeros(t,k);
for i=1:k
    parameterEnd = parameterStart + (1+p(i)+o(i)+q(i)) - 1;
    tarchParameters = parameters(parameterStart:parameterEnd);
    [data(:,i),htMat(:,i)]=tarch_simulate(randomNums(:,i),tarchParameters,p(i),o(i),q(i),'NORMAL',2);
    parameterStart = parameterEnd + 1;
end

%Initialize the covariance
ht=repmat(R,[1 1 t]);
pseudorc=zeros(k,k,t);
for i=1:t
    vol = sqrt(htMat(i,:));
    ht(:,:,i) = ht(:,:,i) .* (vol'*vol);
    r = intraRandomNums((i-1)*m+1:i*m,:).*repmat(vol,m,1);
    pseudorc(:,:,i) = r'*r;
end

%Truncate the data and the covariance
simulatedata=data(burnin+1:t,:);
pseudorc = pseudorc(:,:,burnin+1:t);
ht=ht(:,:,burnin+1:t);