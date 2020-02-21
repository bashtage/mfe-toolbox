function [simulatedata,ht,pseudorc]=matrix_garch_simulate(t,k,parameters,p,o,q,m)
% Simulation of symmetric and asymmetric MATRIX multivariate GARCH models
%
% USAGE:
%   [SIMULATEDATA, HT, PSEUDORC] = matrix_garch_simulate(T, K, PARAMETERS, P, O, Q, M)
%
% INPUTS:
%   T            - Length of the time series to be simulated
%   K            - Cross-sectional dimension
%   PARAMETERS   - A 3-D matrix of K by K matrices where
%                    CC' = PARAMETERS(:,:,1), AA'(j)=PARAMETERS(:,:,1+j),
%                    GG'(j) = PARAMETERS(:,:,1+P+j), BB'(j)=PARAMETERS(:,:,1+P+Q+j)
%                    -OR- a K(K+1)/2*(1+P+O+Q) by 1 vector of parameters of
%                    the form returned my calling matrix_garch
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of asymmetric lags to include
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   M            - [OPTIONAL] Number of ``intradaily'' returns to simulate to pseudo-Realized
%                    Covariance. If omitted, set to 72.
%
% OUTPUTS:
%   SIMULATEDATA - A time series with constant conditional correlation covariance
%   HT           - A [k k t] matrix of simulated conditional covariances
%   PSEUDORC     - A [k k t] matrix of pseudo-Realized Covariances
%
% COMMENTS:
%    The conditional variance, H(t), of a MATRIX GARCH is modeled as follows:
%
%      H(t) = CC' + AA'(1).*r_{t-1}'*r_{t-1} + ... + AA'(P).*r_{t-P}'*r_{t-P}
%                 + GG(1)'.*n_{t-1}'*n_{t-1} + ... + GG(O)'.*n_{t-P}'*n_{t-P}
%                  + BB(1)'.*H(t-1) +...+ BB(Q)'.*H(t-q)
%
%    where n_{t} = r_{t} .* (r_{t}<0).  If using realized measures, the
%    RM_{t-1} replaces r_{t-1}'*r_{t-1}, and the asymmetric version
%    replaces n_{t-1}'*n_{t-1}
%
%   Pseudo Realized Covariances are simulated by generating m-intra daily returns from a N(0,1/m)
%   and computing the Realized Covariance of these. These were used in Patton and Sheppard (2009)
%   when evaluating variance and covariance specifications in a Monte Carlo.  If M=1, then PSEUDORC
%   is just the outer product of the SIMULATEDATA.
%
%   NOTE: This program generates 2000 more than required to minimize any start-up bias
%
% See also TARCH_SIMULATE, CCC_GARCH_SIMULATE, MATRIX_GARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/10/2011


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
if ~isscalar(p) || p<1 || floor(p)~=p
    error('P must be a positive integer if scalar.');
end

% o

if ~isscalar(o) || o<0 || floor(o)~=o
    error('O must be a non-neagative integer if scalar.');
end
% q

if ~isscalar(q) || q<0 || floor(q)~=q
    error('Q must be a non-neagative integer if scalar.');
end
% m
if isempty(m)
    m = 72;
end
if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
% parameters
k2 = k*(k+1)/2;
if ismatrix(parameters)
    parameterCount = k*(k+1)/2 *(1+ p + o + q);
    if size(parameters,2)>size(parameters,1)
        parameters = parameters';
    end
    if length(parameters)~=parameterCount
        error('PARAMETERS must be K(K+1)/28(1+P+O+Q) when using a vector.')
    end
    parameterMatrices= zeros(k,k,1+p+o+q);
    index = 0;
    for i=1:(1+p+o+1)
        temp = vec2chol(parameters(index+1:index+k2));
        parameterMatrices(:,:,i) = (temp*temp'+temp*temp')/2;
        index=index+k2;
    end
    parameters = parameterMatrices;
elseif ndims(parameters)==3
    if any(size(parameters)~=[k k 1+p+o+q])
        error('PARAMETERS must be K by K by (1+P+O+Q) when using a 3-D matrix.')
    end
    parameterMatrices = parameters;
else
    error('The size of PARAMETERS is not compatible with this function.')
end

% Check stationarity
sumParameters = zeros(k);
for i=2:size(parameterMatrices,3)
    if i<=(p+1) || i>(1+p+o)
        w=1;
    else
        w=0.5;
    end
    sumParameters = sumParameters + w*parameterMatrices(:,:,i);
end
sumParameters = (sumParameters+sumParameters')/2;
if max(diag(sumParameters))>1
    warning('MFE::Nonstationary','The parameters do not correspond to a stationary solution.  Check for overflow.')
    uncond = parameters(:,:,1)/.005;
else
    uncond = parameters(:,:,1)./(ones(k)-sumParameters);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the burnin amount, 2000 is probably reasonable
burnin=2000;
%Up t by the burnin amount
t=t+burnin;

%Draw some normal random numbers
intraRandomNums=randn(m*t,k)*sqrt(1/m);
randomNums = cumsum(intraRandomNums);
randomNums = diff([zeros(1,k);randomNums(m:m:m*t,:)]);


%Perform the recursion
backCast = uncond;
backCastAsym = uncond;
T = t;
Ht = zeros(k,k,T);
simulatedData = zeros(T,k);
for t=1:T;
    Ht(:,:,t)=parameterMatrices(:,:,1);
    for j=1:p
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,j+1).*backCast;
        else
            r = simulatedData(t-j,:);
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,j+1).*(r'*r);
        end
    end
    for j=1:o
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+j+1).*backCastAsym;
        else
            n = simulatedData(t-j,:).*(simulatedData(t-j,:)<0);
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+j+1).*(n'*n);
        end
    end    
    for j=1:q
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+o+j+1).*backCast;
        else
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+o+j+1).*Ht(:,:,t-j);
        end
    end
    simulatedData(t,:) = randomNums(t,:)*Ht(:,:,t)^(0.5);
end





%Initialize the covariance
pseudorc=zeros(k,k,t);
for t=1:T
    r = intraRandomNums((t-1)*m+1:t*m,:)*Ht(:,:,t)^(0.5);
    pseudorc(:,:,t) = r'*r;
end

%Truncate the data and the covariance
simulatedata=simulatedData(burnin+1:t,:);
pseudorc = pseudorc(:,:,burnin+1:t);
ht=Ht(:,:,burnin+1:t);