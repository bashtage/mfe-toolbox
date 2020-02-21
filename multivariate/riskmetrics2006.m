function [Ht,weights] = riskmetrics2006(data,tau0,tau1,kmax,rho)
% Computes the Riskmetrics 2006 covariance, which is a weighted average of
% EWMA covariances
%
% USAGE:
%  [HT,WEIGHTS] = riskmetrics2006(DATA,TAU0,TAU1,KMAX,RHO)
%
% INPUTS:
%   DATA     - A T by K matrix of zero mean residuals -OR-
%                K by K by T array of covariance estimators (e.g. realized covariance)
%   TAU0     - [OPTIONAL] Half-life of slowest EWMA
%   TAU1     - [OPTIONAL] Half-life of fastest EWMA
%   KMAX     - [OPTIONAL] Number of EWMA components to use
%   RHO      - [OPTIONAL] Decay factor to use in half-lives.  The
%                half-lives used are TAU1, TAU1*RHO, TAU1*RHO^2, ...
%
% OUTPUTS:
%   HT       - A [K K T] dimension matrix of conditional covariances
%   WEIGHTS  - A T by T matrix which contains the final weights used
%                computing all conditional covariances.  WEIGHT(i,j)
%                contains the weight of the outer-product in time period j
%                in the covariance forecast at period i+1
%
% COMMENTS:
%   The conditional variance, H(t), of a RM2006 covariance model
%
%      H(t) = w(i)*Htilde(t,i)
%
%   where Htilde(t,i) is an EWMA covariance, i=1,2,..., and 
%
%   w(i) = 1-log(TAU1^((i-1)*RHO))/log(TAU0)
%
%   where the weights have been normalized so that they sum to 1.  All
%   EWMAs are initialized using a backward EWMA with the same decay.
%
% EXAMPLES:
%   RiskMetrics 2006 methodology using the suggested reference values
%     tau0   = 1560;
%     tau1   = 4;
%     taumax = 512;
%     rho    = sqrt(2);
%     kmax   = round(log(taumax/tau1)/log(rho)); % 14
%     Ht     = riskmetrics2006(data,tau0,tau1,kmax,rho)
%
%   RiskMetrics 1994 methodology using .94 as a special case of RM2006
%     tau0 = 1560;         % Does not matter
%     rho  = 1;            % Does not matter
%     tau1 = -1/log(.94);
%     kmax = 1;
%     Ht = riskmetrics2006(data,tau0,tau1,kmax,rho)


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 03/10/2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        tau0=[];
        tau1=[];
        kmax=[];
        rho =[];
    case 2
        tau1=[];
        kmax=[];
        rho =[];
    case 3
        kmax=[];
        rho =[];
    case 4
        rho =[];
    case 5
        % Nothing
    otherwise
        error('1 to 5 inputs required.')
end

if isempty(tau0)
    tau0 = 1560;
end
if isempty(tau1)
    tau1 = 4;
end
if isempty(kmax)
    kmax = 14;
end
if isempty(rho)
    rho = sqrt(2);
end

if tau1*rho^(kmax-1)>tau0
    error('The inputs must satisfy: TAU1*RHO^(KMAX-1)<TAU0')
end
if tau1<0
    error('TAU1 must be positive')
end
if tau0<0
    error('TAU0 must be positive')
end
if kmax<1 || floor(kmax)~=kmax
    error('KMAX must be an integer (weakly) larger than 1.')
end
if rho<0
    error('RHO must be positive')
end

if ndims(data)==2
    [T,K] = size(data);
    temp = zeros(K,K,T);
    for t=1:T
        temp(:,:,t) = data(t,:)'*data(t,:);
    end
    data = temp;
else
    K = size(data,1);
    T = size(data,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauks = tau1*rho.^((1:kmax)-1);
w= 1-log(tauks)/log(tau0);
w=w/sum(w); 
Ht = zeros(K,K,T);
Httilde = zeros(K,K,T);


for k=1:kmax
    tauk = tauks(k);
    mu = exp(-1/tauk);
    % back casting
    endPoint = max(min(floor(log(.01)/log(mu)),T),k);
    weights = (1-mu).*mu.^(0:endPoint-1);
    weights = weights/sum(weights);
    backCast = zeros(K);
    for i=1:endPoint
        backCast = backCast + weights(i)*data(:,:,i);
    end
    Httilde(:,:,1) = backCast;

    for t=2:T
        Httilde(:,:,t) = mu*Httilde(:,:,t-1) + (1-mu)*data(:,:,t-1);
    end
    Ht = Ht + w(k) * Httilde;
end


if nargout>1
    weights = zeros(T,T);
    for k=1:kmax
        tauk = tauks(k);
        mu = exp(-1/tauk);
        weightMatrix = [mu.^(0:T-1)' zeros(T,T-1)];
        weightMatrix(T,T) = (1-mu);
        for j=1:(T-2);
            weightMatrix(T,T-j) = mu * weightMatrix(T,T-j+1);
        end
        for t=2:T-1
            weightMatrix(t,2:t) = weightMatrix(T,T+(-t+2:0));
        end
        weights = weights + w(k) * weightMatrix;
    end
end