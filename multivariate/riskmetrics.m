function Ht = riskmetrics(data,lambda,backCast)
% Computes the Riskmetrics or other EWMA covariance
%
% USAGE:
%  [HT] = riskmetrics(DATA,LAMBDA,BACKCAST)
%
% INPUTS:
%   DATA     - A T by K matrix of zero mean residuals -OR-
%                K by K by T array of covariance estimators (e.g. realized covariance)
%   LAMBDA   - EWMA smoothing parameter 0<LAMBDA<1
%   BACKCAST - [OPTIONAL] Covariance matrix to use as initial value.  If
%                omitted, a backward EWMA is used
%
% OUTPUTS:
%   HT       - A [K K T] dimension matrix of conditional covariances
%
% COMMENTS:
%   The conditional variance, H(t), of an EWMA covariance model
%      H(t) = (1-lambda)*r(t-1,:)'*r(t-1,:) + lambda*H(t-1)
%
% EXAMPLES:
%   % The standard RiskMetrics EWMA covariance is computed using
%   Ht = riskmetrics(data,.94)
%   % The standard RiskMetrics EWMA covariance using the unconditional
%   % covariance as the backcast
%   Ht = riskmetrics(data,.94,cov(data))

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 03/10/2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        backCast = [];
    case 3
        % nothing
    otherwise
        error('2 or 3 inputs required.')
end


if ndims(data)==2
    [T,k] = size(data);
    temp = zeros(k,k,T);
    for t=1:T
        temp(:,:,t) = data(t,:)'*data(t,:);
    end
    data = temp;
end

if lambda<=0 || lambda>=1
    error('LAMBDA must be between 0 and 1.')
end

if isempty(backCast)
    endPoint = max(min(floor(log(.01)/log(lambda)),T),k);
    weights = (1-lambda).*lambda.^(0:endPoint-1);
    weights = weights/sum(weights);
    backCast = zeros(k);
    for i=1:endPoint
        backCast = backCast + weights(i)*data(:,:,i);
    end
    
end

backCast = (backCast+backCast)/2;
if min(eig(backCast))<0
    error('BACKCAST must be positive semidefinite if provided.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ht = zeros(k,k,T);
Ht(:,:,1) = backCast;
for i=2:T
    Ht(:,:,i) = (1-lambda)*data(:,:,i-1) + lambda * Ht(:,:,i-1);
end