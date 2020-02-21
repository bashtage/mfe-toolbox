function [parameters,ht,w,pc]=o_mvgarch(data,numfactors,p,o,q,startingVals,options)
% Estimates a multivariate GARCH model using Orthogonal or Factor Garch.  Involves Principal Component
% Analysis to reduce volatility modelling to univariate garches.  See Carrol 2000 (An Inrtoduction to O-Garch)
%
% USAGE:
%   [PARAMETERS,HT,W,PC] = o_mvgarch(DATA,NUMFACTORS,P,O,Q);
%
% INPUTS:
%   DATA       - A T by K matrix of zero mean residuals
%   NUMFACTORS - The number of principal components to include in the MV_GARCH model
%   P          - Positive, scalar integer representing the number of symmetric innovations -OR-
%                   K by 1 vector of individual symmetric innovations order
%   O          - Non-negative, scalar integer representing the number of asymmetric lags -OR-
%                   K by 1 vector of individual asymmetric innovations order
%   Q          - Non-negative, scalar integer representing the number of conditional variance lags -OR-
%                   K by 1 vector of individual conditional variance lags
%
%   OPTIONS    - [OPTIONAL] A fminunc options structure to use in estimating the factor conditional
%                  variances
%
% OUTPUTS:
%   PARAMETERS - a parameter vector of the form
%                [tarch(1)' tarch(2)'  ... tarch(NUMFACTORS)' diag(OMEGA)]
%                where each set of TARCH parameters is
%                tarch(i) =
%                [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                 ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
%                and where OMEGA is the diagonal matrix of idiosyncratic variances.
%                If NUMFACTORS = k diag(OMEGA) is omitted since all values are 0.
%   HT         - A [K K T] dimension matrix of conditional covariances
%   W          - a K by K matrix of componet weights
%   PC         - a T by K matrix of principal componets
%
% COMMENTS:
%   Models the conditional covariance of assets using the NUMFACTORS principal components with the
%   highest contribution to total variation in the data.  The conditional covariance is
%
%   H_t = W*F_t*W' + Omega
%
%   where F_t is a diagonal matrix of the conditional factor variances, W is a K by NUMFACTORS
%   matrix of factor loadings and Omega is a diagonal matrix of (time-invariant) idiosyncratic
%   variances.  If NUMFACTORS = K when Omega = 0.
%
% See also PCA, CCC_MVGARCH, DCC, GOGARCH


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 10/28/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        options = [];
        startingVals = [];
    case 6
        options = [];
    case 7
        % Nothing
    otherwise
        error('5 to 7 arguments required.')
end

[t,k] = size(data);

if min(t,k)<2 || t<k
    error('DATA must be a T by K matrix, T>K>1');
end

%p, o, q much be non-negative scalars

% p
if length(p)==1
    if p<1 || floor(p)~=p
        error('P must be a positive integer if scalar.');
    end
    p = ones(numfactors,1) * p;
else
    if length(p)~=numfactors || min(size(p))~=1 || any(p<1) || any(floor(p)~=p)
        error('P must contain K positive integer elements if a vector.');
    end
end

% o
if length(o)==1
    if o<0 || floor(o)~=o
        error('O must be a non-neagative integer if scalar.');
    end
    o = ones(numfactors,1) * o;
else
    if length(o)~=numfactors || min(size(o))~=1 || any(o<0) || any(floor(o)~=o)
        error('O must contain K non-negative integer elements if a vector.');
    end
end

% q
if length(q)==1
    if q<0 || floor(q)~=q
        error('Q must be a non-neagative integer if scalar.');
    end
    q = ones(numfactors,1) * q;
else
    if length(q)~=numfactors || min(size(q))~=1 || any(q<0) || any(floor(q)~=q)
        error('Q must contain K non-negative integer elements if a vector.');
    end
end

if  ~isempty(startingVals)
    if size(startingVals,2)>size(startingVals,1)
        startingVals = startingVals';
    end
    if length(startingVals)<(k+sum(p)+sum(o)+sum(q))
        error('STARTINGVALS should be a K+sum(P)+sum(O)+sum(Q) by 1 vector');
    end
end

%Make sure options is a valid option structure
if isempty(options)
    options = optimset('fminunc');
    options.Display = 'off';
    options.LargeScale = 'off';
end

try
    optimset(options);
catch ME
    error('OPTIONS is not a valid options structure');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w, pc] = pca(data,'outer');



weights = w(1:numfactors,:);
pcs = pc(:,1:numfactors);

htMat = zeros(t,numfactors);
tarchParameters = cell(numfactors,1);
tarchHt = cell(numfactors,1);
tarchVcv = cell(numfactors,1);
tarchScores = cell(numfactors,1);
for i=1:numfactors
    % Prepare the data
    tarchStartingVals = [];
    if ~isempty(startingVals)
        parameterEnd = i + sum(p(1:i)) + sum(o(1:i)) + sum(q(1:i));
        parameterStart = parameterEnd - 1 - p(i) - o(i) - q(i) + 1;
        tarchStartingVals = startingVals(parameterStart:parameterEnd);
    end
    [tarchParameters{i},~,tarchHt{i},~,tarchVcv{i},tarchScores{i}]=tarch(pcs(:,i),p(i),o(i),q(i),[],2,tarchStartingVals,options);
    htMat(:,i) = tarchHt{i};
end

if numfactors<k
    errors = data - pcs * weights;
    omega = diag(mean(errors.^2));
else
    omega = zeros(k);
end

ht = zeros(k,k,t);

for i=1:t
    ht(:,:,i) = weights' * diag(htMat(i,:)) * weights + omega;
end

parameters=zeros(numfactors+sum(p+o+q),1);
count = 1;
for i=1:numfactors;
    parameters(count:count+p(i)+o(i)+q(i))=tarchParameters{i};
    count = count + 1 + p(i) + o(i) + q(i);
end
