function [ll,lls,Ht] = gogarch_likelihood(parameters,data,p,q,gjrType,P,L,isOgarch,isInference)
% Log-likelihood for use in estimation GOGARCH and OGARCH models
%
% USAGE:
%  [LL,LLS,HT] = gogarch_likelihood(PARAMETERS,DATA,P,Q,GJRTYPE,P,L,ISOGARCH,ISINFERNCE)
%
% INPUTS:
%   PARAMETERS - K*(K-1)/2 + sum(P) + sum(Q) by 1 vector of parameters (only sum(P) + sum(Q) if ISOGARCH)
%   DATA       - K by K by Tarray of data (either realized measures or outer-products of daily data)
%   P          - K by 1 vector of positive, scalar integer representing the number of symmetric innovations
%   Q          - K by 1 vector of non-negative, scalar integer representing the number of conditional covariance lags
%   GJRTYPE    - K by 1 vector containing either 1 (TARCH/AVGARCH) or 2 (GJRGARCH/GARCH/ARCH)
%   P          - Eigenvector of unconditional covariance matrix of the data
%   L          - Diagonal matrix containing the eigenvalues of unconditional covariance matrix of the data
%   ISOGARCH   - Boolean indicating that the model is OGARCH (otherwise GOGARCH)
%   ISINFERNCE - Boolean indicating the likelihood is being used for inference, so that the first
%                  K(K+1)/2 parameters are ivech(S)
%
% OUTPUTS:
%   LL         - The log likelihood computed at PARAMETERS
%   LLS        - T by 1 vector of log-likelihoods
%   HT         - K by K by T vector of conditional covariances
%
% COMMENTS:

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/15/2012

[k,~,T] = size(data);
if size(parameters,1)>size(parameters,2)
    parameters = parameters';
end

offset = 0;
if isInference
    S =  ivech(parameters(1:k*(k+1)/2));
    [P,L] = eig(S);
    P = P';
    offset = offset + k*(k+1)/2;
end

if ~isOgarch
    phi = parameters(offset + (1:k*(k-1)/2));
    U = phi2u(phi);
    offset = offset + k*(k-1)/2;
else
    U = eye(k);
end
Z = P*L^(0.5)*U;
Zinv = U'*L^(-0.5)*P';
stdData = zeros(k,k,T);
for t=1:T
    stdData(:,:,t) = Zinv*data(:,:,t)*Zinv';
end
% Univariate GARCH models
V = zeros(T,k);
w = .06 * .94.^(0:sqrt(T));
w = w/sum(w);
likData = zeros(T,k);
for i=1:k
    count = p(i) + q(i);
    volParameters = parameters(offset + (1:count));
    volParameters = max(volParameters,0);
    volParameters = [1-sum(volParameters) volParameters]; %#ok<AGROW>
    offset = offset + count;
    likData(:,i) = squeeze(stdData(i,i,:));
    volData = sqrt(likData(:,i));
    backCast = w*volData(1:length(w)).^2;
    v = tarch_core_simple(volData,volParameters,backCast,0,p(i),0,q(i),gjrType(i));
    V(:,i) = v;
end

likConst = k*log(2*pi);
logdetZZp = log(det(Z*Z'));
lls = 0.5*(likConst + logdetZZp + sum(log(V),2) + sum(likData./V,2));
ll = sum(lls);

if ~isreal(ll) || isnan(ll) || isinf(ll)
    ll = 1e7;
end

if nargout>2
    Ht = zeros(k,k,T);
    for t=1:T
        Ht(:,:,t) = Z*diag(V(t,:))*Z';
    end
end