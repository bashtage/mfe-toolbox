function [ll, lls] = ccc_mvgarch_likelihood(parameters,data,htMat)
% Log likelihood for CCC_MVGARCH inference.  This function is only used for computing scores.
%
% USAGE:
%   [LL, LLS] = ccc_mvgarch_likelihood(PARAMETERS,DATA,HTMAT)
%
% INPUTS:
%   PARAMETERS - A vector of constant correlation parameters, K*(K-1)/2 by 1 where 
%                  R = corr_ivech(PARAMETERS)
%   DATA       - K x K x T 3D array of positive semi-definite matrices.  
%   HTMAT      - T by K matrix of fit conditional variances.
%
% OUTPUTS:
%   LL         - Minus 1 times the log likelihood
%   LLS        - Time series of log likelihoods (Also multiplied by -1)
%
% COMMENTS:
%   See also CCC_MVGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/28/2009

[k, nothing, t] = size(data);
% Construct the correlation matrix
R = corr_ivech(parameters);
% Compute the likelihood
Rinv = inv(R);
LogDetR = log(det(R));
llconst = k*log(2*pi);
lls = zeros(t,1);
for i=1:t
    stdResid = data(:,:,i)./sqrt(htMat(i,:)'*htMat(i,:));
    lls(i) = 0.5 * (llconst + sum(log(htMat(i,:))) + LogDetR + trace(Rinv*stdResid));
end
ll = sum(lls);
