function [ll, lls, ht] = ccc_mvgarch_joint_likelihood(parameters,data,volData,p,o,q)
% Log likelihood for CCC_MVGARCH inference.  This function is only used for computing the Hessian.
%
% USAGE:
%   [LL, LLS] = ccc_mvgarch_joint_likelihood(PARAMETERS,DATA,VOLDATA,P,O,Q)
%
% INPUTS:
%   PARAMETERS - a parameter vector of the form
%                  [tarch(1)' tarch(2)'  ... tarch(k)' corr_vech(R)]
%                    where each set of TARCH parameters is
%                    tarch(i) =
%                    [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                   ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
%                    and where R is the constant conditional correlation.
%   DATA       - K by K by T 3D array of positive semi-definite matrices.
%   VOLDATA    - T by K matrix of data used in estimating the volatility models.
%   P          - K by 1 vector of symmetric innovation orders for the TARCH models
%   O          - K by 1 vector of asymmetric innovation orders for the TARCH models
%   Q          - K by 1 vector of lagged variance innovation orders for the TARCH models
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

[k,nothing,t]=size(data);
htMat = zeros(t,k);
% Construct the variances
backCastlength = max(floor(t^(1/2)),1);
backCastWeights = .05*(.9.^(0:backCastlength ));
backCastWeights = backCastWeights/sum(backCastWeights);

parameterCount = 1;
for i=1:k
    epsilon = volData(:,i);
    m = max([p(i) o(i) q(i)]);
    fdata =  [mean(epsilon.^2)*ones(m,1) ; epsilon.^2];
    %fIepsilon is fepsilon*(epsilon<0)
    fIdata =  [0.5*mean(epsilon.^2)*ones(m,1) ; epsilon.^2.*(epsilon<0)];
    backCast = backCastWeights*(volData(1:backCastlength+1,i).^2);
    if backCast==0
        backCast=mean(epsilon.^2);
    end
    
    T = length(fdata);
    parameterCountEnd = parameterCount + (1+p(i)+o(i)+q(i)) - 1;
    tarchParameters = parameters(parameterCount:parameterCountEnd);
    variance=tarch_core(fdata,fIdata,tarchParameters,backCast,p(i),o(i),q(i),m,T,2);
    htMat(:,i) = variance(m+1:T);
    parameterCount = parameterCountEnd+1;
end
% Construct the correlation matrix
R = corr_ivech(parameters(parameterCount:parameterCount+(k*(k-1)/2)-1));
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
if nargout>2
    ht = zeros(k,k,t);
    for i=1:t
        h12 = sqrt(htMat(i,:));
        ht(:,:,i) = R .* (h12'*h12);
    end
end
