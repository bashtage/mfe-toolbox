function [LL,LLS,ht] = figarch_likelihood(parameters,p,q,epsilon,epsilon2,truncLag,errorType,estimFlag)
% Log likelihood for FIGARCH(Q,D,P) estimation
%
% USAGE:
%   [LL, LLS, HT] = figarch_likelihood(PARAMETERS, P, Q, EPSILON, EPSILON2, TRUNCLAG, ERRORTYPE, ESTIMFLAG)
%
% INPUTS:
%   PARAMETERS - A vector of FIGARCH process parameters
%                   [omega phi d beta [nu lambda]]'
%   EPSILON    - T by 1 Column vector of mean zero residuals
%   EPSILON2   - TRUNCLAG + T by 1 column vector containing TRUNCLAG backcasts followed by
%                  EPSILON.^2.  That is [zeros(TRUNCLAG,1)+BACKCAST;EPSILON.^2]
%   P          - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q          - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERRORTYPE  - The type of error being assumed, valid types are:
%                     1 if 'NORMAL'
%                     2 if 'STUDENTST'
%                     3 if 'GED'
%                     4 if 'SKEWT'
%   ESTIMFLAG  - [OPTIONAL] Flag (0 or 1) to indicate if the function is being used in estimation.
%                  If it is 1, then the parameters are transformed from unconstrained values to
%                  constrained by standard garch model constraints
%
% OUTPUTS:
%   LL             - Minus 1 times the log likelihood
%   LLS            - Time series of log likelihoods (Also multiplied by -1)
%   HT             - Time series of conditional variances
%
% COMMENTS:
%  See also FIGARCH, FIGARCH_PARAMETER_CHECK, FIGARCH_WEIGHTS
%  FIGARCH_STARTING_VALUES, FIGARCH_TRANSFORM, FIGARCH_ITRANSFORM

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009



if nargin==8 && estimFlag
    [parameters,nu,lambda] = figarch_itransform(parameters,p,q,errorType);
else
    if errorType == 2 || errorType == 3
        nu = parameters(end);
        parameters = parameters(1:end-1);
    elseif errorType==4
        nu = parameters(end-1);
        lambda = parameters(end);
        parameters = parameters(1:end-2);
    end
end

omega = parameters(1);
figarchWeightParameters = parameters(2:2+p+q);
T = size(epsilon,1);
archWeights = figarch_weights(figarchWeightParameters,p,q,truncLag);
tau = truncLag+1:truncLag+T;
ht = zeros(size(epsilon2));
for t = tau;
    ht(t) = omega + archWeights'*epsilon2(t-1:-1:t-truncLag);
end
ht = ht(tau);

%Compute the log likelihoods
switch errorType
    case 1
        [LL, LLS] = normloglik(epsilon,0,ht);
        LLS = -LLS;
        LL = -LL;
    case 2
        [LL, LLS] = stdtloglik(epsilon,0,ht,nu);
        LLS = -LLS;
        LL = -LL;
    case 3
        [LL, LLS] = gedloglik(epsilon,0,ht,nu);
        LLS = -LLS;
        LL = -LL;
    case 4
        [LL, LLS] = skewtloglik(epsilon,0,ht,nu,lambda);
        LLS = -LLS;
        LL = -LL;
end