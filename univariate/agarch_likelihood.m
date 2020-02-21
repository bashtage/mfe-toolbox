function [LL, LLS, ht] = agarch_likelihood(parameters, epsilon_aug, p, q, model_type, error_type, transform_bounds, back_cast, T, estim_flag)
% Log likelihood for AGARCH(P,Q) and NAGARCH(P,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = agarch_likelihood(PARAMETERS, EPSILON_AUG, P, Q, ERROR_TYPE, MODEL_TYPE, TRANSFORM_BOUNDS, BACK_CAST, T, ESTIM_FLAG)
%
% INPUTS:
%   PARAMETERS    - A vector of GARCH process parameters
%                   [omega alpha gamma beta [nu lambda]]
%   EPSILON_AUG   - A column of mean zero data audmented with M backcasts
%   P             - The lag order length for ARCH
%   Q             - The lag order length for GARCH
%   ERROR_TYPE    - The type of error being assumed, valid types are:
%                     1 if 'NORMAL'
%                     2 if 'STUDENTST'
%                     3 if 'GED'
%                     4 if 'SKEWT'
%   MODEL_TYPE    - Variable indicating model type
%                     1 for AGARCH
%                     2 for NAGARCH
%   BACK_CAST     - The value used for variance recursion
%   TRANS_BOUNDS  - 2 by 1 vector containing the .01 and .99 quantiles
%                     of EPSILON for use in parameter transformation
%   T             - Length of data
%   ESTIM_FLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                   is being used in estimation.  If it is 1, then the parameters are
%                   transformed from unconstrained values to constrained by standard
%                   garch model constraints
%
% OUTPUTS:
%   LL            - Minus 1 times the log likelihood
%   LLS           - Time series of log likelihoods (Also multiplied by -1)
%   HT            - Time series of conditional variances
%
% COMMENTS:
%   See also AGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


if nargin==10 && estim_flag
    % If for estimation, transform the parameters
    [parameters,nu,lambda]=agarch_itransform(parameters,p,q,model_type,error_type,transform_bounds);
else
    % Otherwise the parameters simply must be parsed
    if error_type==2 || error_type==3
        % Seperate nu from the remaning parameters
        nu=parameters(p+q+3);
        parameters=parameters(1:2+p+q);
    elseif error_type==4
        lambda=parameters(p+q+4);
        nu=parameters(p+q+3);
        parameters=parameters(1:2+p+q);
    end
end
% Backcast length
m  = max([p q]);
% Compute the conditional variances
ht = agarch_core(epsilon_aug,parameters,back_cast,p,q,m,T,model_type);

% Indices for the relevant opservations
t    = (m + 1):T;
ht   = ht(t);
epsilon = epsilon_aug(t);
% Compute the log likelihoods
switch error_type
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