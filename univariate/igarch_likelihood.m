function [LL, LLS, ht] = igarch_likelihood(parameters, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, estimFlag)
% Log likelihood for IGARCH(P,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = igarch_likelihood(PARAMETERS,EPSILON,FEPSILON,P,Q,ERRORTYPE,IGARCHTYPE,CONSTANT,BACKCAST,T,ESTIMFLAG)
%
% INPUTS:
%   PARAMETERS   - A vector of IGARCH process parameters
%                    [omega alpha beta [nu lambda]]
%   DATA         - Vector of mean zero residuals
%   EPSILON      - A column of mean zero data
%   FEPSILON     - Either abs(EPSILON) or EPSILON.^2, depending on IGARCHTYPE
%   P            - Positive, scalar integer representing the number of
%                    symmetric innovations
%   Q            - Non-negative, scalar integer representing the number
%                    of lags of conditional variance (0 for ARCH)
%   ERRORTYPE    - [OPTIONAL] The error distribution used, valid types are:
%                    'NORMAL'    - Gaussian Innovations [DEFAULT]
%                    'STUDENTST' - T distributed errors
%                    'GED'       - Generalized Error Distribution
%                    'SKEWT'     - Skewed T distribution
%   IGARCHTYPE   - [OPTIONAL] The type of variance process, either
%                    1 - Model evolves in absolute values
%                    2 - Model evolves in squares [DEFAULT]
%   CONSTANT     - [OPTIONAL] Logical value indicating whether model
%                    should include a constant.  Default is true (include).
%   BACKCAST     - The value used for variance recursion
%   T            - Length of data
%   ESTIMFLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                    is being used in estimation.  If it is 1, then the parameters are
%                    transformed from unconstrained values to constrained by standard
%                    garch model constraints
%
% OUTPUTS:
%   LL             - Minus 1 times the log likelihood
%   LLS            - Time series of log likelihoods (Also multiplied by -1)
%   HT             - Time series of conditional variances
%
% COMMENTS:
%   See also IGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009

if nargin==11 && estimFlag
    %If for estimation, transform the parameters
    [parameters,nu,lambda]=igarch_itransform(parameters,p,q,errorType,constant);
else
    %Otherwise the parameters simply must be parsed
    if errorType==2 || errorType==3
        %Seperate nu from the remaning parameters
        nu=parameters(p+q+constant);
        parameters=parameters(1:constant+p+q-1);
    elseif errorType==4
        lambda=parameters(p+q+constant+1);
        nu=parameters(p+q+constant);
        parameters=parameters(1:constant+p+q-1);
    end
end

%Backcast length
m  =  max([p q]);
%Compute the conditional variances
ht=igarch_core(fepsilon,parameters,backCast,p,q,m,T,igarchType,constant);
%Indices for the relevant opservations
t  = (m + 1):T;
ht = ht(t);
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
