function [LL, LLS, ht] = aparch_likelihood(parameters, data_aug, abs_data_aug, p, o, q, errorType, T, deltaIsEstimated, userDelta, estim_flag)
% Log likelihood for APARCH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = aparch_likelihood(PARAMETERS, DATA_AUG, ABS_DATA_AUG, P, O, Q, ERRORTYPE, T, deltaIsEstimated ESTIM_FLAG)
%
% INPUTS:
%   PARAMETERS    - A vector of GARCH process parameters
%                   [omega alpha gamma beta [nu lambda]]
%   DATA_AUG      - Vector of mean zero residuals augmented with zeros
%   ABS_DATA_AUG  - Absolute value of augmented data
%   P             - The lag order length for ARCH
%   O             - The lag order of asymmetric terms
%   Q             - The lag order length for GARCH
%   ERRORTYPE    - The type of error being assumed, valid types are:
%                     1 if 'NORMAL'
%                     2 if 'STUDENTST'
%                     3 if 'GED'
%                     4 if 'SKEWT'
%   T             - Length of data
%   deltaIsEstimated   - 1 if no user delta has been provided, 0 otherwise
%   USERDELTA     - The value of delta to use if specified ot empty otherwise
%   ESTIM_FLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                   is being used in estimation.  If it is 1, then the parameters are
%                   transformed from unconstrained values to constrained by standard
%                   garch model constraints
%
% OUTPUTS:
%   LL             - Minus 1 times the log likelihood
%   LLS            - Time series of log likelihoods (Also multiplied by -1)
%   HT             - Time series of conditional variances
%
% COMMENTS:
%   See also APARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005
if nargin==11 && ~isempty(estim_flag)
    %If for estimation, transform the parameters
    oldp=parameters;
    [parameters,nu,lambda]=aparch_itransform(parameters,p,o,q,errorType,deltaIsEstimated);
    if ~deltaIsEstimated
        parameters=[parameters;userDelta];
    end
else
    %Otherwise the parameters simply must be parsed
    if errorType==2 || errorType==3
        %Seperate nu from the remaning parameters
        nu=parameters(p+o+q+2+deltaIsEstimated);
        parameters=parameters(1:2+p+o+q);
    elseif errorType==4
        nu=parameters(p+o+q+2+deltaIsEstimated);
        lambda=parameters(p+o+q+3+deltaIsEstimated);
        parameters=parameters(1:2+p+o+q);
    end
    if ~deltaIsEstimated
        parameters=[parameters;userDelta];
    end
end
%Backcast length
m  =  max([p o q]);

delta = parameters(1+p+o+q+1);
% Local back casting
data=data_aug(m+1:T);
back_cast_length = max(floor(length(data)^(1/2)),1);
back_cast_weights = .05*(.9.^(0:back_cast_length ));
back_cast_weights = back_cast_weights/sum(back_cast_weights);
back_cast = back_cast_weights*(abs(data(1:back_cast_length+1)).^delta);
if back_cast==0
    back_cast=mean(abs(data_aug(m+1:T)).^delta);
end

%Upper and lower bounds, should be moved out of this function
% Updated to not depend on delta
LB = (cov(data_aug(m+1:T))/100000);
UB = (100*max(data_aug.^2));

%Compute the conditional variances

ht=aparch_core(data_aug,abs_data_aug,parameters,p,o,q,m,T,back_cast,LB,UB);
%ht=aparch_core2(data_aug,abs_data_aug,parameters,p,o,q,m,T,back_cast,LB,UB);

%Indices for the relevant opservations
t = (m + 1):T;
ht=ht(t);
data=data_aug(t);
%Compute the log likelihoods
switch errorType
    case 1
        [LL, LLS] = normloglik(data,0,ht);
        LLS = -LLS;
        LL = -LL;
    case 2
        if ~isreal(nu)
            keyboard
        end
        [LL, LLS] = stdtloglik(data,0,ht,nu);
        LLS = -LLS;
        LL = -LL;
    case 3
        [LL, LLS] = gedloglik(data,0,ht,nu);
        LLS = -LLS;
        LL = -LL;
    case 4
        [LL, LLS] = skewtloglik(data,0,ht,nu,lambda);
        LLS = -LLS;
        LL = -LL;
end