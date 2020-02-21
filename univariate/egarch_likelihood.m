function [LL, LLS, ht] = egarch_likelihood(parameters, data, p, o, q, error_type, back_cast, T, estim_flag)
% Log likelihood for EGARCH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = egarch_likelihood(PARAMETERS, DATA, P, O, Q, ERROR_TYPE, BACK_CAST, T, ESTIM_FLAG)
%
% INPUTS:
%   PARAMETERS    - A vector of GARCH process parameters
%                   [omega alpha gamma beta [nu lambda]]
%   DATA          - Vector of mean zero residuals
%   P             - The lag order length for ARCH
%   O             - The lag order of asymmetric terms
%   Q             - The lag order length for GARCH
%   ERROR_TYPE    - The type of error being assumed, valid types are:
%                     1 if 'NORMAL'
%                     2 if 'STUDENTST'
%                     3 if 'GED'
%                     4 if 'SKEWT'
%   BACK_CAST     - The value used for variance recursion
%   T             - Length of data
%   ESTIM_FLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                   is being used in estimation.  If it is 1, then the parameters are
%                   transformed from unconstrained values to constrained by standard
%                   garch model constraints
%
% OUTPUTS:
%   LL             - Minus 1 times the log likelihood
%   HT             - Time series of conditional variances
%   LLS            - Time series of log likelihoods (Also multiplied by -1)
%
% COMMENTS:
%   See also EGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin==9 && length(estim_flag==1) && all(estim_flag==1); %If for estimation, transform the parameters
    [parameters]=egarch_itransform(parameters,p,o,q,error_type);
elseif nargin==8 || (length(estim_flag)==1 && all(estim_flag==0))
    %Nothing to do
else
    %This is the case where we are simply estimating the intercept for the starting values,
    %estim_flag contains the other parameter values
    parameters=[parameters;estim_flag];
end

if error_type==2 || error_type==3
    %Separate nu from the remaining parameters
    nu=parameters(p+o+q+2);
    parameters=parameters(1:1+p+o+q);
elseif error_type==4
    lambda=parameters(p+o+q+3);
    nu=parameters(p+o+q+2);
    parameters=parameters(1:1+p+o+q);
end

%Backcast length
m  =  max([p o q]);
%Upper is a bound for the variances to keep things from getting crazy
upper=10000*max(data.^2);
%Compute the conditional variances
ht=egarch_core(data,parameters,back_cast,upper,p,o,q,m,T);

%Indices for the relevant observations
t = (m + 1):T;
ht=ht(t);

data=data(t);
%Compute the log likelihoods
switch error_type
    case 1
        [LL, LLS] = normloglik(data,0,ht);
        LLS = -LLS;
        LL = -LL;
    case 2
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

if isnan(LL) || isinf(LL)
   keyboard
end