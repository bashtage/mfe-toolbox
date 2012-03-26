function [LL, LLS, ht] = tarch_likelihood(parameters, data, fdata, fIdata, p, o, q, error_type, tarch_type, back_cast, T, estim_flag)
% Log likelihood for TARCH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = tarch_likelihood(PARAMETERS, DATA, FDATA, FIDATA, P, O, Q, ERROR_TYPE, TARCH_TYPE, BACK_CAST, T, ESTIM_FLAG)
%
% INPUTS:
%   PARAMETERS    - A vector of GARCH process parameters
%                   [omega alpha gamma beta [nu lambda]]
%   DATA          - Vector of mean zero residuals
%   FDATA         - Either abs(data) or data.^2, depending on tarch_type
%   FIDATA        - fdata times an indicator for negative, e.g. fdata.*(data<0)
%   P             - The lag order length for ARCH
%   O             - The lag order of asymmetric terms
%   Q             - The lag order length for GARCH
%   ERROR_TYPE    - The type of error being assumed, valid types are:
%                     1 if 'NORMAL'
%                     2 if 'STUDENTST'
%                     3 if 'GED'
%                     4 if 'SKEWT'
%   TARCH_TYPE    - 1 for absolute vale of return
%                 - 2 for squared returns (standard case)
%   BACK_CAST     - The value used for variance recursion
%   T             - Length of data
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
%   See also TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin==12 && estim_flag
    %If for estimation, transform the parameters
    [parameters,nu,lambda]=tarch_itransform(parameters,p,o,q,error_type);
else
    %Otherwise the parameters simply must be parsed
    if error_type==2 || error_type==3
        %Seperate nu from the remaning parameters
        nu=parameters(p+o+q+2);
        parameters=parameters(1:1+p+o+q);
    elseif error_type==4
        lambda=parameters(p+o+q+3);
        nu=parameters(p+o+q+2);
        parameters=parameters(1:1+p+o+q);
    end
end

%Backcast length
m  =  max([p o q]);
%Compute the conditional variances
ht=tarch_core(fdata,fIdata,parameters,back_cast,p,o,q,m,T,tarch_type);

%Indices for the relevant opservations
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