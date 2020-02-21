function [startingvals,nu,lambda]=agarch_starting_values(startingvals,epsilon,p,q,model_type,error_type)
% Provides starting values for AGARCH(P,Q) or NAGARCH(P,Q) estimation.
% If starting values are user supplied (and thus nonempty), reformats
% starting values depending on error_type.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUT_PARAMETERS] = ...
%        agarch_starting_values(STARTINGVALS,EPSILON,P,Q,MODEL_TYPE,ERROR_TYPE)
%
% INPUTS:
%   STARTINGVALS  - A vector of starting values or empty to perform a grid search
%   EPSILON       - A column of mean zero data
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH-type model)
%   MODEL_TYPE    - [OPTIONAL] The type of variance process, either
%                     'AGARCH'  - Asymmetric GARCH, Engle (1990) [DEFAULT]
%                     'NAGARCH' - Nonlinear Asymmetric GARCH, Engle & Ng (1993)
%   ERROR_TYPE    - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (2+p+q) by 1
%   NU                - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA            - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   See also AGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


%No starting values provided
if isempty(startingvals)
    nu=[];
    lambda=[];
    
    options  =  optimset('fminunc');
    options.TolFun = 1e-005;
    options.TolX = 1e-005;
    options.Display = 'off';
    options.Diagnostics = 'off';
    options.LargeScale = 'off';
    options.MaxFunEvals = 1200;
    % Call tarch to get starting parameters using a standard GARCH model
    startingvals=tarch(epsilon,p,0,q,'NORMAL',[],[],options);
    
    omega = startingvals(1);
    alpha = startingvals(2:p+1);
    beta = startingvals(p+2:p+q+1);
    
    if model_type == 2
        if (sum(alpha) + sum(beta))>.97
            scale = .97 * (sum(alpha) + sum(beta));
            alpha = alpha*scale;
            beta = beta*scale;
        end
    end
    % Put in a value of 0 for gamma
    startingvals = [omega;alpha;0;beta];
    %Use generic values for nu and lambda if needed
    if error_type==2
        nu=8;
    elseif error_type==3
        nu=1.9;
    elseif error_type==4
        nu=8;
        lambda=-.1;
    end
else
    %Values provided, only parse them
    nu=[];
    lambda=[];
    if error_type==2 || error_type==3 || error_type==4
        nu=startingvals(p+q+3);
    end
    if error_type==4
        lambda=startingvals(p+q+4);
    end
    startingvals=startingvals(1:p+q+2);
end