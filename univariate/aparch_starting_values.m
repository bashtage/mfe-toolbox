function [startingvals,nu,lambda,LLs,output_parameters]=aparch_starting_values(startingvals,data,p,o,q,errorType,deltaIsEstimated)
% Perform a grid search to find decent starting values for APARCH(P,O,Q)
% esimtation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on ERRORTYPE.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUT_PARAMETERS] = ...
%        aparch_starting_values(STARTINGVALS,DATA,P,O,Q,ERROR_TYPE);
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   DATA             - Vector of mean zero residuals
%   P                - The lag order length for ARCH
%   O                - The lag order of asymmetric terms
%   Q                - The lag order length for GARCH
%   ERRORTYPE       - The type of error being assumed, valid types are:
%                        1 if 'NORMAL'
%                        2 if 'STUDENTST'
%                        3 if 'GED'
%                        4 if 'SKEWT'
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (1+p+o+q) by 1
%   NU                - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA            - Distribution asymmetry parameter, empty if not applicable
%   LLS               - A vector of log likelihoods corresponding to OUTPUT_PARAMETERS
%   OUTPUT_PARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also APARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005



%No starting values provided
if isempty(startingvals)
    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-003);
    options  =  optimset(options , 'TolX'        , 1e-004);
    options  =  optimset(options , 'Display'     , 'off');
    options  =  optimset(options , 'Diagnostics' , 'off');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxFunEvals' , 200*(2+p+o+q));
    
    [parameters, LLs] = tarch(data,p,0,q,[],1,[],options);
    if deltaIsEstimated
        parameters=[parameters(1:p+1); zeros(o,1) ; parameters(p+2:p+q+1) ;1];
    else
        parameters=[parameters(1:p+1); zeros(o,1) ; parameters(p+2:p+q+1)];
    end
    output_parameters=parameters;
    
    nu=[];
    lambda=[];
    %Possible starting values 
    startingvals=parameters;
    %Use generic values for nu and lambda if needed
    if errorType==2
        nu=8;
    elseif errorType==3
        nu=1.9;
    elseif errorType==4
        nu=8;
        lambda=-.1;
    end
else
    %Values provided, only parse them
    nu=[];
    lambda=[];
    if errorType==2 || errorType==3 || errorType==4
        nu=startingvals(p+o+q+2+deltaIsEstimated);
    end
    if errorType==4
        lambda=startingvals(p+o+q+3+deltaIsEstimated);
    end
    startingvals=startingvals(1:p+o+q+1+deltaIsEstimated);
    LLs=[];
    output_parameters=[];
end


