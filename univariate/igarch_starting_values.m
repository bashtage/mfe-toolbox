function [startingvals,nu,lambda,LLs,outputParameters]=igarch_starting_values(startingvals,epsilon,fepsilon,p,q,T,errorType,igarchType,constant)
% Perform a grid search to find decent starting values for IGARCH(P,Q)
% esimtation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on ERRORTYPE.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUTPARAMETERS] = ...
%        igarch_starting_values(STARTINGVALS,EPSILON,EPSILON2,P,Q,T,ERRORTYPE,IGARCHTYPE,CONSTANT)
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   EPSILON          - A column of mean zero data
%   FEPSILON         - Either abs(EPSILON) or EPSILON.^2, depending on IGARCHTYPE
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   T                - Length of EPSILON
%   ERRORTYPE        - [OPTIONAL] The error distribution used, valid types are:
%                       'NORMAL'    - Gaussian Innovations [DEFAULT]
%                       'STUDENTST' - T distributed errors
%                       'GED'       - Generalized Error Distribution
%                       'SKEWT'     - Skewed T distribution
%   IGARCHTYPE       - [OPTIONAL] The type of variance process, either
%                        1 - Model evolves in absolute values
%                        2 - Model evolves in squares [DEFAULT]
%   CONSTANT         - [OPTIONAL] Logical value indicating whether model
%                       should include a constant.  Default is true (include).
%
% OUTPUTS:
%   STARTINGVALS     - A vector of starting values (CONSTANT+p+q) by 1
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%   LLS              - A vector of log likelihoods corresponding to OUTPUTPARAMETERS
%   OUTPUTPARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also IGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


%Initialize variables
LLs=[];
outputParameters=[];

%No starting values provided
if isempty(startingvals)
    nu=[];
    lambda=[];
    
    %Procedure is to find best starting values, using a grid search find values for normal, then
    
    %Possible starting values based on commonly estimated values
    a=[.05 .1 .2];
    la=length(a);
    ab=1;
    lb=length(ab);
    
    %Many outputParameters and LLs
    outputParameters=zeros(la*lb,1+p+q-1);
    LLs=zeros(la,1);
    
    %Adjustment is needed to the intercept.  Assumes normality
    if igarchType==1
        adjFactor=sqrt(2/pi);
        backCast=mean(abs(epsilon));
    else
        adjFactor=1;
        backCast=cov(epsilon);
    end
    
    %Use an index to count
    index=1;
    
    for i=1:la
        %Loop over a
        alpha=a(i);
        tempAlpha=alpha*ones(p,1)/p;
        %Beta must satisfy the unit root
        beta=1-sum(tempAlpha);
        %Pick omega to match the unconditional
        if constant
            omega=backCast*.01*adjFactor;
        else
            omega = [];
        end
        %Build the parameter vector
        
        if q==0
            parameters=[omega; tempAlpha];
        else
            parameters=[omega; tempAlpha; beta*ones(q-1,1)/q];
        end
        %Set the output parameters
        outputParameters(index,:)=parameters';
        %Set the log likelihoods
        LLs(index)=igarch_likelihood(parameters, epsilon, fepsilon, p, q, 1, igarchType, constant, backCast, T, false);
        %Increment
        index=index+1;
    end
    %Sort the LLs so the best (lowest, since we minimize the -1*LL)
    [LLs,index]=sort(LLs);
    %Order the starting values
    startingvals=outputParameters(index(1),:)';
    %Order the ouputParameters
    outputParameters=outputParameters(index,:);
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
        nu=startingvals(constant+p+q);
    end
    if errorType==4
        lambda=startingvals(constant+p+q+1);
    end
    startingvals=startingvals(1:constant+p+q-1);
end
