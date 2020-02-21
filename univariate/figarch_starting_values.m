function [startingvals,nu,lambda,LLs,outputParameters]=figarch_starting_values(startingvals,epsilon,epsilon2,p,q,errorType,truncLag)
% Perform a grid search to find starting values for FIGARCH(Q,D,P)
% estimation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on ERRORTYPE.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUTPARAMETERS] = ...
%        figarch_starting_values(STARTINGVALS,EPSILON,EPSILON2AUGMENTED,P,Q,ERRORTYPE,TRUNCLAG)
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   EPSILON          - Vector of mean zero residuals
%   EPSILON2         - TRUNCLAG + T by 1 column vector containing TRUNCLAG backcasts followed by
%                        EPSILON.^2.  That is [zeros(TRUNCLAG,1)+BACKCAST;EPSILON.^2]
%   P                - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q                - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERRORTYPE        - The type of error being assumed, valid types are:
%                        1 if 'NORMAL'
%                        2 if 'STUDENTST'
%                        3 if 'GED'
%                        4 if 'SKEWT'
%   TRUNCLAG
%
% OUTPUTS:
%   STARTINGVALS     - A vector of starting values (1+p+o+q) by 1
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%   LLS              - A vector of log likelihoods corresponding to OUTPUTPARAMETERS
%   OUTPUTPARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also FIGARCH
 
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
    ds = [.2 .5 .7];
    phiRatio = [.2 .5 .8];
    betaRatio = [.1 .5 .9];
    
    covar = cov(epsilon);
    % Construct all of the parameters
    ld=length(ds);
    lp=length(phiRatio);
    lq=length(betaRatio);
    outputParameters=zeros(ld,2+p+q);
    index = 1;
    for i=1:lp
        for j=1:lq
            for k=1:ld
                d = ds(k);
                if p
                    phi = (1-d)/2 * phiRatio(i);
                else
                    phi = 0;
                end
                beta = (d + phi) * betaRatio(j);
                if p && q
                    temp = [phi d beta];
                elseif p
                    temp = [phi d];
                elseif q
                    temp = [d beta];
                else
                    temp = d;
                end
                lambda = figarch_weights(temp,p,q,truncLag);
                omega = covar * (1-sum(lambda));
                outputParameters(index,:)  = [omega temp];
                index = index + 1;
            end
        end
    end
    
    outputParameters = unique(outputParameters,'rows');
    %Many outputParameters and LLs
    LLs=zeros(size(outputParameters,1),1);
 
    %Use an index to count
    for i=1:size(outputParameters,1);
        LLs(i) = figarch_likelihood(outputParameters(i,:),p,q,epsilon,epsilon2,truncLag,1,false);
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
        nu=startingvals(p+q+2);
    end
    if errorType==4
        lambda=startingvals(p+q+3);
    end
    startingvals=startingvals(1:p+q+1);
end


