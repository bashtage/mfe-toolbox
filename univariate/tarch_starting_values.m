function [startingvals,nu,lambda,LLs,output_parameters]=tarch_starting_values(startingvals,data,fdata,fIdata,p,o,q,T,error_type,tarch_type)
% Perform a grid search to find decent starting values for TARCH(P,O,Q)
% esimtation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on error_type.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUT_PARAMETERS] = ...
%        tarch_starting_values(STARTINGVALS,DATA,FDATA,FIDATA,P,O,Q,T,ERROR_TYPE,TARCH_TYPE);
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   DATA             - Vector of mean zero residuals
%   FDATA            - Either abs(data) or data.^2, depending on tarch_type
%   FIDATA           - fdata times an indicator for negative, e.g. fdata.*(data<0)
%   P                - The lag order length for ARCH
%   O                - The lag order of asymmetric terms
%   Q                - The lag order length for GARCH
%   T                - Length of data
%   ERROR_TYPE       - The type of error being assumed, valid types are:
%                        1 if 'NORMAL'
%                        2 if 'STUDENTST'
%                        3 if 'GED'
%                        4 if 'SKEWT'
%   TARCH_TYPE        - 1 for absolute vale of return
%                     - 2 for squared returns (standard case)
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (1+p+o+q) by 1
%   NU                - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA            - Distribution asymmetry parameter, empty if not applicable
%   LLS               - A vector of log likelihoods corresponding to OUTPUT_PARAMETERS
%   OUTPUT_PARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%Initialize variables
LLs=[];
output_parameters=[];

%No starting values provided
if isempty(startingvals)
    nu=[];
    lambda=[];

    %Procedure is to find best starting values, using a grid search find values for normal, then

    %Possible starting values based on commonly estimated values
    a=[.05 .1 .2];
    la=length(a);
    g=[.01 .05 .2];
    lg=length(g)+1;
    agb=[.5 .8 .9 .95 .99];
    lb=length(agb);

    %Many output_parameters and LLs
    output_parameters=zeros(la*lb*lg,1+p+o+q);
    LLs=zeros(la*lb*lg,1);

    %Adjustment is needed to the intercept.  Assumes normality
    if tarch_type==1
        adj_factor=sqrt(2/pi);
        back_cast=mean(abs(data));
    else
        adj_factor=1;
        back_cast=cov(data);
    end
    covar=cov(data);

    %Use an index to count
    index=1;

    for i=1:la
        %Loop over a
        alpha=a(i);
        for j=1:lg
            %Loop over g
            if j==lg
                gamma=-alpha/2;
            else
                gamma=g(j);
            end

            for k=1:lb
                %Loop over beta
                temp_gamma=[];
                temp_alpha=alpha*ones(p,1)/p;
                %Make sure gamma satisfies necessary constraints
                if o>0
                    temp_gamma=gamma*ones(o,1)/o;
                    for n=1:o
                        if n<=p
                            temp_gamma(n)=max(temp_gamma(n),-alpha/(2*p));
                        else
                            temp_gamma(n)=0;
                        end
                    end
                end
                %Beta must also satisfy the same constraints
                beta=agb(k)-sum(temp_alpha)-0.5*sum(temp_gamma);
                %Pick omega to match the unconditional
                if tarch_type==1
                    omega=mean(abs(data))*(1-alpha*adj_factor-0.5*gamma-beta);
                else
                    omega=covar*(1-sum(temp_alpha)*adj_factor-0.5*sum(temp_gamma)-sum(beta));
                end
                %Build the parameter vector
                parameters=omega;
                parameters=[parameters; temp_alpha];
                parameters=[parameters; temp_gamma];
                if q>0
                    parameters=[parameters; beta*ones(q,1)/q];
                end
                %Set the output parameters
                output_parameters(index,:)=parameters';
                %Set the log likelihoods
                LLs(index)=tarch_likelihood(parameters, data, fdata, fIdata, p, o, q, 1, tarch_type, back_cast, T);
                %Increment
                index=index+1;
            end
        end
    end
    %Sort the LLs so the best (lowest, since we minimize the -1*LL)
    [LLs,index]=sort(LLs);
    %Order the starting values
    startingvals=output_parameters(index(1),:)';
    %Order the ouput_parameters
    output_parameters=output_parameters(index,:);
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
        nu=startingvals(p+o+q+2);
    end
    if error_type==4
        lambda=startingvals(p+o+q+3);
    end
    startingvals=startingvals(1:p+o+q+1);
end
