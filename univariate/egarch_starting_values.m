function [startingvals,nu,lambda,LLs,output_parameters]=egarch_starting_values(startingvals,data,p,o,q,T,error_type)
% Perform a grid search to find decent starting values for EGARCH(P,O,Q)
% esimtation.  If starting values is user supplied (and thus nonempty), reformats
% starting values depending on error_type.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUT_PARAMETERS] = ...
%        egarch_starting_values(STARTINGVALS,DATA,P,O,Q,T,ERROR_TYPE);
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   DATA             - Vector of mean zero residuals
%   P                - The lag order length for ARCH
%   O                - The lag order of asymmetric terms
%   Q                - The lag order length for GARCH
%   T                - Length of data
%   ERROR_TYPE       - The type of error being assumed, valid types are:
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
%   See also EGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%These are the options for the search.  Should be pretty loose
options=optimset('fminbnd');
options.TolX=1e-3;

LLs=[];
output_parameters=[];
if isempty(startingvals)
    nu=[];
    lambda=[];

    %Proceedure is to find best starting values, usign a grid search find values for normal, then
    %find nu lambda

    a=[.01 .05 .2 .5];
    la=length(a);
    g=[-.2  .01 .1];
    lg=length(g);
    b=[.01 .5 .9 ];
    lb=length(b);
    output_parameters=zeros(la*lb*lg,1+p+o+q);
    LLs=zeros(la*lb*lg,1);

    covar=cov(data);
    back_cast=log(covar);
    index=1;
    for i=1:la
        alpha=a(i);
        for j=1:lg
            gamma=g(j);
            for k=1:lb
                if p>0
                    temp_alpha=alpha*ones(p,1)/p;
                else
                    temp_alpha=[];
                end
                if o>0
                    temp_gamma=gamma*ones(o,1)/o;
                else
                    temp_gamma=[];
                end
                if q>0
                    temp_beta=b(k)*ones(q,1)/q;
                else
                    temp_beta=[];
                end
                %Pick omega in a sensible way
                %Nelson's formulation:
                %ln(ht)=omega + a *abs(r)-a*mean(abs(r)) + g * r + b*ln(ht-1)
                omega_UB=1.1*log(covar);
                omega_LB=log(covar)*.0001;
                if  omega_UB<omega_LB
                    temp=omega_UB;
                    omega_UB=omega_LB;
                    omega_LB=temp;
                end
                const_parameters=temp_alpha;
                const_parameters=[const_parameters; temp_gamma];
                const_parameters=[const_parameters; temp_beta];
                
                [omega,LL]=fminbnd('egarch_likelihood',omega_LB,omega_UB,options,data,p, o, q, 1, back_cast, T,const_parameters);
                parameters=[omega;const_parameters];
                output_parameters(index,:)=parameters';
                LLs(index)=LL;
                index=index+1;
            end
        end
    end
    [LLs,index]=sort(LLs);
    startingvals=output_parameters(index(1),:)';
    output_parameters=output_parameters(index,:);
    %Set up the shape parameters if needed
    if error_type==2
        nu=8;
    elseif error_type==3
        nu=1.9;
    elseif error_type==4
        nu=8;
        lambda=-.1;
    end
else
    %In the case of supplies startig values, simply parse them
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
