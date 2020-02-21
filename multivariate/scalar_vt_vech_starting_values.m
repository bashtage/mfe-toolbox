function [startingvals,lls,output_parameters]=scalar_vt_vech_starting_values(startingvals,data,dataAsym,p,o,q,C,Casym,kappa,useComposite,backCast,backCastAsym)
% Perform a grid search to find decent starting values for SCALAR_VT_VECH(P,Q) esimtation.  If
% starting values is user supplied (and thus nonempty), does nothing.
%
% USAGE:
%   [STARTINGVALS,LLS,OUTPUT_PARAMETERS] = ...
%        scalar_vt_vech_starting_values(STARTINGVALS,DATAAUG,P,Q,T,C);
%
% INPUTS:
%   STARTINGVALS - A vector of starting values or empty to perform a grid search
%   PARAMETERS   - A vector of vech GARCH process parameters: [alpha beta]'
%   DATA         - Augmented (by m back cast values) matrix of mean zero residuals
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   T            - Length of the original data
%   C            - The unconditional covariance of the data (cov(data)
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (p+q) by 1
%   LLS               - A vector of log likelihoods corresponding to OUTPUT_PARAMETERS
%   OUTPUT_PARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%This code performs a quick grid search of potential starting values
t = size(data,3);
a=[.005 .01 .025 .05]; % The innovation term is usually small
g=[.005 .02 .5 .1];
apgpb=[.99 .95 .90 .5]; % The smoothing parameter is usually large
count = 1;
params = zeros(length(a)*length(g)*length(apgpb),p+o+q);
for i=1:length(a)
    alpha = ones(1,p)*a(i)/p;
    for j=1:length(g);
        gamma = ones(1,o)*g(j)/o;
        for k = 1:length(apgpb)
            beta = apgpb(k)-sum(alpha)-sum(gamma)/kappa;
            if beta<=apgpb(k)/2;
                scale = 0.5 * apgpb(k)/(sum(alpha)+sum(gamma)/kappa);
                alpha = alpha*scale;
                gamma = gamma*scale;
                beta = apgpb(k)-sum(alpha)-sum(gamma)/kappa;
                % Reduce alpha and gamma
            end
            beta = beta*ones(1,q)/q;
            params(count,:) = [alpha gamma beta];
            count = count+1;
        end
    end
end
params = unique(params,'rows');
params = params(~all(params==0,2),:);



if isempty(startingvals)
    lls=zeros(length(params),1);
    %Evaluate the likelihood at each parameter
    for i=1:size(params,1)
        parameters=params(i,:);
        lls(i)=scalar_vt_vech_likelihood(parameters,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,false);
    end
    [lls,pos]=sort(lls);
    output_parameters=params(pos,:);
    startingvals=output_parameters(1,:);
else
    lls=scalar_vt_vech_likelihood(startingvals,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,false);
    output_parameters=startingvals;
end
