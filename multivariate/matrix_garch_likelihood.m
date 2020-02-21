function [ll,lls,Ht]=matrix_garch_likelihood(parameters,data,dataAsym,p,o,q,backCast,backCastAsym)
% Log likelihood for MATRIX_GARCH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = matrix_garch_likelihood(PARAMETERS,DATA,DATAASYM,P,O,Q,BACKCAST,BACKCASTASYM)
%
% INPUTS:
%   PARAMETERS   - A vector of vech GARCH process parameters: [alpha beta]' or 
%                    [vech(C)' alpha beta]' if ISJOINT
%   DATA         - K by K by T matrix of covariance innovations
%   DATAASYM     - K by K by T matrix of asymmetric covariance innovations
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of lags of asymmetric process
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   BACKCAST     - Back cast value for starting the recursion
%   BACKCASTASYM - Back cast value (asymetric terms) for starting the recursion
%
% OUTPUTS:
%   LL           - Minus 1 times the log likelihood
%   LLS          - Time series of log likelihoods (Also multiplied by -1)
%   HT           - Time series of conditional covariances
%
% COMMENTS:
%   See also MATRIX_GARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/28/2009

[k,nothing,T]=size(data);
k2 = k*(k+1)/2;
%If nargin~=7, then we don't need to transform the parameters
parameterMatrices= zeros(k,k,1+p+o+q);
index = 0;
for i=1:(1+p+o+q)
    temp = vec2chol(parameters(index+1:index+k2));
    parameterMatrices(:,:,i) = (temp*temp'+temp*temp')/2;
    index=index+k2;
end
%Initialize the covariance
Ht=repmat(backCast,[1 1 T]);

%Initialize the log likelihood
lls=zeros(T,1);

%Compute the likelihood constant.
likconst=k*log(2*pi);

%Perform the recursion
for t=1:T;
    Ht(:,:,t)=parameterMatrices(:,:,1);
    for j=1:p
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,j+1).*backCast;
        else
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,j+1).*data(:,:,t-j);
        end
    end
    for j=1:o
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+j+1).*backCastAsym;
        else
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+j+1).*dataAsym(:,:,t-j);
        end
    end    
    for j=1:q
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+o+j+1).*backCast;
        else
            Ht(:,:,t)=Ht(:,:,t)+parameterMatrices(:,:,p+o+j+1).*Ht(:,:,t-j);
        end
    end
    % Replace these lines to make it work better with poorle conditioned data
    % likelihoods(i)=likconst+(log(det(Ht(:,:,t)))+data(i,:)*Ht(:,:,t)^(-1)*data(i,:)');
    
    %This is a trick to ensure better numerical stability
    Q=sqrt(diag(Ht(:,:,t)));
    R=Ht(:,:,t)./(Q*Q');
    stdresid=data(:,:,t)./(Q*Q');
    lls(t)=0.5*(likconst+2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
end
ll = sum(lls);

% This is a hack since the estimation is unconstrained
if isnan(ll) || isinf(ll) || ll>1e7
    ll = 1e7;
end