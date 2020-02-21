function [ll,lls,Ht] = bekk_likelihood(parameters,data,dataAsym,p,o,q,backCast,backCastAsym,type)
% Likelihood for BEKK(p,q) multivarate volatility model estimation
%
% USAGE:
%  [LL,LLS,HT] = bekk_likelihood(PARAMETERS,DATA,P,O,Q,BACKCAST,TYPE)
%
% INPUTS:
%   PARAMETERS   - Vector of parameters required to compute the (negative) of the log-likelihood
%   DATA         - K by K by T array of data
%   DATAASYM     - K by K by T array of asymmetric data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative, scalar integer representing the number of asymmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   BACKCAST     - K by K matrix to use for back casting
%   TYPE         - Number indicating type: 
%                    1 - Scalar
%                    2 - Diagonal
%                    3 - Full
%
% OUTPUTS:
%   LL           - The log likelihood evaluated at the PARAMETERS
%   LLS          - A T by 1 vector of log-likelihoods
%   HT           - A [K K T] dimension matrix of conditional covariances
%
% COMMENTS:
%
% See also BEKK

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

% Get the parameters together
T = size(data,3);
k = size(data,2);

[C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type);

Ht = zeros(k,k,T);
lls = zeros(T,1);
logLikConst = k*log(2*pi);

for i=1:T
    Ht(:,:,i) = C;
    for j=1:p
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*backCast*A(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*data(:,:,i-j)*A(:,:,j);
        end
    end
    for j=1:o
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*backCastAsym*G(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*dataAsym(:,:,i-j)*G(:,:,j);
        end
    end    
    for j=1:q
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*backCast*B(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*Ht(:,:,i-j)*B(:,:,j);
        end
    end
    lls(i) = 0.5*(logLikConst + log(det(Ht(:,:,i))) + sum(diag(Ht(:,:,i)^(-1)*data(:,:,i))));
end
ll = sum(lls);

if isnan(ll) || isinf(ll) || ~isreal(ll)
    ll = 1e7;
end