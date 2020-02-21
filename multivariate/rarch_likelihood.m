function [ll, lls, Ht] = rarch_likelihood(parameters,data,p,q,C,backCast,type,isJoint,isCChol)
% Likelihood for RARCH(p,q) multivariate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [LL,LLS,HT] = rarch_likelihood(PARAMETERS,DATA,P,Q,C,BACKCAST,TYPE,ISJOINT)
%
% INPUTS:
%   PARAMETERS   - Vector of parameters required to compute the (negative) of the log-likelihood
%   DATA         - K by K by T array of data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   C            - Unconditional covariance of the data
%   BACKCAST     - K by K matrix to use for back casting
%   TYPE         - Number indicating type: 
%                    1 - Scalar
%                    2 - Common Persistence
%                    3 - Diagonal
%   ISJOINT      - Boolean indicating wether the estimation is joint or not
%
% OUTPUTS:
%   LL           - The log likelihood evaluated at the PARAMETERS
%   LLS          - A T by 1 vector of log-likelihoods
%   HT           - A [K K T] dimension matrix of conditional covariances
%
% COMMENTS:
%   When ISJOINT==1, the first K*(K+1)/2 of the parameters are chol2vec(chol(C)').
%
% See also RARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

% Get the parameters together
T = size(data,3);
k = size(data,2);

[C,A,B] = rarch_parameter_transform(parameters,p,q,k,C,type,isJoint,isCChol);
% Fix B in case of CP
B(B<0)= 0;

Gt = zeros(k,k,T);
e = zeros(k,k,T);
lls = zeros(T,1);
C12 = C^(0.5);
Cm12 = C^(-0.5);
intercept = eye(k) - sum(A.^2,3) - sum(B.^2,3);
dint = diag(intercept);
dint(dint<.000001)=.000001;
intercept = diag(dint);
logLikConst = k*log(2*pi);

for i=1:T
    e(:,:,i) = Cm12 * data(:,:,i) * Cm12;
    Gt(:,:,i+1) = intercept;
    for j=1:p
        if (i-j)<=0
            Gt(:,:,i) = Gt(:,:,i) + A(:,:,j)*backCast*A(:,:,j);
        else
            Gt(:,:,i) = Gt(:,:,i) + A(:,:,j)*e(:,:,i-j)*A(:,:,j);
        end
    end
    for j=1:q
        if (i-j)<=0
            Gt(:,:,i) = Gt(:,:,i) + B(:,:,j)*backCast*B(:,:,j);
        else
            Gt(:,:,i) = Gt(:,:,i) + B(:,:,j)*Gt(:,:,i-j)*B(:,:,j);
        end
    end
    V = C12*Gt(:,:,i)*C12;
    lls(i) = 0.5*(logLikConst + log(det(V)) + sum(diag(V^(-1)*data(:,:,i))));
end
ll = sum(lls);

if isnan(ll) || isinf(ll) || ~isreal(ll)
    ll = 1e7;
end

if nargout>2
    Ht = zeros(k,k,T);
    for i=1:T
        Ht(:,:,i) = C12*Gt(:,:,i)*C12;
    end
end