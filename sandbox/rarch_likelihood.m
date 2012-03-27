function [ll, lls, Ht] = rarch_likelihood(parameters,data,p,q,v,C,backCast,type,isJoint)
% Likelihood for RARCH(p,q) multivarate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [LL] = rarch_likelihood(DATA,P,Q)
%  [PARAMETERS,LL,HT,VCV,SCORES] = rarch(PARAMETERS,DATA,P,Q,C,BACKCAST,TYPE,ISJOINT)
%
% INPUTS:
%   PARAMETERS   - Vector of parameters required to compute the (negative) of the log-likelihood
%   DATA         - K by K by T array of data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   V            - Number of components allowed to have time-varying conditional covariance
%   C            - Unconditional covariance of the data
%   BACKCAST     - K by K matrix to use for back casting
%   TYPE         - Number indicating type: 
%                    1 - Scalar
%                    2 - Common Persistence
%                    3 - Diagonal
%   ISJOINT      - Boolean indicating wether the estimation is joint or not
%
% OUTPUTS:
%   PARAMETERS   -
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

% Get the parameters together
T = size(data,3);
k = size(data,2);

[C,A,B] = rarch_parameter_transform(parameters,p,q,v,k,C,type,isJoint,cpTransform);

Gt = zeros(k,k,T);
e = zeros(k,k,T);
lls = zeros(k,k,T);
C12 = C^(0.5);
Cm12 = C^(-0.5);
intercept = eye(k) - sum(A.^2,3) - sum(B.^2,3);
logLikConst = log(2*k*pi);
for i=1:T
    e(:,:,i) = Cm12 * data(:,:,i) * Cm12;
    Gt(:,:,i+1) = intercept;
    for j=1:p
        if i-j<=0
            Gt(:,:,i+1) = Gt(:,:,i+1) + A(:,:,j)*backCast*A(:,:,j);
        else
            Gt(:,:,i+1) = Gt(:,:,i+1) + A(:,:,j)*e(:,:,i-j)*A(:,:,j);
        end
    end
    for j=1:q
        if i-j<=0
            Gt(:,:,i+1) = Gt(:,:,i+1) + B(:,:,j)*backCast*B(:,:,j);
        else
            Gt(:,:,i+1) = Gt(:,:,i+1) + B(:,:,j)*Gt(:,:,i+1)*B(:,:,j);
        end
    end
    V = C12*Gt(:,:,i)*C12;
    lls(i) = logLikConst + log(det(V)) + sum(diag(V^(-1)*data(:,:,i)));
end
ll = sum(lls);

if isnan(ll) || isinf(ll)
    ll = 1e7;
end

if nargout>2
    Ht = zeros(k,k,T);
    for i=1:T
        Ht(:,:,i) = C12*Gt(:,:,i)*C12;
    end
end