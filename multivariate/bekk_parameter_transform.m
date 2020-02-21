function [C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type)
% Parameter transformation for BEKK(p,o,q) multivariate volatility model simulation and estimation
%
% USAGE:
%  [C,A,G,B] = bekk_parameter_transform(PARAMETERS,P,O,Q,K,TYPE)
%
% INPUTS:
%   PARAMETERS - Vector of parameters governing the dynamics.  See BEKK or BEKK_SIMULATE
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   O          - Non-negative, scalar integer representing the number of asymmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance lags
%   K          - Number of assets
%   TYPE       - Integer indicating type
%                  1: Scalar
%                  2: Diagonal
%                  3: Full
%
% OUTPUTS:
%   C - K by K covariance model intercept
%   A - K by K by P matrix of symmetric innovation parameters
%   G - K by K by O matrix of asymmetric innovation parameters
%   B - K by K by Q matrix of smoothing parameters
%
% COMMENTS:
%
% See also BEKK, BEKK_SIMULATE, BEKK_LIKELIHOOD

if type==1
    numParams = 1;
elseif type==2
    numParams = k;
else
    numParams = k*k;
end

k2 = k*(k+1)/2;
C = parameters(1:k2);
C = vec2chol(C);
C = C*C';
offset = k2;
[V,D] = eig(C);
D = diag(D);
if (min(D))<(2*eps*max(D))
    D((D/max(D))<eps) = 2*max(D)*eps;
    C = V*diag(D)*V';
    C=(C+C)/2;
end

m = p+o+q;
temp = zeros(k,k,m);
for j=1:m
    tempP = parameters(offset+(1:numParams));
    offset = offset+numParams;
    if type==1
        temp(:,:,j) = tempP*eye(k);
    elseif type==2
        temp(:,:,j) = diag(tempP);
    else
        temp(:,:,j) = reshape(tempP,k,k);
    end
end
A = temp(:,:,1:p);
G = temp(:,:,p+1:p+o);
B = temp(:,:,p+o+1:p+o+q);