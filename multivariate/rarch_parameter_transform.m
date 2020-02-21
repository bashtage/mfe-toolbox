function [C,A,B] = rarch_parameter_transform(parameters,p,q,k,C,type,isJoint,isCChol)
% Parameter transformation for RARCH(p,q) multivariate volatility model simulation and estimation
%
% USAGE:
%  [C,A,B] = rarch_parameter_transform(PARAMETERS,P,Q,K,C,TYPE,ISJOINT)
%
% INPUTS:
%   PARAMETERS - Vector of parameters governing the dynamics.  See RARCH or RARCH_SIMULATE
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance
%                  lags.  When using 'CP' model, 0<=q<=1
%   K          - Number of assets
%   C          - Unconditional covariance of data. Ignored if ISJOINT==1
%   TYPE       - Integer indicating type
%                  1: Scalar
%                  2: CP
%                  3: Diagonal
%   ISJOINT    - Boolean indicating whether PARAMETER contains only the dynamics (ISJOINT==false) or
%                  both the Cholesky of the unconditional and the dynamic parameters
% OUTPUTS:
%   C - K by K unconditional covariance matrix
%   A - K by K by P matrix of innovation parameters
%   B - K by K by Q matrix of smoothing parameters
%
% COMMENTS:
%   When ISJOINT==1, the first K*(K+1)/2 of the parameters are chol2vec(chol(C)').
%
% See also RARCH, RARCH_SIMULATE, RARCH_LIKELIHOOD

% Scalar defaults
numA = 1;
numB = 1;
if type>=2 % Common Persistence
    numA = k;
end
if type==3 % Diagonal
    numB = k;
end

parameterCount = 0;
if isJoint
    k2 = k*(k+1)/2;
    C = parameters(1:k2);
    if isCChol
        C = vec2chol(C);
        C = C*C';
        % Ensure C is PSD
        [V,D] = eig(C);
        D = diag(D);
        if (min(D))<(2*eps*max(D))
            D((D/max(D))<eps) = 2*max(D)*eps;
            C = V*diag(D)*V';
            C=(C+C)/2;
        end
    else
        C = ivech(C);
    end
    parameterCount = parameterCount + k2;
end
A = zeros(k,k,p);
for i=1:p
    temp = parameters(parameterCount+1:parameterCount+numA);
    parameterCount = parameterCount + numA;
    if numA == 1
        A(:,:,i) = temp * eye(k);
    else
        A(:,:,i) = diag(temp);
    end
end
B = zeros(k,k,q);
for i=1:q
    temp = parameters(parameterCount+1:parameterCount+numB);
    parameterCount = parameterCount + numB;
    if numB == 1
        B(:,:,i) = temp * eye(k);
    else
        B(:,:,i) = diag(temp);
    end
end

% Common Persistence uses a different parameterization
if type==2 && q>0
    theta = B(1,1,1).^2;
    B = theta - diag(sum(A.^2,3));
    B(B<0)=0;
    B = diag(sqrt(B));
end
