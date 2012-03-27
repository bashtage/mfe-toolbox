function [c,ceq] = rarch_constraint(parameters,data,p,q,C,backCast,type,isJoint)

ceq = [];



% Get the parameters together
T = size(data,3);
k = size(data,2);
k2 = k*(k+1)/2;
% Scalar defaults
numA = 1;
numB = 1;
if type>=2 % Common Persistence
    numA = k;
elseif type==3 % Diagonal
    numB = k;
end

parameterCount = 0;
if isJoint
    C = parameters(1:k2);
    C = ivech(C);
    C = C*C';
    % Ensure C is PSD
    [V,D] = eig(C);
    if (min(D)/max(D))<eps
        D((D/max(D))<eps) = 2*eps;
        C = V*D*V';
        C=(C+C)/2;
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
if type==2
    % B is really A+B
    % A is the if B of A
    theta = B(1,1,1);
    A = sqrt(A*theta);
    B = theta - sum(A.^2,3);
end

c = diag(sum(A.^2,3)+sum(B.^2,3) - 1);
