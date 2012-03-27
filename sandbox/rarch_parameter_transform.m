function [C,A,B] = rarch_parameter_transform(parameters,p,q,v,k,C,type,isJoint,cpTransform)

% Scalar defaults
numA = 1;
numB = 1;
if type>=2 % Common Persistence
    numA = v;
elseif type==3 % Diagonal
    numB = v;
end

parameterCount = 0;
if isJoint
    k2 = k*(k+1)/2;
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
        A(1:v,1:v,i) = temp * eye(v);
    else
        A(1:v,1:v,i) = diag(temp);
    end
end
B = zeros(k,k,q);
for i=1:q
    temp = parameters(parameterCount+1:parameterCount+numB);
    parameterCount = parameterCount + numB;
    if numB == 1
        B(1:v,1:v,i) = temp * eye(v);
    else
        B(1:v,1:v,i) = diag(temp);
    end
end

% Common Persistence uses a different parameterization
if type==2 && cpTransform
    % B is really A+B
    % A is the % of B of A
    theta = B(1,1,1).^2;
    A = A*theta;
    B = theta - sum(A.^2,3);
elseif type==2 
    theta = B(1,1,1).^2;
    B = theta - sum(A.^2);
end
