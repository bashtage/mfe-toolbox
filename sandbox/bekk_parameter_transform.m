function [C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type)

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