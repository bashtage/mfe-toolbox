function [O,A,B] = heavy_parameter_transform(parameters,p,q,K)

O = zeros(K,1);
O(:) = parameters(1:K);
offset = K;
A = zeros(K,K,max(max(p)));
B = zeros(K,K,max(max(q)));
for i=1:K
    for j=1:K
        A(i,j,1:p(i,j)) = parameters(offset+(1:p(i,j)));
    end
end
for i=1:K
    for j=1:K
        B(i,j,1:q(i,j)) = parameters(offset+(1:q(i,j)));
    end
end
