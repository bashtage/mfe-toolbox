function [c,ceq] = bekk_constraint(parameters,data,dataAsym,p,o,q,backCast,backCastAsym,type) %#ok<*INUSL>

ceq = [];
[C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type);
k = size(C,1);

switch type
    case 1
        c = sum(A(1,1,:).^2,3) + sum(B(1,1,:).^2,3) + 0.5*sum(G(1,1,:).^2,3) - 1;
    case 2
        c = diag(sum(A.^2,3) + sum(B.^2,3) + 0.5*sum(G.^2,3) - 1);
    case 3
        m = zeros(k*k);
        for i=1:p
            m = m + kron(A(:,:,i),A(:,:,i));
        end
        for i=1:o
            m = m + 0.5*kron(G(:,:,i),G(:,:,i));
        end
        for i=1:q
            m = m + kron(B(:,:,i),B(:,:,i));
        end
        c = eig(m) - 1;
end