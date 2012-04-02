function [data,h] = heavy_simulate(T,K,parameters,p,q,m)

e = zeros(T,K);
data = zeros(T,K);
signs = zeros(T,K);
for i=1:K
    temp = randn(T,m(i));
    signs(:,i) = sign(sum(temp,2));
    e(:,i) = sum(temp.^2,2);
end


[O,A,B] = heavy_parameter_transform(parameters,p,q);

h = ones(K,T);
backCast = ones(1,K);
for t=1:T
    data(i,:) = e(i,:).*sqrt(h(:,t)');
    h(:,t)= O;
    for j=1:pMax
        if (t-j)<=0
            h(:,t) = h(:,t) + A(:,:,j)*data(t-j,:)';
        else
            h(:,t) = h(:,t) + A(:,:,j)*backCast';
        end
    end
    for j=1:qMax
        if (t-j)<=0
            h(:,t) = h(:,t) + B(:,:,j)*h(t-j,:)';
        else
            h(:,t) = h(:,t) + B(:,:,j)*backCast';
        end
    end
end

data(:,m==1) = signs(:,m==1).*sqrt(data(:,m==1));