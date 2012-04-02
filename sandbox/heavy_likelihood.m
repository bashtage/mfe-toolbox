function [ll,lls,h] = heavy_likelihood(parameters,data,p,q,pMax,qMax,backCast)

[T,K] = size(data);
lls = zeros(T,1);
h = zeros(K,T);
[T,K] = size(data);
[O,A,B] = heavy_parameter_transform(parameters,p,q,pMax,qMax);
for t=1:T
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
    lls(t)  = 0.5*(K*log(2*pi) + sum(log(h(:,t))) + sum(data(t,:)'/h(:,t)));
end

ll = sum(lls);
if nargout>2
    h = h';
end