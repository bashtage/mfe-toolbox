function parameters = heavy(data,p,q)


[T,K] = size(data);
isReturn = any(data<0);
data2 = data.^2;
data2(:,~isReturn) = data(:,~isReturn);



