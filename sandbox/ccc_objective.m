function [obj,objs] = ccc_objective(parameters,data,univariate)


[k,~,T] = size(data);
offset = 0;
H = zeros(T,k);

for i=1:k
    u = univariate{i};
    count = u.p+u.o+u.q+1;
    parameters = parameters(offset + (1:count));
    offset = offset+count;
    ht = tarch_core(u.fdata,u.fIdata,parameters,u.back_cast,u.p,u.o,u.q,u.m,u.T,u.tarch_type);
    H(:,i) = ht(m+1:T);
end
stdData = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    stdData(:,:,t) = data./(h*h');
end
z = parameters(offset + (1:k*(k-1)/2));
R = corr_ivech(z);

scales = diag(mean(stdData,3));
objs = zeros(T,1);
for i=2:k
    for j=1:i-1
        scale = sqrt(scales(i)*scales(j));
        errors = squeeze(stdData(i,j,:))/scale - R(i,j);
        objs = objs + 0.5*(errors.^2);
    end
end
obj = sum(objs);

