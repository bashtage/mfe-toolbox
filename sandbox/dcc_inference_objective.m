function [obj,objs] = dcc_inference_objective(parameters,data,dataAsym,m,l,n,univariate)
% "Objective" function which has R and N as solution to problem

[k,~,T] = size(data);
% Parse parameters
count = 0;
for i=1:k
    u = univariate{i};
    count = count + u.p+u.o+u.q+1;
end
garchParameters = parameters(1:count);
offset = count;
% R is next
count = k*(k-1)/2;
R = corr_ivech(parameters(offset + (1:count)));
offset = offset + count;
if l>0
    count = k*(k+1)/2;
    N = ivech(parameters(offset+(1:count)));
end

offset = 0;
H = ones(T,k);
for i=1:k
    u = univariate{i};
    count = u.p+u.o+u.q+1;
    volParameters = garchParameters(offset + (1:count));
    offset = offset+count;
    ht = tarch_core(u.fdata,u.fIdata,volParameters,u.back_cast,u.p,u.o,u.q,u.m,u.T,u.tarch_type);
    H(:,i) = ht(u.m+1:u.T);
end
stdData = zeros(k,k,T);
stdDataAsym = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    stdData(:,:,t) = data(:,:,t)./(h'*h);
    stdDataAsym(:,:,t) = dataAsym(:,:,t)./(h'*h);
end

scales = diag(mean(stdData,3));
objs = zeros(T,1);
for j=1:k-1 % Cols
    for i=j+1:k % Rows
        scale = sqrt(scales(i)*scales(j));
        errors = squeeze(stdData(i,j,:))/scale - R(i,j);
        objs = objs + 0.5*(errors.^2);
    end
end

if l>0
    for j=1:k
        for i= j:k
            errors = squeeze(stdDataAsym(i,j,:)) - N(i,j);
            objs = objs + 0.5*(errors.^2);
        end
    end
end

obj = sum(objs);

