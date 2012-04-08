function [parameters,LL,Ht,VCVrobust,likelihoods,scores,diagnostics]=dcc_mvgarch(data,dataAsym,m,l,n,p,o,q,gjrType,type,composite,startingVals,options)










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data)==2
    data2d = data;
    eta = data.*(data<0);
    data = zeros(k,k,T);
    dataAsym = zeros(k,k,T);
    for t=1:T
        data(:,:,t) = data2d(t,:)'*data2d(t,:);
        dataAsym(:,:,t) = eta(t,:)'*eta(t,:);
    end
else
    data2d = zeros(K,K,T);
    for i=1:K
        data2d(:,i) = squeeze((1-2*(dataAsym(i,i,:)==0)) .* sqrt(data(i,i,:)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate volatility models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = zeros(T,k);
univariate = cell(K,1);
for i=1:K
    [parameters, ~, ht, ~, ~, diagnostics] = tarch(data2d(:,i),p(i),o(i),q(i), [], gjrType(i), [], univariteOptions);
    % Store output for later use
    univariate{i}.p = p(i);
    univariate{i}.o = o(i);
    univariate{i}.q = q(i);
    univariate{i}.fdata = diagnostics.fdata;
    univariate{i}.fIdata = diagnostics.fIdata;
    univariate{i}.back_cast = diagnostics.back_cast;
    univariate{i}.m = diagnostics.m;
    univariate{i}.T = diagnostics.T;
    univariate{i}.tarch_type = gjrType(i);
    univariate{i}.parameters = parameters;
    univariate{i}.ht = ht;
    univariate{i}.A = diagnostics.A;
    H(:,i) = ht;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LB = zeros(length(startingVals));
UB = zeros(length(startingVals));
A = ones(1,length(startingVals));
scale = min(1./eig(R^(-0.5)*N*R^(-0.5)));
A(p+1:p+o) = scale;
b = 1;
parameters = fmincon(@dcc_likelihood,startingVals,A,b,[],[],LB,UB,[],options,stdData,stdDataAsym,m,l,n,R,N,backCase,backCastAsym,3,false,univariate);
if type==2
    % Prepend r2z of the correlation intercept
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


