function [parameters,LL,Ht,VCVrobust,likelihoods,scores,diagnostics]=dcc(data,dataAsym,m,l,n,p,o,q,gjrType,type,composite,startingVals,options)










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        p = [];
        o = [];
        q = [];
        gjrType = [];
        type = [];
        composite = [];
        startingVals = [];
        options = [];
    case 6
        o = [];
        q = [];
        gjrType = [];
        type = [];
        composite = [];
        startingVals = [];
        options = [];
    case 7
        q = [];
        gjrType = [];
        type = [];
        composite = [];
        startingVals = [];
        options = [];
    case 8
        gjrType = [];
        type = [];
        composite = [];
        startingVals = [];
        options = [];
    case 9
        type = [];
        composite = [];
        startingVals = [];
        options = [];
    case 10
        composite = [];
        startingVals = [];
        options = [];
    case 11
        startingVals = [];
        options = [];
    case 12
        options = [];
    case 13
        % Nothing
    otherwise
        error('5 to 13 inputs required.')
end

if ndims(data)==2
    [T,k] = size(data);
    data2d = data;
    eta = data.*(data<0);
    data = zeros(k,k,T);
    dataAsym = zeros(k,k,T);
    for t=1:T
        data(:,:,t) = data2d(t,:)'*data2d(t,:);
        dataAsym(:,:,t) = eta(t,:)'*eta(t,:);
    end
elseif ndims(data)==3
    [k,~,T]=size(data);
    data2d = zeros(K,K,T);
    for i=1:K
        data2d(:,i) = squeeze((1-2*(dataAsym(i,i,:)==0)) .* sqrt(data(i,i,:)));
    end
else
    error('DATA must be either a K by T matrix of a K by K by T 3-dimensional array.')
end

if ~isempty(dataAsym)
    if ndims(dataAsym)~=3
        error('DATAASYM must be a K by K by T array.')
    end
    [k1,k2,T2] = size(dataAsym);
    if k1~=k || k2~=k || T2~=T
        error('DATAASYM must be a K by K by T array.')
    end
end

if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
if isempty(l)
    l = 0;
end
if ~isscalar(l) || floor(l)~=l || l<0
    error('L must be a non-negative integer.')
end
if isempty(n)
    n = 0;
end
if ~isscalar(n) || floor(n)~=n || n<0
    error('N must be a non-negative integer.')
end

if isempty(p)
    p = 1;
end
if isscalar(p)
    p = p * ones(k,1);
end
if isempty(o)
    o = 0;
end
if isscalar(o)
    o = o * ones(k,1);
end
if isempty(q)
    q = 0;
end
if isscalar(q)
    q = q * ones(k,1);
end
if any(floor(p)~=p) || any(p<1) || numel(p)~=k
    error('All elements in P must be positive, and P must be either scalar or K by 1')
end
if any(floor(o)~=o) || any(o<1) || numel(o)~=k
    error('All elements in O must be non-negative, and O must be either scalar or K by 1.')
end
if any(floor(q)~=q) || any(q<1) || numel(q)~=k
    error('All elements in Q must be non-negative, and Q must be either scalar or K by 1.')
end
if isempty(gjrType)
    gjrType = 2;
end
if isscalar(gjrType)
    gjrType = ones(k,1)*gjrType;
end
if any(~ismembber(gjrType,[1 2])) || numel(gjrType)~=k
    error('GJRTYPE must be in {1,2} and mustbe either scalar of K by 1.')
end

if isempty(type)
    type = '3-stage';
end
type = lower(type);
if ~ismember(type,{'3-stage','2-stage'})
    error('TYPE must be either ''3-stage'' or ''2-stage''.')
end
if strcmpi(type,'3-stage')
   type = 3;
else
    type = 2;
end


if isempty(composite)
    composite = 'none';
end
composite = lower(composite);
if ~ismember(composite,{'none','diagonal','full'})
    error('COMPOSITE must be one of ''None'', ''Diagonal'' or ''Full''.')
end
if strcmpi(composite,'none')
   composite = 0;
elseif strcmpi(composite,'diagonal')
    composite = 1;
else
    composite = 2;
end
if type == 2 && composite==1
    warning('oxfordMFE:incorrectOption','When TYPE is ''2-stage'', COMPOSITE must be either ''None'' or ''Full''.')
    composite = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


