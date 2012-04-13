function [parameters, ll ,Ht, VCV, scores, diagnostics]=dcc(data,dataAsym,m,l,n,p,o,q,gjrType,type,composite,startingVals,options)










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
if any(floor(o)~=o) || any(o<0) || numel(o)~=k
    error('All elements in O must be non-negative, and O must be either scalar or K by 1.')
end
if any(floor(q)~=q) || any(q<0) || numel(q)~=k
    error('All elements in Q must be non-negative, and Q must be either scalar or K by 1.')
end
if isempty(gjrType)
    gjrType = 2;
end
if isscalar(gjrType)
    gjrType = ones(k,1)*gjrType;
end
if any(~ismember(gjrType,[1 2])) || numel(gjrType)~=k
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
    stage = 3;
else
    stage = 2;
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
if stage == 2 && composite==1
    warning('oxfordMFE:incorrectOption','When TYPE is ''2-stage'', COMPOSITE must be either ''None'' or ''Full''.')
    composite = 2;
end

if ~isempty(startingVals)
    count = k + sum(p)+ sum(o) + sum(q);
    count = count + m + l + n;
    count = count + k*(k-1)/2;
    if length(startingVals)~=count
        error('STARTINGVALS does not contain the correct number of parameters.')
    end
end

if isempty(options)
    options = optimset('fmincon');
    options.Display = 'iter';
    options.Diagnostics = 'on';
    options.Algorithm = 'interior-point';
    
end
try
    optimset(options);
catch ME
    error('OPTIONS does not appear to be a valid options structure.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate volatility models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = zeros(T,k);
univariate = cell(k,1);
univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';
for i=1:k
    [parameters, ~, ht, ~, ~, scores, diagnostics] = tarch(data2d(:,i),p(i),o(i),q(i), [], gjrType(i), [], univariteOptions);
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
    univariate{i}.scores = scores;
    H(:,i) = ht;
end
stdData = data;
stdDataAsym = dataAsym;
for t=1:T
    h = sqrt(H(t,:));
    hh = h'*h;
    stdData(:,:,t) = stdData(:,:,t)./hh;
    stdDataAsym(:,:,t) = stdDataAsym(:,:,t)./hh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = mean(stdData,3);
r = sqrt(diag(R));
R = R ./ (r*r');
N = mean(stdDataAsym,3);
scale = max(eig(R^(-0.5)*N*R^(-0.5)));
% FIXME
startingVals = [.02 .01/scale .96];

LB = zeros(length(startingVals),1);
UB = ones(length(startingVals),1);
A = ones(1,length(startingVals));
A(m+1:m+l) = scale;
b = .99998;

if (startingVals*A'-b) >= 0
    error('STARTINGVALS for DCC parameters are not comparible with a positive definite intercept.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back casts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = .06*.94.^(0:sqrt(T));
w = w'/sum(w);
backCast = zeros(k);
backCastAsym = zeros(k);
for i=1:length(w)
    backCast = backCast + w(i)*stdData(:,:,i);
    backCastAsym = backCastAsym + w(i)*stdDataAsym(:,:,i);
end
isJoint = false;
isInference = false;
parameters = fmincon(@dcc_likelihood,startingVals,A,b,[],[],LB,UB,[],options,stdData,stdDataAsym,m,l,n,R,N,backCast,backCastAsym,3,composite,isJoint,isInference);

gScale = diag(N);
a = parameters(1:m);
g = parameters(m+1:m+l);
b = parameters(m+l+1:m+l+n);

if stage==2
    intercept = R*(1-sum(a)-sum(b)) - N*sum(g);
    [~,rescaledIntercept] = cov2corr(intercept);
    z = r2z(rescaledIntercept);
    startingVals = [z' parameters];
    
    LB = [-inf*ones(1,k*(k-1)/2) zeros(1,length(parameters))];
    UB = [inf*ones(1,k*(k-1)/2) ones(1,length(parameters))];
    A = [zeros(1,k*(k-1)/2) A];
    b = .99998;
    parameters = fmincon(@dcc_likelihood,startingVals,A,b,[],[],LB,UB,[],options,stdData,stdDataAsym,m,l,n,R,N,backCast,backCastAsym,2,composite,isJoint,isInference,gScale);
    z = parameters(1:k(k-1)/2);
    R = z2r(z);
    parameters = parameters(k(k-1)/2+1:length(parameters));
    parameters  = [corr_vech(R)' parameters];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariances and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
garchParameters = [];
for i=1:k
    garchParameters  = [garchParameters univariate{i}.parameters']; %#ok<AGROW>
end
if stage==3
    parameters = [garchParameters corr_vech(R)' vech(N)' parameters];
elseif stage==2
    parameters = [garchParameters parameters];
end

isJoint = true;
isInference = true;
[ll,~,Rt] = dcc_likelihood(parameters,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
Ht = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    Ht(:,:,t) = Rt(:,:,t).*(h'*h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout<=3
    return
end

v = length(parameters);
A = zeros(v);
scores = zeros(T,v);
offset = 0;

for i=1:k
    u = univariate{i};
    count = 1 + u.p + u.o + u.q;
    ind = offset+(1:count);
    A(ind,ind) = u.A;
    offset = offset + count;
    scores(:,ind) = u.scores;
end

% FIXME : Better gradient function
if stage==2
    % 1. dcc_likelihood
    count = k*(k-1)/2 + m + l + n;
    H = hessian_2sided_nrows(@dcc_likelihood,parameters',count,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    A(offset+(1:count),:) = H/T;
    [g,s]=gradient_2sided(@dcc_likelihood,parameters',data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = cov(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
elseif stage==3
    % 1. dcc_inference_objective
    count = k*(k-1)/2;
    if l>0
        count = count + k*(k+1)/2;
    end
    tempParams = parameters(1:offset+count);
    [g,s]=gradient_2sided(@dcc_inference_objective, tempParams', data,dataAsym,m,l,n,univariate);
    scores(:,offset+(1:count))=s(:,offset+(1:count));
    H = hessian_2sided_nrows(@dcc_inference_objective, tempParams', count, data,dataAsym,m,l,n,univariate);
    A(offset+(1:count),1:(count+offset)) = H/T;
    offset = offset + count;
    % 2. dcc_likelihood
    count = m + l + n;
    H = hessian_2sided_nrows(@dcc_likelihood,parameters',count,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    A(offset+(1:count),:) = H/T;
    [g,s]=gradient_2sided(@dcc_likelihood,parameters',data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = cov(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
end

diagnostics = [];
diagnostics.gScale = gScale;