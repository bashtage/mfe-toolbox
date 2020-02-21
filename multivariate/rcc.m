function [parameters, ll ,Ht, VCV, scores, diagnostics]=rcc(data,dataAsym,m,n,p,o,q,gjrType,type,method,composite,startingVals,options)
% Estimation of scalar RCC(m,n) multivarate volatility model with with TARCH(p,o,q) or 
% GJRGARCH(p,o,q) conditional variances
%
% USAGE:
%  [PARAMETERS] = rcc(DATA,[],M,N)
%  [PARAMETERS] = rcc(DATA,DATAASYM,M,N,P,O,Q,GJRTYPE,TYPE,METHOD,COMPOSITE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
%   M            - Order of symmetric innovations in RCC model
%   N            - Order of lagged correlation in RCC model
%   P            - [OPTIONAL] Positive, scalar integer representing the number of symmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   O            - [OPTIONAL] Non-negative, scalar integer representing the number of asymmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 0.
%   Q            - [OPTIONAL] Non-negative, scalar integer representing the number of conditional covariance lags in
%                    the univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   GJRTYPE      - [OPTIONAL] Either 1 (TARCH/AVGARCH) or 2 (GJR-GARCH/GARCH/ARCH). Can also be a K by 1 vector
%                    containing the model type for each for each series. Default is 2.
%   TYPE         - [OPTIONAL] String, one of 'Scalar' (Default) ,'CP' (Common Persistence) or 'Diagonal'
%   METHOD       - [OPTIONAL] String, one of '3-stage' (Default) or '2-stage'.  Determines whether
%                    the model is estimated using the 3-stage estimator, or if the correlation intercepts
%                    are jointly estimated along with the dynamic parameters.
%   COMPOSITE    - [OPTIONAL] String value, either 'None' (Default), 'Diagonal' or 'Full'.  None
%                    uses standard QMLE.  'Diagonal' and 'Full' both uses composite likelihood where
%                    'Diagonal' uses all pairs of the form i,i+1 while 'Full' uses all pairs.
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)
%
% OUTPUTS:
%   PARAMETERS   - Estimated parameters.  Output depends on METHOD.
%                    Scalar: [VOL(1) ... VOL(K) corr_vech(R)' alpha  beta]
%                    CP: [VOL(1) ... VOL(K) corr_vech(R)' diag(A(:,:,1)) ... diag(A(:,:,m)) theta]
%                    Diagonal: [VOL(1) ... VOL(K) corr_vech(R)' diag(A(:,:,1)) ... diag(A(:,:,m)) diag(B(:,:,1)) ... diag(B(:,:,m))]
%                    where VOL(j) is a (1+P(i)+O(i)+Q(i)) vector containing the parameters from
%                    volatility model i.
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%   DIAGNOSTICS  - A structure containing further outputs.
%
% COMMENTS:
%
% EXAMPLES:
%   % RCC(1,1)
%   parameters = rcc(data,[],1,1)
%   % RCC(1,1), 2-stage
%   parameters = rcc(data,[],1,1,[],[],[],[],[],'2-stage')
%
% See also DCC, CCC_MVGARCH, BEKK, RARCH, SCALAR_VT_VECH, MATRIX_GARCH, TARCH


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 4
        p = [];
        o = [];
        q = [];
        gjrType = [];
        type = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
    case 5
        o = [];
        q = [];
        gjrType = [];
        type = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
    case 6
        q = [];
        gjrType = [];
        type = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
    case 7
        gjrType = [];
        type = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
    case 8
        type = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
    case 9
        method = [];
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
        error('4 to 13 inputs required.')
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
    data2d = zeros(k,k,T);
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
    q = 1;
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
    error('GJRTYPE must be in {1,2} and must be either scalar of K by 1.')
end

if isempty(type)
    type = 'scalar';
end
type = lower(type);
if ~ismember(type,{'scalar','cp','diagonal'})
    error('TYPE must be one of ''Scalar'', ''CP'' or ''Diagonal''');
end
switch type
    case 'scalar'
        type = 1;
    case 'cp'
        type = 2;
    case 'diagonal'
        type = 3;
end

if type==2 && n~=1
    warning('oxfordMFE:IncompatibleParameterization','When using ''CP'', N must always equal 1.  Setting N to 1');
    n = 1;
end

if isempty(method)
    method = '3-stage';
end
method = lower(method);
if ~ismember(method,{'3-stage','2-stage'})
    error('TYPE must be either ''3-stage'' or ''2-stage''.')
end
if strcmpi(method,'3-stage')
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
    if type==1
        count = count + m + n;
    elseif type==2
        count = count + m*k + n;
    else 
        count = count + (m+n)*k;
    end
    count = count + k*(k-1)/2;
    if length(startingVals)~=count
        error('STARTINGVALS does not contain the correct number of parameters.')
    end
    tarchStartingVals = startingVals(1:k+sum(p)+sum(o)+sum(q));
    offset = k+sum(p)+sum(o)+sum(q) + k*(k-1)/2;
    rccStartingVals= startingVals(offset + (1:(m+n)));
else
    tarchStartingVals = [];
    rccStartingVals = [];
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
[H,univariate] = dcc_fit_variance(data2d,p,o,q,gjrType,tarchStartingVals);
stdData = data;
for t=1:T
    h = sqrt(H(t,:));
    hh = h'*h;
    stdData(:,:,t) = stdData(:,:,t)./hh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = mean(stdData,3);
Rm12 = R^(-0.5);
rScale = sqrt(diag(R));
R = R ./ (rScale*rScale');
% TODO : Better starting values
w = .06*.94.^(0:sqrt(T));
w = w'/sum(w);
backCast = zeros(k);
for i=1:length(w)
    backCast = backCast + w(i) * (Rm12*stdData(:,:,i)*Rm12);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(rccStartingVals)
    startingVals = rccStartingVals;
else
    a = [.01 .03 .05 .1];
    theta = [.99 .97 .95];
    [a,theta] = ndgrid(a,theta);
    parameters = sqrt(unique([a(:) theta(:)-a(:)],'rows'));
    minLL= inf;
    isJoint = false;
    isInference = false;
    for i=1:size(parameters,1)
        ll = rcc_likelihood(parameters(i,:),stdData,1,1,R,backCast,3,1,composite,isJoint,isInference,rScale,univariate);
        if ll<minLL
            startingVals = parameters(i,:);
            minLL = ll;
        end
    end
    a = startingVals(1).^2;
    b = startingVals(2).^2;
    startingVals = sqrt([a*ones(1,m)/m b*ones(1,n)/n]);
end
if size(startingVals,1)>size(startingVals,2)
    startingVals = startingVals';
end
UB = ones(size(startingVals));
LB = -UB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back casts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isJoint = false;
isInference = false;
parameters = fmincon(@rcc_likelihood,startingVals,[],[],[],[],LB,UB,@rcc_constraint,options,stdData,m,n,R,backCast,3,1,composite,isJoint,isInference,rScale,univariate);
if type==2 
    A = ones(k,1)*parameters(1:m);
    theta = sqrt(sum(parameters.^2));
    startingVals = [A(:)' theta];
    UB = ones(size(startingVals));
    LB = -UB;
elseif type ==3
    % Use parameters as starting vals
    A = ones(k,1)*parameters(1:m);
    B = ones(k,1)*parameters(m+(1:n));
    startingVals = [A(:)' B(:)'];
    UB = ones(size(startingVals));
    LB = -UB;
end
if type>1
    parameters = fmincon(@rcc_likelihood,startingVals,[],[],[],[],LB,UB,@rcc_constraint,options,stdData,m,n,R,backCast,3,type,composite,isJoint,isInference,rScale,univariate);
end

if stage==2
    z = r2z(R);
    startingVals = [z' parameters];
    
    LB = [-inf*ones(1,k*(k-1)/2) -ones(size(parameters))];
    UB = [inf*ones(1,k*(k-1)/2) ones(size(parameters))];
    parameters = fmincon(@rcc_likelihood,startingVals,[],[],[],[],LB,UB,@rcc_constraint,options,stdData,m,n,R,backCast,2,type,composite,isJoint,isInference,rScale,univariate);
    
    z = parameters(1:k*(k-1)/2);
    R = z2r(z);
    parameters = parameters(k*(k-1)/2+1:length(parameters));
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
    parameters = [garchParameters corr_vech(R)' parameters];
elseif stage==2
    parameters = [garchParameters parameters];
end

isJoint = true;
isInference = true;
[ll,~,Rt] = rcc_likelihood(parameters,data,m,n,R,backCast,2,type,composite,isJoint,isInference,rScale,univariate);
ll = -ll;
Ht = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    Ht(:,:,t) = Rt(:,:,t).*(h'*h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = parameters';
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

% TODO : Better gradient function
if stage==2
    % 1. rcc_likelihood
    count = k*(k-1)/2 + m + n;
    H = hessian_2sided_nrows(@rcc_likelihood,parameters,count,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@rcc_likelihood,parameters,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = cov(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
elseif stage==3
    % 1. dcc_inference_objective
    count = k*(k-1)/2;
    tempParams = parameters(1:(offset+count));
    [~,s]=gradient_2sided(@dcc_inference_objective,tempParams,data,dataAsym,m,0,n,univariate);
    scores(:,offset+(1:count))=s(:,offset+(1:count));
    H = hessian_2sided_nrows(@dcc_inference_objective,tempParams,count,data,dataAsym,m,0,n,univariate);
    A(offset+(1:count),1:(count+offset)) = H/T;
    offset = offset + count;
    % 2. rcc_likelihood
    count = m + n;
    H = hessian_2sided_nrows(@rcc_likelihood,parameters,count,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@rcc_likelihood,parameters,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = covnw(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
end

diagnostics = [];
diagnostics.rScale = rScale;