function [parameters, ll, Ht, VCV, scores] = rarch(data,p,q,type,method,startingVals,options)
% Estimation of RARCH(p,q) multivarate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [PARAMETERS] = rarch(DATA,P,Q)
%  [PARAMETERS,LL,HT,VCV,SCORES] = rarch(DATA,P,Q,TYPE,METHOD,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   TYPE         - [OPTIONAL] String, one of 'Scalar' (Default) ,'CP' (Common Persistence) or 'Diagonal'
%   METHOD       - [OPTIONAL] String, one of '2-stage' (Default) or 'Joint'
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)
%
% OUTPUTS:
%   PARAMETERS   - Estimated parameters in the order: 
%                    Scalar: [vech(C) a(1) ... a(P) b(1) ... b(Q)]
%                    CP: [vech(C) diag(A(:,:,1)) ...  diag(A(:,:,P)) theta]
%                    Diagonal: [vech(C) diag(A(:,:,1)) ...  diag(A(:,:,P)) diag(B(:,:,1)) ...  diag(B(:,:,Q))]
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:
%   The dynamics of a RARCH model are identical to that of a BEKK, except
%   that the model evolves in the rotated space.
%
%   G(:,:,t) = (eye(K) - sum(A.^2,3) - sum(B.^2,3)) +
%       A(:,:,1)*OP(:,:,t-1)*A(:,:,1) + ... A(:,:,p)*OP(:,:,t-1)*A(:,:,p) +
%       B(:,:,1)*G(:,:,t-1)*B(:,:,1) + ... B(:,:,p)*OP(:,:,t-1)*B(:,:,p)
%
%   where in the scalar model A(:,:,i) = a(i)*eye(K), B(:,:,j)=b(j)*eye(K)
%   and in the CP model, B(:,:,j) = theta - sum(A.^2,3).  OP is the outer product of the
%   unconditionally standardized data.
%
%   When ISJOINT=1, the first K*(K+1)/2 element of STARTINGBALS mus tbe vech(C).
%
% EXAMPLES:
%   % Scalar (1,1)
%   parameters = rarch(data,1,1,'Scalar')
%   % Scalar (1,1), jointly estimated
%   parameters = rarch(data,1,1,'Scalar','Joint')
%   % Common persistence (2,1)
%   parameters = rarch(data,2,1,'CP')
%   % Diagons (2,2)
%   parameters = rarch(data,2,2,'Diagonal')

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    error('3 to 8 inputs required.')
end
[T,k] = size(data);
if ndims(data)==3
    [k,~,T] = size(data);
end

switch nargin
    case 3
        type = 'Scalar';
        method = '2-stage';
        startingVals = [];
        options = [];
    case 4
        method = '2-stage';
        startingVals = [];
        options = [];
    case 5
        startingVals = [];
        options = [];
    case 6
        options = [];
    case 7
        % Nothing
    otherwise
        error('3 to 7 inputs required.')
end
if ndims(data)>3
    error('DATA must be either a T by K matrix or a K by K by T array.')
end
if T<=k
    error('DATA must be either a T by K matrix or a K by K by T array, and T must be larger than K.')
end

if p<1 || floor(p)~=p
    error('P must be a positive scalar.')
end
if q<0 || floor(q)~=q
    error('Q must be a non-negative scalar.')
end

if strcmpi(type,'Scalar')
    type = 1;
elseif strcmpi(type,'CP')
    type = 2;
elseif strcmpi(type,'Diagonal')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''CP'' or  ''Diagonal''.')
end

if strcmpi(method,'2-stage')
    isJoint = false;
elseif strcmpi(method,'Joint')
    isJoint = true;
else
    error('METHOD must be either ''2-stage'' or  ''Joint''.')
end

if ~isempty(startingVals)
    switch type
        case 1
            count = p+q;
        case 2
            count = p*k+q;
        case 3
            count = (p+q)*k;
    end
    if isJoint
        count = count + k*(k+1)/2;
    end
    if length(startingVals)~=count
        error('STARTINGVALS does not have the expected number of elements.')
    end
end

if isempty(options)
    options = optimset('fmincon');
    options.Display = 'iter';
    options.Diagnostics = 'on';
    options.Algorithm = 'interior-point';
else
    try
        options = optimset(options);
    catch ME
        error('The user provided options structure is not valid.')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data)==2
    temp = zeros(k,k,T);
    for i=1:T
        temp(:,:,i) = data(i,:)'*data(i,:);
    end
    data = temp;
end
C = mean(data,3);
stdData = zeros(k,k,T);
Cm12 = C^(-0.5);
for i=1:T
    stdData(:,:,i) = Cm12*data(:,:,i)*Cm12;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = .06*.94.^(0:ceil(sqrt(T)));
w = w/sum(w);
backCast = zeros(k);
for i=1:length(w)
    backCast = backCast + w(i)*stdData(:,:,i);
end

k2 = k*(k+1)/2;
if isempty(startingVals)
    as = .02:.03:.11;
    apbs = .9:.03:.99;
    startingLLs = zeros(length(as),length(apbs));
    for i=1:length(as)
        for j=1:length(apbs)
            a = as(i);
            b = apbs(j) - a;
            parameters = sqrt([a b]);
            startingLLs(i,j) = rarch_likelihood(parameters,data,1,1,C,backCast,1,false,false);
        end
    end
else
    % Ignore C in the starting vals
    startingVals = startingVals(k2+1:length(startingVals));
end

[i,j]=find(startingLLs==min(min(startingLLs)));
a = as(i);
b = apbs(j)-a;
switch type
    case 1
        startingVals = sqrt([a/p * ones(1,p) b/q*ones(1,q)]);
    case 2
        startingVals = sqrt([a/p * ones(1,k*p) (a+b)]);
    case 3
        startingVals = sqrt([a/p * ones(1,k*p) b/q*ones(1,k*q)]);
end
startingVals = startingVals';
UB = .99998 * ones(size(startingVals));
LB = -UB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = fmincon(@rarch_likelihood,startingVals,[],[],[],[],LB,UB,@rarch_constraint,options,data,p,q,C,backCast,type,false,false);
[ll,~,Ht] = rarch_likelihood(parameters,data,p,q,C,backCast,type,false,false);
if isJoint
    CChol = chol2vec(chol(C)');
    startingValJoint =[CChol;parameters];
    LBJoint = [-inf*ones(size(CChol));-ones(size(parameters))];
    UBJoint = abs(LBJoint);
    parameters = fmincon(@rarch_likelihood,startingValJoint,[],[],[],[],LBJoint,UBJoint,@rarch_constraint,options,data,p,q,C,backCast,type,true,true);
    [ll,~,Ht] = rarch_likelihood(parameters,data,p,q,C,backCast,type,true,true);
end
ll = -ll;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isJoint
    C = parameters(1:k2);
    C = vec2chol(C);
    C = C*C';
    C = (C+C')/2;
    parameters = [vech(C);parameters(k2+1:length(parameters))];
    [VCV,~,~,scores] = robustvcv(@rarch_likelihood,parameters,0,data,p,q,C,backCast,type,true,false);
else
    scores1 = zeros(T,k2);
    for i=1:T
        scores1(i,:) = vech(data(:,:,i)-C)';
    end
    [~,scores2] = gradient_2sided(@rarch_likelihood,parameters,data,p,q,C,backCast,type,false,false);
    scores = [scores1 scores2];
    B = covnw(scores,ceil(1.2*T^(0.25)));
    m = length(parameters);
    parameters = [vech(C);parameters];
    A1 = -eye(k2);
    A2 = hessian_2sided_nrows(@rarch_likelihood,parameters,m,data,p,q,C,backCast,type,true,false);
    A2 = A2/T;
    A = [[A1 zeros(k2,m)];
        A2];
    Ainv = A\eye(length(A));
    VCV = Ainv*B*Ainv'/T;
end
