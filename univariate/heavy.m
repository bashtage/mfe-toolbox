function [parameters, ll, ht, VCV, scores] = heavy(data,p,q,cons,startingVals,options)
% Estimation of HEAVY volatility model of Shephard and Sheppard.  Also estimates general volatility
% spillover models with 2 or more dimensions.
%  
% USAGE:
%  [PARAMETERS,LL,HT,VCV,SCORES] = heavy(DATA,P,Q,CONS,STARTINGVALS,OPTIONS)
%  
% INPUTS:
%   DATA         - T by K matrix of input data.  Data can be either returns or realized-measure type
%                    data.  Returns are detected by examining a series for negative values, and are
%                    squared for estimation of the model.
%   P            - A K by K matrix containing the lag length of model innovations.  Position (i,j)
%                    indicates the number of lags of series j in the model for series i
%   Q            - A K by K matrix containing the lag length of conditional variances.  Position (i,j)
%                    indicates the number of lags of series j in the model for series i
%   CONS         - [OPTIONAL] String indicating the type of constraint to use on parameters:
%                    'None' - Only constrain the intercepts to be positive
%                    'Positive' - Restrict all parameters to be non-negative
%   STARTINGVALS - [OPTIONAL] A number of parameters by 1 vector of starting values.  See COMMENTS.
%   OPTIONS      - [OPTIONAL] Optimization option structure (fmincon)
%  
% OUTPUTS:
%   PARAMETERS   - A sum(sum(P)) + sum(sum(Q)) + K by 1 vector of estimated parameters. See COMMENTS.
%   LL           - The log likelihood at the optimum.
%   HT           - A T by K matrix of conditional variances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%  
% COMMENTS:
%   Dynamics are given by:
%     h(t,:)' = O + A(:,:,1)*f(data(t-1,:))' + ... + A(:,:,maxP)*f(data(t-maxP,:))' + ...
%                 + B(:,:,1)*h(t-1,:)' + ... + B(:,:,maxQ)*h(t-maxQ,:)'
%  
%   PARAMETERS are ordered:
%   [O' A(1,1,1:p(1,1)) A(1,2,1:p(1,2)) ... A(1,K,1:p(1,K)) A(2,1,1:p(2,1)) ... A(2,K,1:p(2,K)) 
%       ... A(K,1,1:p(K,1)) ... A(K,K,1:p(K,K)) ... B(1,1,1:q(1,1)) ... B(1,K,1:q(1,K)) 
%       ... B(K,1,1:q(K,1)) ... B(K,K,1:q(K,K)) ]
%  
% EXAMPLES:
%   % Standard HEAVY model
%   p = [0 1;0 1]
%   q = eye(2)
%   data = [r rm];
%   parameters = heavy(data,p,q,'None')
%  
% See also HEAVY_SIMULATE, TARCH, EGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/2/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        cons = [];
        startingVals = [];
        options = [];
    case 4
        startingVals = [];
        options = [];
    case 5
        options = [];
    case 6
        % Nothing
    otherwise
        error('3 to 5 inputs required.')
end

if size(data,2)>size(data,1)
    data = data';
end
[T,K] = size(data);

if any(size(p)~=K) || any(any(floor(p)~=p)) || any(any(p<0))
    error('P must be a K by K matrix of non-negative integers.')
end

if any(size(q)~=K) || any(any(floor(q)~=q)) || any(any(q<0))
    error('Q must be a K by K matrix of non-negative integers.')
end

if isempty(cons)
    cons = 'None';
end
switch lower(cons)
    case 'none'
        consMode = 0;
    case 'positive'
        consMode = 1;
    otherwise
        error('Unknown values for CONS')
end

parameterCount = K + sum(sum(p)) + sum(sum(q));
if ~isempty(startingVals)
    if length(startingVals)~=parameterCount
        error('Incorrect number of parameters in STARTINGVALS.')
    end
    if startingVals(1:K)<=0
        error('Initial values for O in STARTINGVALS must be strictly positive.')
    end
end
if ~isempty(options)
    try 
        optimset(options);
    catch ME
        error('OPTIONS does not appear to be a valid options structure.')
    end
else
    options = optimset('fmincon');
    options.Display='iter';
    options.Diagnostics = 'on';
    options.Algorithm = 'interior-point';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isReturn = any(data<0);
data2 = data.^2;
data2(:,~isReturn) = data(:,~isReturn);
scale = mean(data2);
data2 = bsxfun(@times,data2,1./scale);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strategy is very naive
% More sp=ophisticated would need some sort of ARCH-X or GLS regression
if isempty(startingVals)
    startingVals = mean(data2)*.05;
    for i=1:K
        As = sum(p(i,:));
        startingVals  = [startingVals .1*ones(1,As)/As]; %#ok<AGROW>
    end
    for i=1:K
        Bs = sum(p(i,:));
        startingVals  = [startingVals .8*ones(1,Bs)/Bs]; %#ok<AGROW>
    end
    startingVals = startingVals';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation using rescaled data for constrain purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volLB = min(data2)'/10000;
volUB = max(data2)'*100000;
w = .06*.94.^(0:ceil(T^(0.5)));
w = w'/sum(w);
backCast = sum(bsxfun(@times,w,data2(1:length(w),:)));

LB = zeros(parameterCount,1);
UB = inf*ones(parameterCount,1);
UB(K+1:length(LB)) = ones(size(K+1:length(LB)));
if consMode==0
    LB(K+1:length(LB)) = -ones(size(K+1:length(LB)));
end
parameters = fmincon(@heavy_likelihood,startingVals,[],[],[],[],LB,UB,[],options,data2',p,q,backCast,volLB,volUB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and data rescaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2 = bsxfun(@times,data2,scale);
volLB = min(data2)'/10000;
volUB = max(data2)'*100000;
backCast = sum(bsxfun(@times,w,data2(1:length(w),:)));
% Rescale the As and Os
[O,A] = heavy_parameter_transform(parameters,p,q,K);
O = O .* scale';
for i=1:K
    A(:,i,:) = A(:,i,:).*(scale'./scale(i));
end
% Rebuild parameters
temp = O';
for i=1:K
    for j = 1:K
        temp = [temp squeeze(A(i,j,1:p(i,j)))]; %#ok<AGROW>
    end
end
parameters = [temp parameters(length(temp)+1:length(parameters))']';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ll,~,ht] = heavy_likelihood(parameters,data2',p,q,backCast,volLB,volUB);
[VCV,~,~,scores] = robustvcv(@heavy_likelihood,parameters,0,data2',p,q,backCast,volLB,volUB);