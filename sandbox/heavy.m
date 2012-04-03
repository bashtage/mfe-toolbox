function [parameters, ll, ht, VCV, scores] = heavy(data,p,q,cons,startingVals,options)

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