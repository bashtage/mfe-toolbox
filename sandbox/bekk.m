function [parameters, ll, Ht, VCV, scores] = bekk(data,dataAsym,p,o,q,type,startingVals,options)
% Estimation of symmetric and asymmetric BEKK(p,o,q) multivarate volatility models
%
% USAGE:
%  [PARAMETERS] = bekk(DATA,[],P,O,Q)
%  [PARAMETERS,LL,HT,VCV,SCORES] = bekk(DATA,DATAASYM,P,O,Q,TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative, scalar integer representing the number of asymmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   TYPE         - [OPTIONAL] String, one of 'Scalar' (Default) ,'Diagonal' or 'Full'
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)
%
% OUTPUTS:
%   PARAMETERS   -
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:

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
    case 5
        type = 'Scalar';
        startingVals = [];
        options = [];
    case 6
        startingVals = [];
        options = [];
    case 7
        options = [];
    case 8
        % Nothing
    otherwise
        error('5 to 8 inputs required.')
end
if ndims(data)>3
    error('DATA must be either a T by K matrix or a K by K by T array.')
end
if T<=k
    error('DATA must be either a T by K matrix or a K by K by T array, and T must be larger than K.')
end
if ndims(data)==3 && o>0
    if isempty(dataAsym) || ndims(dataAsym)~=3 || any(size(dataAsym)~=size(data))
        error('DATAASYM must be provided when O>0 and DATA is a 3D array.')
    end
end

if p<1 || floor(p)~=p
    error('P must be a positive scalar.')
end
if isempty(o)
    o=0;
end
if o<=0 || floor(q)~=o
    error('O must be a non-negative scalar.')
end
if isempty(q)
    q=0;
end
if q<=0 || floor(q)~=q
    error('Q must be a non-negative scalar.')
end

if strcmpi(type,'Scalar')
    type = 1;
elseif strcmpi(type,'Diagonal')
    type = 2;
elseif strcmpi(type,'Full')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''Diagonal'' or ''Full''.')
end

k2 = k*(k+1)/2;
if ~isempty(startingVals)
    switch type
        case 1
            count = p+o+q;
        case 2
            count = (p+o+q)*k;
        case 3
            count = (p+o+q)*k*k;
    end
    count = count + k2;
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
        eta = (data(i,:).*(data(i,:)<0));
        dataAsym(:,:,i) =  eta'*eta;
    end
    data = temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[startingVals,~,~,~,~,d] = scalar_vt_vech(data,dataAsym,p,o,q,[],options);
C = d.intercept;
C = chol2vec(chol(C)');
switch type
    case 1
        shape = 1;
    case 2
        shape = ones(k,1);
    case 3
        shape = eye(k);
end
sv = [];
for i=1:p+o+q
    temp = sqrt(startingVals(i));
    temp = temp * shape;
    sv = [sv;temp(:)]; %#ok<AGROW>
end
startingVals = [C;sv];
UB = .99998 * ones(size(startingVals));
UB(1:k2) = inf;
LB = -UB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use stdData and set C = eye(K)
%rarch_likelihood(startingVals,data,p,q,C,backCast,type,false,false);
parameters = fmincon(@bekk_likelihood,startingVals,[],[],[],[],LB,UB,@bekk_constraint,options,data,dataAsym,p,o,q,backCast,backCastAsym,type);
[ll,~,Ht] = bekk_likelihood(parameters,data,dataAsym,p,o,q,backCast,backCastAsym,type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VCV,~,~,scores] = robustvcv(@bekk_likelihood,parameters,0,data,dataAsym,p,o,q,backCast,backCastAsym,type);
