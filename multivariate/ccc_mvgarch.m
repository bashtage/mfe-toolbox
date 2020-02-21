function  [parameters, ll, Ht, VCV, scores] = ccc_mvgarch(data,dataAsym,p,o,q,gjrType,startingVals,options)
% Estimation of Constant Conditional Correlation MV GARCH with TARCH(p,o,q) or GJRGARCH(p,o,q)
% conditional variances
%
% USAGE:
%  [PARAMETERS] = ccc_mvgarch(DATA,DATAASYM,P,O,Q)
%  [PARAMETERS,LL,HT,VCV,SCORES] = ccc_mvgarch(DATA,DATAASYM,P,O,Q,GJRTYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL]
%   P            - Positive, scalar integer representing the number of symmetric innovations -OR-
%                     K by 1 vector of individual symmetric innovations order
%   O            - Non-negative, scalar integer representing the number of asymmetric lags -OR-
%                     K by 1 vector of individual asymmetric innovations order
%   Q            - Non-negative, scalar integer representing the number of conditional variance lags -OR-
%                     K by 1 vector of individual conditional variance lags
%   GJRTYPE      - [OPTIONAL] Scalar integer in {1,2} indicating whether to use a TARCH-type model (1) 
%                     or a GJR-GARCH-type model (2, Default)
%   STARTINGVALS - [OPTIONAL] Vector of starting values for the K TARCH models.  It should have the
%                     form  [tarch(1)' tarch(2)'  ... tarch(k)']  
%                     where each set of TARCH parameters is
%                     tarch(i) =
%                     [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                    ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
%   OPTIONS      - [OPTIONAL] Options to use in the TARCH model optimization (fminunc)
%
% OUTPUTS:
%   PARAMETERS   - a parameter vector of the form
%                  [tarch(1)' tarch(2)'  ... tarch(k)' corr_vech(R)]
%                    where each set of TARCH parameters is
%                    tarch(i) =
%                    [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                   ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
%                    and where R is the constant conditional correlation.
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:
%   The conditional variance, H(t), of a constant conditional correlation model is
%      H(t) = Sigma(t) * R * Sigma(t)
%
%       where Sigma(t) is a diagonal matrix with TARCH(P,O,Q) volatilities on its diagonal.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 4/13/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        gjrType = [];
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
       error('Between 5 and 8 arguments required.')
end

%data should be TxK, T>K
if ndims(data)==2
    [T,k]=size(data);
    if ~isempty(dataAsym)
        error('If DATA is a T by K matrix, DATAASYM must be empty.');
    end
    temp = zeros(k,k,T);
    dataAsym = zeros(k,k,T);
    for i=1:T
        temp(:,:,i) = data(i,:)'*data(i,:);
        dataAsym(:,:,i) = (data(i,:).*(data(i,:)<0))'*(data(i,:).*(data(i,:)<0));
    end
    data = temp;
elseif ndims(data)==3
    [k,m,T] = size(data);
    if m~=k
        error('DATA must be K by K by T is a 3D array.');
    end
    if ~isempty(dataAsym)
        if ndims(dataAsym)~=3
            error('DATAASYM must be a 3D array with the same dimensions as DATA');
        end
        [k2,m2,t2]=size(dataAsym);
        if any([k m T]~=[k2 m2 t2])
            error('DATAASYM must be a 3D array with the same dimensions as DATA');
        end
    end
end

if min(T,k)<2 || T<k
    error('DATA must be a T by K matrix or a K by K by T 3D array, T>K>1');
end

%p, o, q much be non-negative scalars
% p
if length(p)==1
    if p<1 || floor(p)~=p
        error('P must be a positive integer if scalar.');
    end
    p = ones(k,1) * p;
else
    if length(p)~=k || min(size(p))~=1 || any(p<1) || any(floor(p)~=p)
        error('P must contain K positive integer elements if a vector.');
    end
end
% o
if length(o)==1
    if o<0 || floor(o)~=o
        error('O must be a non-neagative integer if scalar.');
    end
    o = ones(k,1) * o;
else
    if length(o)~=k || min(size(o))~=1 || any(o<0) || any(floor(o)~=o)
        error('O must contain K non-negative integer elements if a vector.');
    end
end
% q
if length(q)==1
    if q<0 || floor(q)~=q
        error('Q must be a non-neagative integer if scalar.');
    end
    q = ones(k,1) * q;
else
    if length(q)~=k || min(size(q))~=1 || any(q<0) || any(floor(q)~=q)
        error('Q must contain K non-negative integer elements if a vector.');
    end
end
if any(o>0) && isempty(dataAsym)
    error('DATAASYM must be non-empty if O>0.')
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



if  ~isempty(startingVals)
    if size(startingVals,2)>size(startingVals,1)
        startingVals = startingVals';
    end
    if length(startingVals)<(k+sum(p)+sum(o)+sum(q))
        error('STARTINGVALS should be a K+sum(P)+sum(O)+sum(Q) by 1 vector');
    end
end

%Make sure options is a valid option structure
if isempty(options)
    options = optimset('fminunc');
    options.Display = 'off';
    options.LargeScale = 'off';
end

try
    optimset(options);
catch ME
    error('OPTIONS is not a valid options structure');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data)==3
    [k,~,T]=size(data);
    if isempty(dataAsym)
        dataAsym = repmat(eye(k),[1 1 T]);
    end
    data2d = zeros(T,k);
    for t=1:T
        data2d(t,:) = sqrt(diag(data(:,:,t)))';
        data2d(t,:) = data2d(t,:) .* (2*(diag(dataAsym(:,:,t))>0)-1)';
    end
else
    data2d = data;
    [T,k] = size(data);
    temp = zeros(k,k,T);
    for i=1:T
        temp(:,:,i) = data(i,:)'*data(i,:);
        dataAsym(:,:,i) = (data(i,:).*(data(i,:)<0))'*(data(i,:).*(data(i,:)<0));
    end
    data = temp;
end

[H,univariate] = dcc_fit_variance(data2d,p,o,q,gjrType,startingVals);
htArray = zeros(k,k,T);
for i=1:k
    for j=i:k
        htArray(i,j,:) = reshape(sqrt(H(:,i).*H(:,j)),[1 1 T]);
        htArray(j,i,:) = htArray(i,j,:);
    end
end

stdData = data./htArray;
R = mean(stdData,3);
R = R./sqrt(diag(R)*diag(R)');
Ht = bsxfun(@times,R,htArray);
ll = 0;
likConst = k*log(2*pi);
for t=1:T
    ll = ll + 0.5*(likConst + log(det(Ht(:,:,t))) + sum(diag((Ht(:,:,t)\eye(k))*data(:,:,t))));
end

ll = -ll;

% Format parameters
parameters = zeros(sum(p)+sum(o)+sum(q)+k*(k-1)/2,1);
offset = 0;
for i=1:k
    u = univariate{i};
    count = sum(1+p(i)+o(i)+q(i));
    parameters(offset + (1:count)) = u.parameters;
    offset = offset + count;
end
corrIndex = offset + (1:k*(k-1)/2);
parameters(corrIndex) = corr_vech(R);

% Compute the LL and covariance
if nargout>3
    v = length(parameters);
    A = zeros(v);
    scores = zeros(T,v);
    offset = 0;
    for i=1:k
        u = univariate{i};
        count = sum(1+p(i)+o(i)+q(i));
        index = offset + (1:count);
        offset = offset + count;
        A(index,index) = u.A;
        scores(:,index) = u.scores;
    end
    H = hessian_2sided_nrows(@dcc_inference_objective,parameters,k*(k-1)/2,data,dataAsym,0,0,0,univariate);
    [~,s] = gradient_2sided(@dcc_inference_objective,parameters,data,dataAsym,0,0,0,univariate);
    scores(:,corrIndex) = s(:,corrIndex);
    A(corrIndex,:) = H/T;
    B = covnw(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
end