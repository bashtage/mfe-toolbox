function  [parameters, ll, ht, VCV, scores] = ccc_mvgarch(data,dataAsym,p,o,q,startingVals,options)
% Estimation of TARCH(p,o,q) - Constant Conditional Correlation MV GARCH 
%
% USAGE:
%  [PARAMETERS,LL,HT,VCV,SCORES] = ccc_mvgarch(DATA,DATAASYM,P,O,Q,STARTINGVALS,OPTIONS)
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
%   STARTINGVALS - [OPTIONAL] Vector of starting values for the K TARCH models.  It should have the
%                    form  [tarch(1)' tarch(2)'  ... tarch(k)']  
%                    where each set of TARCH parameters is
%                    tarch(i) =
%                    [omega(i) alpha(i,1) ... alpha(i,p(i)) gamma(i,1)
%                                   ... gamma(i,o(i)) beta(i,1) ... beta(i,q(i))]'
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
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)*t^(-1))
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:
%   The conditional variance, H(t), of a constant conditional correlation model is
%      H(t) = Sigma(t) * R * Sigma(t)
%
%       where Sigma(t) is a diagonal matrix with TARCH(P,O,Q) volatilities on its diagonal.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 10/28/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        startingVals = [];
        options = [];
    case 6
        options = [];
    case 7
        % Nothing
    otherwise
       error('Between 5 and 7 arguments required.')
end

%data should be TxK, T>K
if ndims(data)==2
    [t,k]=size(data);
    if ~isempty(dataAsym)
        error('If DATA is a T by K matrix, DATAASYM must be empty.');
    end
    temp = zeros(k,k,t);
    dataAsym = zeros(k,k,t);
    for i=1:t
        temp(:,:,i) = data(i,:)'*data(i,:);
        dataAsym(:,:,i) = (data(i,:).*(data(i,:)<0))'*(data(i,:).*(data(i,:)<0));
    end
    data = temp;
elseif ndims(data)==3
    [k,m,t] = size(data);
    if m~=k
        error('DATA must be K by K by T is a 3D array.');
    end
    if ~isempty(dataAsym)
        if ndims(dataAsym)~=3
            error('DATAASYM must be a 3D array with the same dimensions as DATA');
        end
        [k2,m2,t2]=size(dataAsym);
        if any([k m t]~=[k2 m2 t2])
            error('DATAASYM must be a 3D array with the same dimensions as DATA');
        end
    end
end

if min(t,k)<2 || t<k
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
    [k,nothing,t]=size(data);
else
    [t,k] = size(data);
    temp = zeros(k,k,t);
    for i=1:t
        temp(:,:,i) = data(i,:)'*data(i,:);
        dataAsym(:,:,i) = (data(i,:).*(data(i,:)<0))'*(data(i,:).*(data(i,:)<0));
    end
    data = temp;
end

tarchParameters = cell(k,1);
ht=cell(k,1);
tarchVcv=cell(k,1);
tarchScores = cell(k,1);
htMat = zeros(t,k);
volData = zeros(t,k);

for i=1:k
    % Prepare the data
    tempData = sqrt(data(i,i,:));
    if o(i)>0
        tempData = tempData .* ((dataAsym(i,i,:)<=0) -(dataAsym(i,i,:)>0));
    end
    volData(:,i) = squeeze(tempData);
    tarchStartingVals = [];
    if ~isempty(startingVals)
        parameterEnd = i + sum(p(1:i)) + sum(o(1:i)) + sum(q(1:i));
        parameterStart = parameterEnd - 1 - p(i) - o(i) - q(i) + 1;
        tarchStartingVals = startingVals(parameterStart:parameterEnd);
    end
    [tarchParameters{i},nothing,ht{i},nothing,tarchVcv{i},tarchScores{i}]=tarch(volData(:,i),p(i),o(i),q(i),[],2,tarchStartingVals,options);
    htMat(:,i) = ht{i};
end

htArray = zeros(k,k,t);
for i=1:k
    for j=i:k
        htArray(i,j,:) = reshape(sqrt(htMat(:,i).*htMat(:,j)),[1 1 t]);
        htArray(j,i,:) = htArray(i,j,:);
    end
end

stdData = data./htArray;
R = mean(stdData,3);
R = R./sqrt(diag(R)*diag(R)');
parameters = zeros(sum(p)+sum(o)+sum(q)+k*(k-1)/2,1);
parameterCount = 1;


for i=1:k
    parameterCountEnd = parameterCount + sum(1+p(i)+o(i)+q(i));
    parameters(parameterCount:parameterCountEnd-1) = tarchParameters{i};
    parameterCount = parameterCountEnd;
end
corrIndex = parameterCount:parameterCount+(k*(k-1)/2) - 1;
parameters(corrIndex) = corr_vech(R);

% Compute the LL and covariance
if nargout>1
    [ll, lls, ht] = ccc_mvgarch_joint_likelihood(parameters,data,volData,p,o,q);
    ll = - ll;
end
% Inference
if nargout>3
    A = zeros(length(parameters));
    parameterCount = 1;
    scores = zeros(t,length(parameters));
    for i=1:k
        parameterCountEnd = parameterCount + sum(1+p(i)+o(i)+q(i)) - 1;
        index = parameterCount:parameterCountEnd;
        A(index,index) = tarchVcv{i}^(-1);
        scores(:,index) = tarchScores{i};
        parameterCount = parameterCountEnd + 1;
    end
    corrParameters = parameters(corrIndex);
    A(corrIndex,:) = hessian_2sided_nrows(@ccc_mvgarch_joint_likelihood,parameters,k*(k-1)/2,data,volData,p,o,q);
    [nothing,scores(:,corrIndex)] = gradient_2sided(@ccc_mvgarch_likelihood,corrParameters,data,htMat);
    A = A/t;
    B = covnw(scores,0,0);
    Ainv = A^(-1);
    VCV = Ainv*B*Ainv'/t;
end