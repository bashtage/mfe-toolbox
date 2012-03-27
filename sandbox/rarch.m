function [parameters, ll, ht, VCV, scores] = rarch(data,p,q,v,type,method,startingvals,options)
% Estimation of RARCH(p,q) multivarate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [PARAMETERS] = rarch(DATA,P,Q)
%  [PARAMETERS,LL,HT,VCV,SCORES] = rarch(DATA,P,Q,V,TYPE,METHOD,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   Q            - Non-negative, scalar integer representing the number of conditional covariance lags
%   V            - [OPTIONAL] Number of eigenvalues which are allowed to have dynamics (Default is K)
%   TYPE         - [OPTIONAL] String, one of 'Scalar' (Default) ,'CP' (Common Persistence) or 'Diagonal'
%   METHOD       - [OPTIONAL] String, one of '2-stage' (Default) or 'Joint'
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
% Revision: 3    Date: 3/27/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        v = k;
        type = 'Scalar';
        method = '2-stage';
        startingvals = [];
        options = [];
    case 4
        type = 'Scalar';
        method = '2-stage';
        startingvals = [];
        options = [];
    case 5
        method = '2-stage';
        startingvals = [];
        options = [];
    case 6
        startingvals = [];
        options = [];
    case 7
        options = [];
    case 8
        % Nothing
    otherwise
        error('3 to 8 inputs required.')
end
if ndims(data)>3
    error('DATA must be either a T by K matrix or a K by K by T array.')
end
if T<=K
    error('DATA must be either a T by K matrix or a K by K by T array, and T must be larger than K.')
end

if p<=1 || floor(p)~=p
    error('P must be a positive scalar.')
end
if q<=0 || floor(q)~=q
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
    ifJoint = false;
elseif strcmpi(type,'Joint')
    ifJoint = true;
else
    error('METHOD must be either ''2-stage'' or  ''Joint''.')
end

if isempty(options)
    options = optimset('fmincon');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = mean(data,3);
stdData = zeros(k,k,T);
Cm12 = C^(-0.5);
for i=1:T
    stdData(:,:,i) = Cm12*data(:,:,i)*Cm12;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicParameters = fmincon(@rarch_likelihood,startingVals,[],[],[],[],LB,UB,@rarch_constraint,options,stdData,p,q,eye(K),backCast,type,false);
if isJoint
    allParameters = fmincon(@rarch_likelihood,startingValAll,[],[],[],[],LBall,UBall,@rarch_constraint,options,data,p,q,C,backCast,type,true);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
