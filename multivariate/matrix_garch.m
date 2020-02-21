function [parameters, ll, ht, VCV, scores, diagnostics] = matrix_garch(data,dataAsym,p,o,q,startingvals,options)
% Estimation of symmetric and asymmetric MATRIX multivariate GARCH models 
%
% USAGE:
%  [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = matrix_garch(DATA,DATAASYM,P,O,Q,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance
%                    estimators (e.g. RC scaled by indicator functions)
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of asymmetric lags to include
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   STARTINGVALS - [OPTIONAL] INCOMPLETE
%   OPTIONS      - [OPTIONAL] Options to use in the optimization (fminunc)
%
% OUTPUTS:
%   PARAMETERS   - K(K+1)/2*(1+P+O+Q) by 1 vector of parameters of the form
%                    [vech(C)' vech(A(1))' ... vech(A(P))' vech(G(1))' ...
%                    vech(G(O))' vech(B(1))' ... vech(B(Q))']'
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)*t^(-1))
%   SCORES       - A T by numParams matrix of individual scores
%   DIAGNOSTICS  - Structure containing some diagnostic information
%
% COMMENTS:
%    The conditional variance, H(t), of a MATRIX GARCH is modeled as follows:
%
%      H(t) = CC' + A(1)A(1)'.*r_{t-1}'*r_{t-1} + ... + A(P)A(P)'.*r_{t-P}'*r_{t-P}
%                 + G(1)G(1)'.*n_{t-1}'*n_{t-1} + ... + G(O)G(O)'.*n_{t-P}'*n_{t-P}
%                  + B(1)B(1)'.*H(t-1) +...+ B(Q)B(Q)'.*H(t-q)
%
%    where n_{t} = r_{t} .* (r_{t}<0).  If using realized measures, the
%    RM_{t-1} replaces r_{t-1}'*r_{t-1}, and the asymmetric version
%    replaces n_{t-1}'*n_{t-1}
%
% EXAMPLES:
%     % Estimation of a symmetric MATRIX GARCH(1,0,1) model
%     parameters = matrix_garch(data,[],1,0,1);
%     % Estimation of an asymmetric MATRIX GARCH(1,0,1) model
%     parameters = matrix_garch(data,[],1,1,1);
%     % Estimation of a symmetric MATRIX GARCH(1,0,1) model using realized covariacne
%     data = RC % K by K by T 3D of realized covariacne
%     parameters = matrix_garch(RC,[],1,0,1);

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 10/28/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        o=0;
        q=0;
        startingvals=[];
        options = [];
    case 4
        q=0;
        startingvals=[];
        options = [];
    case 5
        startingvals=[];
        options = [];
    case 6
        options = [];
    case 7
    otherwise
        error('Between 2 and 6 arguments required.')
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
k2 = k*(k+1)/2;
if min(t,k)<2 || t<k
    error('DATA must be a T by K matrix or a K by K by T 3D array, T>K>1');
end

%p, o, q much be non-negative scalars
if length(p)>1 || any(p<1) || floor(p)~=p
    error('P must be a positive scalar');
end
if isempty(o)
    o = 0;
end
if length(o)>1 || any(o<0) || floor(o)~=o
    error('O must be a non-negative scalar');
end
if o>0 && isempty(dataAsym)
    error('DATAASYM must be non-empty if O>0.')
end

if isempty(q)
    q = 0;
end
if length(q)>1 || any(q<0) || floor(q)~=q
    error('Q must be a non-negative scalar');
end

% Startingvals must have (k*(k+1)/2)(1+p+o+q) parameters
if  ~isempty(startingvals)
    if size(startingvals,2)>size(startingvals,1)
        startingvals = startingvals';
    end
    if length(startingvals)<(p+o+q)
        error('STARTINGVALS should be a P+O+Q by 1 vector');
    end
    % Only validate if provided
    % FIXME: Fix Validation
    % TODO: Fix Validation
    kappa = 2;
    A=startingvals(1:p);
    G=startingvals(p+1:p+o);
    B=startingvals(p+o+1:p+o+q);
    if (sum(A)+sum(G)/kappa+sum(B))>=.999998
        error('Weighted sum of STATINGVALUES must be less than 1. See Comments.');
    end
    if any(A<0) || any(B<0) || any(G<0)
        error('STARTINGVALS must all be nonnegative.');
    end
end

%Make sure options is a valid option structure
if isempty(options)
    options=optimset('fminunc');
    options.Display='iter';
    options.Diagnostics='on';
    options.LargeScale='off';
    options.MaxFunEvals = 1000*k2*(1+p+o+q);
end
try
    optimset(options);
catch ME
    error('OPTIONS is not a valid options structure');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the backCast
backCast = zeros(k);
backCastAsym = zeros(k);
tau = max(ceil(sqrt(t)),k);
weights = .06 * .94.^(0:tau);
weights = weights / sum(weights);
for i=1:tau
    backCast = backCast + weights(i) * data(:,:,i);
    if o>0
        backCastAsym = backCastAsym + weights(i) * dataAsym(:,:,i);
    end
end

%Get starting values from scalar_vt_vech
if isempty(startingvals)
    startingOptions=optimset('fminunc');
    startingOptions.Display='off';
    startingOptions.Diagnostics='off';
    startingOptions.LargeScale='off';
    startingOptions.TolX = 1e-4;
    startingOptions.TolFun = 1e-4;
    [scalarVechStartingvals,~,~,~,~,~,diagnostics] = scalar_vt_vech(data,dataAsym,p,o,q,[],[],startingOptions);
    CpC = diagnostics.intercept;
    
    startingvals = zeros(k2 * (1+p+o+q),1);
    startingvals(1:k2) = chol2vec(chol(CpC)');
    index = k2;
    for i=1:(p+o+q)
        matrixParameters = scalarVechStartingvals(i)*(.02*eye(k) + .98*ones(k));
        startingvals(index+1:index+k2) = chol2vec(chol(matrixParameters)');
        index = index + k2;
    end
end


warning('off','MATLAB:illConditionedMatrix')
ll0 = matrix_garch_likelihood(startingvals,data,dataAsym,p,o,q,backCast,backCastAsym);
[parameters,ll,exitflag,output]=fminunc('matrix_garch_likelihood',startingvals,options,data,dataAsym,p,o,q,backCast,backCastAsym);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustification
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exitflag<=0 &&  ll<ll0
    %Did not converge, but function improved
    options.MaxFunEvals=4*100*(p+q);
    options.MaxIter=2*100*(p+q);
    parameters=fminunc('matrix_garch_likelihood',parameters,options,data,dataAsym,p,o,q,backCast,backCastAsym);
end
warning('on','MATLAB:illConditionedMatrix')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustification
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1
    [ll,lls,ht]=matrix_garch_likelihood(parameters,data,dataAsym,p,o,q,backCast,backCastAsym);
    ll=-ll;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the VCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>3
    [VCV,A,B,scores]=robustvcv('matrix_garch_likelihood',parameters,0,data,dataAsym,p,o,q,backCast,backCastAsym);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the VCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagnostics = [];
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
parameterMatrices = zeros(k,k,1+p+o+q);
index = 0;
for i=1:(1+p+o+1)
    temp = vec2chol(parameters(index+1:index+k2));
    parameterMatrices(:,:,i) = temp*temp';
    index=index+k2;
end
diagnostics.C = parameterMatrices(:,:,1);
diagnostics.A = parameterMatrices(:,:,2:p+1);
if o>0
    diagnostics.G = parameterMatrices(:,:,p+2:p+o+1);
else
    diagnostics.G = [];
end

if q>0
    diagnostics.B = parameterMatrices(:,:,p+o+2:p+o+q+1);
else
    diagnostics.B = [];
end
