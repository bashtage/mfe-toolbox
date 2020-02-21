function [parameters, ll, ht, intercept, VCV, scores, diagnostics] = scalar_vt_vech(data,dataAsym,p,o,q,composite,startingvals,options)
% Estimation of symmetric and asymmetric scalar multivariate vech ARCH models using variance
% targeting to reduce the number of parameters needing to be estimated simultaneously
%
% USAGE:
%  [PARAMETERS,LL,HT,INTERCEPT,VCV,SCORES,DIAGNOSTICS] = scalar_vt_vech(DATA,DATAASYM,P,O,Q,COMPOSITE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of asymmetric lags to include
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   COMPOSUTE    - [OPTIONAL] String, one of:
%                    'None' - Standard K-dimensional QMLE (Default)
%                    'Diagonal' - Use only pairwise likelihoode for i,i+1, i=1,2,...,K-1
%                    'Full' - Use all pairwise likelihoods
%   STARTINGVALS - [OPTIONAL] (p+o+q) x 1 vector of starting values
%   OPTIONS      - [OPTIONAL] Options to use in the optimization (fminunc)
%
% OUTPUTS:
%   PARAMETERS   - A p+o+q column vector of parameters.  The intercept is reported in DIAGNOSTICS
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   INTERCEPT    - K by K matrix containing the intercept computed from the unconditional variance
%                    and parameters
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)*t^(-1))
%   SCORES       - A T by numParams matrix of individual scores
%
% COMMENTS:
%    The conditional variance, H(t), of a scalar variance-targeting vech is modeled
%    as follows:
%
%      H(t) = (1-alpha(1)+...+alpha(p)-beta(1)-...-beta(q))*C + ...
%              alpha(1)*r_{t-1}'*r_{t-1} + ... + alpha(p)*r_{t-p}'*r_{t-p}+...
%              gamma(1)*n_{t-1}'*n_{t-1} + ... + gamma(o)*n_{t-p}'*n_{t-p}+...
%              beta(1)*H(t-1) +...+ beta(q)*H(t-q)
%
%    where n_{t} = r_{t} .* (r_{t}<0)
%
% EXAMPLES:
%     Estimation of a Scalar VECH(1,0,1)
%       [simulatedData, Ht, pseudoRC] = scalar_vt_vech_simulate(1000, [.02 .04 .95], .01*eye(2), 1, 1, 1, 72);
%       parameters = scalar_vt_vech(simulatedData,[],1,0,1)
%     Estimation of an asymmetric Scalar VECH(1,1,1)
%       parameters = scalar_vt_vech(simulatedData,[],1,1,1)
%     Estimation of an asymmetric Scalar VECH(1,1,1) using realized-type data
%       asymPseudoRC = zeros(size(pseudoRC));
%       for i = 1:1000
%         asymPseudoRC(:,:,i) = pseudoRC(:,:,i).*(double(simulatedData(i,:)<0)'*double(simulatedData(i,:)<0));
%       end
%       parameters = scalar_vt_vech(pseudoRC,asymPseudoRC,1,1,1);

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
        composite = 'None';
        startingvals=[];
        options = [];
    case 4
        q=0;
        composite = 'None';
        startingvals=[];
        options = [];
    case 5
        composite = 'None';
        startingvals=[];
        options = [];
    case 6
        startingvals=[];
        options = [];
    case 7
        options = [];
    case 8
        % Nothing
    otherwise
        error('Between 3 and 8 arguments required.')
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

if isempty(composite)
    composite = 'None';
end
switch lower(composite)
    case {'none'}
        useComposite = 0;
    case {'diagonal'}
        useComposite = 1;
    case {'full'}
        useComposite = 2;
end

% Startingvals must have p+o+q parameters sum(alpha)+kappa*sum(gamma)+sum(beta)<1
C = mean(data,3);
if o>0
    Casym = mean(dataAsym,3);
    kappa = 1 / (max(eig(C^(-0.5)*Casym*C^(-0.5))) + eps);
else
    Casym = zeros(k);
    kappa = 2;
end
if  ~isempty(startingvals)
    if size(startingvals,2)>size(startingvals,1)
        startingvals = startingvals';
    end
    if length(startingvals)<(p+o+q)
        error('STARTINGVALS should be a P+O+Q by 1 vector');
    end
    %Only validate if provided
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

%Augment the data with backcasts
[startingvals,lls,output_parameters] = scalar_vt_vech_starting_values(startingvals,data,dataAsym,p,o,q,C,Casym,kappa,useComposite,backCast,backCastAsym); %#ok<NASGU>

% finally to transform the parameters to the unrestricted equivalents
startingvals=scalar_vt_vech_transform(startingvals,p,o,q,kappa);

[parameters,ll,exitflag,output]=fminunc('scalar_vt_vech_likelihood',startingvals,options,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustification
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exitflag<=0 &&  ll<lls(1)
    %Did not converge, but function improved
    options.MaxFunEvals=4*100*(p+q);
    options.MaxIter=2*100*(p+q);
    parameters=fminunc('scalar_vt_vech_likelihood',parameters,options,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,true);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustification
%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters=scalar_vt_vech_itransform(parameters,p,o,q,kappa);
if nargout>1
    [ll,~,ht]=scalar_vt_vech_likelihood(parameters,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,false);
    ll=-ll;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the VCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>4
    k2 = k*(k+1)/2;
    CVech=vech(C);
    if o==0
        momentCount = k2;
        CasymVech = [];
    else
        momentCount = 2*k2;
        CasymVech=vech(Casym);
    end
    
    A = zeros(momentCount + p + o + q);
    A(1:momentCount, 1:momentCount) = - t * eye(momentCount);
    jointParameters=[CVech;CasymVech;parameters];
    
    % A will be k2 + p + q square
    scores = zeros(t,momentCount+p+q);
    scoreCount = 1;
    for i=1:k
        for j=i:k
            scores(:,scoreCount) = squeeze(data(j,i,:));
            if o>0
                scores(:,k2+scoreCount) = squeeze(dataAsym(j,i,:));
            end
            scoreCount = scoreCount  + 1;
        end
    end
    [~,gt] = gradient_2sided(@scalar_vt_vech_likelihood,parameters,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,false,useComposite,false);
    scores(:,momentCount+1:momentCount+p+o+q) = gt;
    A(momentCount+1:momentCount+p+o+q,:) = hessian_2sided_nrows(@scalar_vt_vech_likelihood,jointParameters,p+o+q,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,true,useComposite,false);
    
    A = A/t;
    B = covnw(scores);
    Ainv = inv(A);
    VCV = Ainv*B*Ainv'/t; %#ok<MINV>
    VCV = VCV(momentCount+1:momentCount+p+o+q,momentCount+1:momentCount+p+o+q);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the VCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
diagnostics.kappa=kappa;
alpha = sum(parameters(1:p));
gamma = sum(parameters(p+1:p+o));
beta = sum(parameters(p+o+1:p+o+q));
diagnostics.intercept = C*(1-alpha-beta);
diagnostics.composite = composite;
if o>0
    diagnostics.intercept  = diagnostics.intercept  - gamma*Casym;
end
intercept = diagnostics.intercept;