function [parameters, LL, errors, SEregression, diagnostics, VCVrobust, VCV, likelihoods, scores]=armaxfilter(y,constant,p,q,x,startingVals,options,holdBack,sigma2)
% ARMAX(P,Q) estimation
%
% USAGE:
%   [PARAMETERS]=armaxfilter(Y,CONSTANT,P,Q)
%   [PARAMETERS, LL, ERRORS, SEREGRESSION, DIAGNOSTICS, VCVROBUST, VCV, LIKELIHOODS, SCORES]
%               =armaxfilter(Y,CONSTANT,P,Q,X,STARTINGVALS,OPTIONS,HOLDBACK,SIGMA2)
%
% INPUTS:
%   Y            - A column of data
%   CONSTANT     - Scalar variable: 1 to include a constant, 0 to exclude
%   P            - Non-negative integer vector representing the AR orders to include in the model.
%   Q            - Non-negative integer vector representing the MA orders to include in the model.
%   X            - [OPTIONAL]  a T by K  matrix of exogenous variables. These line up exactly with
%                    the Y's and if they are time series, you need to shift them down by 1 place,
%                    i.e. pad the bottom with 1 observation and cut off the top row [ T by K].  For
%                    example, if you want to include X(t-1) as a regressor, Y(t) should line up
%                    with X(t-1)
%   STARTINGVALS - [OPTIONAL] A (CONSTANT+length(P)+length(Q)+K) vector of starting values.
%                   [constant ar(1) ... ar(P) xp(1) ... xp(K) ma(1) ... ma(Q) ]'
%   OPTIONS      - [OPTIONAL] A user provided options structure. Default options are below.
%   HOLDBACK     - [OPTIONAL] Scalar integer indicating the number of observations to withhold at
%                    the start of the sample. Useful when testing models with different lag lengths
%                    to produce comparable likelihoods, AICs and SBICs. Should be set to the highest
%                    lag length (AR or MA) in the models studied.
%   SIGMA2       - [OPTIONAL] T by 1 vector containing the conditional variance of the error.
%                    Allows for GLS estimation.  Default value is ones(T,1).
%
% OUTPUTS:
%   PARAMETERS   - A 1+length(p)+size(X,2)+length(q) column vector of parameters with
%                  [constant ar(1) ... ar(P) xp(1) ... xp(K) ma(1) ... ma(Q) ]'
%   LL           - The log-likelihood of the regression
%   ERRORS       - A T by 1 length vector of errors from the regression
%   SEREGRESSION - The standard error of the regressions
%   DIAGNOSTICS  - A structure of diagnostic information containing:
%                     P          - The AR lags used in estimation
%                     Q          - The MA lags used in estimation
%                     C          - Indicator if constant was included
%                     nX         - Number of X variables in the regression
%                     AIC        - Akaike Information Criteria for the estimated model
%                     HQC        - Hannan-Quinn  Information Criteria for the estimated model
%                     SBIC       - Bayesian (Schwartz) Information Criteria for the
%                                   estimated model
%                     ADJT       - Length of sample used for estimation after HOLDBACK adjustments
%                     T          - Number of observations
%                     ARROOTS    - The characteristic roots of the ARMA
%                                   process evaluated at the estimated parameters
%                     ABSARROOTS - The absolute value (complex modulus if
%                                    complex) of the ARROOTS
%   VCVROBUST    - Robust parameter covariance matrix%
%   VCV          - Non-robust standard errors (inverse Hessian)
%   LIKELIHOODS  - A T by 1 vector of log-likelihoods
%   SCORES       - Matrix of scores (# of params by T)
%
% COMMENTS:
%   The ARMAX(P,Q) model is:
%      y(t) = const + arp(1)*y(t-1) + arp(2)*y(t-2) + ... + arp(P) y(t-P) +
%                   + ma(1)*e(t-1)  + ma(2)*e(t-2)  + ... + ma(Q) e(t-Q)
%                   + xp(1)*x(t,1)  + xp(2)*x(t,2)  + ... + xp(K)x(t,K)
%                   + e(t)
%
%    The main optimization is performed with lsqnonlin with the default options:
%      options =  optimset('lsqnonlin');
%      options.MaxIter = 10*(maxp+maxq+constant+K);
%      options.Display='iter';
%
%   You should use the MEX file (or compile if not using Win64 Matlab) for armaxerrors.c as it
%   provides speed ups of approx 10 times relative to the m file version armaxerrors.m
%
% EXAMPLE:
%   To fit a standard ARMA(1,1), use
%       parameters = armaxfilter(y,1,1,1)
%   To fit a standard ARMA(3,4), use
%       parameters = armaxfilter(y,1,[1:3],[1:4])
%   To fit an ARMA that includes lags 1 and 3 of y and 1 and 4 of the MA term, use
%       parameters = armaxfilter(y,1,[1 3],[1 4])
%
% See also ARMAXFILTER_SIMULATE, HETEROGENEOUSAR, ARMAXERRORS

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/19/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        p=0;
        q=0;
        x=[];
        startingVals=[];
        options=[];
        holdBack = [];
        sigma2 = [];
    case 3
        q=0;
        x=[];
        startingVals=[];
        options=[];
        holdBack = [];
        sigma2 = [];
    case 4
        x=[];
        startingVals=[];
        options=[];
        holdBack = [];
        sigma2 = [];
    case 5
        startingVals=[];
        options=[];
        holdBack = [];
        sigma2 = [];
    case 6
        options=[];
        holdBack = [];
        sigma2 = [];
    case 7
        holdBack = [];
        sigma2 = [];
    case 8
        sigma2 = [];
    case 9
        %Nothing
    otherwise
        error('Between 2 and 9 inputs required');
end
K=size(x,2);
%%%%%%%%%%%%%%%
% y
%%%%%%%%%%%%%%%
if size(y,2) > 1 || length(y)==1
    error('y series must be a column vector.')
elseif isempty(y)
    error('y is empty.')
end

%%%%%%%%%%%%%%%
% P
%%%%%%%%%%%%%%%
if size(p,2)>size(p,1)
    p=p';
end
if isempty(p)
    p=0;
end
if min(size(p))~=1
    error('P must be a column vector of included lags')
end
if  any(p<0) || any(floor(p)~=p)
    error('P must contain non-negative integers only')
end
if max(p)>=(length(y)-max(p))
    error('Too many lags in the AR.  max(P)<T/2')
end
maxp=max(p);
if size(p,1)==1 && p==0
    p=[];
end
if length(unique(p))~=length(p)
    error('P must contain at most one of each lag')
end
nP=length(p);
%%%%%%%%%%%%%%%
% Q
%%%%%%%%%%%%%%%
if size(q,2)>size(q,1)
    q=q';
end
if isempty(q)
    q=0;
end
if min(size(q))~=1
    error('Q must be a column vector of included lags')
end
if  any(q<0) || any(floor(q)~=q)
    error('Q must contain non-negative integers only')
end
if max(q)>=length(y)
    error('Too many lags in the AR.  max(Q)<T')
end
maxq=max(q);
if size(q,1)==1 && q==0
    q=[];
end
if length(unique(q))~=length(q)
    error('Q must contain at most one of each lag')
end
nQ=length(q);
%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%
if ~constant && isempty(p) && isempty(q)
    error('At least one of CONSTANT, P or Q must be nonnegative')
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1')
end

%%%%%%%%%%%%%%%
% startign vals
%%%%%%%%%%%%%%%
if nargin>5
    if ~isempty(startingVals)
        if size(startingVals,2)>1
            startingVals=startingVals';
        end
        %Validate starting vals
        if length(startingVals)~=(nP+nQ+K+constant) || size(startingVals,2)~=1
            error('STARTINGVALS must be a column vector with Constant+length(p)+length(q)+K elements');
        end
    end
else
    startingVals=[];
end

%%%%%%%%%%%%%%%
% options
%%%%%%%%%%%%%%%
if nargin>6 && ~isempty(options)
    try
        optimset(options);
    catch ME
        error('OPTIONS is not a valid minimization option structure');
    end
elseif q>0
    %Setup the options in case of none provided
    options =  optimset('lsqnonlin');
    options.MaxIter = 1000*(maxp+maxq+constant+K);
    options.MaxFunEvals = 1000*(maxp+maxq+constant+K)^2;
    options.Display='iter';
end

%%%%%%%%%%%%%%
% sigma2
%%%%%%%%%%%%%%
if ~isempty(sigma2)
    if size(sigma2,2)>size(sigma2,1)
        sigma2 = sigma2';
    end
    if size(sigma2,1)~=size(y,1) || size(sigma2,2)~=1 || ndims(sigma2)>2 || any(sigma2<0)
        error('SIGMA2 must be a T by 1 vector containing strictly positive numbers.')
    end
else
    sigma2=ones(size(y));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = size(y,1);
T50 = max(ceil(log(T)),2);

np = length(p);
nq = length(q);
zerosToPad = max(maxq-maxp,0);
m = max(maxp,maxq);
if ~isempty(holdBack)
    if holdBack > maxp;
        m = m + (holdBack-maxp);
    end
else
    holdBack = maxp;
end
yAugmentedForMA=[zeros(zerosToPad,1);y];
if ~isempty(x)
    xAugmentedForMA=[zeros(zerosToPad,K);x];
else
    xAugmentedForMA = [];
end
sigmaAugmentedForMA = [ones(zerosToPad,1);sqrt(sigma2)];
tempSigma = ones(size(yAugmentedForMA));
T = size(yAugmentedForMA,1);


startingValsProblem = false;
if ~isempty(startingVals) && maxq > 0
    yTilde = y(holdBack+1:length(y));
    if constant
        yTilde = yTilde - mean(yTilde);
    end
    worstCaseSSE = yTilde'*yTilde;
    % Check for invertable MA
    MAParameters = startingVals(constant+np+K+1:constant+np+K+nq);
    newMAParameters = convert_ma_roots(MAParameters,q);
    if any(newMAParameters~=MAParameters)
        startingVals(constant+np+K+1:constant+np+K+nq) = newMAParameters;
    end
    e = armaxerrors(startingVals,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,tempSigma);
    SSE_startingVals = e'*e;
    if isnan(SSE_startingVals) || (worstCaseSSE <  SSE_startingVals)
        startingValsProblem = true;
        warning('MFEToolbox:StartingValues','User provided STARTINGVALS produce a higher SSE then a naive model and so are being ignored.')
    end
end




% Need to handle important cases
if maxq==0
    % AR and X only, so easy to handle
    [Y,regressors]=newlagmatrix(y,maxp,constant);
    if constant
        regressors = regressors(:,[1 p'+1]);
    else
        regressors = regressors(:,p');
    end
    if K>0
        regressors = [regressors  x(maxp+1:T,:)];
    end
    parameters=regressors\Y;
elseif  maxq>0 && (isempty(startingVals) || startingValsProblem)
    % Use a high order approx to get an estimate of e
    [Y,regressors]=newlagmatrix(y,max(maxp,T50),constant);
    if K>0
        regressors = [regressors  x(max(maxp,T50)+1:size(x,1),:)];
    end
    b=regressors\Y;
    e=Y-regressors*b;
    
    % Set up the dependant variable Y.  This variable should be T - maxp by 1
    [Y,regressors]=newlagmatrix(y,maxp,constant);
    if constant
        regressors = regressors(:,[1 p'+1]);
    else
        regressors = regressors(:,p');
    end
    if K>0
        regressors = [regressors  x(maxp+1:size(x,1),:)];
    end
    startingValRegressors = regressors;
    
    if ~isempty(regressors)
        b = regressors\Y;
        startingVals0 = [b' zeros(1,nq)]';
        e2 = armaxerrors(startingVals0,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,tempSigma);
        SSE0 = e2'*e2;
    else
        e2 = Y;
        SSE0 = e2'*e2;
        startingVals0 = zeros(nq,1);
    end
    
    
    count = 1;
    
    iteratviveStartingVals = zeros(20,constant+np+K+nq);
    iteratviveSSE = zeros(20,1) + (SSE0 + 1);
    
    while count<= 20;
        if count>= 3 && max(abs(startingVals - lastStartingVals))<.01
            break
        end
        % e should have the same size as y.  If not, augment
        if size(e,1)<size(y,1)
            e = [zeros(size(y,1)-size(e,1),1);e]; %#ok<AGROW>
        end
        % If will lose max(maxp,maxq) datapoints here, so need to augment if maxq>maxp
        if maxq>maxp
            e = [zeros(maxq-maxp,1);e]; %#ok<AGROW>
        end
        % Need the double max here to handle the case where maxp>maxq
        [temp,elags]=newlagmatrix(e,max(maxq,maxp),0); %#ok<ASGLU>
        % Select the correct columns of elags
        elags = elags(:,q');
        regressors = [startingValRegressors elags];
        
        if count>1
            lastStartingVals = startingVals;
        end
        startingVals=regressors\Y;
        e2 = armaxerrors(startingVals,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,tempSigma);
        iteratviveStartingVals(count,:) = startingVals';
        iteratviveSSE(count) = e2'*e2;
        e = Y - regressors*startingVals;
        count = count + 1;
    end
    [minSSE, pos] = min(iteratviveSSE);
    startingVals = iteratviveStartingVals(pos,:)';
    if SSE0 < minSSE
        startingVals = startingVals0;
    end
    MAparameters = zeros(1,maxq);
    MAparameters(q) = startingVals(constant+np+K+1:constant+np+K+nq);
    if max(abs(roots([1 MAparameters])))>1
        % Try to invert
        MAParameters = startingVals(constant+np+K+1:constant+np+K+nq);
        newMAParameters = convert_ma_roots(MAParameters, q);
        if any(newMAParameters~=MAParameters)
            startingVals(constant+np+K+1:constant+np+K+nq) = newMAParameters;
        end
    end
end


if maxq>0
    % Only do nonlinear estimation if q>0
    parameters = lsqnonlin('armaxerrors',startingVals,[],[],options,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,sigmaAugmentedForMA);
    
    % Check is MA is invertible
    MAparameters = zeros(1,maxq);
    MAparameters(q) = parameters(constant+np+K+1:constant+np+K+nq);
    
    if max(abs(roots([1 MAparameters])))>1
        % Try to convert
        MAParameters = parameters(constant+np+K+1:constant+np+K+nq);
        newMAParameters = convert_ma_roots(MAParameters, q);
        
        if any(newMAParameters~=MAParameters)
            parameters(constant+np+K+1:constant+np+K+nq) = newMAParameters;
        else
            warning('MFEToolbox:Invertability','MA parameters are not invertible, and it is not possible to invert them with in the selected irregular MA specification.');
        end
        parameters = lsqnonlin('armaxerrors',parameters,[],[],options,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,sigmaAugmentedForMA);
    end
    
    % This is a finalizer.  Probably superfluous
    % optionssearch  =  optimset('fminsearch');
    % optionssearch.MaxIter = 10*(maxp+maxq+constant+K);
    % optionssearch.Display='off';
    % parameters = fminsearch('armaxfilter_likelihood',parameters,optionssearch,p,q,constant,yAugmentedForMA,xAugmentedForMA,m);
end

[LL,likelihoods,errors]=armaxfilter_likelihood(parameters,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,sigmaAugmentedForMA);
likelihoods=-likelihoods;
LL=-LL;
SEregression=sqrt(errors'*errors/(length(errors)-length(parameters)));

if nargout>=4
    diagnostics.P    = p;
    diagnostics.Q    = q;
    diagnostics.adjT = length(y) - holdBack;
    diagnostics.T    = length(y);
    diagnostics.K    = K;
    [aic, hqc, sbic]     = aichqcsbic(errors,constant,p,q,x);
    diagnostics.AIC  = aic;
    diagnostics.HQC  = hqc;
    diagnostics.SBIC = sbic;
    diagnostics.C    = constant;
    diagnostics.nX  = size(x,2);
    [arroots, absarroots]  = armaroots(parameters,constant,p,q,x);
    diagnostics.arroots    = arroots;
    diagnostics.absarroots = absarroots;
    diagnostics.holdBack = holdBack;
end

if nargout>=5
    % Analytical in the cast of AR
    if maxq==0
        [leeds,lags] = newlagmatrix(y,maxp,0);
        lags = lags(:,p);
        tau = max(size(leeds,1),size(lags,1));
        if constant
            lags = [ones(tau,1) lags];
        end
        if ~isempty(x)
            T = size(x,1);
            X = [lags x(T-tau+1:T,:)];
        else
            X = lags;
        end
        
        T = length(errors);
        e = errors(maxp+1:T);
        T = length(e);
        XpXi = (X'*X/T)\eye(size(X,2));
        XeeX = zeros(size(X,2));
        for t=1:T
            XeeX = XeeX + e(t)^2*X(t,:)'*X(t,:);
        end
        XeeX = XeeX/T;
        VCVrobust = XpXi*XeeX*XpXi/T;
        VCV = SEregression^2*XpXi/T;
    else
        [VCVrobust,~,B,scores]=robustvcv('armaxfilter_likelihood',parameters,0,p,q,constant,yAugmentedForMA,xAugmentedForMA,m,sigmaAugmentedForMA);
        VCV=B^(-1)/(T-m);
    end
end
