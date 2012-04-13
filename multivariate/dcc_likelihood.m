function [ll,lls,Rt] = dcc_likelihood(parameters,data,dataAsym,m,l,n,R,N,backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate)
% Likelihood for estimation of scalar DCC(m,n) and ADCC(m,l,n) multivarate volatility models 
% with with TARCH(p,o,q) or GJRGARCH(p,o,q) conditional variances
%
% USAGE:
%  [LL,LLS,RT] = dcc_likelihood(PARAMETERS,DATA,DATAASYM,M,L,N,R,N,BACKCAST,BACKCASTASYM,STAGE,COMPOSITE,ISJOINT,ISINFERENCE,GSCALE,UNIVARIATE)
%
% INPUTS:
%   PARAMETERS   - Vector of ADCC parameters, and possibly volatility and intercept parameters, 
%                    depending on other inputs
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
%   M            - Order of symmetric innovations in DCC model
%   L            - Order of asymmetric innovations in ADCC model
%   N            - Order of lagged correlation in DCC model
%   R
%   N
%   BACKCAST     - K by K  matrix to use for back casting symetric terms
%   BACKCASTASYM - K by K  matrix to use for back casting asymetric terms
%   STAGE        - Integer, either 2 or 3 indicating whether 2-stage ro 3-stage estimator is being used
%   COMPOSITE    - Integer, one of 0 (None, use QMLE), 1 (Use diagonal composite) or 2 (full composite)
%   ISJOINT      - Boolean indicating whether PARAMETERS includes volatility parameters
%   ISINFERENCE  - Boolean indicating whether likelihood is used for making inference, in which case
%                    no transformations are made to parameters.
%   GSCALE       - K by 1 vector used in 2-stage to scale the intercept.  See DCC.
%   UNIVARIATE   - Cell array of structures containing information needed to compute volatilities.
%
% OUTPUTS:
%   LL           - The log likelihood evaluated at the PARAMETERS
%   LLS          - A T by 1 vector of log-likelihoods
%   HT           - A [K K T] dimension matrix of conditional correlations
%
% COMMENTS:
%
% See also DCC

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012


% Modes:
% 1-stage: Univariate, R (corr vech), DCC (N-scale ignored)
% 2-stage (Estimation): R (transformed corr vech), DCC  (N-scale treated as constants)
% 2-stage (Inference): Univariate, R (corr vech), DCC (N-scale ignored)
% 3-stage (Estimation): DCC
% 3-stage (Inferecne): Univariate, R (corr vech), N (vech), DCC

[k,~,T] = size(data);
offset = 0;
% Parse Parameters
if stage==1 || isJoint
    count = 0;
    for i=1:k
        u = univariate{i};
        count = count + u.p+u.o+u.q+1;
    end
    garchParameters = parameters(1:count);
    offset = offset + count;
    computeVol = true;
else
    computeVol = false;
end
if stage<=2 || isJoint
    count = k*(k-1)/2;
    R = parameters(offset + (1:count));
    offset = offset + count;
end
if stage==3 && isJoint && l>0 % N only if 3 stage, joint and asymmetric
    count = k*(k+1)/2;
    N = parameters(offset + (1:count));
    % Transform N to be a matrix
    N = ivech(N);
    offset = offset + count;
end
a = parameters(offset + (1:m));
g = parameters(offset + (m+1:m+l));
b = parameters(offset + (m+l+1:m+l+n));

% Compute volatilities
H = ones(T,k);
if computeVol
    offset = 0;
    for i=1:k
        u = univariate{i};
        count = u.p+u.o+u.q+1;
        volParameters = garchParameters(offset + (1:count));
        offset = offset+count;
        ht = tarch_core(u.fdata,u.fIdata,volParameters,u.back_cast,u.p,u.o,u.q,u.m,u.T,u.tarch_type);
        H(:,i) = ht(u.m+1:u.T);
    end
    stdData = zeros(k,k,T);
    stdDataAsym = zeros(k,k,T);
    for t=1:T
        h = sqrt(H(t,:));
        stdData(:,:,t) = data(:,:,t)./(h'*h);
        stdDataAsym(:,:,t) = dataAsym(:,:,t)./(h'*h);
    end
    logdetH = sum(log(H),2);
else
    stdData = data;
    stdDataAsym = dataAsym;
    logdetH = zeros(T,1);
end

% Transfor R & N, if needed
if stage<=2 || isJoint
    % Transform R
    if isInference
        R = corr_ivech(R);
    else
        R = z2r(R);
    end
end

% Compute intercept
if stage==3
    intercept = R*(1-sum(a)-sum(b))-N*sum(g);
else
    scale = (1-sum(a)-sum(b)) - gScale*sum(g);
    scale = sqrt(scale);
    intercept = R.*(scale*scale');
end
% Check eigenvalues?

if composite == 0
    likconst = k*log(2*pi);
else
    likconst  = 2*log(2*pi);
end
I = eye(k);
Qt = zeros(k,k,T);
Rt = zeros(k,k,T);
lls = zeros(T,1);
for t=1:T
    Qt(:,:,t) = intercept;
    for i = 1:m
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + a(i)*stdData(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + a(i)*backCast;
        end
    end
    for i = 1:l
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + g(i)*stdDataAsym(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + g(i)*backCastAsym;
        end
    end
    for i = 1:n
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + b(i)*Qt(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + b(i)*backCast;
        end
    end
    q = sqrt(diag(Qt(:,:,t)));
    Rt(:,:,t) = Qt(:,:,t)./ (q*q');
    if composite == 0
        lls(t) = 0.5*(likconst + logdetH(t) + log(det(Rt(:,:,t))) + sum(diag((Rt(:,:,t)\I)*stdData(:,:,t))));
    elseif composite == 1
        scale = (k-1);
        for i=1:k-1
            % FIXME : Sclaar optimization of CL badly needed
            j = i + 1;
            lls(t) = lls(t) + 0.5*(likconst + sum(log(H([i j]))) + log(det(Rt([i j],[i j],t))) + sum(diag((Rt([i j],[i j],t)\I)*stdData([i j],[i j],t))))/scale;
        end
    else
        scale = k*(k-1)/2;
        for i=1:k
            for j=i+1:k
                % FIXME : Scalar optimization of CL badly needed
                lls(t) = lls(t) + 0.5*(likconst + sum(log(H([i j]))) + log(det(Rt([i j],[i j],t))) + sum(diag((Rt([i j],[i j],t)\I)*stdData([i j],[i j],t))))/scale;
            end
        end
    end
end
ll = sum(lls);

if isnan(ll) || ~isreal(ll) || isinf(ll)
    ll = 1e7;
end