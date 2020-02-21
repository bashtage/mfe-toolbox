function [parameters,stderr,tstat,pval,const,conststd,r2,errors,s2,paramvec,vcv] = vectorar(y,constant,lags,het,uncorr)
% Estimate a Vector Autoregression and produce the parameter variance-covariance matrix under a
% variety of assumptions on the covariance of the errors:
%   * Conditionally Homoskedastic and Uncorrelated
%   * Conditionally Homoskedastic but Correlated
%   * Heteroskedastic but Conditionally Uncorrelated
%   * Heteroskedastic and Correlated
%
% USAGE:
%   [PARAMETERS]=vectorar(Y,CONSTANT,LAGS)
%   [PARAMETERS,STDERR,TSTAT,PVAL,CONST,CONSTSTD,R2,ERRORS,S2,PARAMVEC,VEC]
%        = vectorar(Y,CONSTANT,LAGS,HET,UNCORR)
%
% INPUTS:
%   Y             - A T by K matrix of data
%   CONSTANT      - Scalar variable: 1 to include a constant, 0 to exclude
%   LAGS          - Non-negative integer vector representing the VAR orders to include in the model.
%   HET           - [OPTIONAL] A scalar integer indicating the type of covariance estimator
%                      0 - Homoskedastic
%                      1 - Heteroskedastic [DEFAULT]
%   UNCORR        - [OPTIONAL] A scalar integer indicating the assumed structure of the error covariance
%                     matrix
%                      0 - Correlated errors  [DEFAULT]
%                      1 - Uncorrelated errors
%
% OUTPUTS:
%   PARAMETERS    - Cell structure containing K by K matrices in the position of the indicated in
%                     LAGS.  For example if LAGS = [1 3], PARAMETERS{1} would be the K by K
%                     parameter matrix for the 1st lag and PARAMETERS{3} would be the K by K matrix
%                     of parameters for the 3rd lag
%   STDERR        - Cell structure with the same form as PARAMETERS containing parameter standard
%                     errors estimated according to UNCORR and HET
%   TSTAT         - Cell structure with the same form as PARAMETERS containing parameter t-stats
%                     computed using STDERR
%   PVAL          - P-values of the parameters
%   CONST         - K by 1 vector of constants
%   CONSTSTD      - K by 1 vector standard errors corresponding to constant
%   R2            - K by 1 vector of R-squares
%   ERRORS        - K by T vector of errors
%   S2            - K by K matrix containing the estimated error variance
%   PARAMVEC      - K*((# lags) + CONSTANT)  by 1 vector of estimated parameters.  The first (# lags
%                     + CONSTANT) correspond to the first row in the usual var form:
%                   [CONST(1) P1(1,1) P1(1,2) ... P1(1,K) P2(1,1) ... P2(1,K) ...]
%                   The next (# lags + CONSTANT) are the 2nd row
%                   [CONST(1) P1(2,1) P1(2,2) ... P1(2,K) P2(2,1) ... P2(2,K) ...]
%                   and so on through the Kth row
%                   [CONST(K) P1(K,1) P1(K,2) ... P1(K,K) P2(K,1) ... P2(K,K) ...]
%   VCV           - A K*((# lags) + CONSTANT) by K*((# lags) + CONSTANT) matrix of estimated
%                     parameter covariances computed using HET and UNCORR
% COMMENTS:
%   Estimates a VAR including any lags.
%   y(:,t)' = CONST + P(1) * y(:,t-1) + P(2)*y(:,t-2) + ... + P(1)*y(:,t-K)'
%
%   where P(j) are K by K parameter matrices and CONST is a K by 1 parameter matrix (if CONSTANT==1)
%
% EXAMPLE:
%   To fit a VAR(1) with a constant
%       parameters = vectorar(y,1,1)
%   To fit a VAR(3) with no constant
%       parameters = armaxfilter(y,0,[1:3])
%   To fit a VAR that includes lags 1 and 3 with a constant
%       parameters = armaxfilter(y,1,[1 3])
%
% See also IMPULSERESPONSE, GRANGERCAUSE,  VECTORARVCV

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.1    Date: 3/1/2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3 || nargin>5
    error('3 to 5 inputs required')
end
if nargin==3
    het=1;
    uncorr=0;
elseif nargin==4
    uncorr=0;
end
% Check Y
if ndims(y)~=2
    error('Y must be T by K')
end
if any(var(y)==0)
    error('At least one component of Y is constant. All components of Y are required to be non-constant (i.e. have some variation).')
end
% Check constant
if ~isscalar(constant)
    error('CONSTANT must be either 0 or 1')
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be either 0 or 1')
end
% Check lags
if ndims(lags)~=2
    error('LAGS must be a vector of positive integers containing lags to include')
end
if size(lags,1)>size(lags,2)
    lags = lags';
end
if ~all(lags>=0)
    error('LAGS must be a vector of nonnegative integers containing lags to include')
end
if ~all(floor(lags)==lags)
    error('LAGS must be a vector of positive integers containing lags to include')
end
if length(lags)~=length(unique(lags))
    error('LAGS must be a vector of unique elements')
end
% Check if lags contains 0 and some positive
if any(lags==0) && any(lags>0)
    lags=setdiff(lags,0);
end
lags=unique(lags);
if lags==0
    if ~constant
        error('If LAGS=0, a constant must be included')
    end
end
% Check het
if ~isscalar(het)
    error('HET must be either 0 or 1')
end
if ~ismember(het,[0 1])
    error('HET must be either 0 or 1')
end
% Check uncorr
if ~isscalar(uncorr)
    error('UNCORR must be either 0 or 1')
end
if ~ismember(uncorr,[0 1])
    error('UNCORR must be either 0 or 1')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First estimate the parameters, the
T=size(y,1);
K=size(y,2);
m=max(lags);
if isempty(m)
    m = 0;
end
ylags=cell(K,1);
X = zeros(T-m,K*sum(lags>0)+constant);
Y = y(m+1:T,:);
index=1;

if constant
    X(:,1)=ones(T-m,1);
    index=index+1;
end
if any(lags>0)
    for k=1:K
        [~,ylags{k}]=newlagmatrix(y(:,k),m,0);
    end
    p=length(lags);
    for i=1:p
        for k=1:K
            X(:,index)=ylags{k}(:,lags(i));
            index=index+1;
        end
    end
end

% Each column contains the parameters for a single y
paramvec = X\Y;
errors = Y-X*paramvec;
s2=errors'*errors/T;
% VCV estimation
vcv=vectorarvcv(X,errors,het,uncorr);


% Need to reshape the parameters into a formatted structure
% The parameters are ordered constant, x(-1) y(-1) x(-2) y(-2)
% First reshape to be a rectangle
Np=numel(paramvec)/K;
tempParam = paramvec';
tempStd   = reshape(sqrt(diag(vcv)),Np,K)';
if constant
    const=tempParam(:,1);
    conststd=tempStd(:,1);
    tempParam=tempParam(:,2:size(tempParam,2));
    tempStd=tempStd(:,2:size(tempStd,2));
else
    const=[];
    conststd=[];
end

parameters=cell(m,1);
stderr=cell(m,1);
tstat=cell(m,1);
pval=cell(m,1);
if any(lags>0)
    for i = 1:length(lags)
        parameters{lags(i)} = tempParam(:,(i-1)*K+1:i*K);
        stderr{lags(i)} = tempStd(:,(i-1)*K+1:i*K);
        tstat{lags(i)} = parameters{lags(i)}./stderr{lags(i)};
        pval{lags(i)} = 2 - 2 * normcdf(abs(tstat{lags(i)}));
    end
end

if constant
    ytil = y-repmat(mean(y),T,1);
    r2=1-sum(errors.^2)./sum(ytil.^2);
else
    r2=1-sum(errors.^2)./sum(y.^2);
end
r2=r2';
paramvec=paramvec(:);