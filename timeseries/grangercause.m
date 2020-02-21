function [stat,pval,statAll,pvalAll]=grangercause(y,constant,lags,het,uncorr,inference)
% Granger causality testing with a variance-covariance matrix estimated under a variety of
% assumptions on the covariance of the errors:
%   * Conditionally Homoskedastic and Uncorrelated
%   * Conditionally Homoskedastic but Correlated
%   * Heteroskedastic but Conditionally Uncorrelated
%   * Heteroskedastic and Correlated
%
% USAGE:
%   [STAT] = grangercause(Y,CONSTANT,LAGS)
%   [STAT,PVAL,STATALL,PVALALL] = grangercause(Y,CONSTANT,LAGS,HET,UNCORR,INFERENCE)
%
% INPUTS:
%   Y             - A T by K matrix of data
%   CONSTANT      - Scalar variable: 1 to include a constant, 0 to exclude
%   LAGS          - Non-negative integer vector representing the VAR orders to include in the model.
%   HET           - [OPTIONAL] A scalar integer indicating the type of covariance estimator
%                      0 - Homoskedastic
%                      1 - Heteroskedastic [DEFAULT]
%   UNCORR        - [OPTIONAL] A scalar integer indicating the assumed structure of the error
%                     covariance matrix
%                      0 - Correlated errors  [DEFAULT]
%                      1 - Uncorrelated errors
%   INFERENCE     - [OPTIONAL] Inference method
%                      1 - Likelihood ratio [DEFAULT]
%                      2 - LM test
%                      3 - Wald test
%
% OUTPUTS:
%   STAT          - K by K matrix of Granger causality statistics computed using the specified
%                     covariance estimator and inference method. STAT(i,j) corresponds to a test that
%                     y(i) is not caused by y(j)
%   PVAL          - K by K matrix of p-values corresponding to STAT
%   STATALL       - K by 1 vector of Granger causality statistics computed using the specified
%                     covariance estimator and inference method. STATALL(i) corresponds to a test that
%                     y(i) is not caused by any y(j), j neq i
%   PVALALL       - K by 1 vector of p-values corresponding to STATALL
%
% COMMENTS:
%   Granger causality tests based on a VAR including any lags.
%
%   y(:,t)' = CONST + P(1) * y(:,y-1) + P(2)*y(:,y-2) + ... + P(1)*y(:,t-K)'
%
%   where P(j) are K by K parameter matrices and CONST is a K by 1 parameter matrix (if CONSTANT==1)
%
% EXAMPLE:
%   Conduct GC testing in a VAR(1) with a constant
%        parameters = grangercause(y,1,1)
%   Conduct GC testing in a VAR(3) with no constant
%        parameters = grangercause(y,0,[1:3])
%   Conduct GC testing in a VAR that includes lags 1 and 3 with a constant
%        parameters = grangercause(y,1,[1 3])
%
% See also VECTORAR, VECTORARVCV

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0    Date: 1/1/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3 || nargin>6
    error('3 to 6 inputs required')
end
if nargin==3
    het=1;
    uncorr=0;
    inference=1;
elseif nargin==4
    uncorr=0;
    inference=1;
elseif nargin==5
    inference=1;
end
% Check Y
if ndims(y)~=2
    error('Y must be T by K')
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
if ~all(lags>0)
    error('LAGS must be a vector of positive integers containing lags to include')
end
if ~all(floor(lags)==lags)
    error('LAGS must be a vector of positive integers containing lags to include')
end
if length(lags)~=length(unique(lags))
    error('LAGS must be a vector of unique elements')
end
lags=sort(lags);
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
% Inference mathod
if ~isscalar(inference)
    error('INFERENCE must be either 0 or 1')
end
if ~ismember(inference,[1 2 3])
    error('INFERENCE must be either 0 or 1')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inference Methods
% 1. LR MLE : Relatively easy, need to iteratively depete the correct
% columns and reestimate.  Use a LR test
% 2. Wald: Just compute the VCV, then figure out how to select the correct
%       rows, and you are done
% 3. LM: Have to re-estimate under the null and use those errors to compute
%       the scores.  Test is based on scores.
% 4. LR Robust: Just like LM above but using errors estimated under the
%       alternative


% LR MLE
% First generate the unrestricted errors, then iterate across series and
% drop regressors

T=size(y,1);
K=size(y,2);
P=length(lags);
m=max(lags);
ylags=cell(K,1);
for k=1:K
    ytemp = [ones(m,1)*mean(y(:,k));
        y(:,k)];
    [nothing,ylags{k}]=newlagmatrix(ytemp,m,0);
end

X=zeros(T,K+constant);
index=1;
if constant
    X(:,1)=ones(T,1);
    index=index+1;
end
p=length(lags);
for i=1:p
    for k=1:K
        X(:,index)=ylags{k}(:,lags(i));
        index=index+1;
    end
end

% Each column contains the parameters for a single y
paramvec = X\y;
errors = y-X*paramvec;
s2=errors'*errors/T;
e=cell(K,K);
for i=1:K
    for j=1:K
        % i will be the LHS variable, j will be the lags
        tempX=X;
        allcols = 1:size(X,2);
        % Drop columns j, K+j, 2K+j,...P(K-1)*j (+constant if present)
        if constant
            drop = (j:K:size(X,2)-1)+constant;
        else
            drop = (j:K:size(X,2));
        end
        remain = setdiff(allcols,drop);
        beta = tempX(:,remain)\y(:,i);
        % These errors are estiamted under the null and can be used for LR and LM testing
        e{i,j} = y(:,i)-tempX(:,remain)*beta;
    end
end
% Construct the regressions needed to test the "all" hypothesis
eall = cell(K,1);
for i=1:K
    tempX = X;
    % Keep columns i, K+i, 2K+i, P(K-1)*i (+constant if present)
    if constant
        remain = [constant (i:K:(size(X,2)-1))+constant];
    else
        remain = i:K:size(X,2);
    end
    beta = tempX(:,remain)\y(:,i);
    eall{i} = y(:,i)-tempX(:,remain)*beta;
end

Np = p*K+constant;
stat = zeros(K);
statAll = zeros(K,1);
% Homoskedastic Likelihood Ratio
if inference==1 && het==0
    for i=1:K
        for j=1:K
            e2=errors;
            e2(:,i)=e{i,j};
            sR=e2'*e2/T;
            stat(i,j)=(T-P*K^2+P)*(log(det(sR))-log(det(s2)));
        end
    end
    for i=1:K
        e2=errors;
        e2(:,i) = eall{i};
        sR=e2'*e2/T;
        statAll(i)=(T-P*K^2+(K-1)*P)*(log(det(sR))-log(det(s2)));
    end
    % Heteroskedasticity robust LR;  key here is score covariance is computed
    % using errors estimated under the alternative
elseif inference==1 && het==1
    for i=1:K
        for j=1:K
            e2=errors;
            e2(:,i)=e{i,j};
            X2=repmat(X,1,K);
            e2=reshape(repmat(e2,Np,1),T,K*Np);
            s=X2.*e2;
            sbar=mean(s);
            % Estimate S using the errors, not e2
            S=vectorarscorecov(errors,X,het,uncorr,K,T,Np);
            stat(i,j)=(T-P*K^2+P)*sbar*S^(-1)*sbar';
        end
    end
    for i=1:K
        e2=errors;
        e2(:,i)=eall{i};
        X2=repmat(X,1,K);
        e2=reshape(repmat(e2,Np,1),T,K*Np);
        s=X2.*e2;
        sbar=mean(s);
        % Estimate S using the errors, not e2
        S=vectorarscorecov(errors,X,het,uncorr,K,T,Np);
        statAll(i)=(T-P*K^2+(K-1)*P)*sbar*S^(-1)*sbar';
    end
    % All LM tests since function will handle;  the key here is that the errors
    % are computed under the null
elseif inference==2
    for i=1:K
        for j=1:K
            e2=errors;
            e2(:,i)=e{i,j};
            X2=repmat(X,1,K);
            e2=reshape(repmat(e2,Np,1),T,K*Np);
            s=X2.*e2;
            sbar=mean(s);
            e2=errors;
            e2(:,i)=e{i,j};
            S=vectorarscorecov(e2,X,het,uncorr,K,T,Np);
            stat(i,j)=(T-P*K^2+P)*sbar*S^(-1)*sbar';
        end
    end
    for i=1:K
        e2=errors;
        e2(:,i)=eall{i};
        X2=repmat(X,1,K);
        e2=reshape(repmat(e2,Np,1),T,K*Np);
        s=X2.*e2;
        sbar=mean(s);
        e2=errors;
        e2(:,i)=eall{i};
        S=vectorarscorecov(e2,X,het,uncorr,K,T,Np);
        statAll(i)=(T-P*K^2+(K-1)*P)*sbar*S^(-1)*sbar';
    end
    % Wald tests
elseif inference==3
    XpXi = ((X'*X)/T)^(-1);
    Ainv = kron(eye(K),XpXi);
    B = vectorarscorecov(errors,X,het,uncorr,K,T,Np);
    V = Ainv*B*Ainv;
    %Vinv=V^(-1);
    % Now I have to cleverly select the parameters to set to 0 and then
    % compute the wald tests
    % Easier to iterate over the j's in the inside loop
    for i=1:K
        for j=1:K
            temp=zeros(size(paramvec))';
            % Select the columns of the parameters to test
            if constant
                pl = (j:K:size(X,2)-1)+constant;
            else
                pl = (j:K:size(X,2));
            end
            
            temp(i,pl)=1;
            temp=temp';
            temp=temp(:);
            pl = find(temp);
            p = paramvec;
            p = p(:);
            stat(i,j)=(T-P*K^2+P)*(p(pl)'*V(pl,pl)^(-1)*p(pl));
        end
    end
    for i=1:K
        temp=zeros(size(paramvec))';
        % Select the columns of the parameters to test
        if constant
            remain = [constant (i:K:size(X,2)-1)+constant];
        else
            remain = (i:K:size(X,2)-1);
        end
        pl = setdiff(1:size(X,2),remain);
        temp(i,pl)=1;
        temp=temp';
        temp=temp(:);
        pl = find(temp);
        p = paramvec;
        p = p(:);
        statAll(i)=(T-P*K^2+(K-1)*P)*(p(pl)'*V(pl,pl)^(-1)*p(pl));
    end
end
% All stats have the same dist
pval=1-chi2cdf(stat,length(lags));
pvalAll = 1-chi2cdf(statAll,(K-1)*length(lags));



function S=vectorarscorecov(errors,X,het,uncorr,K,T,Np)
s2=errors'*errors/T;
XpX=X'*X/T;
if ~het && uncorr
    S=kron(diag(diag(s2)),XpX);
elseif ~het && ~uncorr
    S=kron(s2,XpX);
elseif het && uncorr
    X2=repmat(X,1,K);
    e2=reshape(repmat(errors,Np,1),T,K*Np);
    s=X2.*e2;
    s=s-repmat(mean(s),T,1);
    S=zeros(Np*K);
    for i=1:K
        sel=(i-1)*Np+1:i*Np;
        temp = s(:,sel);
        S(sel,sel)=temp'*temp/T;
    end
elseif het && ~uncorr
    X2=repmat(X,1,K);
    e2=reshape(repmat(errors,Np,1),T,K*Np);
    s=X2.*e2;
    s=s-repmat(mean(s),T,1);
    S=s'*s/T;
end
