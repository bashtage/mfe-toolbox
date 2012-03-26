function [V,lagsused]=covvar(data,maxlag,method)
% Long-run covariance estimation using the VAR-based method of DenHaan and Levin (1996)
%
%   USAGE:
%     V=covvar(DATA)
%     [V,LAGSUSED]=covvar(DATA,MAXLAG,METHOD)
%
%   INPUTS:
%     DATA    - T by K vector of dependent data
%     MAXLAG  - Non-negative integer containing the maximum lag length to use.  If empty or not
%                 included, MAXLAG=min(1.2*floor(T^(1/3)),floor(T/K)) is used
%     METHOD  - An integer value 1, 2, 3, 4 or 5. 2 is the DEFAULT.
%                 1 - Use lags 1 to MAXLAG to compute the VAR based covariance
%                 2 - [DEFAULT] Use up to MAXLAG but select using SIC
%                 3 - Use up to MAXLAG but select using AIC
%                 4 - Use up to MAXLAG but select using SIC and a global search
%                 5 - Use up to MAXLAG but select using AIC and a global search
%
%   OUTPUTS:
%     V        - A K by K covariance matrix estimated using the VAR based esitmator
%     LAGSUSED - A row vector indicating the vlags used in computing the VAR-based long-run
%                  covariance.  If empty, no lags selected, and the estimated covariance is
%                  identical to the usual covariance.
%
%   COMMENTS:
%    If options 4 or 5 are used, all combinations of lags between none and MAXLAG are tried (e.g. [0
%    (constant)],[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3] ...).  These two options are only practically
%    viable for smallish MAXLAG (~<10 to 20)  and small K

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 5/1/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T,K]=size(data);
if nargin==1
    maxlag=min(floor(T/K),floor(1.2*T^(1/3)));
    method=2;
elseif nargin==2
    method=2;
end
if isempty(maxlag)
    maxlag=min(floor(T/K),floor(1.2*T^(1/3)));
end
if isempty(method)
    method=2;
end
if ~ismember(method,1:5)
    error('MATHOD must be a scalar between 1 and 5.')
end
if floor(maxlag)~=maxlag || maxlag<0 || maxlag>floor(T/K)
    error('MAXLAG must be a non-negative integer with MAXLAG<=floor(T/K).')
end
if ndims(data)>2
    error('DATA must be a T by K matrix of data.')
end
switch method
    case {1,2,4}
        ICtype=1;
    case {3,5}
        ICtype=2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the indices to use
if method==1
    indices{1}=1:maxlag;
elseif method==2 || method==3
    indices=cell(maxlag+1,1);
    for i=1:maxlag
        indices{i}=1:i;
    end
    indices{maxlag+1}=[];
else
    indices=cov_VAR_dec2bin(maxlag);
end

% Need to use the same amount of 'y' data for each to make SIC comparable
indep=cell(K,1);
Y=zeros(T-maxlag,K);
X=zeros(T-maxlag,K*maxlag);

for i=1:K;
    [dep,indep{i}]=newlagmatrix(data(:,i),maxlag,0);
    Y(:,i)=dep;
end
index=1;
for j=1:maxlag
    for i=1:K;
        temp=indep{i};
        X(:,index)=temp(:,j);
        index=index+1;
    end
end
X=[ones(size(X,1),1) X];

% Now that the X and Y matrices are set up, I can loop over the indices to
% run the regressions.  I need to save the SIC/AIC
N=length(indices);
IC=zeros(N,1);
T2=(T-maxlag);
for i=1:N
    P=length(indices{i});
    if ~isempty(indices{i});
        cols=repmat(indices{i}-1,K,1)*K+repmat((0:K-1)',1,P)+2;
    else
        cols=[];
    end
    regressors=X(:,[1 cols(:)']);
    B=regressors\Y;
    e=Y-X(:,[1 cols(:)'])*B;
    covE=e'*e/T2;
    if ICtype==1
        IC(i)=log(det(covE))+log(T2)/T2*(P*K^2+K);
    else
        IC(i)=log(det(covE))+2/T2*(P*K^2+K);
    end
end

% Finally pick the best
[temp,ICpos]=min(IC);

% Use the maximum amount of data to estimate the covariance
if ~isempty(indices{ICpos})
    maxlag=max(indices{ICpos});
    Y=zeros(T-maxlag,K);
    X=zeros(T-maxlag,K*maxlag);
    for i=1:K;
        [dep,indep{i}]=newlagmatrix(data(:,i),maxlag,0);
        Y(:,i)=dep;
    end
    index=1;
    for j=1:maxlag
        for i=1:K;
            temp=indep{i};
            X(:,index)=temp(:,j);
            index=index+1;
        end
    end
    X=[ones(size(X,1),1) X];
    T2=(T-maxlag);
    i=ICpos;
    P=length(indices{i});
    if ~isempty(indices{i});
        cols=repmat(indices{i}-1,K,1)*K+repmat((0:K-1)',1,P)+2;
    else
        cols=[];
    end
    regressors=X(:,[1 cols(:)']);
    B=regressors\Y;
    e=Y-X(:,[1 cols(:)'])*B;
    covE=e'*e/T2;
    B=B(2:P*K+1,:);
    B=reshape(B',[K,K,P]);
    A=(eye(K)-sum(B,3))^(-1);
    V=A*covE*A';
else
    % If empty then it's is the usual covariance estimator
    V=cov(Y)*((T-1)/T);
end
lagsused=indices{ICpos};



function indices=cov_VAR_dec2bin(maxlag)
% This is a helper function that creates the indices for the search
indices=cell(2^maxlag-1,1);
default=1:maxlag;
for p=1:(2^maxlag-1)
    rem = p;
    selector=false(1,maxlag);
    for i=maxlag-1:-1:0
        if floor(rem/(2^i))
            selector(i+1)=true;
            rem=rem-2^i;
        end
    end
    indices{p}=default(selector);
end
indices{2^maxlag}=[];
