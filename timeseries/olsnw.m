function [b, tstat, s2, vcvnw, R2, Rbar, yhat] = olsnw(y,x,c,nwlags)
% Linear regression estimation with Newey-West HAC standard errors.
%
% USAGE:
%   [B,TSTAT,S2,VCVNW,R2,RBAR,YHAT] = olsnw(Y,X,C,NWLAGS)
%
% INPUTS:
%   Y       - T by 1 vector of dependent data
%   X       - T by K vector of independent data
%   C       - [OPTIONAL] 1 or 0 to indicate whether a constant should be included (1: include
%               constant). The default value is 1.
%   NWLAGS  - [OPTIONAL] Number of lags to included in the covariance matrix estimator. If omitted
%               or empty, NWLAGS = floor(T^(1/3)). If set to 0 estimates White's Heteroskedasticity
%               Consistent variance-covariance.
%
% OUTPUTS:
%   B       - A K(+1 is C=1) vector of parameters.  If a constant is included, it is the first parameter
%   TSTAT   - A K(+1) vector of t-statistics computed using Newey-West HAC standard errors
%   S2      - Estimated error variance of the regression, estimated using Newey-West with NWLAGS
%   VCVNW   - Variance-covariance matrix of the estimated parameters computed using Newey-West
%   R2      - R-squared of the regression.  Centered if C=1
%   RBAR    - Adjusted R-squared. Centered if C=1
%   YHAT    - Fit values of the dependent variable
%
% COMMENTS:
%   The model estimated is Y = X*B + epsilon where Var(epsilon)=S2.
%
% EXAMPLES:
%   Regression with automatic BW selection
%       b = olsnw(y,x)
%   Regression without a constant
%       b = olsnw(y,x,0)
%   Regression with a pre-specified lag-length of 10
%       b = olsnw(y,x,1,10)
%   Regression with White standard errors
%       b = olsnw(y,x,1,0)
%
% See also OLS


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        c=1;
        nwlags = [];
    case 3
        nwlags = [];
    case 4
        % nothing
    otherwise
        error('2 to 4 inputs only')
end
% Check y
if size(y,1)<size(y,2)
    y=y';
end
if size(y,2)~=1
    error('Y must be a column vector')
end
T = size(y,1);
% Check X
if ~isempty(x)
    if size(x,1)~=T
        error('X must have the same number of rows as Y.')
    end
    if size(x,2)>T
        error('The number of columns of X must be grater than or equal to T')
    end
    if rank(x)<size(x,2)
        error('X is rank deficient')
    end
end
% Check c
if isempty(c)
    c=1;
end
if ~ismember(c,[0 1])
    error('C must be either 0 or 1.')
end
% Check C,X
if c==0 && size(x,2)==0
    error('The model must include a constant or at least one X')
end
if isempty(nwlags) || nargin==3
    nwlags = floor(T^(1/3));
end
%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%
% Figure out the size of X
K=size(x,2);
% Add a constant if needed
if c==1
    x=[ones(T,1) x];
    K=K+1;
elseif ~ismember(c,[0 1]) % Check for an error
    error('C must be either 1 or 0');
end


% Compute beta
b=x\y;
% And the fit values
yhat=x*b;
% And the errors
epsilon=y-yhat;
% And the estimated residual variance
s2=covnw(epsilon,nwlags,0);
% Once we have E, White's VCV is easy
scores = x.*repmat(epsilon,1,K);
XpXi =(x'*x/T)^(-1);
vcvnw = XpXi * covnw(scores,nwlags,0) * XpXi /T;

% Compute t-stats using White
tstat=b./sqrt(diag(vcvnw));

% Finally the R2 and Rbar
if c==1 % Use centered if a constant is included
    ytilde=y-mean(y);
    R2=1 - (epsilon'*epsilon)/(ytilde'*ytilde);
    Rbar=1 - (epsilon'*epsilon)/(ytilde'*ytilde) * (T-1)/(T-K);
else % Use non centered versions
    R2=1 - (epsilon'*epsilon)/(y'*y);
    Rbar=1 - (epsilon'*epsilon)/(y'*y) * (T-1)/(T-K);
end