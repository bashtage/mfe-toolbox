function [b, tstat, s2, vcv, vcvwhite, R2, Rbar, yhat] = ols(y,x,c)
% Linear regression estimation with homoskedasticity and White heteroskedasticity robust standard
% errors.
%
% USAGE:
%   [B,TSTAT,S2,VCV,VCV_WHITE,R2,RBAR,YHAT] = ols(Y,X,C)
%
% INPUTS:
%   Y        - N by 1 vector of dependent data
%   X        - N by K vector of independent data
%   C        - 1 or 0 to indicate whether a constant should be included (1: include constant)
%
% OUTPUTS:
%   B        - A K(+1 is C=1) vector of parameters.  If a constant is included, it is the first
%                parameter.
%   TSTAT    - A K(+1) vector of t-statistics computed using White heteroskedasticity robust
%                standard errors.
%   S2       - Estimated error variance of the regression.
%   VCV      - Variance covariance matrix of the estimated parameters.  (Homoskedasticity assumed)
%   VCVWHITE - Heteroskedasticity robust VCV of the estimated parameters.
%   R2       - R-squared of the regression.  Centered if C=1.
%   RBAR     - Adjusted R-squared. Centered if C=1.
%   YHAT     - Fit values of the dependent variable
%
% COMMENTS:
%   The model estimated is Y = X*B + epsilon where Var(epsilon)=S2
%
% EXAMPLES:
%   Estimate a regression with a constant
%       b = ols(y,x)
%   Estimate a regression without a constant
%       b = ols(y,x,0)
%
% See also OLSNW

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || nargin>3
    error('2 or 3 inputs only')
end
% Check y
if size(y,1)<size(y,2)
    y=y';
end
if size(y,2)~=1
    error('Y must be a column vector')
end
N = size(y,1);
% Check X
if ~isempty(x)
    if size(x,1)~=N
        error('X must have the same number of rows as Y.')
    end
    if size(x,2)>N
        error('The number of columns of X must be grater than or equal to T')
    end
    if rank(x)<size(x,2)
        error('X is rank deficient')
    end
end
% Check c
if nargin==2
    c = 1;
end
if ~ismember(c,[0 1])
    error('C must be either 0 or 1.')
end
% Check C,X
if c==0 && size(x,2)==0
    error('The model must include a constant or at least one X')
end
%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%

% Some basic error checking
K=size(x,2);
% Add a constant if needed
if rank(x) ~= K
    error('X is must be of full rank.')
end
if c==1
    x=[ones(N,1) x];
    if rank(x) ~= K+1
        error('X appears to contains a constant column.  Use C to add a constant.')
    end
    K=K+1;
end

% Compute beta
b=x\y;
% early return
if nargout==1
    return
end
% And the fit values
yhat=x*b;
% And the errors
epsilon=y-yhat;
% And the estimated residual variance
s2=epsilon'*epsilon/(N-K);
% The estimated parameter covariance
vcv=s2*(x'*x)^(-1);
% Compute covariance of scores
scores = x.*repmat(epsilon,1,K);
XeeX = scores'*scores;
% White's VCV is easy
XpXi = (x'*x)^(-1);
vcvwhite= XpXi * XeeX * XpXi;
% Compute t-stats using White
tstat=b./sqrt(diag(vcvwhite));

% Finally the R2 and Rbar
if c==1 % Use centered if a constant is included
    ytilde=y-mean(y);
    R2=1 - (epsilon'*epsilon)/(ytilde'*ytilde);
    Rbar=1 - (epsilon'*epsilon)/(ytilde'*ytilde) * (N-1)/(N-K);
else % Use non centered versions
    R2=1 - (epsilon'*epsilon)/(y'*y);
    Rbar=1 - (epsilon'*epsilon)/(y'*y) * (N-1)/(N-K);
end