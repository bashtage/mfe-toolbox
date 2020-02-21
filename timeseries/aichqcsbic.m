function [aic, hqc, sbic] = aichqcsbic(errors,constant,p,q,X)
% Computes the Akaike, Hannan-Quinn and Schwartz/Bayes Information Criteria for an ARMA(P,Q) as 
% parameterized in ARMAXFILTER 
%
% USAGE:
%   [AIC] = aicsbic(ERRORS,CONSTANT,P,Q)
%   [AIC,HQC,SBIC] = aicsbic(ERRORS,CONSTANT,P,Q,X)
% 
% INPUTS:
%   ERRORS   - A T by 1 length vector of errors from the regression
%   CONSTANT - Scalar variable: 1 to include a constant, 0 to exclude
%   P        - Non-negative integer vector representing the AR orders to include in the model.
%   Q        - Non-negative integer vector representing the MA orders to include in the model.
%   X        - [OPTIONAL]  a T by K  matrix of exogenous variables.
% 
% OUTPUTS:
%   AIC       - The Akaike Information Criteria 
%   HQC       - The Hannan-Quinn Information Criteria
%   SBIC      - The Schwartz/Bayes Information Criteria
% 
% COMMENTS:
%   This is a helper for ARMAXFILTER and uses the same inputs, CONSTANT, P, Q and X.  ERRORS should
%   be the errors returned from a call to ARMAXFILTER with the same values of P, Q, etc. 
%
% EXAMPLES:
%   Compute AIC and SBIC from an ARMA
%       [parameters, LL, errors] = armaxfilter(y, constant, p, q);
%       [aic,sbic] = aicsbic(errors,constant,p,q)
% 
%  See also ARMAXFILTER, HETEROGENEOUSAR

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 10/19/2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3 || nargin>5
    error('3 to 5 inputs required')
end
if nargin==3
    q = [];
    X=[];
elseif nargin==4
    X=[];
end
%%%%%%%%%%%%%%%
% y
%%%%%%%%%%%%%%%
if size(errors,2) > 1 || length(errors)==1
    error('ERRORS series must be a column vector.')
elseif isempty(errors)
    error('ERRORS is empty.')
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
if max(p)>=(length(errors)-max(p))
    error('Too many lags in the AR.  max(P)<T/2')
end
if size(p,1)==1 && p==0
    p=[];
end
if length(unique(p))~=length(p)
    error('P must contain at most one of each lag')
end
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
if max(q)>=length(errors)
    error('Too many lags in the AR.  max(Q)<T')
end
if size(q,1)==1 && q==0
    q=[];
end
if length(unique(q))~=length(q)
    error('Q must contain at most one of each lag')
end
%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = length(errors);
seregression=sqrt(errors'*errors/T);

lp = length(unique(p));
lq = length(unique(q));
K=constant+lp+lq+size(X,2);

aic = log(seregression^2) + 2*K/T;
hqc = log(seregression^2) + 2*K*log(log(T))/T;
sbic = log(seregression^2) + log(T)*K/T;