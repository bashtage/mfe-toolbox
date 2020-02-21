function [yhattph, ytph, forerr, ystd]=arma_forecaster(y, parameters, constant, p, q, r, h, seregression, holdBack)
% Produces h-step ahead forecasts from ARMA(P,Q) models starting at some point in the sample, R, and
% ending at the end of the sample.  Also shifts the data to align y(t+h) with y(t+h|t) in slot t,
% computes the theoretical forecast standard deviation (assuming homoskedasticity) and the forecast
% errors.
%
% USAGE:
%  [YHATTPH] = arma_forecaster(Y,PARAMETERS,CONSTANT,P,Q,R,H)
%  [YHATTPH,YTPH,FORERR,YSTD] = arma_forecaster(Y,PARAMETERS,CONSTANT,P,Q,R,H,SEREGRESSION)
%
% INPUTS:
%   Y            - A column of data
%   CONSTANT     - Scalar variable: 1 if the model includes a constant, 0 to exclude
%   P            - Non-negative integer vector representing the AR orders included in the model.
%   Q            - Non-negative integer vector representing the MA orders included in the model.
%   R            - Length of sample used in estimation.  Sample is split up between R and P, where
%                    the first R (regression) are used for estimating the model and the remainder
%                    are used for prediction (P) so that R+P=T.
%   H            - The forecast horizon
%   SEREGRESSION - [OPTIONAL] The standard error of the regression.  Used to compute confidence
%                    intervals.  If omitted, SEREGRESSION is set to 1.
%
% OUTPUTS:
%  YHATTPH       - h-step ahead forecasts of Y.  The element in position t of YHATTPH is the time t
%                    forecast of Y(t+h).  The first R elements of YHATTPH are NaN.  The next T-R-H
%                    are pseudo in-sample forecasts while the final H are out-of-sample.
%  YTPH          - Value of original data at time t+h shifted to position t.  The first R elements
%                    of YTPH are NaN.  The next T-R-H are the values y(R+H),...,y(T), and the final
%                    H are NaN since there is no data available for comparing to the final H forecasts.
%  FORERR        - The forecast errors, YHATTPH-YTPH
%  YSTD          - The theoretical standard deviation of the h-step ahead forecast (homoskedasticity)
%
% COMMENTS:
% Values not relevant for the forecasting exercise have NaN returned.
%
% See also ARMAXFILTER, HETEROGENEOUSAR

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007

% Steps:
% 1. Produce the fit series and get estimated residuals
% 2. Using the estimated residuals, produce the H-step ahead forecast
% 3. Shift the y's to t+h
% 4. Compute the forecast e
% 5. Compute the forecast std



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<7 || nargin>9
    error('7 to 9 inputs required.')
end
if nargin<9
    holdBack = [];
end
%%%%%%%%%%%%%%%
% y
%%%%%%%%%%%%%%%
if size(y,2) > 1 || length(y)==1
    error('y series must be a column vector.')
elseif isempty(y)
    error('y is empty.')
end
T=length(y);
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
%maxp=max(p);
if length(unique(q))~=length(q)
    error('P must contain at most one of each lag')
end
np = length(p(p~=0));
% lp = length(p);
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
if length(unique(q))~=length(q)
    error('Q must contain at most one of each lag')
end
nq = length(q(q~=0));
%lq = length(q);
%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%
if ~constant && (isempty(p) || max(p)==0) && (isempty(q) || max(q)==0)
    error('At least one of CONSTANT, P or Q must be nonnegative')
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1')
end
%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%
if size(parameters,2)>size(parameters,1)
    parameters=parameters';
end
if size(parameters,2)~=1
    error('PARAMETERS must be a column vector')
end
if  isempty(parameters) || length(parameters)~=(constant + np + nq)
    error('length of PARAMETERS input must compatible with the constant, AR and MA lags')
end
%%%%%%%%%%%%%
% r
%%%%%%%%%%%%%
if length(r)~=1 || floor(r)~=r || r<=0 || r>T
    error('R must be a postive scalar less that or equal to T')
end
%%%%%%%%%%%%%
% h
%%%%%%%%%%%%%
if length(h)~=1 || floor(h)~=h || h<=0
    error('H must be a postive scalar')
end
%%%%%%%%%%%%%
% seregression
%%%%%%%%%%%%%
if nargin==8
    if length(seregression)~=1 || seregression<=0
        error('SEREGRESSION must be a postive number.')
    end
else
    seregression=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Everything must be an ARMA.  If either/o the AR or MA is 0, set it to a
% single lag and set it's coefficient to 0.  Makes coding much easier.
if isempty(p) || all(p==0)
    p=1;
    if constant==1
        parameters=[parameters(1); 0 ;parameters(2:length(parameters))];
    else
        parameters=[parameters(1); 0 ;parameters(1:length(parameters))];
    end
end

if isempty(q) || all(q==0)
    q=1;
    parameters=[parameters;0];
end
% Redefine these here
maxp=max(p);
maxq=max(q);


yorig = y;
% Steps:
% 1. Produce the fit series and get estimated residuals

zerosToPad = max(maxq-maxp,0);
m = max(maxp,maxq);
if ~isempty(holdBack)
    if holdBack > maxp
        m = m + (holdBack-maxp);
    end
end

y=[zeros(zerosToPad,1);y];

np = length(p);
nq = length(q);

errors = armaxerrors(parameters,p,q,constant,y,[],m,ones(size(y)));
if constant
    constantp=parameters(1);
else
    constantp = 0;
end
if np>0
    arparameters = parameters(constant+1:constant+np);
end
if nq>0
    maparameters = parameters(constant+np+1:constant+np+nq);
end

% 2. Using the estimated residuals, produce the H-step ahead forecast
yhattph = NaN*ones(T,1);

m = maxq;
for t=r:T
    ytemp = zeros(T+m+h,1);
    etemp = zeros(T+m+h,1);
    ytemp(m+1:m+t) = y(1+zerosToPad:t+zerosToPad);
    etemp(m+1:m+t) = errors(1+zerosToPad:t+zerosToPad);
    for i=1:h
        % This is r+m+i for the first one
        % The +1 is because we should be using r+m for the first one
        % It may be the case that r-max(q) is <=0, in which case a backcast of
        % should be used, but it never should be that case that r-max(p)
        % should be <=0 Note: This case has not yet been fixed here.
        ytemp(t+m+i) = constantp + arparameters'*ytemp(t+m-p+i) + maparameters'*etemp(t+m-q+i);
    end
    % I want to time t forecast in slot t, but we have m 'new' observations
    % now, so 1 is 1+m, 2 is 2+m an so on
    yhattph(t)=ytemp(t+m+h);
end
% 3. Shift the y's to t+h
ytph = [NaN*ones(r-1,1); yorig(r+h:length(yorig)); NaN*ones(h,1)];
% NaN to remind that there was no forecast then!

% 4. Compute the forecast e
forerr = ytph - yhattph;

% 5. Compute the forecast std
newar = zeros(max(h,maxp),1);
newar(p) = arparameters;
newma = zeros(max(h,maxq),1);
newma(q) = maparameters;

new = [1;newma];
old = zeros(h,1);

for i=1:h
    old(h-i+1) = new(i);
    for j=1:i-1
        old(h-i+1) = old(h-i+1) + newar(j)*old(h-i+1+j);
    end
end
ystd = sqrt(sum(old.^2))*seregression;