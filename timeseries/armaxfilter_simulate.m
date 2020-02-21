function [y,errors]=armaxfilter_simulate(T,const,ar,ARparams,ma,MAparams,X,Xparams)
% ARMAX(P,Q) simulation with normal errors.  Also simulates AR, MA and ARMA models.
%
% USAGE:
%   AR:
%   [Y,ERRORS] = armaxfilter_simulate(T,CONST,AR,ARPARAMS)
%   MA:
%   [Y,ERRORS] = armaxfilter_simulate(T,CONST,0,[],MA,MAPARAMS)
%   ARMA:
%   [Y,ERRORS] = armaxfilter_simulate(T,CONST,AR,ARPARAMS,MA,MAPARAMS);
%   ARMAX:
%   [Y,ERRORS] = armaxfilter_simulate(T,CONST,AR,ARPARAMS,MA,MAPARAMS,X,XPARAMS);
%
% INPUTS:
%   T        - Length of data series to be simulated  OR
%                T by 1 vector of user supplied random numbers (e.g. rand(1000,1)-0.5)
%   CONST    - Value of the constant in the model.  To omit, set to 0.
%   AR       - Order of AR in model.  To include only selected lags, for example t-1 and t-3, use 3
%                and set the coefficient on 2 to 0 
%   ARPARAMS - AR by 1 vector of parameters for the AR portion of the model
%   MA       - Order of MA in model.  To include only selected lags of the error, for example t-1
%                and t-3, use 3 and set the coefficient on 2 to 0 
%   MAPARAMS - MA by 1 vector of parameters for the MA portion of the model
%   X        - T by K matrix of exogenous variables
%   XPARAMS  - K by 1 vector of parameters on the exogenous variables
%
% OUTPUTS:
%   Y        - A T by 1 vector of simulated data
%   ERRORS   - The errors used in the simulation
%
% COMMENTS:
%   The ARMAX(P,Q) model simulated is:
%      y(t) = const + arp(1)*y(t-1) + arp(2)*y(t-2) + ... + arp(P) y(t-P) +
%                   + ma(1)*e(t-1)  + ma(2)*e(t-2)  + ... + ma(Q) e(t-Q)
%                   + xp(1)*x(t,1)  + xp(2)*x(t,2)  + ... + xp(K)x(t,K)
%                   + e(t)
% EXAMPLES:
%   Simulate an AR(1) with a constant
%       y = armaxfilter_simulate(500, .5, 1, .9)
%   Simulate an AR(1) without a constant
%       y = armaxfilter_simulate(500, 0, 1, .9)
%   Simulate an ARMA(1,1) with a constant
%       y = armaxfilter_simulate(500, .5, 1, .95, 1, -.5)
%   Simulate a  MA(1) with a constant
%       y = armaxfilter_simulate(500, .5, [], [], 1, -.5)
%   Simulate a seasonal MA(4) with a constant
%       y = armaxfilter_simulate(500, .5, [], [], 4, [.6 0 0 .2])
%
% See also ARMAXFILTER, HETEROGENEOUSAR

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        ar=0;
        ma=0;
        ARparams=[];
        MAparams=[];
        X=[];
        Xparams=[];
    case 4
        ma=0;
        MAparams=[];
        X=[];
        Xparams=[];
    case 6
        X=[];
        Xparams=[];
    case 8
        % Nothing
    otherwise
        error('The number of inputs must be 2, 4, 6 or 8.')
end

if length(T)~=1 %User supplied errors
    e=T;
    T=length(e);
    if min(size(e))>1
        error('If using user supplied errors, these much be a column vector (e.g. T=randn(100,1))')
    end
end

if ~isempty(ar) && (length(ARparams)<ar || min(size(ARparams))>1)
    error('Incorrect number of AR parameters')
end

if size(ARparams,2)<size(ARparams,1)
    ARparams=ARparams';
end

if ~isempty(ma) && (length(MAparams)<ma || min(size(MAparams))>1)
    error('Incorrect number of MA parameters')
end

if size(MAparams,2)<size(MAparams,1)
    MAparams=MAparams';
end

if length(Xparams)<size(X,2)
    error('Incorrect number of X parameters.  XPARAMS should be K by 1.')
end

if ~isempty(X) && length(X)~=T
    error('X should be length T');
end

if any([ma,ar]<0) || any([length(ma) length(ar)]>1)
    error('MA and AR must all be non negative scalars.')
end

if ~isscalar(const)
    error('CONST must be a scalar')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(X)
    X=0;
    Xparams=0;
else
    if size(Xparams,2)<size(Xparams,1)
        Xparams=Xparams';
    end
end

% Set the startup
B=2000;
%First, take case of the error situation
if exist('e','var')
    e=[e(ceil(rand(B,1)*T));e];
else
    e=randn(T+B,1);
end

%Easy to set up the exogenous variables
meanX=mean(X);
exog=[repmat(meanX,length(e)-length(X),1)*Xparams';
    X*Xparams']+const;

if ma>0
    e=[zeros(ma,1); e];
    [e,elag]=newlagmatrix(e,ma,0);
    exog=exog+elag*MAparams'+e;    
else
    elag=0;
    exog=exog+e;
end


%So all of the exogenous and the innovation are contained in a single
%parameter, we can now loop over the data.  First a  naive starting value for y

if abs(sum(ARparams))<1
    y0=(const+meanX*Xparams')/(1-sum(ARparams));
else % explosive or unit root
    y0=0;
end

y=repmat(y0,length(e),1);

if ar>0
    for i=ar+1:length(e);
        y(i)=exog(i)+ARparams*y(i-1:-1:i-ar);
    end
else
    y=exog;
end

%Fix the size
y=y(B+1:T+B);
errors=e(B+1:T+B);