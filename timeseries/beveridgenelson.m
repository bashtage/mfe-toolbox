function [trend,cyclic] = beveridgenelson(y,parameters,constant,p,q)
% Beveridge-Nelson decomposition for a trending time series.
%
% USAGE:
%   [TREND,CYCLIC] = beveridgenelson(y,parameters,constant,p,q)
%
% INPUTS:
%   Y          - T by 1 column of trending data data
%   PARAMETERS - 1 by CONSTANT + length(P) + length(Q) vector of parameters estimated on diff(Y)
%   CONSTANT   - Scalar variable: 1 the the fit model included a constant
%   P          - Non-negative integer vector representing the AR orders of the fit model
%   Q          - Non-negative integer vector representing the MA orders of the fit model
%
% OUTPUTS:
%   TREND      - T by 1 vector of containing the trend.  First max(P)+1 observations are the same as Y
%   CYCLIC     - T by 1 vector of containing the cyclic component.  First max(P)+1 observations are
%                  the same as Y
%
% COMMENTS:
%   Uses the exact BN decomposition as presented in
%
%   Paul Newbold, "Precise and efficient computation of the Beveridge-Nelson decomposition of
%     economic time series", Journal of Monetary Economics, Volume 26, Issue 3, December 1990, Pages
%     453-457.
%
%   This decomposition is based on a demeaned ARMA model for the short-run dynamics.  If CONSTANT = 1,
%   this function will demeand the data using the model implied long-run mean (the constant
%   parameter divided by the sum of the AR coefficients).
%
% EXAMPLES:
%   BN decomposition for Real US GDP using an AR(1) for the SR component
%       load GDP
%       parameters = armaxfilter(diff(lnGDP),1,1);
%       [trend,cyclic] = beveridgenelson(lnGDP,parameters,1,1);
%   BN decomposition for Real US GDP using an ARMA(1,4) for the SR component
%       parameters = armaxfilter(diff(lnGDP),1,1,1:4);
%       [trend,cyclic] = beveridgenelson(lnGDP,parameters,1,1,1:4);
%
% See also ARMAXFILTER

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 11/1/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        p=[];
        q=[];
    case 4
        q=[];
    case 5
        % Nothing
    otherwise
        error('3 to 5 inputs required.');
end
if size(y,2)>size(y,1)
    y=y';
end
if size(y,2)~=1
    error('Y must be a T by 1 vector.');
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1.')
end
np = length(p);
nq = length(q);
if np>0
    maxp = max(p);
else
    maxp = 0;
end
if nq>0
    maxq = max(q);
else
    maxq = 0;
end
if length(parameters)~=(constant+np+nq)
    error('PARAMETERS must have CONSTANT + length(P) + length(Q) parameters.');
end
% Remove the constant, if any
if constant
    constant = parameters(1);
    parameters = parameters(2:length(parameters));
    if np>0
        longRunMean = constant/(1-sum(parameters(1:np)));
    else
        longRunMean = constant;
    end
end
deltaY = diff(y) - longRunMean;
if maxp>0
    A = zeros(maxp);
    A(1,p) = parameters(1:np);
    for i=2:maxp
        A(i,i-1) = 1;
    end
    Aeigs = abs(eig(A));
    if  max(Aeigs)>=1
        error('The model for the short-run dynamics is not stationary. Stationarity is required for the BN decomposition to be well defined.');
    end
else
    A = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




T = size(deltaY,1);
% Construct the errors
m = max(maxp,maxq);
zerosToPad = max(maxq-maxp,0);
deltaYAugmentedForMA = [zeros(zerosToPad,1);deltaY];
errors = armaxerrors(parameters,p,q,0,deltaYAugmentedForMA,[],m,ones(size(deltaYAugmentedForMA)));
errors = errors(zerosToPad+1:length(errors));

% Initialize the trend and cyclic component
trend = zeros(T+1,1);
trend(1:maxp+1) = y(1:maxp+1);
cyclic = zeros(T+1,1);
deltaYForecast = zeros(T+m,1);
errorsForecast = zeros(1,T+m);
% Initialize the coefficient needed
e = [1 zeros(1,maxp-1)];
longRunScale = (eye(maxp) - A)^(-1)*A;
longRunScale = e*longRunScale;
for t = maxp + 1 : T
    % Construct maxq forecasts, then use companion matrix to compute sum
    deltaYForecast(:) = 0;
    errorsForecast(:) = 0;
    deltaYForecast(1:t) = deltaY(1:t);
    errorsForecast(1:t) = errors(1:t);
    for h = 1 : maxq
        for j=1:np
            deltaYForecast(t+h) = deltaYForecast(t+h) + parameters(j) * deltaYForecast(t+h-p(j));
        end
        for j=1:nq
            if (t+h-q(j))>0
                deltaYForecast(t+h) = deltaYForecast(t+h) + parameters(np+j) * errorsForecast(t+h-q(j));
            end
        end
    end
    % Select the correct forecasts
    cumForecast = sum(deltaYForecast(t+1:t+maxq));
    if np>0
        deltaYHat=deltaYForecast(t+maxq-maxp+1:t+maxq);
        cumForecast = cumForecast + longRunScale*deltaYHat;
    end
    trend(t+1) = y(t+1) + cumForecast;
    cyclic(t+1) = -cumForecast;
end
