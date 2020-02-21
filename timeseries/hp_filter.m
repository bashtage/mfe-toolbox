function [trend,cyclic] = hp_filter(y,lambda)
% Hodrick-Prescott filtering of multiple time series
%
% USAGE:
%   [TREND,CYCLIC] = hp_filter(Y,LAMBDA)
%
% INPUTS:
%   Y      - A T by K matrix of data to be filtered.
%   LAMBDA - Positive, scalar integer containing the smoothing parameter of the HP filter.
%
% OUTPUTS:
%   TREND  - A T by K matrix containing the filtered trend
%   CYCLIC - A T by K matrix containing the filtered cyclic component
%
% COMMENTS:
%   The cyclic component is simply the original data minus the trend, CYCLIC = Y - TREND.  1600 is
%   the recommended value of LAMBDA for Quarterly Data while 14400 is the recommended value of LAMBDA
%   for monthly data.
%
% EXAMPLES:
%   Load US GDP data
%       load GDP
%   Standard HP Filter with lambda = 1600
%       [trend, cyclic] = hp_filter(log(GDP),1600)
%
%  See also BKFILTER, BEVERIDGENELSON

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 19/10/2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(y)>2 || size(y,1)<1
    error('Y must be a T by K matrix with T>=4');
end
if ~isscalar(lambda) || lambda<0
    error('LAMBDA must be a positive scalar.');
end
lambdaUpperBound = 1e10;
if lambda>lambdaUpperBound
    warning('MFEToolbox:Stability',['HP_FILTER may not be accurate for very large values of LAMBDA.  Values above ' num2str(lambdaUpperBound) ' force a pure linear trend.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=size(y,1);
if lambda<lambdaUpperBound
    % Handly small T cases
    switch T
        case {1,2}
            Gamma = eye(T);
        case 3
            Gamma(1,1) = (1+lambda);
            Gamma(1,2) = -2*lambda;
            Gamma(1,3) = lambda;
            Gamma(2,1) = -2*lambda;
            Gamma(2,2) = 1 + 4 * lambda;
            Gamma(2,3) = -2*lambda;
            Gamma(3,3) = Gamma(1,1);
            Gamma(3,2) = Gamma(1,2);
            Gamma(3,1) = Gamma(1,3);
        otherwise
            Gamma = spalloc(T,T,5*(T-4)+2*4+2*3);
            weigths = [lambda -4*lambda (1+6*lambda) -4*lambda lambda];
            
            for i=3:T-2;
                Gamma(i,i-2:i+2) = weigths; %#ok<SPRIX>
            end
            rowOneWeights = [1+lambda -2*lambda lambda];
            rowTwoWeights = [-2*lambda 1+5*lambda -4*lambda lambda];
            Gamma(1,1:3) = rowOneWeights;
            Gamma(2,1:4) = rowTwoWeights;
            Gamma(T,T-2:T) = fliplr(rowOneWeights);
            Gamma(T-1,T-3:T) = fliplr(rowTwoWeights);
    end
    trend = Gamma\y;
    cyclic = y-trend;
else
    x = [ones(T,1) (1:T)'];
    trend = x*(x\y);
    cyclic = y - trend;
end
