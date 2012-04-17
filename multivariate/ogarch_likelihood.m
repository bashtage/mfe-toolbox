function [ll,lls] = ogarch_likelihood(parameters,data,p,q,gjrType,backCast)
% Targetted log-likelihood for use in estimation of volatility models in OGARCH models
%
% USAGE:
%  [LL,LLS] = ogarch_likelihood(PARAMETERS,DATA,P,Q,GJRTYPE,BACKCAST)
%
% INPUTS:
%   PARAMETERS - P + Q by 1 vector of parameters
%   DATA       - A T by K matrix of zero mean residuals 
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance lags
%   GJRTYPE    - Scalar, either 1 (TARCH/AVGARCH) or 2 (GJRGARCH/GARCH/ARCH)
%   BACKCAST   - Value to use for back cating
%
% OUTPUTS:
%   LL         - The log likelihood computed at PARAMETERS
%
% COMMENTS:
%   Uses 1-sum(PARAMETERS) as the intercept, which is an identifying assumption in OGARCH models

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/15/2012

volParameters = [1-sum(parameters) parameters];
volParameters(volParameters<0) = 0;
v = tarch_core_simple(data,volParameters,backCast,0,p,0,q,gjrType);
lls = 0.5 * (log(2*pi) + log(v) + data.^2./v);
ll = sum(lls);