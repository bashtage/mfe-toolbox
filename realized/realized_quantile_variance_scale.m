function [scales,weights,covar]=realized_quantile_variance_scale(samplesperbin,quantiles,simulations,symmetric)
% Computes the scales needed for estimating the integrated variance using Realized Quantile
% Variance.  Also computes the weights of the optimal combination and non-scaled covariance from
% which the weights are derived.
%
% USAGE:
%   [SCALES,WEIGHTS,COVAR]=realized_quantile_variance_scale(SAMPLESPERBIN,QUANTILES,SIMULATIONS)
%
% INPUTS:
%   SAMPLESPERBIN    - Number of returns to use in each bin when computing the quantiles. NOTE: The
%                        number of returns produced by filtering according to SAMPLINGTYPE and
%                        SAMPLINGINTERVAL must be an integer multiple of SAMPLESPERBIN.
%   QUANTILES        - k by 1 vector of quantile values to use when computing RQ, must satisfy 0.5<QUANTILES<=1.
%                        Quantiles must produce values such that QUANTILES*SAMPLESPERBIN is an
%                        integer.  The simplest method to accomplish this is to specify QUANTILES as
%                        the ratio of the index number of the return to SAMPLESPERBIN (e.g.
%                        QUANTILES = [13 15 19]/20)
%   SIMULATIONS       - [OPTIONAL] Scalar integer.  Should usually be at least 1,000,000.  Default
%                         value is 10,000,000.
%   SYMMETRIC         - [OPTIONAL] Logical value indicating whether the symmetric estimator is being
%                         used. Default is false.
%
% OUTPUTS:
%   SCALES            - k by 1 vector of scales for standardized quantile based estimates of variance
%   WEIGHTS           - k by 1 vector of weights that produce the optimal combination (lowest variance)
%   COVAR             - k by k non-scaled covariance matrix of the quantile realized variance - one
%                         for each quantile used.  Forms the actual (scaled) covariance when
%                         multiplied by the integrated quarticity.
%
% COMMENTS:
%   Uses Monte Carlo integration with 10,000,000 simulations.  Trade offs between accuracy and time
%   can be made by altering SIMULATIONS.
%
%  See also REALIZED_QUANTILE_VARIANCE
%

%
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008


if nargin<4
    symmetric = false;
end
if nargin ==2
    simulations = 10000000;
end

% Save the state and set it to a specific value so that the simulated results will be reproducable
state0=randn('state');
randn('state',datenum('MAR-26-1974'))

k = length(quantiles);

if symmetric
    indices = round(samplesperbin*quantiles);
else
    indicesHigh = round( samplesperbin*quantiles);
    indicesLow = round(samplesperbin*(1-quantiles)+1 );
end



% Need to block eventually incase someone tries to compute this with a
% large number of samplesperbin
rows = min(floor(2^13/samplesperbin),simulations);
iter = ceil(simulations/rows);
remaining = simulations;

scales = zeros(1,k);
opMatrix  = zeros(k);
% The MC integration
for i=1:iter
    numThisIter = min(remaining,rows);
    remaining = remaining - rows;
    
    x=randn(numThisIter,samplesperbin);
    if symmetric
        x=sort(x.^2,2);
        x2=x(:,indices);
    else
        x=sort(x,2);
        x2=x(:,indicesHigh).^2 + x(:,indicesLow).^2;
    end
    scales = scales + sum(x2)/simulations;
    opMatrix = opMatrix + x2'*x2/simulations;
end

% Compute the covariance
covar = opMatrix - scales'*scales;
covar = samplesperbin*diag(1./scales)*covar*diag(1./scales);
% And the weights according to the GMVP
weights=covar^(-1)*ones(k,1)/(ones(1,k)*covar^(-1)*ones(k,1));

% Reset the state
randn('state',state0);