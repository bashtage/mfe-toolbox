function [rpav]=realized_preaveraged_bipower_variation(price,time,timeType,samplingType,samplingInterval,options)
% Estimated integrated variance using Preaveraged Bipower Variation
%
% USAGE:
%   [RPBV] = realized_preaveraged_bipower_variation(PRICE)
%   [RPBV] = realized_preaveraged_bipower_variation(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%   [RPBV] = realized_preaveraged_bipower_variation(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,OPTIONS)
%
% INPUTS:
%   PRICE            - m by 1 vector of high frequency prices
%   TIME             - [OPTIONAL] m by 1 vector of times where TIME(i) corresponds to PRICE(i).
%   TIMETYPE         - [OPTIONAL] String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS.mmm, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGTYPE     - [OPTIONAL] String describing the type of sampling to use when
%                        filtering PRICE
%                        'CalendarTime' - Sample in calendar time using observations separated by
%                          SAMPLINGINTERVAL seconds
%                        'CalendarUniform' - Sample in calendar time using SAMPLINGINTERVAL
%                          observations spread uniformly between TIME(1) and TIME(m)
%                        'BusinessTime' - Sample in business (tick) time using observation separated
%                          by SAMPLINGINTERVAL ticks
%                        'BusinessUniform' - Sample in business (tick) time using observations
%                          uniformly spaced in business time.
%                        'Fixed' - Sample at specific points in time. When using fixed,
%                          SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                          as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - [OPTIONAL] Scalar integer or n by 1 vector whose meaning depends on the
%                         selected SAMPLINGTYPE
%   OPTIONS           - [OPTIONAL] Preaveraged Realized Variance option structure initialized by calling
%                         realized_options('Preaveraging'). See help realized_options for a description of
%                         available options.
%
% OUTPUTS:
%   RPBV - Preaveraged bipower variation estimate
%
% COMMENTS:
%  Follows Christensen, Oomen and Podolski (2014) most closely, with the
%  noise variance estimator used in Hautsch and Podolski (2013)
%
% EXAMPLES:
%
%  See also PREALIZED_PREAVERAGED_VARIANCE, REALIZED_OPTIONS, REALIZED_KERNEL, 
%  REALIZED_NOISE_ESTIMATE, REALIZED_VARIANCE, REALIZED_VARIANCE_OPTIMAL_SAMPLING, 
%  REALIZED_RANGE, REALIZED_QUANTILE_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 2/27/2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        m = length(price);
        time = linspace(9.5*3600,16*3600,m)';
        timeType = 'seconds';
        samplingType = 'businesstime';
        samplingInterval = 1;
        options = realized_options('preaveraging');
    case 5
        options = realized_options('preaveraging');
    case 6
        % Nothing
    otherwise
        error('One, five or six inputs required.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need method to estimate omega, and options.

% Inserted to protect against inputing integer times
time = double(time);

% 1. Filter the rpice
filteredPrice = realized_price_filter(price,time,timeType,samplingType,samplingInterval);
returns = diff(log(filteredPrice));
% 2. Compute rpreaveraged returns
theta = options.theta;
m = length(returns);
K = ceil(theta * sqrt(m));
g = @(x) min(x,1-x);
w = g((1:(K-1))/K);

preav_returns = nan(m-K+2,1);
for i=1:m-K+2
    preav_returns(i) = w*returns(i:i+K-2);
end

% 3. Compute constants
psi_1 = K*sum((g((1:K)/K) - g((0:(K-1))/K)).^2);
psi_2 = 1/K*sum((g((1:(K-1))/K)).^2);
% 4. Estimate noise
[noiseVariance, ~, ~, noiseEstimateOomen] = realized_noise_estimate(price, time, timeType, options);

omega = noiseEstimateOomen;
if omega<0
    omega = noiseVariance;
end

% 5. Compute Bipower Variation
mu_1 = sqrt(2/pi);
const1 = m/(m-2*K+2);
const2 = 1/(K*psi_2*mu_1^2);

bias = psi_1/(theta.^2 * psi_2) * omega^2;

lead = preav_returns(K:end);
lag = preav_returns(1:end-K+1);

rpav = const1 * const2 * sum(abs(lead.*lag)) - bias;
