function [noiseVariance, debiasedNoiseVariance, IQEstimate, noiseEstimateOomen] = realized_noise_estimate(price, time, timeType, options)
% Estimation of the optimal bandwidth to use when estimating the quadratic variation using a Realized Kernel
%
% USAGE:
%   [NOISEVARIANCE,DEBIASEDNOISEVARIANCE,IQESIQESTIMATETIMATE,NOISEVARIANCEOOMEN] 
%                        = realized_kernel_select_bandwidth(PRICE,TIME,TIMETYPE,OPTIONS)
%
% INPUTS:
%   PRICE       - m by 1 vector of prices
%   TIME        - m by 1 vector of times where TIME(i) corresponds to PRICE(i)
%   TIMETYPE    - String describing the way times are measured
%               'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%               'seconds' Time measured in seconds past midnight on the first day
%               'unit'  Unit normalized date format, e.g. .1, .234, .9. Unit normalized times are
%                 more general than the other types and can be applied to data from more than one
%                 calendar day
%   OPTIONS     - Realized kernel options structure.  See help realized_kernel_options
%
% OUTPUTS:
%   NOISEVARIANCE         - Bandi-Russel noise variance estimate 
%   DEBIASEDNOISEVARIANCE - Debiased Bandi-Russel noise variance estimate using the adjustment of BNHLS
%   IQESTIMATE            - An estimate of the lower bound of the IQ based on low-freuqency returns
%   NOISEVARIANCEOOMEN    - Estimate using Oomen (2006) alternative AC(1) estimator
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL.  See Barndorf-Nielsen, Hansen, Lunde and Shephard
%   (2008) for details about the optimal selection of bandwidth
%
%  See also REALIZED_KERNEL, REALIZED_PRICE_FILTER, REALIZED_KERNEL_WEIGHTS
%  REALIZED_KERNEL_SELECT_LAG_LENGTH, REALIZED_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=4
    error('Seven inputs required.')
end
if size(price,2)>size(price,1)
    price=price';
end
if size(price,2)>1
    error('PRICE must be a m by 1 vector.')
end
if size(time,2)>size(time,1)
    time=time';
end
if any(diff(time)<0)
    error('TIME must be sorted and increasing')
end
if size(time,2)>1 || length(time)~=length(price)
    error('TIME must be a m by 1 vector.')
end
time = double(time);
timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','matlab'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''Matlab''.');
end

% List of flat top kernels
flatTopKernelList = {'bartlett','twoscale','2ndorder','epanechnikov',...
    'cubic','multiscale','5thorder','6thorder','7thorder','8thorder','parzen',...
    'th1','th2','th5','th16'};

% List of non flat top kernels
nonFlatTopKernelList = {'nonflatparzen','qs','fejer','thinf','bnhls'};

% Combined kernel list
kernelList = [flatTopKernelList nonFlatTopKernelList];

if ~isfield(options,'medFrequencyKernel') || ~ismember(options.medFrequencyKernel,kernelList)
    error('KERNEL must be a field of OPTIONS and one of the listed types.')
end

medFrequencySamplingType = options.medFrequencySamplingType;
medFrequencySamplingInterval = options.medFrequencySamplingInterval;
medFrequencyKernel = options.medFrequencyKernel;
medFrequencyBandwidth = options.medFrequencyBandwidth;
if ~isscalar(medFrequencySamplingInterval) || floor(medFrequencySamplingInterval)~=medFrequencySamplingInterval || medFrequencySamplingInterval<=0
    error('MEDIUMFREQUENCYTIME must be a postive integer scalar value.');
end

IQEstimationSamplingType = options.IQEstimationSamplingType;
IQEstimationSamplingInterval = options.IQEstimationSamplingInterval;
if ~isscalar(IQEstimationSamplingInterval) || floor(IQEstimationSamplingInterval)~=IQEstimationSamplingInterval || IQEstimationSamplingInterval<=0
    error('LOWFREQUENCYTIME must be a postive integer scalar value.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To compute the "optimal" bandwidth an estimate of the noise and the QV using low frequency data
lowFrequencyRealizedVariance = realized_variance(price,time, timeType , IQEstimationSamplingType, IQEstimationSamplingInterval);

noiseVarianceSamplingType = options.noiseVarianceSamplingType;
noiseVarianceSamplingInterval = options.noiseVarianceSamplingInterval;

% Compute the noiseVariance using the debiased version of Bandi-Russell, using all prices
noiseVariance = realized_variance(price, time, timeType , noiseVarianceSamplingType , noiseVarianceSamplingInterval);
% Compute the effective 'n' to use in the Bandi-Russel estimator
noiseFilteredPrice = realized_price_filter(price, time, timeType , noiseVarianceSamplingType , noiseVarianceSamplingInterval);
if options.useAdjustedNoiseCount
    % Only count the number of non-zero returns
    n = sum(diff(noiseFilteredPrice)~=0);
else
    n = length(noiseFilteredPrice) - 1;
end
noiseVariance = noiseVariance/(2*n);

noiseReturns = diff(log(noiseFilteredPrice));
n = length(noiseReturns);
noiseEstimateOomen = -1/(n-1) * noiseReturns(1:end-1)'*noiseReturns(2:end);

% Require a realized kernel estimate to adjust the variance
medFrequencyOptions = realized_options('kernel');
medFrequencyOptions.kernel = medFrequencyKernel;
medFrequencyOptions.bandwidth = medFrequencyBandwidth;
medFrequencyOptions.endTreatment = 'stagger';

medFreuqencyRealizedKernel = realized_kernel(price, time, timeType, medFrequencySamplingType, medFrequencySamplingInterval, medFrequencyOptions);
medFrequencyRealizedVariance = realized_variance(price,time, timeType , medFrequencySamplingType, medFrequencySamplingInterval);
% Use the BNHLS corrected noise estimate
debiasedNoiseVariance = exp( log(noiseVariance) - medFreuqencyRealizedKernel/medFrequencyRealizedVariance);
IQEstimate = lowFrequencyRealizedVariance^2;