function [rv,rvDebiased,rvSS,rvDebiasedSS,diagnostics] = realized_variance_optimal_sampling(price,time,timeType,samplingType,samplingInterval,subsamples,options)
% Estimates realized variance using Bandi-Russell optimal sampling selection
%
% USAGE:
%   [RV,RVD,RVSS,RVSSD,DIAGNOSTICS] = realized_variance_optimal_sampling(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINVERVAL)
%   [RV,RVD,RVSS,RVSSD,DIAGNOSTICS] = realized_variance_optimal_sampling(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINVERVAL,SUBSAMPLES,OPTIONS)
%
% INPUTS:
%   PRICE            - m by 1 vector of high frequency prices
%   TIME             - m by 1 vector of times where TIME(i) corresponds to PRICE(i)
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGTYPE     - String describing the type of sampling to use when
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
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on SAMPLINGTYPE
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                         variance estimators to average with the original realized variance.
%                         Subsample realized variances are based on prices uniformly spaced between
%                         the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES=1
%                         will compute a subsample realized variance using the mid-point of the
%                         price sample points, 2 will use 1/3 and 2/3, and so on. In general this
%                         number should be small so the subsample estimators will be "sparse". If  
%                         the computed optimal sampling frequency is smaller than SUBSAMPLES, then 
%                         SUBSAMPLES is set to the compute optimal sampling frequency.
%   OPTIONS           - [OPTIONAL] Option structure initialized by realized_options.
%                         See help realized_options for a description of fields.
%
% OUTPUTS:
%   RV                - Realized variance estimated using the Bandi-Russell estimator of the optimal
%                         sampling frequency
%   RVD               - Debiased realized variance estimated using the Bandi-Russell estimator of
%                         the optimal sampling frequency for a debiased estimator
%   RVSS              - Subsample RV using the computed optimal samplig frequency
%   RVSSD             - Debiased version of subsample RV using the computed optimal samplig frequency
%   DIAGNOSTICS       - Structure with fields
%                         OPTIMALSAMPLES          - Optimal number of samples estimated using the
%                                                   Bandi-Russell methodology
%                         OPTIMALSAMPLESDEBIASED  - Optimal number of samples estimated using the
%                                                   Bandi-Russell methodology for use with a
%                                                   debiased estimator.
%                         SAMPLES                 - Number of sampled used.  Should be equal to
%                                                   OPTIMALNUMBEROFSAMPLES unless larger than the
%                                                   number of prices or smaller than 2
%                         SAMPLESDEBIASED         - Number of sampled used in debiased estimator.
%                                                   Should be equal to
%                                                   OPTIMALNUMBEROFSAMPLESDEBIASED unless larger
%                                                   than the number of prices or smaller than 2
%                         NOISEVARIANCE           - Estimated noise variance
%                         DEBIASEDNOISEVARIANCE   - Estimated noise variance bias adjusted
%                         IQESTIMATE              - Estimated IQ
%
% COMMENTS:
%   This function estimates the optimal number of samples to use when computing realized variance
%   using the method of Bandi & Russell (2008).  The role of SAMPLINGTYPE and SAMPLINGINTERVAL are
%   to set the MAXIMUM frequency of returns.  For example, if the maximum frequency was to be 15
%   seconds, set SAMPLINGTYPE='CalendarTime' and SAMPLINGINTERVAL=15.  To use all of the data,
%   set SAMPLINGTYPE='BusinessTime' and SAMPLINGINTERVAL=1.
%
% EXAMPLE:
%  % Optimal sampling in using all data
%  rvos = realized_variance_optimal_sampling(PRICE)
%  % Optimal sampling in business time with subsampling
%  [rvos,rvosSS] = realized_variance_optimal_sampling(PRICE,TIME,'Wall','BusinessTime',1,10)
%  % Optimal sampling in business time, using every 15th trade
%  rvos = realized_variance_optimal_sampling(PRICE,TIME,'Wall','BusinessTime',15)
%  % Optimal sampling where the maximum sampling frequency is 5 seonds
%  rvos = realized_variance_optimal_sampling(PRICE,TIME,'Wall','CalendarTime',5)
%
%  See also REALIZED_VARIANCE, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE,
%  REALIZED_PRICE_FILTER, REALIZED_THRESHOLD_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 1
        time = linspace(wall2seconds(93000),wall2seconds(160000),length(price));
        timeType = 'seconds';
        samplingType = 'businessuniform';
        samplingInterval = 1;
        subsamples = 1;
        options = realized_options('optimal sampling');        
    case 5
        subsamples = 1;
        options = realized_options('optimal sampling');        
    case 6
        options = realized_options('optimal sampling');
    case 7
        if isempty(subsamples)
            subsamples = 1;
        end
        % Nothing
    otherwise
        error('5 to 7 inputs required.')
end

errorMessage = realized_variance_optimal_sampling_parameter_check(price,time,timeType,samplingType,samplingInterval,subsamples,options);
if ~isempty(errorMessage )
    error(errorMessage)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inserted to protect against inputing integer times
time = double(time);
% Filter the price
[filteredPrice,filteredTime] = realized_price_filter(price,time,timeType,samplingType,samplingInterval);

% Compute the number of prices
m = size(filteredPrice,1);

% Next need to estimate IQ using RV squared
[noiseVariance, debiasedNoiseVariance, IQEstimate] = realized_noise_estimate(filteredPrice, filteredTime, timeType, options);
% Need to double these since these are the variance of the error to returns, not to prices which is
% what is returns from realized_noise_esitmate
noiseVariance = 2*noiseVariance;
debiasedNoiseVariance = 2*debiasedNoiseVariance;

% Need to estimate the IQ,  Variance of Noise, and 4th moment of noise (Eq. (16))
optimalNumberOfSamples = round((IQEstimate/(noiseVariance)^2)^(1/3));
if optimalNumberOfSamples>m
    warning('oxfordRealized:excessiveLags','WARNING')
    SamplesUsed = m;
elseif optimalNumberOfSamples<2
    warning('oxfordRealized:tooFewLags','WARNING')
    SamplesUsed = 2;
else
    SamplesUsed = optimalNumberOfSamples;
end
% Compute the realized variance using BT since the prices have been
% filtered
actualSubsamples = min(subsamples,SamplesUsed);
[rv,rvSS] = realized_variance(filteredPrice, filteredTime, timeType,'BusinessTime',SamplesUsed,actualSubsamples);

filteredPriceNoise = realized_price_filter(filteredPrice, filteredTime, timeType, options.noiseVarianceSamplingType, options.noiseVarianceSamplingInterval);
returns = diff(log(filteredPriceNoise));
if options.useAdjustedNoiseCount
    n = sum(returns~=0);
else
    n = length(returns);
end
noise4thPower = sum(returns.^4)/n;
noiseVariance = sum(returns.^2)/n;
% FIXME: Need to change way denominator is estimated
% Eq. 18
optimalNumberOfSamplesDebiased = round((2*IQEstimate / (2*noise4thPower - 3*noiseVariance^2))^(1/2));
if optimalNumberOfSamplesDebiased >m
    warning('oxfordRealized:excessiveLags','WARNING')
    samplesDebiased = m;
elseif optimalNumberOfSamplesDebiased <2
    warning('oxfordRealized:tooFewLags','WARNING')
    samplesDebiased = 2;
else
    samplesDebiased = optimalNumberOfSamplesDebiased;
end

% TODO See if it is possible to do a debiased using 1-m*/m type change
actualSubsamples = min(subsamples,samplesDebiased);
[rvDebiased,rvDebiasedSS] = realized_variance(filteredPrice,filteredTime,timeType,samplingType,samplesDebiased,actualSubsamples);
rvDebiased = rvDebiased - optimalNumberOfSamplesDebiased * noiseVariance;
rvDebiasedSS = rvDebiasedSS - optimalNumberOfSamplesDebiased * noiseVariance;
rvDebiased = rvDebiased/(1-optimalNumberOfSamplesDebiased/m);
rvDebiasedSS = rvDebiasedSS/(1-optimalNumberOfSamplesDebiased/m);

diagnostics.m = length(filteredPrice) - 1;
diagnostics.optimalSamples = optimalNumberOfSamples;
diagnostics.samples = SamplesUsed;
diagnostics.optimalNumberOfSamplesDebiased = optimalNumberOfSamplesDebiased;
diagnostics.samplesDebiased = samplesDebiased;
diagnostics.noiseVariance = noiseVariance;
diagnostics.debiasedNoiseVariance = debiasedNoiseVariance;
diagnostics.IQEstimate = IQEstimate;

function errorMessage = realized_variance_optimal_sampling_parameter_check(price,time,timeType,samplingType,samplingInterval,subsamples,options)
% Support function for realized_variance_optimal_sampling that does input validation
%
% USAGE:
%   [ERRORMESSAGE,OPTIONS] = realized_variance_optimal_sampling_parameter_check(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINVERVAL,OPTIONS)
%
% INPUTS:
%   See realized_variance_optimal_sampling
%
% OUTPUT:
%   ERRORMESSAGE - String containing a description of the error if one is detected.  Empty if no error.
%
% COMMENTS:
%   See realized_options for a description of the other fields in OPTIONS
%
%  See also REALIZED_VARIANCE_OPTIMAL_SAMPLING


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
errorMessage = [];

if size(price,2)>size(price,1)
    price=price';
end
if size(price,2)>1
    errorMessage = 'PRICE must be a m by 1 vector.';
    return
end
if size(time,2)>size(time,1)
    time=time';
end
if any(diff(time)<0)
    errorMessage = 'TIME must be sorted and increasing';
    return
end
if size(time,2)>1 || length(time)~=length(price)
    errorMessage = 'TIME must be a m by 1 vector.';
    return
end

timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    errorMessage = 'TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.';
    return;
end

samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    errorMessage = ('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.');
    return;
end


if ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        errorMessage = 'SUBSAMPLES must be a non-negative scalar.';
        return
    end
end


% SAMPLINGINTERVAL
m=size(price,1);
t0=time(1);
tT=time(m);
if ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform'})
    % Must be a scalar integer if timeType is seconds or wall
    if ismember(timeType,{'wall','seconds'})
        if ~isscalar(samplingInterval) || floor(samplingInterval)~=samplingInterval || samplingInterval<1
            error('SAMPLINGINTERVAL must be a positive integer for the SAMPLINGTYPE selected when using ''wall'' or ''seconds'' as TIMETYPE.')
        end
    else
        if ~isscalar(samplingInterval) || samplingInterval<0
            error('SAMPLINGINTERVAL must be a positive value for the SAMPLINGTYPE selected when using ''unit'' as TIMETYPE.')
        end
    end
else
    if size(samplingInterval,2)>size(samplingInterval,1)
        samplingInterval=samplingInterval';
    end
    if ~(any(samplingInterval>=t0) && any(samplingInterval<=tT))
        error('At least one sampling interval must be between min(TIME) and max(TIME) when using ''Fixed'' as SAMPLINGTYPE.')
    end
    if any(diff(samplingInterval)<=0)
        error('When using ''Fixed'' as SAMPLINGTYPE the vector of sampling times in SAMPLINGINTERVAL must be sorted and strictly increasing.')
    end
end


% Options

 
% List of flat top kernels
flatTopKernelList = {'bartlett','twoscale','2ndorder','epanechnikov',...
    'cubic','multiscale','5thorder','6thorder','7thorder','8thorder','parzen',...
    'th1','th2','th5','th16'};
 
% List of non flat top kernels
nonFlatTopKernelList = {'nonflatparzen','qs','fejer','thinf','bnhls'};
 
% Combined kernel list
kernelList = [flatTopKernelList nonFlatTopKernelList];
 
% Check fields for valid values
optionsFieldNames = fieldnames(options);
% Insert any missing fiedls
defaultOptions = realized_options('Optimal Sampling');
defaultFieldNames = fieldnames(defaultOptions);
missingFields = setdiff(defaultFieldNames,optionsFieldNames);
if ~isempty(missingFields)
    for i = 1:length(missingFields)
        options.(missingFields{i}) = defaultOptions.(missingFields{i});
    end
end
 
for i=1:length(optionsFieldNames)
    fieldName = optionsFieldNames{i};
    fieldValue = options.(fieldName);
    if ischar(fieldValue)
        fieldValue = lower(fieldValue);
        options.(fieldName) = fieldValue;
    end
 
    switch fieldName
        case {'medFrequencyKernel'}
            % Member of kernelList
            if ~ismember(fieldValue,kernelList)
                errorMessage = ['OPTIONS.' fieldName ' must be one of the listed types.'];
                return
            end
        case {'medFrequencyBandwidth'}
            % Non-negative scalar
            if ~isempty(fieldValue) && ~isnonnegativescalar(fieldValue)
                errorMessage = ['OPTIONS.' fieldName ' must a non-negative scalar.'];
                return
            end
        case {'useDebiasedNoise','useAdjustedNoiseCount'}
            % Logical or scalar
            if ~islogical(fieldValue) && ~ismember(fieldValue,[0 1])
                errorMessage = 'OPTIONS.useDebiasedNoise must be a logical value.';
                return
            end
        case {'IQEstimationSamplingType','medFrequencySamplingType','noiseVarianceSamplingType'}
            % One of the sampling types
            if ~ismember(fieldValue,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
                errorMessage = ['OPTIONS.' fieldName ' must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.'];
                return
            end
        case {'medFrequencySamplingInterval','noiseVarianceSamplingInterval','IQEstimationSamplingInterval'}
            % Non-negative scalar, less that 1 if timeType is unit
            if isempty(fieldValue) || ~isnonnegativescalar(fieldValue)
                errorMessage = ['OPTIONS.' fieldName ' must be a non-negative scalar between 0 and 1.'];
                return
            end
            if strcmp(timeType,'unit') && fieldValue>1
                errorMessage = ['OPTIONS.' fieldName ' must be less than 1 if TIMETYPE when ''unit''.'];
                return
            end
    end
end
 
function condition = isnonnegativescalar(x)
% Function that returns logical true if that input is a non-empty scalar >=0
condition = ~isempty(x) && isscalar(x) && x>=0;