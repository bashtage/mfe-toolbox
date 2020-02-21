function [rk,rkAdjusted,diagnostics]=realized_kernel(price,time,timeType,samplingType,samplingInterval,options)
% Estimated quadratic variation using Realized Kernels
%
% USAGE:
%   [RK,RKADJUSTED,DIAGNOSTICS] = realized_kernel(PRICE)
%   [RK,RKADJUSTED,DIAGNOSTICS] = realized_kernel(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%   [RK,RKADJUSTED,DIAGNOSTICS] = realized_kernel(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,OPTIONS)
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
%   OPTIONS           - [OPTIONAL] Realized Kernel option structure initialized by calling
%                         realized_options('Kernel'). See help realized_options for a description of
%                         available options.
%
% OUTPUTS:
%   RK          - Realized kernel estimate
%   RKADJUSTED  - Realized kernel estimate adjusted to recognize that less than the entire sample
%                   was used to estimate the quadratic variation
%   DIAGNOSTICS - Structure of useful diagnostic information.  Fields are
%                   KERNEL                - The kernel used in the RK
%                   BANDWIDTH             - Number of lags used in the RK except in the case of an
%                                           infinite lag kernel.
%                   ADJUSTMENT            - De-biasing adjustment for the RK to recognize that less than
%                                           the entire sample was used to compute the 0-lag term
%                   FILTEREDPRICE         - Prices used in the RK
%                   FILTEREDTIME          - Time stamps of filtered prices used in the RK
%                   NOISEVARIANCE         - Bandi-Russell noise variance estimate. Empty if not needed to
%                                           compute the optimal bandwidth or number of points to jitter
%                   DEBIASEDNOISEVARIANCE - Bias adjusted estimate of the noise variance
%                   JITTERLAGS            - Present only if the end points are jittered, contains the
%                                           number of points that used to implement the pre-averaging.
%
% COMMENTS:
%   For best results:
%     - Use prices sampled at close to their highest frequency.  This improves the estimate of noise.
%     - Use business time sampling with at most a small number of trades between samples (preferably
%       1, but usually <15)
%     - If using calendar time sampling, sample close to the limit of the data availability (e.g. 5-15
%       seconds for a stock that trades a few thousand times per day)
%     - When sampling less frequently (e.g. 5 - 30 minutes) the standard realized variance is
%       probably more appropriate
%
%   This function is essentially a wrapper around six other functions. If interested in performance, functions
%   not used (especially realized_kernel_parameter_check) can be omitted
%     1. realized_kernel_parameter_check - Input validation.
%     2. realized_noise_estimate - Estimate the noise variance and a lower bound for IQ for use in
%        selecting the bandwidth or number of lags to use when jittering end points.
%     3. realized_kernel_bandwidth - Compute the optimal bandwidth using the outputs of
%        realized_noise_variance
%     4. realized_kernel_jitter_lag_length - Compute the optimal number of points to jitter
%     5. realized_price_filter - Filter the original prices to compute the returns used for
%     computing realized autocovariances
%     6. realized_kernel_weights computes the weights for a kernel and bandwidth
%     7. realized_kernel_core computes the estimate using the weights and returns.
%
% EXAMPLES:
%  % Default usage with 'nonflatparzen', 'BusinessTime' sampling with interval 1
%  RK = realized_kernel(PRICE)
%
%  % 10-tick realized kernel with 'nonflatparzen'
%  RK = realized_kernel(PRICE,TIME,'wall','BusinessTime',10)
%
%  % Realized Kernel using a Cubic kernel with other options at their default
%  options = realized_options('Kernel');
%  options.kernel = 'cubic';
%  RK = realized_kernel(PRICE,TIME,'wall','BusinessTime',1,options)
%
%  See also REALIZED_OPTIONS, REALIZED_NOISE_ESTIMATE, REALIZED_VARIANCE,
%  REALIZED_VARIANCE_OPTIMAL_SAMPLING, REALIZED_RANGE, REALIZED_QUANTILE_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

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
        options = realized_options('kernel');
    case 5
        options = realized_options('kernel');
    case 6
        % Nothing
    otherwise
        error('One, five or six inputs required.')
end

[errorMessage, options] = realized_kernel_parameter_check(price,time,timeType,samplingType,samplingInterval,options);

if ~isempty(errorMessage )
    error(errorMessage)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isJittered = strcmpi(options.endTreatment,'jitter');

% Inserted to protect against inputing integer times
time = double(time);

% 2. [Optional] If options.bandWidth is not provided, estimate the Noise Variance. Also needed to estimate the number
% of lags or the amount of jitter
if isempty(options.bandwidth) || (isJittered && isempty(options.jitterLags))
    [noiseVariance, debiasedNoiseVariance, IQEstimate] = realized_noise_estimate(price, time, timeType, options);
    diagnostics.noiseVariance=noiseVariance;
    diagnostics.debiasedNoiseVariance=debiasedNoiseVariance;
    % Select the correct noise variance estimate
    if options.useDebiasedNoise
        selectedNoiseVariance = debiasedNoiseVariance;
    else
        selectedNoiseVariance = noiseVariance;
    end
end


% 3. Compute bandwidth is needed
if isempty(options.bandwidth)
    % Compute the bandwidth
    options.bandwidth = realized_kernel_bandwidth(selectedNoiseVariance, IQEstimate, options);
    % Check that the bandwidth is not too large
end

% 4. If jittering is required, perform the optimal jittering prior to filtering prices.
if isJittered
    if isempty(options.jitterLags)
        jitterLags = realized_kernel_jitter_lag_length(selectedNoiseVariance, IQEstimate, options.kernel, options.filteredN);
        options.jitterLags = jitterLags;
    end
    % Jitter the price
    m = length(price);
    jl = options.jitterLags;
    p0 = mean(price(1:jl));
    p1 = mean(price(m-jl+1:m));
    t0 = time(ceil(mean(1:jl)));
    t1 = time(floor(mean(m-jl+1:m)));
    price = [p0;price(jl+1:m-jl);p1];
    time = [t0;time(jl+1:m-jl);t1];
end

% 5. Filter the price
[filteredPrice,filteredTime] = realized_price_filter(price,time,timeType,samplingType,samplingInterval);
returns = diff(log(filteredPrice));

% 6. Compute the weights
weights=realized_kernel_weights(options);

% 7. Compute the realized kernel
rk = realized_kernel_core(returns, weights, options);


% Clean up
n=size(filteredPrice,1);
m=size(price,1);

% Construct the diagnostics vector
diagnostics.kernel=options.kernel;
diagnostics.bandwidth=options.bandwidth;
diagnostics.filteredPrice=filteredPrice;
diagnostics.filteredTime=filteredTime;
diagnostics.weights = weights;
if isJittered
    diagnostics.jitterLags = options.jitterLags;
end

% Compute the debiasing adjustment
if isJittered
    diagnostics.adjustment = 1;
    % FIXME: Should be a debiasing adjustment here
else
    diagnostics.adjustment = (filteredTime(n-round(options.bandwidth))-filteredTime(1+round(options.bandwidth)))/(time(m)-time(1));
end

% Compute the bias adjusted RK
rkAdjusted = rk / diagnostics.adjustment;

















function [errorMessage, options] = realized_kernel_parameter_check(price,time,timeType,samplingType,samplingInterval,options)
% Support function for realized_kernel that does input validation
%
% USAGE:
%   [ERRORMESSAGE,OPTIONS] = realized_kernel_parameter_check(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINVERVAL,OPTIONS)
%
% INPUTS:
%   See realized_kernel
%
% OUTPUT:
%   ERRORMESSAGE - String containing a description of the error if one is detected.  Empty if no
%                    error.
%   OPTIONS - A realized_kernel options structure with the additional field
%               filteredN - Number of observations that are used to compute the kernel.  This is
%               needed for computing the optimal bandwidth.
%
% COMMENTS:
%   See realized_options for a description of the other fields in OPTIONS
%
%  See also REALIZED_KERNEL, REALIZED_OPTIONS
 
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
 
% initialized errorMessage
errorMessage = [];
 
% Price
if size(price,2)>size(price,1)
    price=price';
end
if size(price,2)>1
    errorMessage = 'PRICE must be a m by 1 vector.';
    return
end
 
% Time
if size(time,2)>size(time,1)
    time=time';
end
if any(diff(time)<0)
    errorMessage = 'TIME must be sorted and increasing';
    return
end
if size(time,2)>1 || length(time)~=length(price)
    errorMessage = 'TIME must be a m by 1 vector the same size as PRICE.';
    return
end
% TimeType
timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    errorMessage = 'TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.';
    return
end
% Sampling Type
samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    errorMessage = 'SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.';
    return
end
 
m=size(price,1);
t0=time(1);
tT=time(m);
if ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform'})
    % Must be a scalar integer
    if ~isscalar(samplingInterval) || floor(samplingInterval)~=samplingInterval || samplingInterval<1
        errorMessage = 'SAMPLINGINTERVAL must be a positive integer for the SAMPLINGTYPE selected.';
        return
    end
else
    if size(samplingInterval,2)>size(samplingInterval,1)
        samplingInterval=samplingInterval';
    end
    if ~(any(samplingInterval>=t0) && any(samplingInterval<=tT))
        errorMessage = 'At least one sampling interval must be between min(TIME) and max(TIME) when using ''Fixed'' as SAMPLINGTYPE.';
        return
    end
    if any(diff(samplingInterval)<=0)
        errorMessage = 'When using ''Fixed'' as SAMPLINGTYPE the vector of sampling times in SAMPLINGINTERVAL must be sorted and strictly increasing.';
        return
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
defaultOptions = realized_options('Kernel');
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
        case {'kernel','medFrequencyKernel'}
            % Member of kernelList
            if ~ismember(fieldValue,kernelList)
                errorMessage = ['OPTIONS.' fieldName ' must be one of the listed types.'];
                return
            end
        case {'bandwidth','medFrequencyBandwidth'}
            % Non-negative scalar
            if ~isempty(fieldValue) && ~isnonnegativescalar(fieldValue)
                errorMessage = ['OPTIONS.' fieldName ' must a non-negative scalar.'];
                return
            end
        case 'endTreatment'
            % String
            if ~ismember(fieldValue,{'jitter','stagger'})
                errorMessage = 'OPTIONS.endTreatment must be either ''Jitter'' or ''Stagger''.';
                return
            end
        case {'jitterLags','maxBandwidth'}
            % Non-negative scalar integer
            if ~isempty(fieldValue) && ~isnonnegativescalarinteger(fieldValue)
                errorMessage = ['OPTIONS.' fieldName ' must be a non-negative scalar integer.'];
                return
            end
        case {'useDebiasedNoise','useAdjustedNoiseCount'}
            % Logical or scalar
            if ~islogical(fieldValue) && ~ismember(fieldValue,[0 1])
                errorMessage = 'OPTIONS.useDebiasedNoise must be a logical value.';
                return
            end
        case 'maxBandwidthPerc'
            % Scalar 0<=fieldvalue<=1
            if ~isempty(fieldValue) && ~isnonnegativescalar(fieldValue) && fieldValue>1
                errorMessage = 'OPTIONS.maxBandwidthPerc must be a non-negative scalar between 0 and 1.';
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
 
% Compute the number of prices available in the filtered series
if ~strcmpi(samplingType,'businesstime') || samplingInterval~=1
    tempFilteredPrice = realized_price_filter(price,time,timeType,samplingType,samplingInterval);
    options.filteredN = length(tempFilteredPrice);
else
    options.filteredN = length(price);
end
 
 
function condition = isnonnegativescalarinteger(x)
% Function that returns logical true if that input is a non-empty, scalar integer >=0
condition = ~isempty(x) && isscalar(x) && x>=0 && x==floor(x);
 
 
function condition = isnonnegativescalar(x)
% Function that returns logical true if that input is a non-empty scalar >=0
condition = ~isempty(x) && isscalar(x) && x>=0;