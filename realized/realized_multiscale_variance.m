function [rvms,rvmsSS,diagnostics] = realized_multiscale_variance(price,time,timeType,samplingType,samplingInterval,subsamples,options)
% Estimated quadratic variation using the Multi-Scale estimator of Zhang
%
% USAGE:
%   RVMS = realized_multiscale_variance(PRICE)
%   [RVMS,RVMSSS,DIAGNOSTICS] = realized_multiscale_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SUBSAMPLES,OPTIONS)
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
%                        filtering PRICE for estimating the fast scale estimator.
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
%   OPTIONS           - [OPTIONAL] Realized Two-Scale option structure initialized by calling
%                         realized_options('Multiscale'). See help realized_options for a description of
%                         available options.
%
% OUTPUTS:
%   RVMS        - Realized multi-scale variance estimate
%   RVMSSS      - Realized multi-scale variance estimate constructed by averaging RVMS across multiple
%                 initial observations.  If SAMPLINGTYPE is 'BusinessTime' and SAMPLINGINTERAL is 1
%                 then subsampling is not possible.
%   DIAGNOSTICS - Structure of useful diagnostic information.  Fields are
%                   BANDWIDTH             - Ratio of slow scale to fast scale where the fast scale
%                                           sampling is determined by SAMPLINGTYPE and SAMPLINGINTERVAL.
%                   NOISEVARIANCE         - Bandi-Russell noise variance estimate. Empty if user
%                                           supplies bandwidth.
%                   DEBIASEDNOISEVARIANCE - Bias adjusted estimate of the noise variance.
%                   IQESTIMATE            - Estimate of IQ used. Empty if not needed. Empty if user
%                                           supplies bandwidth.
%
% COMMENTS:
%   For best results:
%     - Use prices sampled close to the highest frequency available, if not using 'BusinessTime' and
%     - If using calendar time sampling, sample close to the limit of the data availability (e.g. 5-15
%       seconds for a stock that trades a few thousand times per day)
%     - When sampling less frequently (e.g. 5 - 30 minutes) standard realized variance is
%       probably more appropriate
%
% EXAMPLE:
%  % Default usage with 'BusinessTime' sampling with interval 1 and automatic bandwidth selection
%  RVMS = realized_multiscale_variance(PRICE)
%
%  % 5-tick multi-scale variance with automatic bandwidth selection
%  RVMS = realized_multiscale_variance(PRICE,TIME,'wall','BusinessTime',5)
%
%  % Two scale variance with a bandwidth of 30 sampling at every tick
%  options = realized_options('Multiscale');
%  options.bandwidth = 30;
%  RVMS = realized_multiscale_variance(PRICE,TIME,'wall','BusinessTime',1,0,options)
%
%  % 5-tick multi-scale variance with automatic bandwidth selection and dense subsampling
%  RVMS = realized_multiscale_variance(PRICE,TIME,'wall','BusinessTime',5,4)
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
        subsamples = 1;
        options = realized_options('multiscale');
    case 5
        subsamples = 1;
        options = realized_options('multiscale');
    case 6
        options = realized_options('multiscale');
    case 7
        % Nothing
    otherwise
        error('One, five, six or seven inputs required.')
end
samplingType=lower(samplingType);
timeType=lower(timeType);
[errorMessage, options] = realized_multiscale_variance_parameter_check(price,time,timeType,samplingType,samplingInterval,subsamples,options);

if ~isempty(errorMessage )
    error(errorMessage)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logPrice = log(price);
[filteredLogPrice,filteredTimes] = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
% FIXME: Is m length of returns or of price?
m = length(filteredLogPrice) - 1;


if isempty(options.bandwidth)
    % Estimate the optimal bandwidth if not provided
    [noiseVariance, debiasedNoiseVariance, IQEstimate] = realized_noise_estimate(price, time, timeType, options);
    diagnostics.noiseVariance=noiseVariance;
    diagnostics.debiasedNoiseVariance=debiasedNoiseVariance;
    diagnostics.IQEstimate=IQEstimate;
    % Select the correct noise variance estimate
    if options.useDebiasedNoise
        selectedNoiseVariance = debiasedNoiseVariance;
    else
        selectedNoiseVariance = noiseVariance;
    end
    % Compute the bandwidth
    m = length(filteredLogPrice) - 1;
    % FIXME Bandwidth is wrong
    % deriv = 2*(52/35)*IQEstimate*c.^4 - 48/5*c.^(2)*(selectedNoiseVariance*(sqrt(IQEstimate)+selectedNoiseVariance/2))-144*selectedNoiseVariance^2;
    % Never consider < 1/sqrt(m) so that BW is always >=1
    if IQEstimate>0 % makes no sense if IQ == 0
        LB = 2/sqrt(m);
        % NEver use more than sqrt m
        % Comes form Eq. 37 of Zhang assuming that Var(eps^2) = 2 omega^4
        UB = m^(1/2);
        c = UB;
        UBderiv = 2*(52/35)*IQEstimate*c.^4 - 48/5*c.^(2)*(selectedNoiseVariance*(sqrt(IQEstimate)+selectedNoiseVariance/2))-144*selectedNoiseVariance^2;
        for j=1:32
            MP = (UB+LB)/2;
            c = MP;
            MPderiv = 2*(52/35)*IQEstimate*c.^4 - 48/5*c.^(2)*(selectedNoiseVariance*(sqrt(IQEstimate)+selectedNoiseVariance/2))-144*selectedNoiseVariance^2;

            if sign(MPderiv) == sign(UBderiv)
                UB = MP;
                UBderiv = MPderiv;
            else
                LB = MP;
            end
        end
        c = (UB+LB)/2;
        options.bandwidth = c*sqrt(m);
    else
        warning('oxfordRealized:smallBandwidth','IQ Estimate is 0 so optimal bandwidth selection is not possible.')
        options.bandwidth = 5;
    end
end

bandwidth = ceil(options.bandwidth);
if bandwidth<2
    warning('oxfordRealized:smallBandwidth','The selected bandwidth is less then 2, and RVMS is not well defined in this case.  Setting bandwidth to 2 and proceeding.')
    bandwidth = 2;
end
diagnostics.bandwidth = bandwidth;



RVq = zeros(bandwidth,1);
for skip=1:bandwidth
    [RVq(skip),count,naturalCount] = overlap_realized_variance(filteredLogPrice,skip);
    RVq(skip) = RVq(skip) * (naturalCount / count);
end
% Weights come from Eq. 20
M = bandwidth;
weights = zeros(M,1);
for i=1:M
    weights(i) = 12 * (i / (M^2)) * (i/M - 1/2 - 1/(2*M)) / (1- (1 / M^2));
end
rvms = weights'*RVq;

% Subsampled version computed AC0 for multiple starting values, which then leads to different values
% for RVq and the statistic, which can then be averaged

% Subsampled RVq and RV1
subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
RVqs = zeros(subsamples,M);
count = zeros(subsamples,M);
naturalCount = zeros(subsamples,M);
for i=1:subsamples
    filteredLogPrice =  subsampledLogPrices{i};
    for skip=1:M
        [RVqs(i,skip),count(i,skip),naturalCount(i,skip)] = overlap_realized_variance(filteredLogPrice,skip);
    end
end
count = sum(count,1);
scale = naturalCount(1,:) ./ count;
RVqs = sum(RVqs,1);
RVq = RVqs .* scale;
RVq = RVq';
rvmsSS = weights'*RVq;

end


function [rv,count,naturalCount] = overlap_realized_variance(price,skip)

naturalCount = (length(price)-1)/skip;
m = length(price);
returns = price(2+skip:m) - price(1:m-skip-1);
count = length(returns);
rv = returns' * returns;

end



function [errorMessage,options] = realized_multiscale_variance_parameter_check(price,time,timeType,samplingType,samplingInterval,subsamples,options)
% Support function for realized_multiscale_variance that does input validation
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
flatTopKernelList = {'bartlett','multiscale','2ndorder','epanechnikov',...
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

if ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        errorMessage = 'SUBSAMPLES must be a non-negative scalar.';
    end
end
end

function condition = isnonnegativescalar(x)
% Function that returns logical true if that input is a non-empty scalar >=0
condition = ~isempty(x) && isscalar(x) && x>=0;
end