function [rmk,diagnostics]=realized_multivariate_kernel(varargin)
% Computes the multivariate version of the realized kernel of BNHLS using the non-flat-top kernels
% and end-point jittering
%
% USAGE:
%   [RMK,DIAGNOSTICS] = realized_multivariate_kernel(PRICE01,TIME01,PRICE02,TIME02,...,PRICEN,TIMEN,..
%                                                    TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
%   [RMK,DIAGNOSTICS] = realized_multivariate_kernel(PRICES,TIMES,TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
% INPUTS:
%   PRICE01          - m01 by 1 vector of high frequency prices
%   TIME01           - m01 by 1 vector of times where TIME01(i) corresponds to PRICE01(i)
%   PRICE02          - m02 by 1 vector of high frequency prices
%   TIME02           - m02 by 1 vector of times where TIME02(i) corresponds to PRICE02(i)
%   [N price-time pairing as required, where N can be any integer]
%   PRICEN           - mN by 1 vector of high frequency prices
%   TIMEN            - mN by 1 vector of times where TIMEN(i) corresponds to PRICEN(i)
%     -or-
%   PRICES           - Cell array containing prices, e.g. PRICES{1}=PRICE01,
%                        PRICES{2}=PRICE02, etc.
%   TIMES            - Cell array containing observation times, e.g. TIMES{1}=TIME01,
%                        TIMES{2}=TIME02, etc.
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGINTERVAL - Scalar integer or n by 1 vector which indicates the sampling interval of the
%                        refresh-time synchronized prices.  This value is the multivariate extension
%                        of the SAMPLINGINTERVAL for 'BusinessTime' sampling in the univariate kernel.
%                        1 indicates to use every price.  2 indicates to use every other price, and
%                        so on. 1 is usually the correct choice.
%   OPTIONS          - [OPTIONAL] Realized Kernel option structure initialized by calling
%                        realized_options('Kernel'). See help realized_options for a description of
%                        available options.
%
% OUTPUTS:
%   RMK              - Realized covariance estimate
%   DIAGNOSTICS      - Structure of useful diagnostic information.  Fields are
%                     KERNEL                - The kernel used in the RMK
%                     BANDWIDTH             - Number of lags used in the RMK except in the case of an
%                                             infinite lag kernel.
%                     ADJUSTMENT            - De-biasing adjustment for the RK to recognize that less than
%                                             the entire sample was used to compute the 0-lag term
%                     RTPRICE               - Prices used in the RK
%                     RTTIME                - Time stamps of filtered prices used in the RK
%                     NOISEVARIANCE         - Bandi-Russell noise variance estimate. Empty if not needed to
%                                             compute the optimal bandwidth or number of points to jitter
%                     DEBIASEDNOISEVARIANCE - Bias adjusted estimate of the noise variance
%                     JITTERLAGS            - Present only if the end points are jittered, contains the
%                                             number of points that used to implement the pre-averaging.
%
% COMMENTS:
%   The only samplng scheme available is refresh time.  Hence the input SAMPLINGTYPE is not defined.
%   The only end point treatment is jittering (pre-averaging).  This is required to ensure the
%   multivariate realized kernels are positive semi-definite.
%
%   Calling realized_multivariate_kernel with a single price series:
%   [RMK] = realized_multivariate_kernel(PRICE1,TIME1,TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
%   
%   is the same as calling realized_kernel:
%   [RK] = realized_kernel(PRICE1,TIME1,TIMETYPE,'BusinessTime',SAMPLINGINTERVAL,OPTIONS)
%
% EXAMPLES:
%
%  See also REALIZED_COVARIANCE, REALIZED_HAYASHI_YOSHIDA, REALIZED_KERNEL, REALIZED_VARIANCE,
%  REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking and Parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    error('At least 4 inputs required')
end

% Handle case where cell array is used to hold prices and times
if iscell(varargin{1})
    k = length(varargin{1});
    m = length(varargin)-2;
    temp = cell(1,2*k+m);
    for i=1:k
        temp{2*i-1} = varargin{1}{i};
        temp{2*i} = varargin{2}{i};
    end
    for i=1:m
        temp{2*k+i} = varargin{2+i};
    end
    varargin = temp;
end
% Check inputs
[errorMessage, options, numPrices, timeType, samplingInterval] = realized_multivariate_kernel_parameter_check(varargin{:});

if ~isempty(errorMessage )
    error(errorMessage)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 0. Quick filter using refresh time to estimate the number of data points that will be available
rtPrice = realized_refresh_time(timeType,varargin{1:2*numPrices});
if samplingInterval~=1
    rtPrice=rtPrice(1:samplingInterval:size(rtPrice,1),:);
end
nMax = size(rtPrice(1:samplingInterval:size(rtPrice,1),:),1);

% 1. Determine optimal jitter and jitter each series [if not provided]
% and
% 2. Determine optimal number of lags for each series [if not provided]
options.filteredN = nMax;

if isempty(options.jitterLags) && isempty(options.bandwidth)
    % Need both
    bandwidth = zeros(numPrices ,1);
    jitterLags = zeros(numPrices ,1);
    noiseVariance = zeros(numPrices ,1);
    IQEstimate = zeros(numPrices ,1);
    for i=1:numPrices
        price = varargin{2*(i-1)+1};
        time = varargin{2*i};
        
        [noiseVariance(i), debiasedNoiseVariance, IQEstimate(i)] = realized_noise_estimate(price, time, timeType, options);
        bandwidth(i) = realized_kernel_bandwidth(noiseVariance(i), IQEstimate(i), options);

        if bandwidth>(.25*nMax)
            bandwidth = round(.25*nMax);
            warning('oxfordRealized:realizedKernelLength','The estimated bandwidth requires a lag length larger than 25%% of the available data.  Bandwidth has been truncated to 25%% of data.');
        end
        jitterLags(i) = realized_kernel_jitter_lag_length(noiseVariance(i), IQEstimate(i), kernel, nMax);
        jitterLags(i) = max(jitterLags(i),1);
    end
    options.bandwidth = round(mean(bandwidth));
elseif isempty(options.jitterLags) && ~isempty(options.bandwidth)
    % Need jitterLags
    jitterLags = zeros(numPrices ,1);
    noiseVariance = zeros(numPrices ,1);
    IQEstimate = zeros(numPrices ,1);
    for i=1:numPrices
        price = varargin{2*(i-1)+1};
        time = varargin{2*i};
        
        [noiseVariance(i), debiasedNoiseVariance, IQEstimate(i)] = realized_noise_estimate(price, time, timeType, options);

        jitterLags(i) = realized_kernel_jitter_lag_length(noiseVariance(i), IQEstimate(i), kernel, nMax);
        jitterLags(i) = max(jitterLags(i),1);
    end
elseif ~isempty(options.jitterLags) && isempty(options.bandwidth)
    bandwidth = zeros(numPrices ,1);
    jitterLags = zeros(numPrices ,1);
    noiseVariance = zeros(numPrices ,1);
    IQEstimate = zeros(numPrices ,1);
    for i=1:numPrices
        price = varargin{2*(i-1)+1};
        time = varargin{2*i};
        
        [noiseVariance(i), debiasedNoiseVariance, IQEstimate(i)] = realized_noise_estimate(price, time, timeType, options);
        bandwidth(i) = realized_kernel_bandwidth(noiseVariance(i), IQEstimate(i), options);

        if bandwidth>(.25*nMax)
            bandwidth = round(.25*nMax);
            warning('oxfordRealized:realizedKernelLength','The estimated bandwidth requires a lag length larger than 25%% of the available data.  Bandwidth has been truncated to 25%% of data.');
        end
    end
    options.bandwidth = round(mean(bandwidth));
    jitterLags = ones(numPrices,1)*options.jitterLags;
end

% Do the jittering
for i=1:numPrices
    price = varargin{2*(i-1)+1};
    time = varargin{2*(i-1)+2};
    m = size(price,1);
    p0 = mean(price(1:jitterLags(i)));
    p1 = mean(price(m-jitterLags(i)+1:m));
    t0 = time(ceil(mean(1:jitterLags(i))));
    t1 = time(floor(mean(m-jitterLags(i)+1:m)));
    price = [p0;price(jitterLags(i)+1:m-jitterLags(i));p1];
    time = [t0;time(jitterLags(i)+1:m-jitterLags(i));t1];
    varargin{2*(i-1)+1} = price;
    varargin{2*(i-1)+2} = time;
end
% 3. Put into refresh time
[rtPrice,rtTime] = realized_refresh_time(timeType,varargin{1:2*numPrices});

if samplingInterval~=1
    rtPrice=rtPrice(1:samplingInterval:size(rtPrice,1),:);
end
% 4. Compute the weights
weights = realized_kernel_weights(options);
% 5. Compute the kernel
rtReturns = diff(log(rtPrice));
rmk = realized_multivariate_kernel_core(rtReturns,weights);

diagnostics.kernel     = options.kernel;
diagnostics.bandwidth  = options.bandwidth;
diagnostics.jitterLags = jitterLags;
diagnostics.weights    = weights;
diagnostics.rtPrice = rtPrice;
diagnostics.rtTime = rtTime;





function rmk = realized_multivariate_kernel_core(returns, weights)
% Realized Multivariate Kernel core routine that computes the value of a multivariate realized kernel
% given a set of returns and weights.
%
% USAGE:
%   [RMK] = realized_kernel_core(RETURNS,WEIGHTS)
%
% INPUTS:
%   RETURNS   - m by 1 column vector of returns
%   WEIGHTS   - H by 1 column vector of kernel weights corresponding to
%                    lags 1, 2, ..., H.  H should be much smaller than m
%
% OUTPUTS:
%   RMK        - Realized multivariate kernel value
%
% COMMENTS:
%   This is a helper function for REALIZED_MULTIVARIATE_KERNEL
%
%  See also REALIZED_MULTIVARIATE_KERNEL, REALIZED_KERNEL, REALIZED_PRICE_FILTER,
%  REALIZED_KERNEL_WEIGHTS, REALIZED_KERNEL_SELECT_LAG_LENGTH, REALIZED_COVARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=2
    error('Two inputs required.')
end

if size(weights,2)>1
    error('WEIGHTS must be a H by 1 vector.')
end

% Get the number of returns
[m,k]=size(returns);
if length(weights)>=m
    warning('oxfordRealized:realizedKernelLength','The length of WEIGHTS is longer than the length of RETURNS.\n  The weights are being truncated at N-1 where N is the number of returns')
    weights = weights(1:m-1);
elseif length(weights)> .1*length(returns)
    warning('oxfordRealized:realizedKernelLength','The number of WEIGHTS may be excessive.  Consider using a smaller kernel.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the size of the kernel
H=length(weights);


% Use all returns when computing autocovariances
gammaH = zeros(k,k,H);
% Loop and compute the values for gamma
for i=1:H
    returnsMinus  = returns(1:m-i,:);
    returnsPlus   = returns(i+1:m,:);
    gammaH(:,:,i) = returnsMinus' * returnsPlus;
end
% Compute the 0-lag gamma using all returns
gamma0 = returns' * returns;

% Construct the kernel
rmk = gamma0;
if ~isempty(weights)
    for j=1:length(weights)
        rmk = rmk + weights(j)*(gammaH(:,:,j)+gammaH(:,:,j)');
    end
else
    warning('oxfordRealized:realizedKernelLength','The number of lags used was 0.');
end
























function [errorMessage, options, numPrices, timeType, samplingInterval] = realized_multivariate_kernel_parameter_check(varargin)
% Support function for realized_multivariate_kernel that does input validation
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
 
% Determine the number of prices and times
k = length(varargin);
if k/2~=floor(k/2)
    options = varargin{k};
    timeType = varargin{k-2};
    samplingInterval = varargin{k-1};
    numPrices = (k-3)/2;
else
    options = realized_options('Multivariate Kernel');
    timeType = varargin{k-1};
    samplingInterval = varargin{k};
    numPrices = (k-2)/2;
end





% Prices and times
for i=1:numPrices
    price = varargin{2*i-1};
    time = varargin{2*i};
    if size(price,2)>size(price,1)
        price=price';
    end
    if size(price,2)>1
        errorMessage = ['PRICE series ' num2str(i) ' must be a m(i) by 1 vector.'];
        return
    end
    
    % Time
    if size(time,2)>size(time,1)
        time=time';
    end
    if any(diff(time)<0)
        errorMessage = ['TIME for series ' num2str(i) ' must be sorted and increasing'];
        return
    end
    if size(time,2)>1 || length(time)~=length(price)
        errorMessage = ['TIME for series ' num2str(i) ' must be a m(i) by 1 vector the same size as PRICE series (i).'];
        return
    end
end

% TimeType
timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    errorMessage = 'TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.';
    return
end


if samplingInterval<1 || ~isscalar(samplingInterval)
        errorMessage = 'SAMPLINGINTERVAL must be a positive integer (usually 1).';
        return
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
        case 'kernel'
            % Member of kernelList
            if ~ismember(fieldValue,nonFlatTopKernelList)
                errorMessage = ['OPTIONS.' fieldName ' must be one of the listed types.'];
                return
            end
        case 'medFrequencyKernel'
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
            if ~ismember(fieldValue,{'jitter'})
                errorMessage = 'OPTIONS.endTreatment must be ''Jitter''.';
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

 
 
function condition = isnonnegativescalarinteger(x)
% Function that returns logical true if that input is a non-empty, scalar integer >=0
condition = ~isempty(x) && isscalar(x) && x>=0 && x==floor(x);
 
 
function condition = isnonnegativescalar(x)
% Function that returns logical true if that input is a non-empty scalar >=0
condition = ~isempty(x) && isscalar(x) && x>=0;
