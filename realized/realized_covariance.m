function [rc,rcSS]=realized_covariance(varargin)
% Estimates Realized Covariance and subsampled Realized Variance
%
% USAGE:
%   [RC,RCSS] = realized_covariance(PRICE01,TIME01,PRICE02,TIME02,...,PRICEN,TIMEN,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SUBSAMPLES)
%
% INPUTS:
%   PRICE01          - m01 by 1 vector of high frequency prices
%   TIME01           - m01 by 1 vector of times where TIME01(i) corresponds to PRICE01(i)
%   PRICE02          - m02 by 1 vector of high frequency prices
%   TIME02           - m02 by 1 vector of times where TIME02(i) corresponds to PRICE02(i)
%   [N - 3 price-time pairing as required, where N can be any integer]
%   PRICEN           - mN by 1 vector of high frequency prices
%   TIMEN            - mN by 1 vector of times where TIMEN(i) corresponds to PRICEN(i)
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGTYPE     - String describing the type of sampling to use when filtering PRICE01. Once
%                        PRICE01 is filtered, all other prices are filtered using the same
%                        time-stamps as the filtering proceedure produces for PRICE01.
%                        'CalendarTime' - Sample in calendar time using
%                          observations separated by SAMPLINGINTERVAL seconds
%                        'CalendarUniform' - Sample in calendar time using
%                          SAMPLINGINTERVAL observations spread uniformly
%                          between TIME(1) and TIME(m)
%                        'Fixed' - Sample at specific points in time. When
%                          using fixed, SAMPLINGINTERVAL must be a n by 1 vector
%                          of times with the same TIMETYPE as TIME (i.e.
%                          seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on the
%                         selected SAMPLINGTYPE
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                         variance estimators to average with the original realized variance.
%                         Subsample realized variances are based on prices uniformly spaced between
%                         the times.  SUBSAMPLES=1 will compute a subsample realized variance using
%                         the mid-point of the price sample points, 2 will use 1/3 and 2/3, and so
%                         on. In general this number should be small so the subsample estimators
%                         will be "sparse".
%
% OUTPUTS:
%   RC                - Realized covariance estimate
%   RKSS              - Realized covariance including subsample averaging.  If SUBSAMPLE = 0 or is
%                         omitted, RCSUBSAMPLE = RC
%
% COMMENTS:
%
% EXAMPLES:
%
%  See also REALIZED_KERNEL, REALIZED_HAYASHI_YOSHIDA, REALIZED_VARIANCE,
%  REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to parse teh varargin
Q = length(varargin);
if Q<5
    error('At least 5 minutes required.')
end
% Determin number of price, either Q-2 or Q-3 must be a time type
if ismember(varargin{Q-2},{'wall','seconds','unit'})
    % Then the number of inputs is (Q-3)/2
    k = (Q-3)/2;
elseif ismember(varargin{Q-3},{'wall','seconds','unit'}) % This is case where subsample is asked for
    k = (Q-4)/2;
else
    error('Inputs are not valid. The final arguments must be either TIMETYPE, SAMPLINGTYPE and SAMPLINGINTERVAL or TIMETYPE, SAMPLINGTYPE, SAMPLINGINTERVAL and SUBSAMPLES if using subsampling.')
end

% Parse price and time, make sure they are compatible with one another
price = cell(k,1);
time = cell(k,1);
% Gether the times and prices into cell arrays
for i=1:k
    price{i} = varargin{(i-1)*2+1};
    time{i} = double(varargin{i*2});
    if size(price{i},2)>size(price{i},1)
        price{i}=price{i}';
    end
    if size(time{i},2)>size(time{i},1)
        time{i}=time{i}';
    end
    if size(price{i},1)~=size(time{i},1)
        error(['Problem with PRICE' num2str(k) ' and TIME' num2str(k) '.  Their lengths are not compatible or they are not column vectors.'])
    end
    if any(diff(time{i})<0)
        error(['TIME' num2str(k) ' must be sorted and increasing.'])
    end
    time{i} = double(time{i});
end
timeType = varargin{2*k+1};
samplingType = varargin{2*k+2};
samplingInterval = varargin{2*k+3};

if Q==(2*k+4)
    subsamples = varargin{2*k+4};
else
    subsamples  = 1;
end

timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.');
end

samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','fixed'})
    error('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'' or ''Fixed''.');
end


if ismember(samplingType,{'calendartime','calendaruniform'})
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
    for i = 1:k
        m=size(price{i},1);
        t0=time{i}(1);
        tT=time{i}(m);
        if ~(any(samplingInterval>=t0) && any(samplingInterval<=tT))
            error(['At least one sampling interval must be between min(TIME' num2str(k) ' ) and max(TIME' num2str(k) ') when using ''Fixed'' as SAMPLINGTYPE.'])
        end
    end
    if any(diff(samplingInterval)<=0)
        error('When using ''Fixed'' as SAMPLINGTYPE the vector of sampling times in SAMPLINGINTERVAL must be sorted and strictly increasing.')
    end
end

if ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Convert times to unit if not already unit
if ~strcmp(timeType,'unit')
    [time{1},time0,time1,samplingInterval]=realized_convert2unit(time{1},timeType,samplingType,samplingInterval);
    for j=2:k
        % This is custom function to map to the unit interval according to the first price
        time{j}=realized_multivariate_convert2unit(time{j},timeType,time0,time1);
    end
end

% Filter prices and compute the RC
[filteredPrice(:,1),filteredTimes] = realized_price_filter(price{1},time{1},'unit',samplingType,samplingInterval);
% filteredTimes turns a problem with uniform times into a problem with
% fixed times
n  = size(filteredPrice,1);
filteredPrice = [filteredPrice zeros(n,k-1)];
for j=2:k
    filteredPrice(:,j) = realized_price_filter(price{j},time{j},'unit','fixed',filteredTimes);
end

returns = diff(log(filteredPrice));
rc=returns' * returns;
adjFactor = length(returns);


subsampledLogPrices = cell(k,1);
for j=1:k
    logPrice = log(price{j});
    subsampledLogPrices{j} = realized_subsample(logPrice,time{j},timeType,'fixed',filteredTimes,subsamples);
end

rcs = zeros(k,k,subsamples);
totalCount = 0;
for i=1:subsamples
    m = size(subsampledLogPrices{1}{i},1);
    allLogPrices = zeros(m,k);
    for j=1:k
        allLogPrices(:,j) = subsampledLogPrices{j}{i};
    end
    
    returns = diff(allLogPrices);
    if i==1
        baseCount = size(returns,1);
    end
    totalCount = totalCount + size(returns,1);
    rcs(:,:,i) = (returns' * returns);
end
rcSS = sum(rcs,3) * (baseCount / totalCount);



function time=realized_multivariate_convert2unit(time,timeType,time0,time1)
% Quick conversion function
if strcmp(timeType,'wall')
    time = wall2unit(time,time0,time1);
elseif strcmp(timeType,'seconds')
    time = seconds2unit(time,time0,time1);
end

