function [rsvn,rsvp,rsvnSS,rsvpSS]=realized_semivariance(price,time,timeType,samplingType,samplingInterval,subsamples)
% Estimates Realized Semivariance and subsampled Realized Semivariance
%
% USAGE:
%   [RSVN,RSVP,RSVNSS,RSVPSS] = realized_semivariance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SUBSAMPLES)
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
%                        'CalendarTime' - Sample in calendar time using observations separated by SAMPLINGINTERVAL seconds
%                        'CalendarUniform' - Sample in calendar time using SAMPLINGINTERVAL
%                          observations spread uniformly between TIME(1) and TIME(m)
%                        'BusinessTime' - Sample in business (tick) time using observation separated
%                          by SAMPLINGINTERVAL ticks
%                        'BusinessUniform' - Sample in business (tick) time using observations
%                          uniformly spaced in business time.
%                        'Fixed' - Sample at specific points in time. When using fixed,
%                          SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                          as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on the selected SAMPLINGTYPE
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                         semivariance estimators to average with the original realized semivariance.
%                         Subsample realized semivariances are based on prices uniformly spaced between
%                         the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES=1
%                         will compute a subsample realized semivariance using the mid-point of the
%                         price sample points, 2 will use 1/3 and 2/3, and so on. In general this
%                         number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   RSVN              - Negative-part realized semivariance estimate
%   RSVP              - Positive-part realized semivariance estimate
%   RSVNSS            - Subsample averaged negative-part realized semivariance estimate
%   RSVPSS            - Subsample averaged positive-part realized semivariance estimate
%
% COMMENTS:
%
% EXAMPLES:
%  % 5-minute realized semivariance using wall time returns and fixed increments
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [RSVN,RSVP] = realized_semivariance(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % 10-tick realized semivariance
%  [RSVN,RSVP] = realized_semivariance(PRICE,TIME,'wall','BusinessTime',10)
%
%  % 5-minute realized semivariance with subsampling every minute
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [RSVN,RSVP,RSVNSS,RSVPSS] = realized_semivariance(PRICE,TIME,'wall','Fixed',fixedInterval,5)
%
%  See also REALIZED_VARIANCE, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE,
%  REALIZED_PRICE_FILTER

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5 || nargin>6
    error('Five or Six inputs required.')
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
% Inserted to protect against inputing integer times
time = double(time);

timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.');
end
samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    error('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.');
end

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

if nargin<6 || isempty(subsamples)
    subsamples = 1;
end    
if nargin==6 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert times to unit if not already unit
if ~strcmp(timeType,'unit')
    [time,~,~,samplingInterval]=realized_convert2unit(time,timeType,samplingType,samplingInterval);
end
% Since everything is converted to unit, set timeType to unit
timeType='unit';

% Filter prices and compute the RV
logPrice = log(price);
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
returns = diff(filteredLogPrice);
rsvn=sum((returns.^2).*(returns<0));
rsvp=sum((returns.^2).*(returns>0));

subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
rsvns = zeros(subsamples,1);
rsvps = zeros(subsamples,1);
totalCount = 0;
for i=1:subsamples
    returns = diff(subsampledLogPrices{i});
    rsvns(i) = sum((returns.^2).*(returns<0));
    rsvps(i) = sum((returns.^2).*(returns>0));
     if i==1
        baseCount = length(returns);
    end
    totalCount = totalCount + length(returns);
end
rsvnSS = sum(rsvns)  * (baseCount / totalCount);
rsvpSS = sum(rsvps)  * (baseCount / totalCount);
