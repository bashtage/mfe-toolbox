function [medRV,minRV,medRVSS,minRVSS]=realized_min_med_variance(price,time,timeType,samplingType,samplingInterval,subsamples)
% Estimates truncated realized variance using the median and minimum method of Andersen, Dobrev and Schaumburg
%
% USAGE:
%   [MEDRV,MINRV] = realized_medrv(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SAMPLESPERBIN)
%   [MEDRV,MINRV,MEDRVSS,MINRVSS] = realized_medrv(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SAMPLESPERBIN,SUBSAMPLES)
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
%                         the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES should 
%                         be greater than 1, since 1 uses a single subsample and so is equivalent 
%                         to the original RV will compute a subsample realized variance using the 
%                         mid-point of the price sample points, 2 will use the original sample and 
%                         the midpoint between the original sampling times, and so on.  
%                         In general this number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   MEDRV             - Median RV estimate
%   MINRV             - Minimum RV estimate
%   MEDRVSS           - Subsampled version of median RV estimate
%   MINRVSS           - Subsampled version of minimum RV estimate
%
% COMMENTS:
%
% EXAMPLE:
%  % Using all 1-minute returns between 9:30:00 and 16:00:00
%  fixedInterval = seconds2wall(wall2seconds(93000):60:wall2seconds(160000));
%  [MEDRV,MINRV] = realized_medrv(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % Using all 5-minute returns between 9:30:00 and 16:00:00, but subsampling every minute
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [MEDRV,MINRV,MEDRVSS,MINRVSS] = realized_medrv(PRICE,TIME,'wall','Fixed',fixedInterval,5)
%
%  See also REALIZED_BIPOWER_VARIATION, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, 
%  REALIZED_RANGE, REALIZED_PRICE_FILTER

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 10/24/2011

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


logPrice = log(price);
% Convert times to unit if not already unit
if ~strcmp(timeType,'unit')
    [time,~,~,samplingInterval]=realized_convert2unit(time,timeType,samplingType,samplingInterval);
end
% Since everything is converted to unit, set timeType to unit
timeType='unit';
% m is the length of price

% Filter prices and compute the RV
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
% Construct the aggregated return series.  Remember will need to skip
returns = diff(filteredLogPrice);
m = length(returns);

returns2 = (sqrt(m)*returns).^2;
r2 = [returns2(1:m-2) returns2(2:m-1) returns2(3:m)];
medRV = mean(median(r2,2));
medRVScale = pi/(6-4*sqrt(3)+pi);
medRV = medRV*medRVScale ;

r2 = [returns2(1:m-1) returns2(2:m)];
minRV = mean(min(r2,[],2));
minRVscale = pi/(pi-2);
minRV = minRV*minRVscale;



subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
medRVs = nan(m*subsamples,1);
minRVs = nan(m*subsamples,1);
medRVcount = 0;
minRVcount = 0;
for i=1:subsamples
    filteredLogPrice = subsampledLogPrices{i};
    returns = diff(filteredLogPrice);
    m = size(returns,1);
    returns2 = m*returns.^2;
    % MedRV
    r2 = [returns2(1:m-2) returns2(2:m-1) returns2(3:m)];
    n = size(r2,1);
    medRVs(medRVcount+(1:n)) = median(r2,2);
    medRVcount = medRVcount  + n;
    % MinRV    
    r2 = [returns2(1:m-1) returns2(2:m)];
    n = size(r2,1);
    minRVs(minRVcount+(1:n)) = min(r2,[],2);
    minRVcount = minRVcount + n;
end
medRVSS = mean(medRVs(1:medRVcount));
medRVSS = medRVSS * medRVScale;
minRVSS = mean(minRVs(1:minRVcount));
minRVSS = minRVSS * minRVscale;