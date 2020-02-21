function [filteredPrice,filteredTime,actualTime] = realized_price_filter(price,time,timeType,samplingType,samplingInterval)
% Price filtering for computing realized variances and other
% high-frequency price based returns
%
% USAGE:
%   [FILTEREDPRICE,FILTEREDTIME,ACTUALTIME] =  realized_price_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%
% INPUTS:
%   PRICE             - m by 1 column vector of prices
%   TIME              - m by 1 column vector of times corresponding to PRICE.
%                         Must be sorted and increasing.
%   TIMETYPE          - String describing the way times are measured
%                        'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                        'seconds' Time measures in seconds past midnight.
%                          TIME must satisfy 0<=TIME<86400
%                        'unit'  Unit normalized date format, e.g. .1, .234, .9
%                          Unit normalized times are more general than the
%                          other types and can be applied to data from more
%                          than one calendar day
%   SAMPLINGTYPE      - String describing the type of sampling to use when
%                         filtering PRICE
%                         'CalendarTime' - Sample in calendar time using
%                           observations separated by SAMPLINGINTERVAL
%                           seconds. If TIMETYPE is 'unit',
%                           SAMPLINGINTERVAL must be between 0 and 1 and
%                           represents the fraction of the sample to skip
%                           when sampling.
%                         'CalendarUniform' - Sample in calendar time using
%                           SAMPLINGINTERVAL observations spread uniformly
%                           between TIME(1) and TIME(m)
%                         'BusinessTime' - Sample in business (tick) time
%                           using observation separated by SAMPLINGINTERVAL
%                           ticks
%                         'BusinessUniform' - Sample in business (tick)
%                           time using observations uniformly spaced in
%                           business time.
%                         'Fixed' - Sample at specific points in time. When
%                           using fixed, SAMPLINGINTERVAL must be a n by 1 vector
%                           of times with the same TIMETYPE as TIME (i.e.
%                           seconds if TIME is in seconds)
%   SAMPLINGINTERVAL   - Scalar integer or n by 1 vector whose meaning depends on the
%                          selected SAMPLINGTYPE
%
% OUTPUTS:
%   FILTEREDPRICE      - n by 1 vector of filtered prices.  n depends on
%                          SAMPLINGTYPE as well as the data
%   FILTEREDTIME       - n by 1 vector of sampling times corresponding to
%                          the chosen sampling scheme
%   ACTUALTIME         - n by 1 vector. ACTUALTIME(i) provides the actual
%                          observation time that was used when sampling at FILTEREDTIME(i)
%
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL.  Last price
%   interpolation is always used.  If SAMPLINGTYPE is 'Fixed' then any
%   requested samples before the first observation will be back filled
%   with PRICE(1). An indicator for back filled values can be constructed
%   using ACTUALTIME>FILTEREDTIME
%
% EXAMPLES:
%   % 5-minute prices
%   FILTEREDPRICE =  realized_price_filter(PRICE,TIME,'wall','CalendarTime',300)
%
%   % 79 prices spread uniformly in calendar time
%   FILTEREDPRICE =  realized_price_filter(PRICE,TIME,'wall','CalendarUniform',79)
%
%   % 20-tick returns
%   FILTEREDPRICE = realized_price_filter(PRICE,TIME,'wall','CalendarUniform',79)
%
%   % 391 prices spread uniformly in business time
%   FILTEREDPRICE = realized_price_filter(PRICE,TIME,'wall','CalendarUniform',79)
%
%   % 5-minute prices coresponding to 9:30:00, 9:35:00, ...
%   sampleTimes =  seconds2wall(wall2seconds(93000):300:wall2seconds(160000))';
%   FILTEREDPRICE = realized_price_filter(PRICE,TIME,'wall','Fixed',sampleTimes)
%
%   % .05-calendar time sampling when using unit times
%   TIME = wall2unit(TIME,min(TIME),max(TIME));
%   FILTEREDPRICE = realized_price_filter(PRICE,TIME,'unit','CalendarTime',.05)
%
%
%  See also REALIZED_KERNEL, REALIZED_KERNEL_CORE, REALIZED_KERNEL_WEIGHTS
%  REALIZED_KERNEL_SELECT_LAG_LENGTH, REALIZED_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=5
    error('Five inputs required.')
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
elseif any(diff(time)==0)
    warning('oxfordRealized:realizedPriceFilter','TIME contains multiple entries with the same value. This creates an ambiguity and FILTEREDPRICE will contain the last price if TIME does not only unique elements.')
end
if size(time,2)>1 || length(time)~=length(price)
    error('TIME must be a m by 1 vector.')
end
timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.');
end
% Inserted to protect against inputing integer times
time = double(time);

samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    error('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.');
end

m=size(price,1);
t0Original=time(1);
tTOriginal=time(m);
if ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform'})
    % Must be a scalar integer, unless using unit times
    if (~isscalar(samplingInterval) || floor(samplingInterval)~=samplingInterval || samplingInterval<1) && ~strcmp(timeType,'unit')
        error('SAMPLINGINTERVAL must be a positive integer for the SAMPLINGTYPE selected.')
    end
else
    if size(samplingInterval,2)>size(samplingInterval,1)
        samplingInterval=samplingInterval';
    end
    if ~(any(samplingInterval>=t0Original) && any(samplingInterval<=tTOriginal))
        error('At least one sampling interval must be between min(TIME) and max(TIME) when using ''Fixed'' as SAMPLINGTYPE.')
    end
    if any(diff(samplingInterval)<=0)
        error('When using ''Fixed'' as SAMPLINGTYPE the vector of sampling times in SAMPLINGINTERVAL must be sorted and strictly increasing.')
    end
end

if strcmp(timeType,'unit') && strcmp(samplingType,'calendartime')
    % samplingInterval must be between 0 and 1
    if samplingInterval>1
        error('When TIMETYPE is ''unit'' and SAMPLINGTYPE is ''CalendarTime'', SAMPLINGINTERVAL must also be in ''unit'' terms, and so must be between 0 and 1')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the size of the price
m=size(price,1);

% Convert times
if strcmp(samplingType,'fixed')
    time0 = min(samplingInterval);
    time1 = max(samplingInterval);
else
    time0 = min(time);
    time1 = max(time);
end
if strcmp(timeType,'wall')
    time = wall2unit(time,time0,time1);
elseif strcmp(timeType,'seconds')
    time = seconds2unit(time,time0,time1);
end
% If fixed convert these times as well
if strcmp(samplingType,'fixed')
    if strcmp(timeType,'wall')
        samplingInterval = wall2unit(samplingInterval,time0,time1);
    elseif strcmp(timeType,'seconds')
        samplingInterval = seconds2unit(samplingInterval,time0,time1);
    end
end


% Get the first and last time
t0=time(1);
tT=time(m);
% Switch based on the cases
switch samplingType
    case 'calendartime'
        % If calendar time
        % Convert sampling interval to be a fraction of a period
        % Unit strides
        if strcmp(timeType,'wall')
            timeBase = wall2seconds([time0 time1]);
            samplingInterval = samplingInterval/diff(timeBase);
        elseif strcmp(timeType,'seconds')
            samplingInterval = samplingInterval/(time1-time0);
        end
        filteredTime=(t0:samplingInterval:tT)';
        % Append the final observation
        if ~ismember(tT,filteredTime)
            filteredTime=[filteredTime; tT];
        end
        % Compute the filtered prices
        [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime);
    case 'calendaruniform'
        % Sampling interval contains the number of samples
        % Uniform spacing
        filteredTime=linspace(t0,tT,samplingInterval)';
        % Compute the filtered prices
        [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime);
    case 'businesstime'
        indices=(1:samplingInterval:m)';
        if ~ismember(m,indices)
            indices=[indices;m];
        end
        filteredPrice=price(indices);
        filteredTime=time(indices);
        actualTime = filteredTime;
    case 'businessuniform'
        % Sampling interval contains the number of samples
        indices=floor(linspace(1,m,samplingInterval));
        filteredPrice=price(indices);
        filteredTime=time(indices);
        actualTime = filteredTime;
    case 'fixed'
        filteredTime = samplingInterval;
        % Compute the filtered prices
        [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime);
    otherwise
        error('Unrecogniced SAMPLINGTYPE.')
end

% Finally convert filteredTime back to the original time format
if strcmp(timeType,'wall')
    filteredTime = unit2wall(filteredTime,time0,time1);
    actualTime = unit2wall(actualTime,time0,time1);
elseif strcmp(timeType,'seconds')
    filteredTime = unit2seconds(filteredTime,time0,time1);
    actualTime = unit2seconds(actualTime,time0,time1);
end



function [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime)

% Get the size of the inputs
m=size(price,1);
n=size(filteredTime,1);
% This uses a 2-index algorithm so that it is fast for any realistic size.
timeIndex=1;
filteredTimeIndex=1;
% Pl hold the index values.  Initializing to 1 makes the back fill easy
pl=ones(n,1);
while timeIndex<=m && filteredTimeIndex<=n
    if time(timeIndex)<=filteredTime(filteredTimeIndex)
        pl(filteredTimeIndex)=timeIndex;
        timeIndex=timeIndex+1;
    elseif filteredTimeIndex<n
        % Increment filteredTimeIndex until filteredTime>=time(timeIndex)
        while filteredTimeIndex<n && filteredTime(filteredTimeIndex)<time(timeIndex)
            filteredTimeIndex = filteredTimeIndex + 1;
            % Last price interploation for these since there is no new
            % information
            pl(filteredTimeIndex) = pl(filteredTimeIndex-1);
        end
        % Only assign if you aren't on the last one
        if filteredTimeIndex<n
            pl(filteredTimeIndex)=timeIndex;
        end
    elseif time(timeIndex)>filteredTime(n)
        % No more to assign, so break
        break
    end
end

% Clean up any trailing ones
starter = find(diff(pl)<0) + 1;
if ~isempty(starter)
    pl(starter:length(pl)) = pl(starter-1);
end
% Use pl the index the filteredPrice and filterTimeActual
filteredPrice = price(pl);
actualTime = time(pl);

