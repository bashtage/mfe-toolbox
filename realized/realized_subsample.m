function [subsampledPrice,subsampledTimes] = realized_subsample(price,time,timeType,samplingType,samplingInterval, subSamples)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=6
    error('Six inputs required.')
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

if isempty(subSamples) || ~isscalar(subSamples) || subSamples<1 || floor(subSamples)~=subSamples
    error('SUBSAMPLES must be a positive scalar.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subsampledPrice = cell(1,length(subSamples+1));
subsampledTimes = cell(1,length(subSamples+1));

% Convert times to unit if not already unit
if ~strcmp(timeType,'unit')
    [time,time0,time1,samplingInterval]=realized_convert2unit(time,timeType,samplingType,samplingInterval);
end
% Since everything is converted to unit, set timeType to unit
originalTimeType = timeType;
timeType='unit';
% m is the length of price
m=size(price,1);



[filteredPrice,filteredTimes] = realized_price_filter(price,time,timeType,samplingType,samplingInterval);
fullFilteredTimes = filteredTimes;
subsampledPrice{1} = filteredPrice;
subsampledTimes{1} = filteredTimes;

% Get the number of filtered times
n = size(filteredTimes,1);

if subSamples==1
    return
end

switch lower(samplingType)
    case {'calendartime','calendaruniform','fixed'}
        % Uniform sampling in CT so the gap is constant
        gap=diff(filteredTimes);
        % The stepsize is constant
        step = gap /(subSamples);
        for i=2:subSamples
            thisSampleTime = fullFilteredTimes;
            % Need to compute times for this subsample
            thisSampleTime = thisSampleTime(1:n-1)+(i-1)*step;
            % Use these times to filter the price
            [filteredPrice, filteredTimes] = realized_price_filter(price,time,'unit','fixed',thisSampleTime);
            subsampledPrice{i} = filteredPrice;
            subsampledTimes{i} = filteredTimes;
        end
    case {'businesstime', 'businessuniform'}
        % Uniform sampling in business time
        if strcmp(samplingType,'businesstime')
            % If BT then samplingInterval contains the number of ticks to skip
            gap=samplingInterval;
        else
            % If businessuniform then gap will be average
            gap = m/samplingInterval;
        end
        % The step size is constant
        step = floor(gap/subSamples);
        for i=2:subSamples
            % The indices indicate which points to use since the
            % sampling is linear in ticks
            indices = floor(step*(i-1):gap:m);
            indices = indices(indices>0 & indices<=m);
            filteredPrice = price(indices);
            filteredTimes = time(indices);
            subsampledPrice{i} = filteredPrice;
            subsampledTimes{i} = filteredTimes;
        end
end


timeType = originalTimeType;
for i=1:subSamples
    if strcmp(timeType,'wall')
        subsampledTimes{i} = unit2wall(subsampledTimes{i},time0,time1);
    elseif strcmp(timeType,'seconds')
        subsampledTimes{i} = unit2seconds(subsampledTimes{i},time0,time1);
    end
    if ~strcmp(timeType,'unit')
        digits = -log10(max(subsampledTimes{i}) * eps) - 1;
        digits = max(digits,0);
        subsampledTimes{i} = floor(subsampledTimes{i} * 10^digits)/(10^digits);
    end
end