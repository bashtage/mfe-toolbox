function [time,time0,time1,samplingInterval]=realized_convert2unit(time,timeType,samplingType,samplingInterval)
% Helper function that converts wall and seconds times to unit times, assigning 
% 0 to the smallest time and 1 to the largest.
%
% USAGE:
%   [TIME,TIME0,TIME1,SAMPLINGINTERVAL]=realized_convert2unit(TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%
% INPUTS:
%   TIME             - m by 1 vector of times where TIME(i) corresponds to
%                       PRICE(i).  Must be either wall time or measurements in seconds
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on
%                       the first day.         
%   SAMPLINGTYPE     - String describing the type of sampling to use when
%                        filtering PRICE
%                        'CalendarTime' - Sample in calendar time using
%                          observations separated by SAMPLINGINTERVAL
%                          seconds. 
%                        'CalendarUniform' - Sample in calendar time using
%                          SAMPLINGINTERVAL observations spread uniformly
%                          between TIME(1) and TIME(m)
%                        'BusinessTime' - Sample in business (tick) time
%                          using observation separated by SAMPLINGINTERVAL
%                          ticks
%                        'BusinessUniform' - Sample in business (tick)
%                          time using observations uniformly spaced in
%                          business time.
%                        'Fixed' - Sample at specific points in time. When
%                          using fixed, SAMPLINGINTERVAL must be a n by 1 vector
%                          of times with the same TIMETYPE as TIME (i.e.
%                          seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on the
%                         selected SAMPLINGTYPE
%
% OUTPUTS:
%   TIME             - The original times mapped to 0-1
%   TIME0            - The base time in original format. Required to invert the transformation.
%   TIME0            - The final time in original format. Required to invert the transformation.
%   SAMPLINGINTERVAL - The sampling interval converted to the unit scale.
%                       If SAMPLINGYPE is not 'CalendarTime' or 'Fixed', SAMPLINGINTERVAL
%                       is not altered 
%
% COMMENTS:
%
%  See also REALIZED_KERNEL, REALIZED_PRICE_FILTER, REALIZED_VARIANCE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=4
    error('Four inputs required.')
end
if size(time,2)>size(time,1)
    time=time';
end
if any(diff(time)<0)
    error('TIME must be sorted and increasing')
elseif any(diff(time)==0)
    warning('oxfordRealized:realizedPriceFilter','TIME contains multiple entries with the same value. This creates an ambiguity and FILTEREDPRICE will contain the last price if TIME does not only unique elements.')
end
if size(time,2)>1 || length(time)~=length(time)
    error('TIME must be a m by 1 vector.')
end
% Inserted to protect against inputing integer times
    time = double(time);

timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds'})
    error('TIMETYPE must be one of ''wall'' or ''seconds''.');
end
samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    error('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.');
end
 
m=size(time,1);
t0Original=time(1);
tTOriginal=time(m);
if ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform'})
    if (~isscalar(samplingInterval) || floor(samplingInterval)~=samplingInterval || samplingInterval<1) 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

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
if strcmp(samplingType , 'calendartime')
    if strcmp(timeType,'wall')
        timeBase = wall2seconds([time0 time1]);
        samplingInterval = samplingInterval/diff(timeBase);
    elseif strcmp(timeType,'seconds')
        samplingInterval = samplingInterval/(time1-time0);
    end
end

