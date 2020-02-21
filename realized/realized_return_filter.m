function [returns,interval] = realized_return_filter(price,time,timeType,samplingType,samplingInterval,subsamples)
% THESE COMMENTS ARE WRONG!!!
% Computes the Quantile Realized Variance of Christensen, Oomen and Podolskij (2008),
% the MinRV and MedRV estimator of Andersen, Dobrev and Shaumberg and the
% general symmetrized version suggested by Sheppard in a discussion of COP
%
% USAGE:
%   [RETURN,INTERVAL] =  realized_return_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SUBSAMPLES)
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
%                          by SAMPLINGINTERVALticks
%                        'BusinessUniform' - Sample in business (tick) time using observations
%                          uniformly spaced in business time.
%                        'Fixed' - Sample at specific points in time. When using fixed,
%                          SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                          as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL - Scalar integer or n by 1 vector whose meaning depends on the selected SAMPLINGTYPE
%   SUBSAMPLES       - [OPTIONAL] Integer value containing the number of subsamples to use when
%                        computing RRSS. It must be the case that SUBSAMPLE<SAMPLESPERBIN. Setting
%                        SUBSAMPLES>0 will produce an estimator based on jittering the starting
%                        point over floor(linspace(0,(SUBSAMPLES)/(SUBSAMPLES+1),SUBSAMPLES+1)*SAMPLESPERBIN)+1
%                        For example, if SUBSAMPLES = 3 and SAMPLESPERBIN = 20, then the starting
%                        points for the four RR estimators will be observations 1, 6, 11 and 16.  If
%                        omitted SUBSAMPLE = 0 and RRSS=RR.
%
% OUTPUTS:
%   RETURNS          - Quantile realized variance vector (k by 1) corresponding to the values in QUANTILES
%   INTERVAL         - Subsample based version of RQ.  RQSS = RQ if SUBSAMPLES = 0
%
% COMMENTS:
%
% EXAMPLES:
%
%

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InputChecking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5 || nargin>6
    error('Seven or eight inputs required.')
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
    % Must be a scalar integer
    if ~isscalar(samplingInterval) || floor(samplingInterval)~=samplingInterval || samplingInterval<1
        error('SAMPLINGINTERVAL must be a positive integer for the SAMPLINGTYPE selected.')
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

% SUBSAMPLES
if nargin==6 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative integer.')
    end
else
    subsamples=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert if wall to compute the interval
if strcmp(timeType,'wall')
    time=wall2seconds(time);
    timeType = 'seconds';
end

% Initialize the return cell arrays and the interval vector
returns=cell(subsamples+1,1);
interval=zeros(subsamples+1,1);
% 1. Filter the price
logPrice = log(price);
[filteredLogPrice,filteredTimes] = realized_price_filter(logPrice ,time,timeType,samplingType,samplingInterval);
returns{1} = diff(filteredLogPrice);
interval(1) = 1;


if subsamples>0
    % Initialize the place holder for rvSS and rvAdjustment
    % Get the number of filtered times
    n = size(filteredTimes,1);
    % Compute the length of time used in computing the realized variance to
    % bias adjust the subsample RVs
    baseDifference = filteredTimes(n) - filteredTimes(1);
    switch lower(samplingType)
        case {'calendartime','calendaruniform','fixed'}
            % Uniform sampling in CT so the gap is constant
            gap=diff(filteredTimes);
            % The stepsize is constant
            step = gap /(subsamples+1);
            if ismember(samplingType,{'calendartime','calendaruniform'})
                step = [step;mean(step)];
            else
                step = [step;step(length(step))];
            end
            for i=1:subsamples
                thisSampleTime = filteredTimes;
                % Need to compute times for this subsample
                thisSampleTime = thisSampleTime(1:n)+i*step;
                % Use these times to filter the price
                filteredLogPrice = realized_price_filter(logPrice,time,'unit','fixed',thisSampleTime);
                % Compute returns
                returns{i+1} = diff(filteredLogPrice);
                % Compute the ith subsample RV
                % The adjustment depends on the amount of time used in this
                % subsample relative to the full sample RV
                interval(i+1) = (thisSampleTime(n-1)-thisSampleTime(1))/baseDifference;
            end
        case {'businesstime', 'businessuniform'}
            % Uniform sampling in business time
            if strcmp(samplingType,'businesstime')
                % If BT then samplingInterval contains the number of ticks to skip
                originalIndices = 1:samplingInterval:m;
                if ~ismember(m,originalIndices)
                    originalIndices = [originalIndices m];
                end
                gap=samplingInterval;
            else
                % If businessuniform then gap will be average
                originalIndices = floor(linspace(1,m,samplingInterval));
                if ~ismember(m,originalIndices)
                    originalIndices = [originalIndices m];
                end
                gap = m/samplingInterval;
            end
            % The step size is constant
            step = gap/(subsamples+1);
            for i=1:subsamples
                % The indices indicate which points to use since the
                % sampling is linear in ticks
                indices = floor(originalIndices + step*i);
                indices(indices>m) = m;
                returns{i+1} = diff(logPrice(indices));
                interval(i+1) = (time(max(indices))-time(min(indices)))/baseDifference;
            end
    end
end