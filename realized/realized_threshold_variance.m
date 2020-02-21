function [rtv,rtvD,rtvSS,rtvSSD,diagnostics]=realized_threshold_variance(price,time,timeType,samplingType,samplingInterval,thresholdScale,localVar,subsamples)
% Estimates Realized Variance using thresholding to remove jumps and a subsampled version
%
% USAGE:
%   [RTV,RTVD] = realized_threshold_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%   [RTV,RTVD,RTVSS,RTVSSD] = realized_threshold_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,THRESHOLDSCALE,LOCALVAR,SUBSAMPLES)
%
% INPUTS:
%   PRICE           - m by 1 vector of high frequency prices
%   TIME            - m by 1 vector of times where TIME(i) corresponds to PRICE(i)
%   TIMETYPE        - String describing the way times are measured
%                      'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                      'seconds' Time measured in seconds past midnight on the first day.
%                      'unit'  Unit normalized date format, e.g. .1, .234, .9
%                        Unit normalized times are more general than the other types and can be
%                        applied to data from more than one calendar day
%   SAMPLINGTYPE    - String describing the type of sampling to use when
%                       filtering PRICE
%                       'CalendarTime' - Sample in calendar time using observations separated by
%                         SAMPLINGINTERVAL seconds
%                       'CalendarUniform' - Sample in calendar time using SAMPLINGINTERVAL
%                         observations spread uniformly between TIME(1) and TIME(m)
%                       'BusinessTime' - Sample in business (tick) time using observation separated
%                         by SAMPLINGINTERVAL ticks
%                       'BusinessUniform' - Sample in business (tick) time using observations
%                         uniformly spaced in business time.
%                       'Fixed' - Sample at specific points in time. When using fixed,
%                         SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                         as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL - Scalar integer or n by 1 vector whose meaning depends on SAMPLINGTYPE
%   THRESHOLDSCALE   - [OPTIONAL] Scale used when applying the threshold.  The interpretation of
%                        THRESHOLDSCALE depends on THRESHOLDTYPE.  If omitted, the default value is 3.
%   LOCALVAR         - [OPTIONAL] M by 1 vector of local variances for use in thresholding, where M is
%                        the number of returns after filtering according to SAMPLINGTYPE and SAMPLINGINTERVAL.
%                        The local variance is estimated using a local kernel weighted average if omitted.
%   SUBSAMPLES       - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                        variance estimators to average with the original realized variance.
%                        Subsample realized variances are based on prices uniformly spaced between
%                        the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES should
%                        be greater than 1, since 1 uses a single subsample and so is equivalent
%                        to the original RV will compute a subsample realized variance using the
%                        mid-point of the price sample points, 2 will use the original sample and
%                        the midpoint between the original sampling times, and so on.
%                        In general this number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   RTV                - Realized threshold variance estimate
%   RTVD               - Debiased version of realized threshold variance. The debiasing method
%                          depends on THRESHOLDTYPE. See comments.
%   RTVSS              - Realized threshold variance including subsample averaging
%   RTVSSD             - Debiased version of realized threshold variance including subsample averaging
%
% COMMENTS:
%   The debiasing method follows Corsi, Pirino,and Reno.
%   The number of returns is length(filteredPrice)-1 where
%   filteredPrice = realized_price_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%
% EXAMPLE:
%  % 5-minute realized threshold variance using wall time returns and fixed increments, and default  values
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  RV = realized_threshold_variance(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % 5-minute realized threshold variance using wall time returns and fixed increments, and default  values
%  % with subsampling
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [RTV,RTVD,RTVSS,RTVSSD] = realized_variance(PRICE,TIME,'wall','Fixed',fixedInterval,[],[],5)
%
%  See also REALIZED_VARIANCE, REALIZED_BIPOWER_VARIATION, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        thresholdScale = 3;
        localVar = [];
        subsamples = 1;
    case 6
        localVar = [];
        subsamples = 1;
    case 7
        subsamples = 1;
    case 8
        % Nothing
    otherwise
        error('5 to 8 inputs required.')
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


% THRESHOLDSCALE
if isempty(thresholdScale)
    thresholdScale=3;
end
if thresholdScale<=0 || ~isscalar(thresholdScale)
    error('THRESHOLDSCALE must be a positive scalar')
end

% LOCALVAR
if size(localVar,2)>size(localVar,1)
    localVar = localVar';
end
if ~isempty(localVar) && size(localVar,1)~=1
    error('LOCALVAR must be M by 1 where M is the number of returns after filtering.')
end

% SUBSAMPLES
if isempty(subsamples)
    subsamples = 1;
end
if nargin==8 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


logPrice =log(price);
% Filter prices and compute the RV
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
returns = diff(filteredLogPrice);
m = size(returns,1);
c = thresholdScale;
if isempty(localVar)
    L = 25;
    V = inf*ones(size(returns));
    K = -L:L;
    K = 1/sqrt(2*pi)*exp(-(K/L).^2/2);
    K(L+(-1:1)) = 0;
    returns2 = returns.^2;
    
    finished = false;
    while ~finished
        ind = returns2<(c^2.*V);
        Vold  = V;
        for i=1:m
            pl = i-L:i+L;
            valid = pl>1 & pl<=m;
            tempInd = ind(pl(valid));
            w = K(valid);
            V(i) = w*(returns2(pl(valid)).*tempInd)/(w*tempInd);
        end
        if all(V==Vold)
            finished = true;
        end
    end
else
    V = localVar;
    if length(V)~=m
        error('LOCALVAR must have the same number of returns as the filtered series.')
    end
end
expectedValue = 1./(2*normcdf(-c)*sqrt(pi))*(2/c^2) * gamma(3/2).*gammainc(c^2/2,3/2,'upper') * c^2 *V;

threshold = returns2>(c^2.*V);
mOrig = m;
rtv = returns(~threshold)'*returns(~threshold);
rtvD = returns(~threshold)'*returns(~threshold) + sum(expectedValue(threshold));

diagnostics.localVar = V;
diagnostics.returns = returns;


subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);

returns = NaN(size(returns,1),subsamples);
for i=1:subsamples
    temp = diff(subsampledLogPrices{i});
    returns(1:length(temp),i) = temp;
end
returns = returns';
returns = returns(~isnan(returns));
if size(returns,2)>size(returns,1)
    returns = returns';
end

if isempty(localVar)
    L = subsamples *  25;
    V = inf*ones(size(returns));
    K = -L:L;
    K = 1/sqrt(2*pi)*exp(-(K/L).^2/2);
    K(L+(-subsamples:subsamples)) = 0;
    m = size(returns,1);
    returns2 = returns.^2;
    
    finished = false;
    while ~finished
        ind = returns2<(c^2.*V);
        Vold  = V;
        for i=1:m
            pl = i-L:i+L;
            valid = pl>1 & pl<=m;
            tempInd = ind(pl(valid));
            w = K(valid);
            V(i) = w*(returns2(pl(valid)).*tempInd)/(w*tempInd);
        end
        if all(V==Vold)
            finished = true;
        end
    end
else
    V = localVar;
    V = repmat(V,1,subsamples);
    for i=2:size(V,2)
        w = (i-1)/subsamples;
        V(1:m-1,i) = (1-w)*V(1:m-1,1) + w*V(2:m,1);
    end
    V=V';
    V=V(:);
    V=V(1:length(returns));
end
diagnostics.localVarSS = V;
diagnostics.returnsSS = returns;
expectedValue = 1./(2*normcdf(-c)*sqrt(pi))*(2/c^2) * gamma(3/2).*gammainc(c^2/2,3/2,'upper') * c^2 *V;

threshold = returns2>(c^2.*V);
rtvSS = mOrig/m * returns(~threshold)'*returns(~threshold);
rtvSSD = mOrig/m * (returns(~threshold)'*returns(~threshold) + sum(expectedValue(threshold)));
