function [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,QTtype,skip,subsamples)
% Estimates skip-k realized quarticity using either tripower quarticity or quad power quarticity,
% and subsample verisons of these
%
% USAGE:
%   [QT,QTSS] = realized_quarticity(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,QTTYPE,SKIP,SUBSAMPLES)
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
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on the selected SAMPLINGTYPE
%   QTTYPE            - [OPTIONAL] String value, either 'BNS', 'Tripower' or 'Quadpower'. Default value is
%                         'Tripower'.
%   SKIP              - [OPTIONAL] Non-negative integer indicating number of return to skip with computing BPV.
%                         The skip-k TQ is |r_m|^(4/3)|r_{m-(k+1)}|^(4/3)|r_{m-2*(k+1)}|^(4/3) so
%                         when k=0 the usual TQ estimator is produced. If BNS is selected as QTTYPE
%                         then skip is ignored. Default value is 0.
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample TQ estimators
%                         to average with the original TQ.  Subsample TQ are based on prices
%                         uniformly spaced between the times (Calendar sampling) or ticks (Business
%                         sampling).  SUBSAMPLES=1 will compute a subsample TQ using the mid-point
%                         of the price sample points, 2 will use 1/3 and 2/3, and so on. In general
%                         this number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   QT                - Estimate of quarticity
%   QTSS              - Estimate of quarticity using subsample averaging
%   QTD               - Debiased estimate of quarticity equal to QT times m / (m- k * skip - k) where k=3 
%                         for 'Tripower' and 4 for 'Quadpower' (0 for 'BNS')
%   QTSSD             - Debiased estimate of quarticity using subsample averaging
%
% COMMENTS:
%   Barndorf-Nielsen and Shephard's estimator is not robust to jumps.
%
%   Skip-k quarticity will have a small downward bias of the order (N-k)/N where N is the number of returns
%   used in computing the QT.  If k is large relative to N this bias may be worth correcting.
%
% EXAMPLE:
%  % 5-minute quarticity estimation using tripower, wall time returns and fixed increments
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  QT = realized_quarticity(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % 5-minute quarticity estimation using quadpower, wall time returns and fixed increments
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  QT = realized_quarticity(PRICE,TIME,'wall','Fixed',fixedInterval,'Quadpower')
%
%  % 10-tick quarticity estimation
%  QT = realized_quarticity(PRICE,TIME,'wall','BusinessTime',10)
%
%  % 5-minute skip-1 tripower varation
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  QT = realized_quarticity(PRICE,TIME,'wall','Fixed',fixedInterval,[],1)
%
%  % 5-minute quarticity estimation with subsampling every minute
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [QT,QTSS] = realized_quarticity(PRICE,TIME,'wall','Fixed',fixedInterval,[],0,5)
%
%
%  See also REALIZED_VARIANCE, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE,
%  REALIZED_BIPOWER_VARIATION, REALIZED_SEMIVARIANCE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5 || nargin>8
    error('Five to Eight inputs required.')
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
% Inserted to protect against imputing integer times
timeCase = whos('time');
if ~strcmp(timeCase.class,'double')
    time = double(time);
end
 
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

if nargin>=6 && ~isempty(QTtype)
    QTtype = lower(QTtype);
    if ~ismember(QTtype,{'bns','tripower','quadpower'})
        error('QTTYPE must be ''BNS'', ''Tripower'' or ''Quadpower''.')
    end
else
    QTtype = 'tripower';
 
end
 
if nargin>=7 && ~isempty(skip)
    if ~isscalar(skip) || skip<0 || floor(skip)~=skip
        error('SKIP must be a non-negative scalar.')
    end
else
    skip = 0;
end
 
if nargin==8 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
else
    subsamples = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Filter prices and compute the quarticity
logPrice = log(price);
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
returns = diff(filteredLogPrice);
m = length(returns);
if (m-1-skip)<1 || (2+skip)>m
    error('The value of SKIP is too large for the sampling method used. SKIP should be small relative to the number of prices computed using the chosen sampling method.')
elseif (skip/m)>.1
    warning('oxfordRealized:skipTooLarge',['The value chosen for skip is large relative to the number of prices returned for the chosen \n'...
        'SAMPLINGTYPE and SAMPLINGINTERVAL - (skip / #returns) > .1 indicating more than 10%%  of the returns \n'...
        'are being dropped.  Interpret the results with caution.']);
end
 
[qt,count] = realized_quarticity_core(returns,QTtype,skip);
biasScale = count / length(returns);
qtDebiased = qt / biasScale;

subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
qts = zeros(subsamples,1);
totalCount = 0;
for i=1:subsamples
    returns = diff(subsampledLogPrices{i});
    [qts(i),count] = realized_quarticity_core(returns,QTtype,skip);
    if i==1
        baseCount = count;
    end
    totalCount = totalCount + count;
end
% Adjustment factor accounts for the removal of 1 return in the subsamples
qtSS = sum(qts) * baseCount / totalCount;
qtSSDebiased = qtSS / biasScale;
 
 
 
function [qt,count] = realized_quarticity_core(returns,QTtype,skip)
m = length(returns);
if ((m-2-2*skip)<1 && strcmp(QTType,'quadpower')) || ((m-3-3*skip)<1 && strcmp(QTType,'quadpower'))
    error('The value of SKIP is too large for the sampling method used. SKIP should be small relative to the number of prices computed using the chosen sampling method.')
end
switch QTtype
    case 'bns'
        constant = 3;
        qt = constant^(-1) * m * sum(returns.^4);
        count = length(returns);
    case 'tripower'
        constant = 2^(2/3) * gamma(7/6)/gamma(1/2);
        qt =  constant^(-3) * m * ...
            sum(abs(returns(1:m-2-2*skip)).^(4/3) ...
            .* abs(returns(2+skip:m-1-skip)).^(4/3) ...
            .* abs(returns(3+2*skip:m)).^(4/3));
        count = m-2-2*skip;
    case 'quadpower'
        constant = sqrt(2/ pi);
        qt =  constant^(-4) * m * ...
            sum(abs(returns(1:m-3-3*skip)) ...
            .* abs(returns(2+skip:m-2-2*skip)) ...
            .* abs(returns(3+2*skip:m-1-skip)) ...
            .* abs(returns(4+3*skip:m)));
        count = m-3-3*skip;
end