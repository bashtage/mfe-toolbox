function [bv,bvSS,bvDebiased,bvSSDebiased]=realized_bipower_variation(price,time,timeType,samplingType,samplingInterval,skip,subsamples)
% Computes bipower variation (BPV), skip-k bipower variation and subsample verisons of these
%
% USAGE:
%   [BV,BVSS,BVD,BVSSD] = realized_bipower_variation(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SKIP,SUBSAMPLES)
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
%   SKIP              - [OPTIONAL] Non-negative integer indicating number of return to skip with computing BPV.
%                         The skip-k BPV is |r_m||r_{m-1-k}| so when k=0 the usual BPV estimator is
%                         produced. Defulat value is 0.
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample BPV estimators
%                         to average with the original BPV.  Subsample BPV are based on prices
%                         uniformly spaced between the times (Calendar sampling) or ticks (Business
%                         sampling).  SUBSAMPLES=1 will compute asubsample BPV using the mid-point
%                         of the price sample points, 2 will use 1/3 and 2/3, and so on. In general
%                         this number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   BV                - Bipower variation estimate
%   BVSS              - Bipower variation including subsample averaging
%   BVD               - Debiased version of bipower variation equal to BV times m/(m-SKIP-1) where m 
%                         is the number of returns used to compute BV
%   BVSSD             - Debiased version of subsampled bipower variation 
%
% COMMENTS:
%   Skip-k BPV will have a small downward bias of the order (N-k)/N where N is the number of returns
%   used in computing the BPV when there are no jumps in the price.  If k is large relative to N
%   this bias may be worth correcting.
%
% EXAMPLE:
%  % 5-minute bipower variation using wall time returns and fixed increments
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  BV = realized_bipower_variation(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % 10-tick bipower variation
%  BV = realized_bipower_variation(PRICE,TIME,'wall','BusinessTime',10)
%
%  % 5-minute skip-1 bipower varation
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [BV,BVSS] = realized_bipower_variation(PRICE,TIME,'wall','Fixed',fixedInterval,1)
%
%  % 5-minute realized variance with subsampling every minute
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  BV = realized_bipower_variation(PRICE,TIME,'wall','Fixed',fixedInterval,0,5)
%
%
%  See also REALIZED_VARIANCE, REALIZED_QUARTICITY, REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, 
%  REALIZED_RANGE, REALIZED_PRICE_FILTER, REALIZED_KERNEL_WEIGHTS

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5 || nargin>7
    error('Five to Seven inputs required.')
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

if nargin>=6 && ~isempty(skip)
    if ~isscalar(skip) || skip<0 || floor(skip)~=skip
        error('SKIP must be a non-negative scalar.')
    end
else
    skip = 0;
end

if nargin==7 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
else
    subsamples = 1;
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
m = length(returns);
if (m-1-skip)<1 || (2+skip)>m
    error('The value of SKIP is too large for the sampling method used. SKIP should be small relative to the number of prices computed using the chosen sampling method.')
elseif (skip/m)>.1
    warning('oxfordRealized:skipTooLarge',['The value chosen for skip is large relative to the number of prices returned for the chosen \n'...
        'SAMPLINGTYPE and SAMPLINGINTERVAL - (skip / #returns) > .1 indicating more than 10%%  of the returns \n'...
        'are being dropped.  Interpret the results with caution.']);
end

bv= (pi/2) * abs(returns(1:m-1-skip))'*abs(returns(2+skip:m));
biasScale = (m-1-skip) / m ;
bvDebiased = bv / biasScale;

subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
bvs = zeros(subsamples,1);
totalCount  = 0;
for i=1:subsamples
    returns = diff(subsampledLogPrices{i});
    m = length(returns);
    bvs(i) = (pi/2) * abs(returns(1:m-1-skip))'*abs(returns(2+skip:m));
    if i==1
        baseCount = m-1-skip;
    end
    totalCount = totalCount + m-1-skip;
end
bvSS = sum(bvs) * (baseCount / totalCount);
bvSSDebiased = bvSS / biasScale;
