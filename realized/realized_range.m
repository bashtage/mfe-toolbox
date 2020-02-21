function [rr,rrSS]=realized_range(price,time,timeType,samplingType,samplingInterval,samplesperbin,overlap,subsamples)
% Computes Realized Range estimate of quadratic variation for a set of prices
%
% USAGE:
%   [RR,RRSS] = realized_range(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SAMPLESPERBIN)
%   [RR,RRSS] = realized_range(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SAMPLESPERBIN,OVERLAP,SUBSAMPLE)
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
%                         'Fixed' - Sample at specific points in time. When using fixed,
%                           SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                           as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on the selected SAMPLINGTYPE
%   SAMPLESPERBIN     - Number of samples to use in each window when computing the local high and
%                         low prices. NOTE:  It must be the case that (n-1)/(SAMPLESPERBIN-1) is an
%                         integer where n is the number of prices returned from calling
%                         realized_price_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%   OVERLAP           - [OPTIONAL] Boolean indicating whether to use overlapping blocks.  Default is TRUE.
%   SUBSAMPLES        - [OPTIONAL] Integer value containing the number of subsamples to use when
%                         computing RRSS. It must be the case that SUBSAMPLE<SAMPLESPERBIN. Setting
%                         SUBSAMPLES>0 will produce an estimator based on jittering the starting
%                         point over
%                         floor(linspace(0,(SUBSAMPLES)/(SUBSAMPLES+1),SUBSAMPLES+1)*SAMPLESPERBIN)+
%                         1 For example, if SUBSAMPLES = 3 and SAMPLESPERBIN = 20, then the starting
%                         points for the four RR estimators will be observations 1, 6, 11 and 16.
%                         If omitted SUBSAMPLE = 0 and RRSS=RR.
% OUTPUTS:
%   RR          - Realized Range estimate
%   RRSS        - Subsampled Realized Range estimate. If SUBSAMPLES=0, RRSS=RR
%
% COMMENTS:
%   Required scales computed from Monte Carlo integration and 1,000,000 simulations.  These values
%   were then filtered through a concave regression to enforce the concavity of the actual values.
%
% EXAMPLES:
%   % Estimate the RR from prices available from 9:30 to 16:00 using 1 minute returns and 31 samples
%   % per bin.  The corresponds to taking high and low values over 30 minute windows/
%   % NOTE: (391 - 1)/(30 - 1) is an integer.
%   sampleTimes = seconds2wall(wall2seconds(93000):60:wall2seconds(160000))
%   RR=realized_range(price,time,'wall','Fixed',sampleTimes,31);
%
%   % The same but using 4 subsamples
%   [RR,RRSS]=realized_range(price,time,'wall','Fixed',sampleTimes,31,4);
%
%   % Estimate the RR from prices previously filtered using 21 samples per bin
%   % 'BusinessTime' and 1 will pass the data in PRICE through realized_price_filter without altering it
%   time = linspace(0,1,length(price));
%   [RR,RRSS]=realized_range(price,time,'unit','BusinessTime',1,21)
%
%  See also REALIZED_KERNEL, REALIZED_VARIANCE, REALIZED_QUANTILE_VARIANCE, REALIZED_PRICE_FILTER

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InputChecking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<6 || nargin>8
    error('Six tp eight inputs required.')
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
% SAMPLESPERBIN >=2, integer
if ~isscalar(samplesperbin) || samplesperbin<1 || floor(samplesperbin)~=samplesperbin
    error('SAMPLESPERBIN must be a positive integer >= 2.')
end

% Overlap
if nargin>=7
    if isempty(overlap)
        overlap = true;
    end
    if ~isscalar(overlap)
        error('OVERLAP must be a scalar boolean value.')
    end
else
    overlap = true;
end
    
% SUBSAMPLES <SAMPLESPERBIN, integer
if nargin==8 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples || subsamples>=samplesperbin
        error('SUBSAMPLES must be a non-negative integer less than SAMPLESPERBIN.')
    end
end
if nargin<8 || isempty(subsamples)
    subsamples = 1;
end

if subsamples>1 && ~overlap
    warning('oxfordRealized:invalidArgument','SUBSAMPLE can only be used if OVERLAP = true')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Log the price
logPrice = log(price);
% 1. Filter the price
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
% Check the the number of prices are compatible with SAMPLESPERBIN
n = length(filteredLogPrice);
if ~overlap
    if floor((n-1)/(samplesperbin-1))~=((n-1)/(samplesperbin-1))
        error('The number of prices returned from realized_price_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL) minus 1 must be an integer multiple of SAMPLESPERBIN-1.');
    end
end

% Compute the number of bins and initialize hi and lo
expectedBins = (n-1)/(samplesperbin-1);
if overlap
    bins = n - samplesperbin;
    binStart = 1:n-samplesperbin;
else
    bins = expectedBins;
    binStart = 1:(samplesperbin-1):(n-1);
end
overlapScale = bins / expectedBins;
hi = zeros(bins,1);
lo = zeros(bins,1);
% Get the scale by interpolation from pre-computed values
scale = realized_range_scale(samplesperbin);
% Loop over bins
count = 1;
for i=binStart
    pl = i:i+samplesperbin-1;
    hi(count) = max(filteredLogPrice(pl));
    lo(count) = min(filteredLogPrice(pl));
    count = count+1;
end
% Compute the range
range = hi-lo;
% Sum and scale for RR
rr = sum(range.^2/scale);
rr = rr / overlapScale;

if ~overlap
    rrSS = rr;
    return;
end

% Only subsampling if overlap
subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
rrs = zeros(subsamples,1);
totalCount  = 0;
for i=1:subsamples
    filteredLogPrice = subsampledLogPrices{i};
    n = length(filteredLogPrice);
    bins = n - samplesperbin;
    binStart = 1:n-samplesperbin;
    count = 1;
    for j=binStart
        pl = j:j+samplesperbin-1;
        hi(count) = max(filteredLogPrice(pl));
        lo(count) = min(filteredLogPrice(pl));
        count = count+1;
    end
    % Compute the range
    range = hi-lo;
    % Sum and scale for RR
    rrs(i) = sum(range.^2/scale);

    if i==1
        baseCount = bins;
    end
    totalCount = totalCount + bins;
end
rrSS = sum(rrs) * baseCount / totalCount;
rrSS = rrSS / overlapScale;


function scale = realized_range_scale(samplesperbin)

% This function returns the values from realized_range_simulation
m = [2 3 4 5 6 7 9 10  11 13 14 16 19 21 25 26  27 31 37 40 41 46 51 53
    61 66 73 76 79 91 101 105  118 121 131 151 157 181 196 201 226 235 261 301 313 326 361 391
    451 469 521 586 601 651 781 901 937 976 1171 1301 1561 1801 1951 2341 2601 2926 3901 4681 5851 7801 11701 23401]';
m = m(:);

% Original non smoothed scale
%  scale =  [1 1.22738793950581 1.38129316222994 1.49586325323039 1.58382026938316 1.65411741489763 1.76274703793664 1.80567095995311
%      1.84343091968443 1.90653760080108 1.93322755041988 1.97972207275679 2.03597552868968 2.06756054767511 2.11894835858662 2.12983316154042
%      2.14087862754787 2.17817135967418 2.22351292496946 2.24214459170628 2.24788882641980 2.27460218702568 2.29708856820153 2.30549703101574
%      2.33455132668516 2.35008755352443 2.36910007903299 2.37641132661795 2.38334301198733 2.40797148527175 2.42506529311421 2.43121854869946
%      2.44940750892417 2.45313695996941 2.46482117342484 2.48449458835561 2.48968108618955 2.50802029234284 2.51779941215272 2.52073146871522
%      2.53422764214392 2.53865696918306 2.54995123666773 2.56449302119240 2.5682995321247  2.57232562423207 2.58166543320398 2.58895454926453
%      2.60087684266229 2.60410912035026 2.61225354306209 2.62107413464386 2.62284403677135 2.62850877861585 2.6404619875275  2.64916691784836
%      2.65147342238523 2.65381274756804 2.66371302258563 2.66905057072249 2.67772167526478 2.68397033526092 2.68732490655565 2.69437726049380
%      2.69819116664869 2.70225895240861 2.71123872218922 2.71634024135096 2.72194531107221 2.72832937211447 2.73595458781774 2.74588988225517]';

% Smoothed scale using concave regression
scale = [  1.000000 1.228112 1.382860 1.494147 1.581668 1.669189 1.756710 1.804702 1.851427
    1.898152 1.942636 1.987119 2.031603 2.070209 2.108815 2.134815 2.160815 2.186815
    2.212815 2.238815 2.257754 2.276692 2.295631 2.314570 2.333509 2.352447 2.366642
    2.380837 2.395032 2.409227 2.423421 2.435522 2.447622 2.459722 2.471822 2.483922
    2.496022 2.508122 2.517765 2.526830 2.535895 2.544959 2.554024 2.563089 2.570665
    2.578241 2.585817 2.593393 2.600969 2.607880 2.614791 2.621570 2.628348 2.635126
    2.641904 2.648682 2.654688 2.660694 2.666700 2.672706 2.678712 2.684718 2.690724
    2.696730 2.702736 2.708742 2.714748 2.720754 2.726760 2.732766 2.738772 2.744779]';

scale=scale(:);

% Interpolation is in Monte Carlo, otherwise interpolate.  If SAMPLESPERBIN is very large return
% asymptotic value
m = [m;1e6];
scale = [scale;4*log(2)];
if samplesperbin<23401
    scale = interp1(m,scale,samplesperbin);
elseif samplesperbin<100000
    interp1([23401 100000],[scale(length(scale)) 4*log(2)],samplesperbin)
else
    scale = 4*log(2);
end