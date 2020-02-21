function rchy = realized_hayashi_yoshida(priceA,timeA,priceB,timeB,timeType,samplingType,samplingInterval,overlap)
% Computed the Hayashi-Yoshida estimator of quadratic covariation, and allows for
% the empirical-performance motivated K-lead-and-lag version similar to Drost and Nijman (1997).
%
% USAGE:
%   [RCHY] = realized_hayashi_yoshida_(PRICEA,TIMEA,PRICEB,TIMEB,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,OVERLAP)
%
% INPUTS:
%   PRICEA            - mA by 1 vector of high frequency prices
%   TIMEA             - mA by 1 vector of times where TIMEB(i) corresponds to PRICEA(i)
%   PRICEA            - mB by 1 vector of high frequency prices
%   TIMEA             - mB by 1 vector of times where TIMEB(i) corresponds to PRICEB(i)
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
%   OVERLAP            - Number of ticks to overlap when computing the HY estimator.  The original HY
%                          estimator uses 0, which corresponds to the maximum likelihood estimator for
%                          price processes whose observation times are driven by independent Poisson
%                          processes.  Empirically this estimator performs poorly because prices are
%                          not a vector semi-martingale and using a larger number of lags can alleviate
%                          this problem.
%
% OUTPUTS:
%   RCHY             - The K-lead-and-lag Hayashi-Yoshida covariance estimator
%
% COMMENTS:
%   Filtering the price in calendar time allow the creation of a Hayashi-Yoshida corrected
%   calendar-time sampled realized covariance.  The value in SAMPLINGTYPE is applied to price A and
%   then the realized HY estimator is computed using the filtered price of A and the filtered times
%   of A, and all observations of B.  Price filtering is done using realized_price_filter
%
% EXAMPLES:
%   % Standard use with all prices with wall time prices
%   RCHY = realized_hayashi_yoshida(PRICEA,TIMEA,PRICEB,TIMEB,'wall','BusinessTime',1,0)
%
%   % 10 lead and lag RCHY with all prices with wall time prices
%   RCHY = realized_hayashi_yoshida(PRICEA,TIMEA,PRICEB,TIMEB,'wall','BusinessTime',1,10)
%
%   % 1-minute realized covariance with a HY correction
%   RCHY = realized_hayashi_yoshida(PRICEA,TIMEA,PRICEB,TIMEB,'wall','CalendarTime',60,0)
%
%  See also REALIZED_MULTIVARIATE_KERNEL, REALIZED_COVARIANCE, REALIZED_KERNEL, REALIZED_VARIANCE,
%  REALIZED_RANGE, REALIZED_QUANTILE_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<7 || nargin>8
    error('Seven or eight inputs required.')
end
if size(priceA,2)>size(priceA,1)
    priceA=priceA';
end
if size(priceA,2)>1
    error('PRICEA must be a m by 1 vector.')
end
if size(timeA,2)>size(timeA,1)
    timeA=timeA';
end
if any(diff(timeA)<0)
    error('TIMEA must be sorted and increasing')
elseif any(diff(timeA)==0)
    warning('oxfordRealized:realizedPriceFilter','TIMEA contains multiple entries with the same value. This creates an ambiguity and FILTEREDPRICE will contain the last price if TIMEA does not only unique elements.')
end
if size(timeA,2)>1 || length(timeA)~=length(priceA)
    error('TIMEA must be a m by 1 vector.')
end
% Inserted to protect against inputing integer times
timeA = double(timeA);
if size(priceB,2)>size(priceB,1)
    priceB=priceB';
end
if size(priceB,2)>1
    error('PRICEB must be a m by 1 vector.')
end
if size(timeB,2)>size(timeB,1)
    timeB=timeB';
end
if any(diff(timeB)<0)
    error('TIMEB must be sorted and increasing')
elseif any(diff(timeB)==0)
    warning('oxfordRealized:realizedPriceFilter','TIMEB contains multiple entries with the same value. This creates an ambiguity and FILTEREDPRICE will contain the last price if TIMEB does not only unique elements.')
end
if size(timeB,2)>1 || length(timeB)~=length(priceB)
    error('TIMEB must be a m by 1 vector.')
end
% Inserted to protect against inputing integer times
timeB = double(timeB);

% make sure the intersection is non-empty
if ~(min(timeB)<max(timeA) || min(timeA)<max(timeB))
    warning('oxfordRealized:realizedHayashiYoshida','The intersection of TIMESA and TIMESB is empty.  The RCHY will necessarily be 0.');
    intersectionIsEmpty = true;
else
    intersectionIsEmpty = false;
end


timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.');
end

samplingType=lower(samplingType);
if ~ismember(samplingType,{'calendartime','calendaruniform','businesstime','businessuniform','fixed'})
    error('SAMPLINGTYPE must be one of ''CalendarTime'', ''CalendarUniform'', ''BusinessTime'', ''BusinessUniform'' or ''Fixed''.');
end

m=size(priceA,1);
t0Original=timeA(1);
tTOriginal=timeA(m);
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

if nargin==7
    overlap = 0;
elseif overlap<0 || floor(overlap)~=overlap || max(size(overlap))>1
    error('OVERLAP must be a non-negative integer.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Price A is the base price, price B is the other price.  First sample priceA according to
% samplingType and samplingInterval, then use the actual times of these to estimate the HY
% respecting the value in overlap

if ~intersectionIsEmpty
    % First filter the price
    [filteredPriceA,filteredTimeA,actualTimeA] = realized_price_filter(priceA,timeA,timeType,samplingType,samplingInterval);

    % Then call the core routine
    rchy = realized_hayashi_yoshida_core(filteredPriceA,actualTimeA,priceB,timeB,overlap);
else
    rchy = 0;
end





function [rchy,times] = realized_hayashi_yoshida_core(priceA,timeA,priceB,timeB,overlap)
% Core routine for computing the Hayashi-Yoshida estimator of quadratic covariation, and allows for
% the empirical-performance motivated K-lead-and-lag version similar to Drost and Nijman (1997).
%
% USAGE:
%   [RCHY,TIMES] = realized_hayashi_yoshida_core(PRICEA,TIMEA,PRICEB,TIMEB,OVERLAP)
%
% INPUTS:
%   PRICEA           - mA by 1 vector of high frequency prices
%   TIMEA            - mA by 1 vector of times where TIMEB(i) corresponds to PRICEA(i)
%   PRICEA           - mB by 1 vector of high frequency prices
%   TIMEA            - mB by 1 vector of times where TIMEB(i) corresponds to PRICEB(i)
%   OVERLAP          - [OPTIONA] Number of ticks to overlap when computing the HY estimator.  The
%                        original HY estimator uses 0, which corresponds to the maximum likelihood
%                        estimator for price processes whose observation times are driven by
%                        independent Poisson processes.  Empirically this estimator performs poorly
%                        because prices are not a vector semi-martingale and using a larger number
%                        of lags can alleviate this problem. If omitted OVERLAP = 0.
%
% OUTPUTS:
%   RCHY             - The K-lead-and-lag Hayashi-Yoshida covariance estimator
%   TIMES            - A mA by 1 matrix of time stamps contains the times where the prices were
%                        sampled for computing the returns in the HY estimator.  This is mostly for
%                        diagnostic purposes.
%
% COMMENTS:
%   This is a helper function of realized_hayashi_yoshida and does no input checking.  In general
%   it should not be directly called.
%
%  See also REALIZED_MULTIVARIATE_KERNEL, REALIZED_COVARIANCE, REALIZED_KERNEL, REALIZED_VARIANCE,
%  REALIZED_RANGE, REALIZED_QUANTILE_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

priceA = log(priceA);
priceB = log(priceB);
timeA = double(timeA);
timeB = double(timeB);
% Price A is the base price, price B is the one which will move

% First find the first time that B is available before timeA(1).  If it is empty then find the first
% timeA(1) which is weakly after the timeB(1)
nA = size(priceA,1);
nB = size(priceB,1);

if timeA(1)<timeB(1)
    % A starts before B, so the easy solution is to project backward the price of the first B to the
    % time of the first A.  This will generate a 0 return but makes the algorithm easier.
    timeB=[timeA(1);timeB];
    priceB=[priceB(1);priceB];
    nB = size(priceB,1);
end

if timeA(nA)>timeB(nB)
    % A ends after B, so the easy solution is to project forward the last price B to the
    % time of the last A.  This will generate a 0 return but makes the algorithm easier.
    timeB=[timeB;timeA(nA)];
    priceB=[priceB;priceB(nB)];
    nB = size(priceB,1);
end




% Initialize the indices for A and B
indexA = 1;
indexB = find(timeB<=timeA(1), 1,'last');
if overlap>0
    indexB = max(indexB-overlap,1); % Make sure that indexB is >= 1
end

% Initialize rchy and the times
rchy = 0;
times = zeros(nA-1,4);
while indexA<nA
    % Get the two prices at indexA and indexA+1
    pA1 = priceA(indexA);
    pA2 = priceA(indexA+1);
    % get the times
    times(indexA-1,1) = timeA(indexA);
    times(indexA-1,2) = timeA(indexA+1);
    % Increment the index to A
    indexA = indexA+1;

    % Compute the index, including overlap if needed, make sure it doesn't end up before the first
    % index
    indexB_minus_overlap = max(indexB - overlap,1);
    % Get the first price of B
    pB1 = priceB(indexB_minus_overlap);
    % Get the time
    times(indexA-1,3) = timeB(indexB_minus_overlap);
    while timeB(indexB)<timeA(indexA)
        % The time stamp of B is strictly less than the time stamp of the second price in the return
        % to asset A so increment the counter.  Since we do the pre- and post- pend trick above,
        % there is no need to worry about uneven ending.
        indexB = indexB + 1;
    end
    % Compute the index plus the overlap, make sure it doesn't go outsie the number of obs
    indexB_plus_overlap = min(indexB + overlap,nB);
    % Get the second price
    pB2 = priceB(indexB_plus_overlap);
    % Get the final time
    times(indexA-1,4) = timeB(indexB_plus_overlap);
    % Add to rchy
    rchy = rchy + (pA2-pA1)*(pB2-pB1);

    % Rewind indexB by 1 since it is now after A
    if timeB(indexB)>timeA(indexA)
        indexB = indexB - 1;
    end
end
