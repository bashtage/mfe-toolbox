function [rq,rqSS,diagnostics] = realized_quantile_variance(price,time,timeType,samplingType,samplingInterval,quantiles,blockSize,symmetric,overlap,subsamples)
% Computes the Quantile Realized Variance of Christensen, Oomen and Podolskij (2008), the MinRV and
% MedRV estimator of Andersen, Dobrev and Shaumberg and the general symmetrized version suggested by
% Sheppard in a discussion of COP 
%
% USAGE:
%   [RQ] = realized_quantile_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,QUANTILES,BLOCKSIZE)
%   [RQ,RQSS,DIAGNOSTICS] = realized_quantile_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,QUANTILES,BLOCKSIZE,SYMMETRIC,OVERLAP,SUBSAMPLES)
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
%   QUANTILES        - k by 1 vector of quantile values to use when computing RQ.
%                        When SYMMETRIC = false, 0.5<QUANTILES<=1.
%                        When SYMMETRIC = true, 0<QUANTILES<=1
%                        In either case, QUANTILES * BLOCKSIZE must be an integer.
%                        The simplest method to set the quantiles is to determine which of the
%                        ordered returns should be used in the QRV, and then to set
%                        QUANTILES = ORDER/BLOCKSIZE.  For example, when SYMETRIC = FALSE, and
%                        BLOCKSIZE = 20, QUANTILES = [13 15 18]/20 will use the returns in position
%                        13, 15 and 18 (as well as 8 6 and 3).  When SYMMETRIC = TRUE, QUANTILES =
%                        [5 8 10 13 15 18]/20 will use the absolute value of returns in positions 5,
%                        8, 10, 13, 15 and 18.  For a fixed block size, you should generally more
%                        quantiles when SYMMETRIC = TRUE.
%   BLOCKSIZE        - Number of returns to use in each block when computing the quantiles. NOTE: If
%                        OVERLAP = FALSE, then the number of returns produced by filtering according
%                        to SAMPLINGTYPE and SAMPLINGINTERVAL must be an integer multiple of
%                        BLOCKSIZE. Otherwise there is no constraint.
%   SYMMETRIC        - [OPTIONAL] Boolean indicating whether the base the estimator on the absolute value
%                        of returns (symmetric) of on the non-adjusted returns. Default value is
%                        FALSE (as used by COP).
%   OVERLAP          - [OPTIONAL] Boolean indicating whether to use all overlapping blocks (TRUE),
%                        or only non-overlapping blocks. Default value is TRUE.
%   SUBSAMPLES       - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                        variance estimators to average with the original realized variance.
%                        Subsample realized variances are based on prices uniformly spaced between
%                        the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES=1
%                        will compute a subsample realized variance using the mid-point of the
%                        price sample points, 2 will use 1/3 and 2/3, and so on. In general this
%                        number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   RQ               - Quantile realized variance vector (k by 1) corresponding to the values in QUANTILES
%   RQSS             - Subsample based version of RQ.  RQSS = RQ if SUBSAMPLES = 0 or omitted.
%   DIAGNOSTICS      - Structure with fields:
%                        RQINDIV   - k by 1 vector of RQ estimates
%                        RQSSINDIV - k by 1 vector of subsample RQ estimates
%                        RQWEIGHTS - Weights used in combining the individual RQ estimates.  RQ
%                                      weights are computed from RQCOV as the minimum variance
%                                      combination
%                        RQCOV     - k by k non-scaled (integrated quarticity) covariance matrix.
%
% COMMENTS:
%   The cases where SAMPLESPERBIN is in {4 5 6 10 13 15 18 20 25 26 30 36 39 50 60 65 72 75 100 144}
%   (non-symmetric) or {2 3 4 5 6 8 9 10 12 13 15 18 20 24 25 26 30 36 39 40 50 60 65 72 75 100 144}
%   (symmetric) the scales and non-scaled covariances have been pre-computed using Monte Carlo
%   integration using 100,000,000 simulations. Choosing another number for SAMPLESPERBIN requires a
%   run of realized_quantile_variance_scale which can be slow if SAMPLESPERBIN is large. If using
%   another value, 1,000,000 simulations will be used in computing the weights. If repeatedly using
%   another value, pre-computing the appropriate scales and non-scaled covariance using
%   realzed_quantile_weight_simulation is recommended.  
%
% EXAMPLES:
%   % Estimate the RQ from prices available from 9:30 to 16:00 using 1 minute returns, 30 samples per
%   % bin and quantiles [18 22 27]/30 ([0.6 0.733 0.9])
%   sampleTimes = seconds2wall(wall2seconds(93000):60:wall2seconds(160000))
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,[18 22 27]/30,30);
%
%   % The same but using 4 subsamples
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,[18 22 27]/30,30,[],[],4);
%
%   % Estimate the RQ from prices available from 9:30 to 16:00 using 1 minute returns, 30 samples per
%   % bin, the symmetric verion and quantiles [5 10 15 18 22 27]/30 ([0.15 0.33 0.5 0.6 0.733 0.9])
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,[5 10 15 18 22 27]/30,30,true);
%
%   % Estimate the MedRV of ABS
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,2/3,3,true);
%
%   % Estimate the MinRV of ABS
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,1/2,2,true);
%
%   % Estimate a MedRV-like RQV only using blocks of size 5
%   RQ=realized_quantile_variance(price,time,'wall','Fixed',sampleTimes,3/5,5,true);
%
%  See also REALIZED_MIN_MED_VARIANCE, REALIZED_KERNEL, REALIZED_VARIANCE,
%  REALIZED_BIPOWER_VARIATION, REALIZED_THRESHOLD_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InputChecking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<7 || nargin>10
    error('Seven to ten inputs required.')
end
switch nargin
    case 7
        symmetric = false;
        overlap = true;
        subsamples = 1;
    case 8
        overlap = true;
        subsamples = 1;
    case 9
        subsamples = 1;
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

% SAMPLERPERBIN, positive integer >=2
if blockSize<1 || floor(blockSize)~=blockSize
    error('SAMPLESPERBIN must be a positive integer greater than or equal to 2 (and usually at least 20)')
end

% SYMMETRIC
if isempty(symmetric)
    symmetric = false;
end
if ~isscalar(symmetric)
    error('SYMMETRIC must be a scalar logical value.')
end
symmetric = logical(symmetric);

% QUANTILES
if symmetric
    if any(quantiles<=0) || any(quantiles>1) || length(quantiles)>blockSize || length(unique(quantiles*blockSize)) ~= length(quantiles)
        error('QUANTILES must be unique, greater than 0 and less than or equal to 1, and the number of quantiles must be smaller than SAMPLESPERBIN when SYMMETRIC = true.')
    end
else
    if any(quantiles<=0.5) || any(quantiles>1) || length(quantiles)>(blockSize/2) || length(unique(quantiles*blockSize)) ~= length(quantiles)
        error('QUANTILES must be unique, greater than 0.5 and less than or equal to 1, and the number of quantiles must be smaller than 0.5*SAMPLESPERBIN when SYMMETRIC = false.')
    end
end

% OVERLAP
if isempty(overlap)
    overlap = true;
end
if ~isscalar(overlap)
    error('OVERLAP must be a scalar logical value.')
end
overlap = logical(overlap);

% SUBSAMPLES
if isempty(subsamples)
    subsamples = 1;
end
if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
    error('SUBSAMPLES must be a non-negative integer.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Compute returns
logPrice = log(price);
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
returns = diff(filteredLogPrice);
% Check that the number of returns is compatible with SAMPLESPERBIN
n=length(returns);
if ~overlap
    if floor(n/blockSize)~=(n/blockSize)
        error(['The number of returns computed from the prices returned from realized_price_filter(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL) must be an integer multiple of SAMPLESPERBIN.  The number of return produced is ' num2str(n)]);
    end
end

if overlap
    binStart = 1:n-blockSize+1;
else
    binStart = 1:blockSize:n;
end
returns = returns * sqrt(n);
if symmetric
    returns = abs(returns);
end
q = length(quantiles);
rqindiv = zeros(length(binStart),q);



% Compute the indices to use
if symmetric
    indices = round(quantiles*blockSize);
else
    upperIndices = round(quantiles*blockSize);
    lowerIndices = round((1-quantiles)*blockSize+1);
end
% Loop over the bins



count = 1;
for j=binStart
    binreturns = returns(j:j+blockSize-1);
    binreturns = sort(binreturns)';
    if symmetric
        rqindiv(count,:)=binreturns(indices).^2;
    else
        rqindiv(count,:)=binreturns(lowerIndices).^2+binreturns(upperIndices).^2;
    end
    count = count + 1;
end
rqindiv = mean(rqindiv);

[rq,rqindiv,rqcov,rqweights] = realized_quantile_variance_core(rqindiv,blockSize,quantiles,symmetric);

diagnostics.rqindiv = rqindiv;
diagnostics.rqcov = rqcov;
diagnostics.rqweights = rqweights;


subsampledLogPrice = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
rqindiv = nan(subsamples * length(binStart),q);
count = 1;
for i=1:subsamples
    filteredLogPrice = subsampledLogPrice{i};
    returns = diff(filteredLogPrice);
    n = length(returns);
    returns = returns * sqrt(n);
    if symmetric
        returns = abs(returns);
    end
    binStart = 1:n-blockSize+1;
    for j=binStart
        binreturns = returns(j:j+blockSize-1);
        binreturns = sort(binreturns)';
         if symmetric
            rqindiv(count,:)=binreturns(indices).^2;
         else
            rqindiv(count,:)=binreturns(lowerIndices).^2+binreturns(upperIndices).^2;
        end
        count = count + 1;
    end
end
rqindiv = rqindiv(1:count-1,:);
rqindiv = mean(rqindiv);

[rqSS,rqindivSS] = realized_quantile_variance_core(rqindiv,blockSize,quantiles,symmetric);
diagnostics.rqindivSS = rqindivSS;

end


function [rq,rqindiv,rqcov,rqweights] = realized_quantile_variance_core(rqindiv,blockSize,quantiles,symmetric)

scales = load('realized_quantile_scales');
if symmetric
    simulationBlockSize = scales.symmetricSimulationSamplerperbin;
    simulationCovariance = scales.symmetricExpectedCovariance;
    simulationExpectedQuantiles = scales.symmetricExpectedQuantiles;
    simulationQuantile = scales.symmetricSimulationQuantile;
else
    simulationBlockSize = scales.asymmetricSimulationSamplerperbin;
    simulationCovariance = scales.asymmetricExpectedCovariance;
    simulationExpectedQuantiles = scales.asymmetricExpectedQuantiles;
    simulationQuantile = scales.asymmetricSimulationQuantile;
    
end
simulationBlockSize = cell2mat(simulationBlockSize);
if ismember(blockSize,simulationBlockSize )
    % If available load and use the pre-computed value
    
    
    [~,pl]=ismember(blockSize,simulationBlockSize);
    thisScale = simulationExpectedQuantiles{pl};
    thisQuantile = simulationQuantile{pl};
    thisCovariance = simulationCovariance{pl};
    q = length(quantiles);
    indicator = zeros(1,q);
    for i = 1:q
        [~,indicator(i)] = min(abs(thisQuantile-quantiles(i)));
    end
    rqExpectedSquaredQuantileValue = thisScale(indicator);
    % Find index
    rqcov     = thisCovariance(indicator,indicator);
    rqCovInv = rqcov \ eye(q);
    rqweights = rqCovInv*ones(q,1)/(ones(1,q)*rqCovInv*ones(q,1));
else
    % Need to simulate using 1,000,000 simulations
    if blockSize>100
        warning('oxfordRealized:realizedQuantileVariance','Computing the scales needed.  Since SAMPLESPERBIN is very large this may take a long time.  \nConsider pre-computing this value, especially if using this value of SAMPLESPERBIN many times.')
    else
        warning('oxfordRealized:realizedQuantileVariance','Computing the scales needed.  \nConsider pre-computing this value, especially if using this value of SAMPLESPERBIN many times.')
    end
    [rqExpectedSquaredQuantileValue,rqweights,rqcov]=realized_quantile_variance_scale(adjSamplesPerBin,quantiles,1000000);
end


rqindiv  = rqindiv ./ rqExpectedSquaredQuantileValue;
rq = rqweights'*rqindiv';
end