function [rqv,diagnostics]=realized_qmle_variance(price,time,timeType,samplingType,samplingInterval,options)
% Estimated quadratic variation using QMLE
%
% USAGE:
%   [RQV] = realized_qmle_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL)
%   [RQV,DIAGNOSTICS] = realized_qmle_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,OPTIONS)
%
% INPUTS:
%   PRICE            - m by 1 vector of high frequency prices
%   TIME             - [OPTIONAL] m by 1 vector of times where TIME(i) corresponds to PRICE(i).
%   TIMETYPE         - [OPTIONAL] String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS.mmm, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGTYPE     - [OPTIONAL] String describing the type of sampling to use when
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
%   SAMPLINGINTERVAL  - [OPTIONAL] Scalar integer or n by 1 vector whose meaning depends on the
%                         selected SAMPLINGTYPE
%   OPTIONS           - [OPTIONAL] Realized option structure initialized by calling
%                         realized_options('QMLE'). See help realized_options for a description of
%                         available options.
%
% OUTPUTS:
%   RQV               - Realized variace estimate as computed by QMLE
%   DIAGNOSTICS       - Structure of useful diagnostic information.  Fields are
%                         NOISEVARIANCE - Estimated noise variance as computed using QMLE
%                         ITERATIONS    - Number of iterations needed to converge
%
% COMMENTS:
%
% EXAMPLES:
%   % Realized QMLE variance using all data
%   rqv = realized_qmle_variance(price)
%   % Realized QMLE variance using every 10 observations, to avoid dependence in the MMN
%   rqv = realized_qmle_variance(price,[],[],'BusinessTime',10)

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 3/10/2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        time=[];
        timeType=[];
        samplingType=[];
        samplingInterval=[];
        options=[];
    case 2
        timeType=[];
        samplingType=[];
        samplingInterval=[];
        options=[];
    case 3
        samplingType=[];
        samplingInterval=[];
        options=[];
    case 4
        samplingInterval=[];
        options=[];
    case 5
        options=[];
    case 6
        % Nothing
    otherwise
         error('One to six inputs required.')
end

if isempty(timeType) && ~isempty(time)
    error('TIMETYPE must be provided if TIME is user supplied.')
end
    
if isempty(time)
    time = (1:length(price))'/price;
end

if isempty(timeType)
    timeType = 'unit';
end

if isempty(samplingType)
    timeType = 'BusinessTime';
end
if isempty(samplingInterval)
    samplingInterval = 1;
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

if nargin==6 && ~isempty(subsamples)
    if ~isscalar(subsamples) || subsamples<0 || floor(subsamples)~=subsamples
        error('SUBSAMPLES must be a non-negative scalar.')
    end
end

if nargin<7 || isempty(options)
    options = realized_options('QMLE');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logPrice = log(price);
[filteredLogPrice,filteredTimes] = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);
% Need starting values for sigma2 and omega2
noisePrice = realized_price_filter(logPrice,time,timeType,options.noiseVarianceSamplingType,options.noiseVarianceSamplingInterval);
noiseReturn = diff(noisePrice);
omega2 = noiseReturn'*noiseReturn/length(noiseReturn);

r = diff(filteredLogPrice);
n = size(r,1);
sigma2 = realized_variance(logPrice,time,timeType,options.medFrequencySamplingType,options.medFrequencySamplingInterval);
delta = 1/n;
L = spalloc(n,n,4+3*(n-2));
diagInd = (0:n-1)*n+(1:n);
uDiagInd = (1:n-1)*n+(1:n-1);
lDiagInd = (0:n-2)*n+(2:n);
L(diagInd)=2;
L(uDiagInd)=-1;
L(lDiagInd)=-1;
O = spalloc(n,n,4+3*(n-2));

theta = [sigma2;omega2];
thetaOld = [-1;-1];
iter = 0;
while max(abs(theta-thetaOld)>1e-3) && iter<20
    thetaOld = theta;
    O(diagInd) = sigma2 * delta + 2*omega2;
    O(lDiagInd) = -omega2;
    O(uDiagInd) = -omega2;
    Om1 = inv(O);
    Om1 = (Om1+Om1')/2;
    Om2= inv(O*O);
    Om2 = (Om2+Om2')/2;
    Om1L = Om1*L;
    Om2L = Om2*L;
    Om2L2 = Om2L*L;
    Om1LOm1 = Om1L*Om1;
    Om1LOm1 = (Om1LOm1+Om1LOm1')/2;
    
    trOm2L = diag(Om2L)'*ones(n,1);
    trOm2L2 = diag(Om2L2)'*ones(n,1);
    trOm2 = diag(Om2)'*ones(n,1);
    denom = (trOm2L)^2 - (trOm2) * (trOm2L2);
    W1 = n * trOm2L * Om1LOm1 - n * trOm2L2 * Om2;
    W1 = W1 / denom;
    W1 = (W1+W1')/2;
    W2 = trOm2L * Om2 - trOm2 * Om1LOm1;
    W2 = W2/denom;
    W2 = (W2+W2')/2;
    
    sigma2 = r'*W1*r;
    omega2 = r'*W2*r;
    omega2 = max(0,omega2);
    theta = [sigma2;omega2];
    iter = iter + 1;
    disp(iter)
    disp(theta)
end

if max(abs(theta-thetaOld)>1e-3)  && iter>=20
    warning('MFE::Convergence','The realized QMLE estimator did not converge.');
end


rqv = sigma2;
diagnostics.iterations = iter;
diagnostics.noiseVariance = omega2;