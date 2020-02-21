function [prices, refreshTimes, actualTimes] = realized_refresh_time(timeType,varargin)
% Computes refresh time synchronized returns.  See BNHLS (2008) Multivariate Realized Kernels.
%
% USAGE:
%   [REFRESHPRICES] = realized_refresh_time(TIMETYPE,PRICE01,TIME01,PRICE02,TIME02,...)
%   [REFRESHPRICES,REFRESHTIMES,ACTUALTIMES] = realized_refresh_time(TIMETYPE,PRICE01,TIME01,PRICE02,TIME02,...)
%
% INPUTS:
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on
%                       the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the
%                         other types and can be applied to data from more
%                         than one calendar day
%   PRICEXX          - m by 1 vector of prices
%   TIMEXX           - m by 1 vector of times, ascending and must be unique
%                        NOTE: An arbitrarily large number of prices and times can be input, but must be at least 2.
%
% OUTPUTS:
%   REFRESHPRICES    - n by K matrix of refresh time synchronized prices
%   REFRESHTIMES     - n by 1 vector of refresh times
%   ACTUALTIMES      - n by K matrix of actual observation times corresponding to REFRESHPRICES
%
% COMMENTS:
%   The price series inputs should generally be in tick time.  It is possible to compute a calendar
%   time-type refresh time price by sampling the most liquid asset in calendar time using
%   realized_price_filter and then running this function on the filtered prices of the most liquid
%   asset and the original prices of the others.
%
% EXAMPLES:
%   % Wall time prices for 3 assets
%   REFRESHPRICES = realized_refresh_prices('wall',price1,time1,price2,time2,price3,time3)
%
%   % Calendar time-type sampling using 1 minute prices
%   [FILTEREDPRICE1,FILTEREDTIME1,ACTUALTIME1] = realized_price_filter(price1,time1,'wall','CalendarTime',300)
%   REFRESHPRICES = realized_refresh_prices('wall',FILTEREDPRICE1,ACTUALTIME1,price2,time2,price3,time3)
%
%  See also REALIZED_PRICE_FILTER, REALIZED_MULTIVARIATE_KERNEL, REALIZED_VARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3 || floor(nargin/2)==(nargin/2)
    error('The number of inputs must be odd and at least 3.')
end
timeType=lower(timeType);
if ~ismember(timeType,{'wall','seconds','unit'})
    error('TIMETYPE must be one of ''wall'', ''seconds'' or ''unit''.');
end

n = size(varargin,2)/2;
price = cell(n,1);
time = cell(n,1);
for i=1:n
    tempPrice = varargin{2*(i-1)+1};
    tempTime = varargin{2*(i-1)+2};

    % Inserted to protect against inputing integer times
    tempTime = double(tempTime);

    if size(tempTime,2)>size(tempTime,1)
        tempTime = tempTime';
    end
    if size(tempPrice,2)>size(tempPrice,1)
        tempPrice = tempPrice';
    end
    if size(tempTime,2)>1 || size(tempPrice,2)>1 || length(tempPrice)~=length(tempTime)
        error(['Problem reading price and time pair #' num2str(i) '.  Both PRICE' num2str(i) ' and TIME' num2str(i) ' must be column vectors with the same number of elements.'])
    end
    price{i} = tempPrice;
    time{i}  = tempTime;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the union of all the times
utimes = time{1};
for i=2:n
    utimes = union(utimes,time{i});
end

m = size(utimes,1);
% Construct an indicator indicating which of the union of times is present in each time vector
timeIndicator = false(m,n);
for i=1:n
    timeIndicator(:,i) = ismember(utimes,time{i});
end


% This row vector will indicate if prices have been refreshed
refreshIndicator = false(1,n);
refreshTimes = false(m,1);
for i=1:m
    % Use an or.  When all are 1 all prices have been refreshed
    refreshIndicator = refreshIndicator | timeIndicator(i,:);
    if all(refreshIndicator)
        % Then this is a refresh time;
        refreshTimes(i) = true;
        % Reser refresh indicator
        refreshIndicator = false(1,n);
    end
end
% Map the refreshtimes
refreshTimes = utimes(refreshTimes);


% And finally the prices and actual times
prices = zeros(length(refreshTimes),n);
actualTimes = zeros(length(refreshTimes),n);
% Use realized_price_filter to do tha real work using the refresh times and fixed sampling
for i=1:n
    [prices(:,i),temp,actualTimes(:,i)] = realized_price_filter(price{i},time{i},timeType,'Fixed',refreshTimes);
end
