function [prices, refreshTimes, actualTimes] = realized_refresh_time_bivariate(timeType,price1,time1,price2,time2)
% Computes refresh time synchronized returns.  See BNHLS (2008) Multivariate Realized Kernels.
%
% USAGE:
%   [REFRESHPRICES] = realized_refresh_time_bivariate(TIMETYPE,PRICE1,TIME1,PRICE2,TIME2)
%   [REFRESHPRICES,REFRESHTIMES,ACTUALTIMES] = realized_refresh_time_bivariate(TIMETYPE,PRICE1,TIME1,PRICE2,TIME2)
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
%   PRICE1           - m1 by 1 vector of prices
%   TIME1            - m1 by 1 vector of times, ascending and must be unique
%   PRICE2           - m2 by 1 vector of prices
%   TIME2            - m2 by 1 vector of times, ascending and must be unique
%
% OUTPUTS:
%   REFRESHPRICES    - n by 2 matrix of refresh time synchronized prices
%   REFRESHTIMES     - n by 1 vector of refresh times
%   ACTUALTIMES      - n by 2 matrix of actual observation times corresponding to REFRESHPRICES
%
% COMMENTS:
%   The price series inputs should generally be in tick time.  It is possible to compute a calendar
%   time-type refresh time price by sampling the most liquid asset in calendar time using
%   realized_price_filter and then running this function on the filtered prices of the most liquid
%   asset and the original prices of the others.
%
% EXAMPLES:
%   % Wall time prices
%   REFRESHPRICES = realized_refresh_prices('wall',price1,time1,price2,time2)
%
%  See also REALIZED_REFRESH_TIME, REALIZED_MULTIVARIATE_KERNEL, REALIZED_VARIANCE

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

m1 = size(price1,1);
m2 = size(price2,1);

%Find the first time
prices = zeros(min(m1,m2),2);
actualTimes = zeros(min(m1,m2),2);
% Just looping over price 1, in principle this shoudl be the shorted series
count = 1;
ind1 = 1;
ind2 = 1;
while ind1<=m1 && ind2<=m2
    if time1(ind1)<time2(ind2)
        % Update ind1 to find the largest time1<=time2
        while ind1<=m1 && time1(ind1)<time2(ind2) 
            ind1=ind1+1;
        end
        if ind1>m1
            ind1=m1;
        end
        if time1(ind1)>time2(ind2)
            ind1=ind1-1;
        end
    elseif  time2(ind2)<time1(ind1)
        while  ind2<=m2 && time2(ind2)<time1(ind1)
            ind2=ind2+1;
        end
        if ind2>m2
            ind2=m2;
        end
        if time2(ind2)>time1(ind1)
            ind2=ind2-1;
        end
    end
    actualTimes(count,1) = time1(ind1);
    actualTimes(count,2) = time2(ind2);
    prices(count,1) = price1(ind1);
    prices(count,2) = price2(ind2);
    ind1=ind1+1;
    ind2=ind2+1;
    count = count + 1;
end
% Need to append last refresh time
prices = prices(1:count-1,:);
actualTimes = actualTimes(1:count-1,:);
refreshTimes = max(actualTimes,[],2);