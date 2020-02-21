function [medianPrice,medianTime,totalVol,nObs] = realized_compute_median(price,time,volume)
% Computes the median price at each time stamp for a vector of prices which may have multiple
% observations with the same time stamp.
%
% USAGE:
%   [MEDIANPRICE,MEDIANTIME,TOTALVOL,NOBS] = realized_compute_median(PRICE,TIME,VOLUME)
%
% INPUTS:
%   PRICE            - m by 1 vector of high frequency prices
%   TIME             - m by 1 vector of times where TIME(i) corresponds to PRICE(i), sorted in
%                        ascending order
%   VOLUME           - [OPTIONAL] m by 1 vector of transaction volumes or bid/ask size.  If omitted
%                        TOTALVOLUME will be a vector of 0's.
%
% OUTPUTS:
%   MEDIANPRICE      - n by 1 vector of median prices where n is the number of unique elements of TIME
%   MEDIANTIME       - n by 1 vector of time stamps corresponding to MEDIANPRICE
%   TOTALVOLUME      - n by 1 vector capturing the total volume at each time stamp
%   NOBS             - n by 1 vector indicating number of prices at each element of MEDIANTIME
%
% COMMENTS:
%   This is a helper function for most of the Realized toolkit.  Most function in the Realized
%   toolkit expect time stamps to be unique.  Median price is a reasonable and fairly robust (to
%   noise) method of computing a unique price for each time stamp.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 6/12/2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InputChecking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || nargin>3
    error('Six or seven inputs required.')
end
if size(price,2)>size(price,1)
    price=price';
end
if size(price,2)>1
    error('PRICE must be a m by 1 vector.')
end
m = size(price,1);

if size(time,2)>size(time,1)
    time=time';
end
if any(diff(time)<0)
    error('TIME must be sorted and increasing')
end
if size(time,2)>1 || length(time)~=m
    error('TIME must be a m by 1 vector.')
end
time = double(time);

if nargin == 3
    if size(volume,2)>size(volume,1)
        volume=volume';
    end
    if size(volume,2)>1 || length(volume)~=m
        error('VOLUME must be a m by 1 vector.')
    end
else
    volume = zeros(m,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure time and volume are double since I store it was uint32
time = double(time);
volume = double(volume);

% Quickly find the time change points
pl = find(diff(time));
pl = [1;pl+1;length(time)+1];

% Get the number of observations
N = length(pl);

% Pre-allocatte
medianPrice = zeros(N-1,1);
medianTime = time(pl(1:N-1));
totalVol = zeros(N-1,1);
nObs = zeros(N-1,1);

% Loop over unique sizes, with special cases for 1 and 2
[sizes,~,ind]=unique(diff(pl));
for i=1:length(sizes)
    j=find(ind==i);
    loc = bsxfun(@plus,pl(j),0:sizes(i)-1)';
    switch sizes(i)
        case 1
            medianPrice(j) = price(pl(j));
            totalVol(j) = volume(pl(j));
        case 2
            medianPrice(j) = (price(pl(j)) + price(pl(j)+1))/2;
        otherwise
            sortedPrice = sort(price(loc));
            if mod(sizes(i),2)==0
                medianPrice(j) = mean(sortedPrice([sizes(i)/2;sizes(i)/2+1],:))';
            else
                medianPrice(j) = sortedPrice(ceil(sizes(i)/2),:)';
            end
    end
    totalVol(j) = sum(volume(loc))';
    nObs(j) = sizes(i);
end

% Old, slower method
% Loop
% singleton = pldiff==1;
% medianPrice(singleton) = price(pl(singleton));
% totalVol(singleton) = volume(pl(singleton));
% nObs(singleton) = 1;
% notSingleton = find(~singleton);
% 
% for i=1:length(notSingleton)
%     j = notSingleton(i);
%     medianPrice(j) = median(price(pl(j):pl(j+1)-1));
%     totalVol(j) = sum(volume(pl(j):pl(j+1)-1));
%     nObs(j) = pl(j+1)-pl(j)+1;
% end
