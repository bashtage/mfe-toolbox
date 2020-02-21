function [trend, cyclic, noise] = bkfilter(y,p,q,k)
% Baxter-King filtering of multiple time series
%
% USAGE:
%   [TREND,CYCLIC,NOISE] = bkfilter(Y,P,Q,K)
%
% INPUTS:
%   Y      - A T by K matrix of data to be filtered.
%   P      - Number of periods to use in the higher frequency filter (e.g. 6 for quarterly data).
%              Must be at least 1.
%   Q      - Number of periods to use in the lower frequency filter (e.g. 32 for quarterly data). Q
%              can be inf, in which case the low pass filter is a 2K+1 moving average.
%   K      - [OPTIONAL] Number of points to use in the finite approximation bandpass filter.  The
%              default value is 12.  The filter throws away the first  and last K points.
%
% OUTPUTS:
%   TREND  - A T by K matrix containing the filtered trend.  The first and last K points equal Y.
%   CYCLIC - A T by K matrix containing the filtered cyclic component. The first and last K points are 0.
%   NOISE  - A T by K matrix containing the filtered noise component.  The first and last K points are 0.
%
% COMMENTS:
%   The noise component is simply the original data minus the trend and cyclic component, NOISE = Y
%   - TREND - CYCLIC where the trend is produced by the low pasf filter and the cyclic component is
%   produced by the difference of the high pass filter and the low pass filter.  The recommended
%   values of P and Q are 6 and 32 or 40 for quarterly data, or 18 and 96 or 120 for monthly data.
%   Setting Q=P produces a single bandpass filer and the cyclic component will be 0.
%
% EXAMPLES:
%   Load US GDP data
%       load GDP
%   Standard BK Filter with periods of 6 and 32
%       [trend, cyclic] = bkfilter(log(GDP),6,32)
%   BK Filter for low pass filtering only at 40 period, CYCLIC will be 0
%       [trend, cyclic] = bkfilter(log(GDP),40,40)
%   BK Filter using a 2-sided 20 point approximation
%       trend = bkfilter(log(GDP),6,32,20) 
%
% See also HP_FILTER, BEVERIDGENELSON

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 19/10/2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        k = 12;
    case 4
        % Nothing
    otherwise
end
if ndims(y)>2 || size(y,1)<k/2
    error('Y must be a T by K matrix with T>=4');
end
if ~isscalar(p) || p<=0 || p<1 
    error('P must be a positive scalar.');
end 
if ~isscalar(q) || q<=0 || q<p
    error('Q must be a positive scalar with Q>=P.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[T,n] = size(y);
% Compute angular frequency
fp = (2*pi)/p;
fq = (2*pi)/q;

bp=zeros(2*k+1,1);
bp(k+1) = fp/pi;
weights = (sin((1:k)*fp)./((1:k).*pi))';
bp(1:k)=flipud(weights);
bp(k+2:2*k+1)=weights;
thetap = (1-sum(bp))/(2*k+1);
bp = bp + thetap;

bq=zeros(2*k+1,1);
bq(k+1) = fq/pi;
weights = (sin((1:k)*fq)./((1:k).*pi))';
bq(1:k)=flipud(weights);
bq(k+2:2*k+1)=weights;
thetaq = (1-sum(bq))/(2*k+1);
bq = bq + thetaq;

b = bp - bq;

trend = y;
cyclic = zeros(T,n);
noise = zeros(T,n);

for t=k+1:T-k
    % Trend is bq(y)
    trend(t,:) = bq'*y(t-k:t+k,:);
    % Cyclic is bp(y) - bq(y)
    cyclic(t,:) = b'*y(t-k:t+k,:);
    % Noise is y - bp(y)
    noise(t,:) = y(t,:) - bp'*y(t-k:t+k,:);
end