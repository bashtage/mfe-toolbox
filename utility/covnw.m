function V=covnw(data,nlag,demean)
% Long-run covariance estimation using Newey-West (Bartlett) weights 
%  
% USAGE:
%   [V] = covnw(DATA)
%   [V] = covnw(DATA,NLAG,DEMEAN)
%
% INPUTS:
%   DATA   - T by K vector of dependent data
%   NLAG   - Non-negative integer containing the lag length to use.  If empty or not included,
%              NLAG=min(floor(1.2*T^(1/3)),T) is used 
%   DEMEAN - Logical true or false (0 or 1) indicating whether the mean should be subtracted when
%              computing the covariance 
%
% OUTPUTS:
%   V      - A K by K covariance matrix estimated using Newey-West (Bartlett) weights
%   
% COMMENTS:
%
% EXAMPLES:
%   Simulate an AR(1)
%       y = armaxfilter_simulate(1000,0,1,.9);
%   Newey-West covariance with automatic BW selection
%       lrcov = covnw(y)
%   Newey-West covariance with 10 lags
%       lrcov = covnw(y, 10)
%   Newey-West covariance with 10 lags and no demeaning
%       lrcov = covnw(y, 10, 0)
%
% See also COVVAR
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 5/1/2007
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=size(data,1);
if nargin==1
    nlag=min(floor(1.2*T^(1/3)),T);
    demean=true;
elseif nargin==2
    demean=true;    
end    
if isempty(nlag)
    nlag=min(floor(1.2*T^(1/3)),T);
end
if isempty(demean)
    demean=true;
end
if ~ismember(demean,[0 1]) 
    error('DEMEAN must be either logical true or false.')
end
if floor(nlag)~=nlag || nlag<0 
    error('NLAG must be a non-negative integer.')
end
if ndims(data)>2
    error('DATA must be a T by K matrix of data.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if demean
    data=data-repmat(mean(data),T,1);
end
 
% NW weights
w=(nlag+1-(0:nlag))./(nlag+1);
% Start the covariance
V=data'*data/T;
for i=1:nlag
    Gammai=(data((i+1):T,:)'*data(1:T-i,:))/T;
    GplusGprime=Gammai+Gammai';
    V=V+w(i+1)*GplusGprime;
end


