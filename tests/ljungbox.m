function [q, pval] = ljungbox(data, lags)
% Ljung-Box tests for the presence of serial correlation in up to q lags.  Returns LAGS Ljung-Box
% statistics tests, one for tests for each lag between 1 and LAGS. Under the null of no serial
% correlation and assuming homoskedasticity, the Ljung-Box test statistic is asymptotically
% distributed X2(q)   
% 
% USAGE:
%  [Q,PVAL] = ljungbox(DATA,LAGS)
% 
% INPUTS:
%  DATA      - A T by 1 vector of data
%  LAGS      - The maximum number of lags to compute the LB.  The statistic and pval will be
%                returned for all sets of lags up to and including LAGS 
% 
% OUTPUTS:
%  Q         - A LAGS by 1 vector of Q statistics
%  PVAL      - A LAGS by 1 set of appropriate pvals
% 
% COMMENTS:
%  This test statistic is common but often inappropriate since it assumes homoskedasticity.  For a
%  heteroskedasticity consistent serial correlation test, see lmtest1 
%
% SEE ALSO:
%  LMTEST1, SACF, SPACF
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=2
    error('2 inputs required.')
end
T=length(data);
if T<=lags
    error('At least LAGS observations requires')
end
if size(data,1)~=T,
    data=data';
end
if size(data,2)~=1
    error('DATA must be a column vector')
end
if ~isscalar(lags) && lags>0 && floor(lags)==lags
    error('LAGS must be a positive integer')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ac = sacf(data,lags,0,0);

q = zeros(lags,1);
T = length(data);
for L=1:lags
    q(L)=T*(T+2)*sum(ac(1:L).^2./(T-(1:L)'));
end
pval = 1 - chi2cdf(q,(1:lags)');