function [lm, pval] = lmtest1(data,q,robust)
% LM tests for the presence of serial correlation in q lags, with or without heteroskedasticity.
% Returns Q LM tests, one for tests for each lag between 1 and Q. Under the null of no serial
% correlation, the LM-test is asymptotically  distributed X2(q)  
% 
% USAGE:
%  [LM,PVAL] = lmtest1(DATA,Q)
%  [LM,PVAL] = lmtest1(DATA,Q,ROBUST)
% 
% INPUTS:
%  DATA      - A T by 1 vector of data
%  Q         - The maximum number of lags to regress on.  The statistic and pval will be returned
%                for all sets of lags up to and including q 
%  ROBUST    - [OPTIONAL] Logical variable (0 (non-robust) or 1 (robust)) to  indicate whether
%                heteroskedasticity robust standard errors  should be used. Default is to use robust
%                standard errors (ROBUST=1).  
% 
% OUTPUTS:
%  LM        - A Qx1 vector of statistics
%  PVAL      - A Qx1 set of appropriate pvals
% 
% 
% COMMENTS:
% The variance estimator is computed under the alternative (to increase the power of the test).  As a
% result, this test is an LR-class test.  It is otherwise identical to the usual LM test for serial
% correlation.  
%
% SEE ALSO:
%  LJUNGBOX
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>3 || nargin<2
    error('2 or 3 inputs required.')
elseif nargin==2
    robust = true;
end
T=length(data);
if T<q
    error('At least Q observations requires')
end
if size(data,1)~=T,
    data=data';
end
if size(data,2)~=1
    error('DATA must be a column vector')
end
if ~isscalar(q) && q>0 && floor(q)==q
    error('Q must be a positive integer')
end
if nargin==3
    if ~ismember(robust,[0,1])
        error('ROBUST must be either 0 or 1')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lm=zeros(q,1);
for Q=1:q
   data=data-mean(data);
   [y,x]=newlagmatrix(data,Q,0);
   e = y;
   s = x(:,1:Q).*repmat(e,1,Q);
   sbar = mean(s);
   % Transformation to LR-like test
   if robust
       s = bsxfun(@minus, s, sbar);
       S=s'*s/T;
   else
       e = e - mean(e);
       S=(e'*e/T)*(x'*x)/T;
   end
   lm(Q)=T*sbar*S^(-1)*sbar';
end
pval=1-chi2cdf(lm,(1:q)');