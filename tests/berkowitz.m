function [stat,pval,H]=berkowitz(x,type,alpha,dist,varargin)
% Performs a Kolmogorov-Smirnov-like test using the Berkowitz transformation to a univariate normal
% that the data are from a specified distribution.
%
% USAGE:
%   [STAT,PVAL,H] = berkowitz(X,TYPE,ALPHA,DIST,VARARGIN)
%
% INPUTS:
%   X        -  A set of random variable to be tested for distributional correctness
%   TYPE     -  [OPTIONAL] A char string, either 'CS' if the data are cross-sectional or 'TS' for
%                 time series.  The TS checks for autocorrelation in the prob integral transforms
%                 while the CS does not.  'TS' is the default value.  
%   ALPHA    -  [OPTIONAL] The size for the test or use for computing H. 0.05 if not entered or
%                 empty.
%   DIST     -  [OPTIONAL] A char string of the name of the CDF of X, i.e. 'normcdf' for the normal,
%                 'stdtcdf' for standardized Studnet's T, etc.  If not provided or empty, data are
%                 assumed to have a uniform distribution (i.e. that data have already been fed
%                 through a probability integral transform)   
%   VARARGIN -  [OPTIONAL] Arguments passed to the CDF, such as the mean and variance for a normal
%                 or a d.f. for T.  The VARARGIN should be such that DIST(X,VARARGIN) is a valid
%                 function with the correct inputs.  
%
% OUTPUTS:
%   STAT     - The Berkowitz statistic computed as a likelihood ratio of normals
%   PVAL     - The asymptotic probability of significance
%   H        - 1 for reject the null that the distribution is correct using the size provided (or
%                .05 if not), 0 otherwise 
%
% EXAMPLES:
% Test uniform data from a TS model
%     stat = berkowitz(x);
% Test standard normal data from a TS model
%     [stat,pval] = berkowitz(x,'TS',[],'normcdf');
% Test normal mean 1, standard deviation 2 data from a TS model
%     [stat,pval] = berkowitz(x,'TS',[],'normcdf',1,2);
%
% COMMENTS:
%
% See also KOLMOGOROV
 
% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 3/1/2007
 
 
%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECK
%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        type='TS';
        alpha=.05;
        dist=[];
        varargin=[];
    case 2
        alpha=.05;
        dist=[];
        varargin=[];
    case 3
        dist=[];
        varargin=[];
    case 4
        varargin=[];
    case 0
        error('1 or more inputs required.')
end
if min(size(x))~=1 || ndims(x)~=2
    error('X must be a column vector of data')
end
if size(x,2)~=1
    x=x';    
end
if isempty(type)
    type='TS';
end
if isempty(dist)
    if any(x>=1) || any(x<=0)
        error('If DIST is not input, X must satisfy 0<X<1')
    end
else
    if ~ischar(dist)
        error('DIST much be a character string')
    end
end
if ~isscalar(alpha) || alpha<=0 || alpha>=1
    error('TESTSIZE must be a scalar between 0 and 1')
end
if ~ismember(type,{'TS','CS'})
    error('TYPE must be either ''TS'' or ''CS''')
end
%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECK
%%%%%%%%%%%%%%%%%%%
 
if ~isempty(dist)
    try
        if isempty(varargin)
            cdfvals=feval(dist,x);
        else
            cdfvals=feval(dist,x,varargin{:});
        end
    catch
        error('There was an error calling the DIST function with the VARARGIN provided.')
    end
else
    cdfvals=x;
end
workingdata=norminv(cdfvals);
 
if strcmp(type,'CS')
    mu=mean(workingdata);
    e=workingdata-mu;
    sigma2=e'*e/length(e);
    restrictedLL=normloglik(workingdata,0,1);
    unrestrictedLL=normloglik(e,0,sigma2);
    stat=-2*(restrictedLL-unrestrictedLL);
    pval=1-chi2cdf(stat,2);
    H=pval<alpha;
else
    [y,lags]=newlagmatrix(workingdata,1,1);
    p=lags\y;
    e=y-lags*p;
    sigma2=e'*e/length(e);
    restrictedLL=normloglik(y,0,1);
    unrestrictedLL=normloglik(e,0,sigma2);
    stat=-2*(restrictedLL-unrestrictedLL);
    pval=1-chi2cdf(stat,3);
    H=pval<alpha;
end;