function [statistic,pval,H] = jarquebera(data, K, alpha)
% Computes the Jarque-Bera test for normality using the skewness and kurtosis to determine if a
% distribution is normal. 
%
% USAGE:
%   [STATISTIC] = jarquebera(DATA)
%   [STATISTIC,PVAL,H] = jarquebera(DATA,K,ALPHA)
%
% INPUTS:
%   DATA        - A set of data to be tested for deviations from normality
%   K           - [OPTIONAL] The number of dependant variables if any used in constructing the errors
%                   (if omitted K=2) 
%   ALPHA       - [OPTIONAL] The level of the test used for the null of normality.  Default is .05
%
% OUTPUTS:
%   STATISTIC - A scalar representing the statistic
%   PVAL      - A scalar pval of the null
%   H         - A hypothesis dummy (0 for fail to reject the null of normality, 1 otherwise)
%
% COMMENTS:
%   The data entered can be mean 0 or not.  In either case the sample mean is subtracted and the
%   data are standardized by the sample standard deviation before computing the statistic .
%
% EXAMPLES:
%   J-B test on normal data
%       x = randn(100,1);
%       [statistic, pval] = jarquebera(x);
%   J-B test on regression errors where there were 4 regressors (4 mean parameters + 1 variance)
%       x = randn(100,1);
%       [statistic, pval] = jarquebera(x, 5)


% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/1/2008


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        K=2;
        alpha=.05;
    case 2
        alpha=.05;
    case 3
    otherwise
        error('1 to 3 input arguments required')
end


[T,C]=size(data);
if ~any([T C]==1)
    error('DATA must be a vector')
end

if ~isscalar(alpha) || any(alpha>1) || any(alpha<0)
    error('ALPHA must be a scalar value between 1 and 0')
end


if ~isscalar(K) || any(K<0) || K~=floor(K)
    error('K must be an integer scalar value greater than or equal to 0.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = data - mean(data);
z = z/std(z);
statistic = (T-K)*(mean(z.^3)^2/6+(mean(z.^4)-3)^2/24);
pval = 1-chi2cdf(statistic,2);
H= pval < alpha;