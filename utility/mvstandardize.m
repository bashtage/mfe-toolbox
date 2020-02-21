function st=mvstandardize(x,sigma,demean)
% Standardization of a matrix of data where each column returned will be
% mean zero and have variance 1, and will be uncorrelated
%
% USAGE:
%   ST = standardize(X)
%   ST = standardize(X,SIGMA,DEMEAN)
%
% INPUTS:
%   X  - A T by K matrix of data to be standardized
%
% OUTPUTS:
%   ST     - A T by K matrix of standardized data
%   SIGMA  - [OPTIONAL] Either a K by K covariance matrix or a K by K by T matrix of time-varying
%             covariances.  If ommitted, the sample covariance is used
%   DEMEAN - [OPTIONAL] A logical value indicating whether to demean the data (1) or not (0).
%              Default is 1
%
% COMMENTS:
%
% See also deman, standardize

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2005

switch nargin
    case 1
        sigma = cov(x);
        demean = true;
    case 2
        demean = true;
    case 3
    otherwise
end

if ndims(demean)~=2 || max(size(demean))~=1
    error('DEMEAN must be a logical scalar.')
end
[T,K]=size(x);
if ndims(x)>2
    error('X must be a T by K matrix with T>K')
end
if ndims(sigma)==2
    if any(size(sigma)~=K) || min(eig(sigma<=0))
        error('SIGMA by be a K by K positive definite matrix.')
    end
elseif ndims(sigma)==3
    if any(size(sigma)~=[K K T])
        error('If SIGMA is a 3-D matrix, it must be K by K by T.');
    end
    for i=1:T
        if min(eig(sigma(:,:,i))<=0)
            error('If SIGMA is K by K by T, each covariance must be positive definite')
        end
    end
else
    error('SIGMA by be either a K by K matrix, or a K by K by T matrix.')
end

if demean
    x=demean(x);
end

if ndims(sigma)==2
    st = x * sigma^(-0.5);
else
    st = zeros(T,K);
    for i=1:T
        st(i,:) = x(i,:) * sigma(:,:,i)^(-0.5);
    end
end
