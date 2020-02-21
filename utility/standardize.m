function st=standardize(x,demean)
% Standardization of a matrix of data where each column returned will be
% mean zero and have variance 1
%
% USAGE:
%   ST = standardize(X)
%   ST = standardize(X,DEMEAN)
%
% INPUTS:
%   X      - A T by K matrix of data to be standardized
%   DEMEAN - [OPTIONAL] A logical value indicating whether to demean the data (1) or not (0).
%              Default is 1
%
% OUTPUTS:
%   ST - A T by K matrix of standardized data
%
% COMMENTS:
%   For each column j, st(:,j)=x(:,j)-mean(x(:,j))
%                               ----------------
%                                  std(x(:,j))

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2005

switch nargin
    case 1
        demean = true;
    case 2
    otherwise
        error('1 or 2 inputs only.')
end

if ndims(demean)~=2 || max(size(demean))~=1
    error('DEMEAN must be a logical scalar.')
end

[T,K]=size(x);
mu=mean(x);
stdev=std(x);
if demean
    st=(x-repmat(mu,T,1))./repmat(stdev,T,1);
else
    st=x./repmat(stdev,T,1);
end