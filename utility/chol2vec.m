function stackedData = chol2vec(matrixData)
% Transform a lower triangular matrix its half-vec representation for use in covariance applications
% 
% USAGE:
%   STACKEDDATA = chol2vec(MATRIXDATA)
% 
% INPUTS:
%   MATRIXDATA   - a K by K lower triangular matrix
%
% OUTPUTS:
%   STACKEDDATA - A K(K+1)/2 vector of stacked data
%
% COMMENTS:
%   The data is stacked according to 
%     [ data(1) 0          0           ...    ...               0
%       data(2) data(K+1)  0           ...    ...               0
%       data(3) data(K+2)  data(K+2)   0      ...               0
%       ...     ....       ...         ...    ...               0
%       data(K) data(2K-1) ...         ...    data(K(K+1)/2-1)  data(K(K+1)/2) ]
%
% See also vec2chol

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k,l] = size(matrixData);
pl = ~tril(true(k));
if k~=l || any(matrixData(pl)~=0)
    error('MATRIXDATA must be a lower triangular matrix');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sel = tril(true(k));
stackedData = matrixData(sel);
