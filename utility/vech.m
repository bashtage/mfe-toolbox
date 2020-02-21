function stackedData = vech(matrixData)
% Transform a symmetric matrix to its half-vec representation for use in covariance applications
% 
% USAGE:
%   STACKEDDATA = vech(MATRIXDATA)
% 
% INPUTS:
%   MATRIXDATA   - a K by K symmetric matrix
%
% OUTPUTS:
%   STACKEDDATA - A K(K+1)/2 vector of stacked data
%
% COMMENTS:
%   The data is stacked according to 
%     [ data(1) data(2)    data(3)     ...               data(K)
%       data(2) data(K+1)  data(K+2)   ...               ...
%       data(3) data(K+2)  data(2K)    ...               ...
%       ...     ....       ...         ...               data(K(K+1)/2-1)
%       data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]
%
% See also IVECH

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k,l] = size(matrixData);
if k~=l || any(any(matrixData~=matrixData'))
    error('MATRIXDATA must be a symmetric matrix');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sel = tril(true(k));
stackedData = matrixData(sel);