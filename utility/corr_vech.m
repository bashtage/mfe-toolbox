function stackedData = corr_vech(matrixData)
% Transform a correlation matrix into a vector for use in correlation applications
%
% USAGE:
%   STACKEDDATA = corr_vech(MATRIXDATA)
%
% INPUTS:
%   MATRIXDATA   - a K by K symmetric matrix 
%
% OUTPUTS:
%   STACKEDDATA:   A K(K-1)/2 vector of stacked data 
%
% COMMENTS:
%   The data is stacked according to 
%     [ 1         data(1)    data(2)       ...   data(K-1)
%       data(1)   1          data(K)       ...   ...
%       data(2)   data(K)    data(2K-4)    ...   ...
%       ...       ....       ...           ...   data(K(K-1)/2)
%       data(K-1) data(2K-3) ...           ...   1              ]
%
% See also VECH, CORR_IVECH

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
sel = ~triu(true(k));
stackedData = matrixData(sel);
 