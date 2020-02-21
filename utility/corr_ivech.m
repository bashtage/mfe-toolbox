function matrixData=corr_ivech(stackedData)
% Transform a vector into a symmetric correlation matrix for use in correlation applications
% 
% USAGE:
%   MATRIXDATA = corr_ivech(STACKEDDATA)
% 
% INPUTS:
%   STACKEDDATA:   A K(K-1)/2 vector of data to be transformed 
%
% OUTPUTS:
%   MATRIXDATA   - a K by K symmetric matrix
%
% COMMENTS:
%   The data is stacked according to 
%     [ 1         data(1)    data(2)       ...   data(K-1)
%       data(1)   1          data(K)       ...   ...
%       data(2)   data(K)    data(2K-4)    ...   ...
%       ...       ....       ...           ...   data(K(K-1)/2)
%       data(K-1) data(2K-3) ...           ...   1              ]
%
% See also IVECH, CORR_VECH

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(stackedData,2)>size(stackedData,1)
    stackedData=stackedData';
end

if size(stackedData,2)~=1
    error('STACKED_DATA must be a column vector.')
end

K2=size(stackedData,1);
K=(-1+sqrt(1+8*K2))/2 + 1;

if floor(K)~=K
    error(['The number of elemeents in STACKED_DATA must be conformable to' ...
    'the inverse vech operation.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrixData = zeros(K);
loc = ~triu(true(K));
matrixData(loc) = stackedData;
matrixData = matrixData + matrixData' + eye(K);