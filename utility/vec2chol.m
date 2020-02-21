function matrixData= vec2chol(stackedData)
% Transform a lower triangular matrix its half-vec representation for use in covariance applications
% 
% USAGE:
%   MATRIXDATA = chol2vec(STACKEDDATA)
% 
% INPUTS:
%   MATRIXDATA  - A K(K+1)/2 vector of stacked data
%
% OUTPUTS:
%   STACKEDDATA - a K by K lower triangular matrix
%
% COMMENTS:
%   The output data is stacked according to 
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
if size(stackedData,2)>size(stackedData,1)
    stackedData=stackedData';
end

if size(stackedData,2)~=1
    error('STACKED_DATA must be a column vector.')
end

K2=size(stackedData,1);
K=(-1+sqrt(1+8*K2))/2;

if floor(K)~=K
    error(['The number of elemeents in STACKED_DATA must be conformable to' ...
    'the inverse chol2vec operation.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the output data
matrixData=zeros(K);

% Use a logical trick to inverse the vech
pl=tril(true(K));
matrixData(pl)=stackedData;
