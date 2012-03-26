function matrixData=ivech(stackedData)
% Transform a vector in to a symmetric matrix for use in covariance applications
% 
% USAGE:
%   MATRIXDATA = ivech(STACKEDDATA)
% 
% INPUTS:
%   STACKEDDATA:   A K(K+1)/2 vector of data to be transformed 
%
% OUTPUTS:
%   MATRIXDATA   - a K by K symmetric matrix of the form 
%                  [ data(1) data(2)    data(3)     ...               data(K)
%                    data(2) data(K+1)  data(K+2)   ...               ...
%                    data(3) data(K+2)  data(2K)    ...               ...
%                    ...     ....       ...         ...               data(K(K+1)/2-1)
%                    data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]
%
% COMMENTS:
%
% See also VECH

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
    'the inverse vech operation.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the output data
matrixData=zeros(K);

% Use a logical trick to inverse the vech
pl=tril(true(K));
matrixData(pl)=stackedData;
diag_matrixData=diag(diag(matrixData));
matrixData=matrixData+matrixData'-diag_matrixData;

