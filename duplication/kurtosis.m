function k=kurtosis(x, flag, dim)
% Estimate the 4th standardized moment from a vector or matix of data.
%
%   K=kurtosis(X,FLAG,DIM)
%
% INPUTS
%   X     - Data to be used in skewness calculation
%   FLAG  - Not implemented
%   DIM   - Dimension to compute the kurtosis on if x is not a vector

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/25/2006

if nargin==1
    dim=1;
end

ncm1 = mean(x,dim);
ncm2 = mean(x.^2,dim);
ncm3 = mean(x.^3,dim);
ncm4 = mean(x.^4,dim);

cm2 = ncm2-ncm1.^2;
cm4 = ncm4 - 4*ncm3.*ncm1 + 6*ncm1.^2.*ncm2 - 4.*ncm1.^4 + ncm1.^4;

k = cm4./cm2.^2;
