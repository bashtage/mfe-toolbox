function sk=skewness(x, flag, dim)
% Computes the standardized 3rd moment
%
% USAGE 
%   SK=skewness(X,FLAG,DIM)
%
% INPUTS
%   X     - Data to be used in skewness calculation
%   FLAG  - Not implemented
%   DIM   - Dimension to compute the skewness on if x is not a vector

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/25/2006

if nargin==1
    dim=1;
end

ncm1 = mean(x,dim);
ncm2 = mean(x.^2,dim);
ncm3 = mean(x.^3,dim);

cm2 = ncm2-ncm1.^2;
cm3 = ncm3 - 3*ncm1.*ncm2 + 2.*ncm1.^3;

sk = cm3./cm2.^(1.5);
