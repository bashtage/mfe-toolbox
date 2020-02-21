function x = demean(y)
% Demeans an N-dimensional vector
% 
% USAGE:
%   DEMEANED = vech(DATA)
% 
% INPUTS:
%   DATA     - An N-dimensional matrix
%
% OUTPUTS:
%   DEMEANED - A mean 0 matrix with the same size as DATA  (column-by-column)
%
% COMMENTS:
%
% See also STANDARDIZE

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2010

mu = mean(y,1);
x = bsxfun(@minus,y,mu);

