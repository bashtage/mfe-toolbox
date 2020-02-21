function r = skewtrnd(v,lambda,varargin)
% Generate random variables from Hansen's Skewed T distribution
% with V degrees of freedom and skewness parameter LAMBDA
%
% USAGE:
%   R=skewtrnd(V,LAMBDA)
%   R=skewtrnd(V,LAMBDA,S1)
%   R=skewtrnd(V,LAMBDA,S1,S2,S3,...,SN)
%   R=skewtrnd(V,LAMBDA,[S1 S2 S3 ... SN])
%
% INPUTS:
%   V      - Degree of freedom parameter.  Either scalar or matrix.
%   LAMBDA - Skewness parameter
%   Sx     - [OPTIONAL] Size of output.
%            R will be:
%               1 by 1 if V is scalar and no S are entered
%               size(V) if V is non-scalar and no S are entered
%               S1 by S1 if V is scalar and only S1 is entered
%               [S1 S2 ... SN] otherwise
%               If V is non-scalar and S are provided, size(V) must equal [S1 S2 ... SN]
%
% OUTPUTS:
%   R     - Skewed T distributed random variables
%
% COMMENTS:
%   V>2
%   -.99<LAMBDA<.99
%   Uses SKEWTINV
%
% REFERENCES:
%   [1] Hansen (1994), Intl.Econ.Rev. (35)
%
% See also SKEWTPDF, SKEWTCDF, SKEWTINV, SKEWTLOGLIK

% Copyright: Andrew Patton
% a.patton@lse.ac.uk
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin < 2
    error('Requires at least two input arguments.');
end

[err, errtext, sizeOut] = iscompatible(2,v,lambda,varargin{:});

if err
    error(errtext)
end
u = rand(sizeOut);
r = skewtinv(u,v,lambda);