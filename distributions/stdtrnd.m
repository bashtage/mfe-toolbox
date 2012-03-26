function r=stdtrnd(v,varargin)
% Generate random variables from a Standardized T distribution with V degrees of freedom
%
% USAGE:
%   R=stdtrnd(V)
%   R=stdtrnd(V,S1)
%   R=stdtrnd(V,S1,S2,S3,...,SN)
%   R=stdtrnd(V,[S1 S2 S3 ... SN])
%
% INPUTS:
%   V     - Degree of freedom parameter.  Either scalar or matrix.
%   Sx    - [OPTIONAL] Size of output.
%           R will be:
%               1 by 1 if V is scalar and no S are entered
%               size(V) if V is non-scalar and no S are entered
%               S1 by S1 if V is scalar and only S1 is entered
%               [S1 S2 ... SN] otherwise
%               If V is non-scalar and S are provided, size(V) must equal [S1 S2 ... SN]
%
% OUTPUTS:
%   R     - Standardized T distributed random variables
%
% COMMENTS:
%   V>2
%   Uses TRND
%
% REFERENCES:
%   [1] Cassella and Berger (1990) 'Statistical Interence'
%
% See also STDTPDF, STDTCDF, STDTINV, STDTLL, TINV

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%v must either be scalar or the same size as that in varargin if present

if nargin<1
    error('At least one input required.')
end

[err, errtext, sizeOut, v] = iscompatible(1,v,varargin{:});

if err
    error(errtext)
end

r=trnd(v,sizeOut);
stdev=sqrt(v./(v-2));
stdev(v<=2)=NaN;
r=r./stdev;