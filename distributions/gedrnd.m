function r = gedrnd(v,varargin)
% Generate random variables from a Generalized Error Distribution (GED) with V degrees of freedom
%
% USAGE:
%   R=gedrnd(V)
%   R=gedrnd(V,S1)
%   R=gedrnd(V,S1,S2,S3,...,SN)
%   R=gedrnd(V,[S1 S2 S3 ... SN])
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
%   R     - GED distributed random variables
%
% COMMENTS:
%   V>=1
%   A scalar GED r.v. with variance normalized to 1 has probability
%   density given by:
%       f(x,v) = [v/(lda*2^(1+1/v)*gamma(1/v))]*exp(-0.5*|x/lda|^v)
%       lda = [2^(-2/v)*gamma(1/v)/gamma(3/v)]^0.5
%   If X~GED(v) then abs(X)^v = Y~Gamma(1/v) (cf. Tadikamalla 1980).
%   GAMRND does the computational work.
%
% REFERENCES:
%   [1] Tadikamalla (1980), J.Am.Stat.Assoc. (75)
%   [2] Nelson (1991), Econometrica
%
% See also GEDPDF, GEDCDF, GEDINV, GEDLL, GAMRND

% Copyright: Ivana Komunjer 
% komunjer@hss.caltech.edu
% Modifications Copyright:
% Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if nargin < 1,
    error('Requires at least one input argument.');
end

v(v<1)=NaN;

[err, errtext, sizeOut] = iscompatible(1,v,varargin{:});

if err
    error(errtext)
end

%   Return NaN if the argument V is outside its limit.
r = gamrnd(1./v,1,sizeOut);
r = r.^(1./v);
rndsgn = 2*((rand(sizeOut)>0.5)-0.5);  % generates a random sign -1 or +1
r = r.*rndsgn;

scalex = (gamma(3./v)./gamma(1./v)).^0.5; % standardizes the obtained values
r = r./scalex;