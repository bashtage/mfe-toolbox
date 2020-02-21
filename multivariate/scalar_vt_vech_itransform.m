function parameters=scalar_vt_vech_itransform(tparameters,p,o,q,kappa)
% SCALAR_VT_VECH(P,Q) inverse parameter transformation.  Used to map parameters
% from the real line to those appropriate for a scalar MVGARCH process. 
% Used in the estimation of SCALAR_VT_VECH.
%
% USAGE:
%   [PARAMETERS]=scalar_vt_vech_itransform(TPARAMETERS,P,O,Q,KAPPA)
%
% INPUTS:
%   TPARAMETERS      - Column vector of transformed parameters (-inf,inf)
%   P                - Positive, scalar integer representing the number of symmetric innovations
%   Q                - Non-negative, scalar integer representing the number of lags of conditional variance 
%
% OUTPUTS:
%   PARAMETERS       - A 1+p+q column vector of parameters corresponding to
%                      [alpha(1),...,alpha(p), beta1 ... beta(q)]'
%
% COMMENTS:
%   Output parameters will satisfy:
%    (1) alpha(i) >= 0 for i = 1,2,...,p
%    (2) beta(i)  >= 0 for i = 1,2,...,q
%    (3) sum(alpha) + sum(beta) < 1
%
% See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%Upper bound to keep it a bit away from 1
UB=.999998;
%Make sure there are not inf causing values
tparameters(tparameters>50)=50;
%Map them to [0,1]
parameters=exp(tparameters)./(1+exp(tparameters));
%Normalize them, starting with the second
scale=UB;
for i=1:p 
    parameters(i)=parameters(i)*scale;
    scale=scale-parameters(i);
end
for i=p+1:p+o
    parameters(i)=parameters(i)*scale*kappa;
    scale=scale-parameters(i)/kappa;
end
for i=p+o+1:p+o+q
     parameters(i)=parameters(i)*scale;
    scale=scale-parameters(i);
end
