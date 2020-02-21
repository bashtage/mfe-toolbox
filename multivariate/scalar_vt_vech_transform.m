function tparameters=scalar_vt_vech_transform(parameters,p,o,q,kappa)
% SCALAR_VT_VECH(P,Q) parameter transformation.  Used to map parameters
% from a scalar MVGARCH process to the real line. Used in the estimation of SCALAR_VT_VECH.
%
% USAGE:
%   [TPARAMETERS]=scalar_vt_vech_transform(PARAMETERS,P,O,Q,KAPPA)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of symmetric innovations
%   Q                - Non-negative, scalar integer representing the number of lags of conditional variance 
%
% OUTPUTS:
%   TPARAMETERS      - A 1+p+q column vector of transformed parameters corresponding to
%                      [alpha(1),...,alpha(p), beta1 ... beta(q)]'
%
% COMMENTS:
%   Input parameters must satisfy:
%    (1) alpha(i) >= 0 for i = 1,2,...,p
%    (2) beta(i)  >= 0 for i = 1,2,...,q
%    (3) sum(alpha) + sum(beta) < 1
%
% See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


if size(parameters,2)>size(parameters,1)
   parameters = parameters';
end
%Upper bound to keep it a bit away from 1
UB=.999998;

alpha=parameters(1:p);
gamma = parameters(p+1:p+o);
beta=parameters(p+o+1:p+o+q);
%Check that the parameters satisfy the necessary constraints
if  any(alpha<0) || any(beta<0) || any(gamma<0) || (sum(alpha)+sum(gamma)/kappa+sum(beta))>=UB
    error('These do not conform to the necessary set of restrictions to be transformed.')
end

%Alpha, beta cannot be exactly zero or there will be problems with log()
alpha(alpha<1e-8)=1e-8;
gamma(gamma<1e-8)=1e-8;
beta(beta<1e-8)=1e-8;
%Up the upper bound a small amount to make sure it is satisfied
UB=UB+1e-8*(p+o+q);

%Set the scale
scale=UB;
%Initialize the transformed parameters
parameters=[alpha;gamma;beta];
tparameters=[alpha;gamma;beta];
for i=1:(p)
    %Scale the parameters
    tparameters(i)=tparameters(i)/scale;
    %Use an inverse logistic
    tparameters(i)=log(tparameters(i)/(1-tparameters(i)));
    %Update the scale
    scale=scale-parameters(i);
end
for i=p+1:p+o
    %Scale the parameters
    tparameters(i)=tparameters(i)/scale/kappa;
    %Use an inverse logistic
    tparameters(i)=log(tparameters(i)/(1-tparameters(i)));
    %Update the scale
    scale=scale-parameters(i)/kappa;
end

for i=p+o+1:p+o+q
    %Scale the parameters
    tparameters(i)=tparameters(i)/scale;
    %Use an inverse logistic
    tparameters(i)=log(tparameters(i)/(1-tparameters(i)));
    %Update the scale
    scale=scale-parameters(i);
end