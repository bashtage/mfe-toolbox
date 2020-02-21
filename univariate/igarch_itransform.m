function [transParameters,nu,lambda]=igarch_itransform(parameters,p,q,errorType,constant)
% IGARCH(P,Q) inverse parameter transformation.  Used to map parameters from
% the real line to a set of parameters appropriate for a IGARCH model.
% Used in the estimation of IGARCH.
%
% USAGE:
%   [TRANSPARAMETERS,NU,LAMBDA]=igarch_itransform(PARAMETERS,P,Q,ERRORTYPE,CONSTANT)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   ERRORTYPE       - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%   CONSTANT         - 1 if model includes a constant, 0 otherwise
%
% OUTPUTS:
%   TRANSPARAMETERS - A CONSTANT+p+q column vector of parameters with
%                      [omega,alpha(1),...,alpha(p) beta(1) ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Output parameters satisfy:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) sum(alpha(i) + beta(j)) = 1 for i = 1,2,...p and j = 1,2,...q
%    (5) nu>2 of Students T and nu>1 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
% See also IGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009

nu=[];
lambda=[];

%Upper constraint to make sure that there is no overflow
parameters(parameters>100)=100;

%Parse the parameters
if constant
    omega=parameters(1);
end
alpha=parameters(constant+1:p+constant);
beta=parameters(p+constant+1:p+q+constant-1);

%Handle the transformation of nu and lambda
if errorType==2 || errorType==4
    %Square here
    nu=parameters(p+q+constant);
    nu=2.01+nu^2;
elseif errorType==3
    %Logistic here, have to be careful about overflow
    nu=parameters(p+q+constant);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
end
%If skewt, use a logistic to map to -.99,.99
if errorType==4
    lambda=parameters(p+q+constant+1);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
end

%Upper bound of transform
UB=.999998;

%Simple transform of omega
if constant
    tomega = exp(omega);
else
    tomega = [];
end

%Initialize the transformed parameters
talpha=alpha;
tbeta=beta;

%Set the scale
scale=UB;

for i=1:p
    %Alpha is between 0 and scale
    talpha(i)=(exp(talpha(i))/(1+exp(talpha(i))))*scale;
    %Update the scale
    scale=scale-talpha(i);
end
if q>1
    for i=1:q-1
        %Beta is between 0 and scale
        tbeta(i)=(exp(tbeta(i))/(1+exp(tbeta(i))))*scale;
        %Update the scale
        scale=scale-tbeta(i);
    end
else
    tbeta = [];
end
%Regroup the transformed parameters.
transParameters=[tomega;talpha;tbeta];
