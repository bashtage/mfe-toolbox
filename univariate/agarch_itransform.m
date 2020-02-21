function [trans_parameters,nu,lambda]=agarch_itransform(parameters,p,q,model_type,error_type,transform_bounds)
% AGARCH(P,Q) inverse parameter transformation.  Used to map parameters from
% the real line to a set of parameters appropriate for a AGRCH model.
% Used in the estimation of AGARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=agarch_itransform(PARAMETERS,P,Q,MODEL_TYPE,ERROR_TYPE,TRANS_BOUNDS)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   MODEL_TYPE       - The type of variance process, either
%                        1 - AGARCH
%                        2 - NAGARCH
%   ERROR_TYPE       - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%   TRANS_BOUNDS     - 2 by 1 vector containing the .01 and .99 quantiles
%                        of EPSILON for use in parameter transformation
%
% OUTPUTS:
%   TRANS_PARAMETERS - A 2+p+q column vector of parameters with
%                      [omega,alpha(1),...,alpha(p),gamma,beta1 ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Output parameters satisfy:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) -q(.01,EPSILON)<gamma<q(.99,EPSILON) for AGARCH 
%    (5) sum(alpha(i) + beta(k)) < 1 for i = 1,2,...p and k=1,2,...,q for
%    AGARCH and sum(alpha(i)*(1+gamma^2) + beta(k)) < 1 for NAGARCH
%    (6) nu>2 of Students T and nu>1 for GED
%    (7) -.99<lambda<.99 for Skewed T
%
% See also AGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


nu=[];
lambda=[];

%Upper constraint to make sure that there is no overflow
parameters(parameters>100)=100;

%Parse the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
gamma=parameters(p+2);
beta=parameters(p+3:p+q+2);

%Handle the transformation of nu and lambda
if error_type==2 || error_type==4
    %Square here
    nu=parameters(p+q+3);
    nu=2.01+nu^2;
elseif error_type==3
    %Logistic here, have to be careful about overflow
    nu=parameters(p+q+3);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
end
%If skewt, use a logistic to map to -.99,.99
if error_type==4
    lambda=parameters(p+q+4);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
end

%Simple transform of omega
tomega=exp(omega);
% gamma is also relatively simple if AGARCH
tgamma=exp(gamma)/(1+exp(gamma));
if model_type == 1
    LowerQuantile = transform_bounds(1);
    UpperQuantile = transform_bounds(2);
    tgamma=tgamma * (UpperQuantile-LowerQuantile) + LowerQuantile;
else
    tgamma=(tgamma - 0.5)*sqrt(10);
end


%Initialize the transformed parameters
talpha=alpha;
tbeta=beta;

%Upper bound of transform
UB=.9998;
scale = UB;
%Set the scale, note that scale depends on gamma if NAGARCH
for i=1:p
    %Alpha is between 0 and scale
    talpha(i)=(exp(talpha(i))/(1+exp(talpha(i))))*scale;
    %Update the scale
    scale=scale-talpha(i);
end
if model_type == 2
    talpha = talpha / (1+tgamma^2);
end
for i=1:q
    %Beta is between 0 and scale
    tbeta(i)=(exp(tbeta(i))/(1+exp(tbeta(i))))*scale;
    %Update the scale
    scale=scale-tbeta(i);
end
%Regroup the transformed parameters.

trans_parameters=[tomega;talpha;tgamma;tbeta];