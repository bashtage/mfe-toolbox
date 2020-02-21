function [trans_parameters,nu,lambda]=aparch_itransform(parameters,p,o,q,errorType,deltaIsEstimated)
% APARCH(P,O,Q) inverse parameter transformation.  Used to map parameters from
% the real line to a set of parameters appropriate for a APARCH model.
% Used in the estimation of APARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=aparch_itransform(PARAMETERS,P,O,Q,ERRORTYPE,DELTAISESTIMATED)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   O                - Non-negative scalar integer representing the number
%                      of asymmetric innovations (0 for symmetric processes)
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   ERRORTYPE        - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%   DELTAISESTIMATED - 1 or 0 to indicate whether a user supplied value of
%                        delta has been provided (0) or delta is jointly
%                        estimated (1)
%
% OUTPUTS:
%   TRANS_PARAMETERS - A 1+p+o+q+DELTAISESTIMATED column vector of parameters with
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o),beta1 ... beta(q) delta]'
%                      delta is only included if no user provided value is supplied
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) -1<gamma(i)<1  for i=1,...,o
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) sum(alpha(i) + beta(k)) < 1 for i = 1,2,...p and k=1,2,...,q
%    (5) 0.3<delta<4
%    (6) nu>2 of Students T and 1<nu<49 for GED
%    (7) -.99<lambda<.99 for Skewed T
%
% See also APARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

nu=[];
lambda=[];

%Upper constraint to make sure that there is no overflow
parameters(parameters>100)=100;

%Parse the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
gamma=parameters(p+2:p+o+1);
beta=parameters(p+o+2:p+o+q+1);
if deltaIsEstimated
    delta=parameters(p+o+q+2);
end

%Handle the transformation of nu and lambda
if errorType==2 || errorType==4
    %Square here
    nu=parameters(p+o+q+2+deltaIsEstimated);
    nu=2.01+nu^2;
elseif errorType==3
    %Logistic here, have to be careful about overflow
    nu=parameters(p+o+q+2+deltaIsEstimated);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
end
%If skewt, use a logistic to map to -.99,.99
if errorType==4
    lambda=parameters(p+o+q+3+deltaIsEstimated);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
end

%Upper bound of transform
UB=.9998;

%Simple transform of omega
tomega=exp(omega);

%Initialize the transformed parameters
talpha=alpha;
tbeta=beta;
tgamma = 1.999*(exp(gamma)./(1+exp(gamma)))-.9995;
if deltaIsEstimated
    tdelta = 0.3+3.7*exp(delta)/(1+exp(delta));
else
    tdelta=[];
end
%Set the scale
scale=UB;

for i=1:p
    %Alpha is between 0 and scale
    talpha(i)=(exp(talpha(i))/(1+exp(talpha(i))))*scale;
    %Update the scale
    scale=scale-talpha(i);
end
for i=1:q
    %Beta is between 0 and scale
    tbeta(i)=(exp(tbeta(i))/(1+exp(tbeta(i))))*scale;
    %Update the scale
    scale=scale-tbeta(i);
end
%Regroup the transformed parameters.
trans_parameters=[tomega;talpha;tgamma;tbeta;tdelta];
if any(~isreal(trans_parameters))
    keyboard
end
