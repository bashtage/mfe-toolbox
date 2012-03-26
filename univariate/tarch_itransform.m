function [trans_parameters,nu,lambda]=tarch_itransform(parameters,p,o,q,error_type)
% TARCH(P,O,Q) inverse parameter transformation.  Used to map parameters from
% the real line to a set of parameters appropriate for a TARCH model.
% Used in the estimation of TARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=tarch_itransform(parameters,p,o,q,error_type)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   O                - Non-negative scalar integer representing the number
%                      of asymmetric innovations (0 for symmetric processes)
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   ERROR_TYPE       - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%
% OUTPUTS:
%   TRANS_PARAMETERS - A 1+p+o+q column vector of parameters with
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o),beta1 ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Output parameters satisfy:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) gamma(i) + alpha(i) > 0 for i=1,...,o
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) sum(alpha(i) + 0.5*gamma(j) + beta(k)) < 1 for i = 1,2,...p and
%    j = 1,2,...o, k=1,2,...,q
%    (5) nu>2 of Students T and nu>1 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
% See also TARCH

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

%Handle the transformation of nu and lambda
if error_type==2 || error_type==4
    %Square here
    nu=parameters(p+o+q+2);
    nu=2.01+nu^2;
elseif error_type==3
    %Logistic here, have to be careful about overflow
    nu=parameters(p+o+q+2);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
end
%If skewt, use a logistic to map to -.99,.99
if error_type==4
    lambda=parameters(p+o+q+3);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
end

%Upper bound of transform
UB=.9998;

%Simple transform of omega
tomega=exp(omega);

%Initialize the transformed parameters
talpha=alpha;
tgamma=gamma;
tbeta=beta;

%Set the scale
scale=UB;

for i=1:p
    %Alpha is between 0 and scale
    talpha(i)=(exp(talpha(i))/(1+exp(talpha(i))))*scale;
    %Update the scale
    scale=scale-talpha(i);
end
for i=1:o;
    %Gamma is between -Alpha(i) and 2*scale
    %first map it into 0,1
    tgamma(i)=(exp(tgamma(i))/(1+exp(tgamma(i))));
    if p>=i %alpha exists
        %Then map it into 0, 2*scale+alpha(i)
        tgamma(i)=tgamma(i)*(2*scale+talpha(i));
        %Finally map it into -alpha(i),2*scale
        tgamma(i)=tgamma(i)-talpha(i);
    else
        %No alpha, map it into 0, 2*scale
        tgamma(i)=tgamma(i)*(2*scale);
    end
    %Update the scale
    scale=scale-0.5*tgamma(i);
end
for i=1:q
    %Beta is between 0 and scale
    tbeta(i)=(exp(tbeta(i))/(1+exp(tbeta(i))))*scale;
    %Update the scale
    scale=scale-tbeta(i);
end
%Regroup the transformed parameters.
trans_parameters=[tomega;talpha;tgamma;tbeta];