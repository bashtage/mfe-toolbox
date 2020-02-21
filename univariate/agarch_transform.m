function [trans_parameters,nu,lambda]=agarch_transform(parameters,p,q,model_type,error_type,transform_bounds)
% AGARCH(P,Q) parameter transformation.  Used to map parameters from a AGARCH
% process to the real line. Used in the estimation of AGARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=agarch_itransform(PARAMETERS,DATA,P,Q,MODEL_TYPE,ERROR_TYPE)
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
%   TRANS_PARAMETERS - A 2+p+q column vector of transformed parameters corresponding to
%                      [omega,alpha(1),...,alpha(p),gamma, beta1 ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
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

%Handle nu, >2.01 in the case of a T, 1.01<nu<50 in the case of a GED
%T uses square, ged uses a logistic
if error_type==2 || error_type==4
    nu=parameters(p+q+3);
    nu=sqrt(nu-2.01);
elseif error_type==3
    nu=parameters(p+q+3);
    temp=(nu-1)/49;
    nu=log(temp/(1-temp));
end

%Lambda must be between -.99 and .99.  Use a logistic
if error_type==4
    lambda=parameters(p+q+4);
    temp=(lambda+0.995)/1.99;
    lambda=log(temp/(1-temp));
end

%Parse the parameters
omega=parameters(1);
alpha=parameters(2:p+1);
beta=parameters(p+3:p+q+2);
gamma = parameters(p+2);

%Use log transform for omega
tomega=log(omega);

% Gamma
if model_type == 1
    LowerQuantile = transform_bounds(1);
    UpperQuantile = transform_bounds(2);    
    tgamma = (gamma - LowerQuantile)/(UpperQuantile-LowerQuantile);
else
    tgamma = (gamma + sqrt(10))/(2*sqrt(10));
end
tgamma = log(tgamma/(1-tgamma));




%Upper bound to keep it a bit away from 1
UB=.999998;

%Check that the parameters satisfy the necessary constraints
if  any(alpha<0) || any(beta<0) || (sum(alpha)*(1+gamma^2)+sum(beta))>=UB
    error('These do not conform to the necessary set of restrictions to be transformed.')
end

%Finally, must be certain that none of the parameters are exactly zero, and
%that the alpha2+gamma2>0
alpha(alpha==0)=1e-8;
beta(beta==0)=1e-8;

%Set the scale
scale=UB;
%Initialze the transformed alpha
% Note that alpha is enlarged here since it is shrunk in the inverse
% transform
talpha=alpha * (1+gamma^2);
for i=1:p
    % Note that 
    %Scale the alpha
    talpha(i)=alpha(i)./scale;
    %Use an inverse logistic
    talpha(i)=log(talpha(i)./(1-talpha(i)));
    %Update the scale
    scale=scale-alpha(i);
end

%Initialize the beta
tbeta=beta;
%Iterate over betas
for i=1:q
    %Scale the betas
    tbeta(i)=tbeta(i)./scale;
    %Use an inverse logistic
    tbeta(i)=log(tbeta(i)./(1-tbeta(i)));
    %Update the scale
    scale=scale-beta(i);
end

%Regroup the parameters
trans_parameters=[tomega;talpha;tgamma;tbeta];
