function [trans_parameters,nu,lambda]=aparch_transform(parameters,p,o,q,errorType,deltaIsEstimated)
% APARCH(P,O,Q) parameter transformation.  Used to map parameters from a APARCH
% process to the real line. Used in the estimation of APARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=aparch_transform(PARAMETERS,P,O,Q,ERRORTYPE)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - Positive, scalar integer representing the number of
%                      symmetric innovations
%   O                - Non-negative scalar integer representing the number
%                      of asymmetric innovations (0 for symmetric processes)
%   Q                - Non-negative, scalar integer representing the number
%                      of lags of conditional variance (0 for ARCH)
%   ERRORTYPE       - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%
% OUTPUTS:
%   TRANS_PARAMETERS - A 1+p+o+q+1 column vector of transformed parameters corresponding to
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o), beta1 ... beta(q) delta]'
%                       delta is only included if no user provided delta is supplied        
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

%Constraints are omega>0
%nu>2.01 or 1.01
%-.99<lambda<.99
% alpha>0
% gamma>-alpha(i)
% beta>0
% sum(alpha)+0.5sum(gamma)+sum(beta)<1

nu=[];
lambda=[];
%Handle nu, >2.01 in the case of a T, 1.01<nu<50 in the case of a GED
%T uses square, ged uses a logistic
if errorType==2 || errorType==4
    nu=parameters(p+o+q+2+deltaIsEstimated);
    nu=sqrt(nu-2.01);
    parameters(p+o+q+2+deltaIsEstimated)=nu;
elseif errorType==3
    nu=parameters(p+o+q+2+deltaIsEstimated);
    temp=(nu-1)/49;
    nu=log(temp/(1-temp));
    parameters(p+o+q+2+deltaIsEstimated)=nu;
end

%Lambda must be between -.99 and .99.  Use a logistic
if errorType==4
    lambda=parameters(p+o+q+3+deltaIsEstimated);    
    temp=(lambda+0.995)/1.99;
    lambda=log(temp/(1-temp));
    parameters(p+o+q+3+deltaIsEstimated)=lambda;
end


%Use log transform for omega
omega=parameters(1);
tomega=log(omega);

%Parse the parameters
alpha=parameters(2:p+1);
gamma=parameters(p+2:p+o+1);
beta=parameters(p+o+2:p+o+q+1);
if deltaIsEstimated
    delta=parameters(p+o+q+2);
end
%Upper bound to keep it a bit away from 1
UB=.999998;

%Check that the parameters satisfy the necessary constraints
if  any(alpha<0) || any(beta<0) || any(gamma<-1) ||any(gamma>1) || (sum(alpha)+sum(beta))>=UB
    error('These do not conform to the necessary set of restrictions to be transformed.')
end

%Finally, must be certain that none of the parameters are exactly zero, and
%that the alpha2+gamma2>0
alpha(alpha==0)=1e-8;
beta(beta==0)=1e-8;

%Finally, up the upper bound a tiny bit
UB=UB+1e-8*(p+o+q);

%Set the scale
scale=UB;
%Initialze the transformed alpha
talpha=alpha;
for i=1:p
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


%Transform gamma
tgamma=(gamma+1)/2;
tgamma=log(tgamma./(1-tgamma));

%Transform delta if needed
if deltaIsEstimated
    tdelta=(delta-.3)/3.7;
    tdelta=log(tdelta./(1-tdelta));
else
    tdelta =[];
end


%Regroup the parameters
trans_parameters=[tomega;talpha;tgamma;tbeta;tdelta];
if any(~isreal(trans_parameters))
    keyboard
end