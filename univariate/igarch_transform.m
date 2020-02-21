function [transParameters,nu,lambda]=igarch_transform(parameters,p,q,errorType,constant)
% IGARCH(P,Q) parameter transformation.  Used to map parameters from a IGARCH
% process to the positive unit simplex. Used in the estimation of IGARCH.
%
% USAGE:
%   [TRANSPARAMETERS,NU,LAMBDA]=igarch_transform(PARAMETERS,P,Q,ERRORTYPE,CONSTANT)
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
%   TRANSPARAMETERS - A CONSTANT+p+q-1 column vector of transformed parameters corresponding to
%                      [omega,alpha(1),...,alpha(p) beta1 ... beta(q-1)]'
%                      where the final beta has been excluded
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
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
    nu=parameters(p+q+constant);
    nu=sqrt(nu-2.01);
elseif errorType==3
    nu=parameters(p+q+constant);
    temp=(nu-1)/49;
    nu=log(temp/(1-temp));
end

%Lambda must be between -.99 and .99.  Use a logistic
if errorType==4
    lambda=parameters(p+q+constant+1);
    temp=(lambda+0.995)/1.99;
    lambda=log(temp/(1-temp));
end


%Use log transform for omega
if constant
    omega=parameters(1);
    tomega=log(omega);
else
    tomega = [];
end

%Parse the parameters
alpha=parameters(constant+1:p+constant);
beta=parameters(constant+p+1:constant+p+q-1);

%Upper bound to keep it a bit away from 1
UB=.999998;

%Check that the parameters satisfy the necessary constraints
if isempty(beta)
    sumbeta = 0;
else
    sumbeta = sum(beta);
end
if  any(alpha<0) || any(beta<0) || (sum(alpha)+sumbeta)>=UB
    error('These do not conform to the necessary set of restrictions to be transformed.')
end

%Finally, must be certain that none of the parameters are exactly zero, and
%that the alpha2+gamma2>0
alpha(alpha==0)=1e-8;
beta(beta==0)=1e-8;

%Finally, up the upper bound a tiny bit
UB=UB+1e-8*(p+q);

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

if q>1
    %Initialize the beta
    tbeta=beta;
    %Iterate over betas
    for i=1:(q-1)
        %Scale the betas
        tbeta(i)=tbeta(i)./scale;
        %Use an inverse logistic
        tbeta(i)=log(tbeta(i)./(1-tbeta(i)));
        %Update the scale
        scale=scale-beta(i);
    end
else
    tbeta = [];
end

%Regroup the parameters
transParameters=[tomega;talpha;tbeta];
