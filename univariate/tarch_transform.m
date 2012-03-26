function [trans_parameters,nu,lambda]=tarch_transform(parameters,p,o,q,error_type)
% TARCH(P,O,Q) parameter transformation.  Used to map parameters from a TARCH
% process to the positive unit simplex. Used in the estimation of TARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=tarch_transform(PARAMETERS,P,O,Q,ERROR_TYPE)
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
%   TRANS_PARAMETERS - A 1+p+o+q column vector of transformed parameters corresponding to
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o), beta1 ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) gamma(i) + alpha(i) > 0 for i=1,...,o
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) sum(alpha(i) + 0.5*gamma(j) + beta(k)) < 1 for i = 1,2,...p 
%        j = 1,2,...o and k=1,2,...,q
%    (5) nu>2 of Students T and 1<nu<49 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
% See also TARCH

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
if error_type==2 || error_type==4
    nu=parameters(p+o+q+2);
    nu=sqrt(nu-2.01);
    parameters(p+o+q+2)=nu;
elseif error_type==3
    nu=parameters(p+o+q+2);
    temp=(nu-1)/49;
    nu=log(temp/(1-temp));
    parameters(p+o+q+2)=nu;
end

%Lambda must be between -.99 and .99.  Use a logistic
if error_type==4
    lambda=parameters(p+o+q+3);
    temp=(lambda+0.995)/1.99;
    lambda=log(temp/(1-temp));
    parameters(p+o+q+3)=lambda;
end


%Use log transform for omega
omega=parameters(1);
tomega=log(omega);

%Parse the parameters
alpha=parameters(2:p+1);
gamma=parameters(p+2:p+o+1);
beta=parameters(p+o+2:p+o+q+1);

%Upper bound to keep it a bit away from 1
UB=.999998;

%Must make sure alpha and gamma sum to >=0 where they both exist
gamma2=[gamma;zeros(max(0,p-o),1)];
alpha2=[alpha;zeros(max(0,o-p),1)];
%Check that the parameters satisfy the necessary constraints
if  any(alpha<0) || any(beta<0) || any(gamma(min(p,o)+1:o)<0) || any((alpha2+gamma2)<0)  || (sum(alpha)+0.5*sum(gamma)+sum(beta))>=UB
    error('These do not conform to the necessary set of restrictions to be transformed.')
end
if p>0 && o>0
    if any( (gamma(1:min(p,o))+alpha(1:min(p,o)))<0 )
        error('These do not conform to the necessary set of restrictions to be transformed.')
    end
end

%Finally, must be certain that none of the parameters are exactly zero, and
%that the alpha2+gamma2>0
alpha(alpha==0)=1e-8;
beta(beta==0)=1e-8;

gamma2=[gamma;zeros(max(0,p-o),1)];
alpha2=[alpha;zeros(max(0,o-p),1)];

pl=find((alpha2+gamma2)==0);
if p>=o
    %If p is larger than o (weakly), just up alpha bit
    alpha(pl)=alpha(pl)+1e-8;
else
    %Otherwise, up the problematic gammas a bit
    gamma(pl)=gamma(pl)+1e-8;
end
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

%Initialize the transformed gamma
tgamma=gamma;
for i=1:o
    %Range of gamma is from -alpha(i) to 2*scale
    tgamma(i)=gamma(i);
    if i<=p %If there is an alpha, then add in the value
        tgamma(i)=tgamma(i)+alpha(i);
        %Now tgamma is between 0 and 2*scale+alpha
        tgamma(i)=tgamma(i)/(2*scale+alpha(i));
    else
        %No alpha, then between 0 and 2*scale
        tgamma(i)=tgamma(i)/(2*scale);
    end
    %Make sure tgamma is not exactly 1 or 0!
    if tgamma(i)==1
        tgamma(i)=1-eps;
    elseif tgamma(i)==0
        tgamma(i)=eps;
    end

    %Tgamma is now between 0 and 1, inverse logistic
    tgamma(i)=log(tgamma(i)/(1-tgamma(i)));
    %Update the scale
    scale=scale-gamma(i)/2;
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
