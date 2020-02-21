function parameters=egarch_transform(parameters,p,o,q,error_type)
% EGARCH(P,O,Q) parameter transformation.  Used to map parameters from a EGARCH
% process to the real line. Used in the estimation of EGARCH.
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=egarch_transform(PARAMETERS,P,O,Q,ERROR_TYPE);
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
%   TRANS_PARAMETERS - A 1+p+o+q column vector of parameters corresponding to
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o), beta1 ... beta(q)]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
%    (1) nu>2 of Students T and 1<nu<49 for GED
%    (2) -.99<lambda<.99 for Skewed T
%    (3) Other parameters are not transformed
%   
% See also EGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


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
