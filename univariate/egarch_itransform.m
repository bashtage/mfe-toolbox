function parameters=egarch_itransform(parameters,p,o,q,error_type)
% EGARCH(P,O,Q) parameter transformation.  Used to map parameters from a EGARCH
% process to the real line. Used in the estimation of EGARCH.
%
% USAGE:
%   [PARAMETERS]=egarch_itransform(PARAMETERS,P,O,Q,ERROR_TYPE)
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
%   PARAMETERS      - A 1+p+o+q column vector of parameters corresponding to
%                      [omega,alpha(1),...,alpha(p),gamma(1),...,gamma(o), beta1 ... beta(q) [nu lambda]]'
%
% COMMENTS:
%   Output parameters will satisfy:
%    (1) nu>2 of Students T and 1<nu<49 for GED
%    (2) -.99<lambda<.99 for Skewed T
%    (3) Other parameters are not transformed
%   
% See also EGARCH

%If t nushoudl be > 2.01
%if GED, 1.01<nu<50
%if Skew-t, nu>2.01, -.99<lambda<.99

if error_type==2 || error_type==4
    nu=parameters(p+o+q+2);
    nu=2.01+nu^2;
    parameters(p+o+q+2)=nu;
elseif error_type==3
    nu=parameters(p+o+q+2);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
    parameters(p+o+q+2)=nu;
end

if error_type==4
    lambda=parameters(p+o+q+3);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
    parameters(p+o+q+3)=lambda;
end
