function [trans_parameters,nu,lambda] = figarch_transform(parameters,p,q,errorType)
% FIGARCH(Q,D,P) parameter transformation.  Used to map parameters from a FIGARCH process to the
% real line. Used in the estimation of FIGARCH. 
%
% USAGE:
%   [TRANS_PARAMETERS,NU,LAMBDA]=figarch_transform(PARAMETERS,P,Q,ERROR_TYPE)
%
% INPUTS:
%   PARAMETERS       - Column parameter vector
%   P                - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q                - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERROR_TYPE       - One of:
%                        1 - Gaussian Innovations
%                        2 - T-distributed errors
%                        3 - Generalized Error Distribution
%                        4 - Skewed T distribution
%
% OUTPUTS:
%   TRANS_PARAMETERS - A 2+p+q column vector of transformed parameters corresponding to
%                      [omega phi d beta]'
%   NU               - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA           - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Input parameters must satisfy:
%    (1) omega > 0
%    (2) 0<= d <= 1
%    (3) 0 <= phi <= (1-d)/2 
%    (3) 0 <= beta <= d + phi 
%    (5) nu>2 of Students T and nu>1 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
% See also FIGARCH
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009
 
 
 
%Handle nu, >2.01 in the case of a T, 1.01<nu<50 in the case of a GED
%T uses square, ged uses a logistic
if errorType==2 || errorType==4
    nu=parameters(p+q+3);
    nu=sqrt(nu-2.01);
elseif errorType==3
    nu=parameters(p+q+3);
    temp=(nu-1)/49;
    nu=log(temp/(1-temp));
else
    nu = [];
end
 
%Lambda must be between -.99 and .99.  Use a logistic
if errorType==4
    lambda=parameters(p+q+4);
    temp=(lambda+0.995)/1.99;
    lambda=log(temp/(1-temp));
else
    lambda = [];
end
 
 
% 0<d<1-2*phi
% 0<beta<phi+d
% omega >0
omegaTrans  = log(parameters(1));
% Find d
if p
    d = parameters(3);
else
    d = parameters(2);
end
% 0<d<1
dTrans = log(d/(1-d));
% phi < (1-d)/2
if p
    phi = parameters(2);
    phiTrans = phi/((1-d)/2);
    phiTrans = log(phiTrans/(1-phiTrans));
    phiplusd = phi + d;
else
    phiTrans = [];
    phiplusd = d;
end
if q
    beta = parameters(3+p);
    betaTrans = beta / phiplusd;
    %0 <  beta < phi + d
    betaTrans = log(betaTrans/(1-betaTrans));
else
    betaTrans = [];
end
trans_parameters= [omegaTrans;phiTrans;dTrans;betaTrans];