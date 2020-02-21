function [parameters,nu,lambda] = figarch_itransform(parameters,p,q,errorType)
% FIGARCH(Q,D,P) parameter inverse transformation.  Used to map parameters from the real line to
% values which satisfy non-negativity constraints. Used in the estimation of FIGARCH. 
%
% USAGE:
%   [PARAMETERS,NU,LAMBDA]=figarch_itransform(PARAMETERS,P,Q,ERROR_TYPE)
%
% INPUTS:
%   PARAMETERS    - Column parameter vector
%   P             - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q             - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERROR_TYPE    - One of:
%                     1 - Gaussian Innovations
%                     2 - T-distributed errors
%                     3 - Generalized Error Distribution
%                     4 - Skewed T distribution
%
% OUTPUTS:
%   PARAMETERS    - A 2+p+q column vector of parameters corresponding to
%                     [omega phi d beta]'
%   NU            - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA        - Distribution asymmetry parameter, empty if not applicable
%
% COMMENTS:
%   Output parameters will satisfy:
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
 
 
 
 
%Handle the transformation of nu and lambda
if errorType==2 || errorType==4
    %Square here
    nu=parameters(p+q+3);
    nu=2.01+nu^2;
elseif errorType==3
    %Logistic here, have to be careful about overflow
    nu=parameters(p+q+3);
    if nu>100
        nu=100;
    end
    nu=exp(nu)/(1+exp(nu));
    nu=49*nu+1.01;
else
    nu = [];
end
%If skewt, use a logistic to map to -.99,.99
if errorType==4
    lambda=parameters(p+q+4);
    lambda=exp(lambda)/(1+exp(lambda));
    lambda=1.98*lambda-.99;
else
    lambda = [];
end
 
 
 
 
 
% 0<d<1-2*phi
% 0<beta<phi+d
% omega >0
omega  = exp(parameters(1));
% Find d
if p
    d = parameters(3);
else
    d = parameters(2);
end
% 0<d<1
d = exp(d)/(1+exp(d));
% phi < (1-d)/2
if p
    phi = parameters(2);
    phi = exp(phi)/(1+exp(phi));
    phi = (1-d)/2 * phi;
    phiplusd = phi + d;
else
    phi = [];
    phiplusd = d;
end
if q
    beta = parameters(3+p);
    %0 <  beta < phi + d
    beta = exp(beta)/(1+exp(beta));
    beta = beta * phiplusd;
else
    beta = [];
end
parameters = [omega;phi;d;beta];
