function [LLF, likelihoods, errors] = armaxfilter_likelihood(parameters,p,q,constant,y,x,m,sigma)
% PURPOSE:
%   Likelihood function for armaxfilter
%
% USAGE:
%   [LLF, LIKELIHOODS, ERRORS] = armaxfilter_likelihood(PARAMETERS,P,Q,CONSTANT,Y,X,M)
%   [LLF, LIKELIHOODS, ERRORS] = armaxfilter_likelihood(PARAMETERS,P,Q,CONSTANT,Y,X,M,SIGMA)
%
% INPUTS:
%   PARAMETERS  - A vector of ARMAX process params of the form [constant, AR, eXog, MA]
%   P           - Vector containing lag indices of the AR component
%   Q           - Vector containing lag indices of the MA component
%   CONSTANT    - Value indicating whether the model contains a constant (1) or not (0)
%   Y           - Dependent variable
%   X           - Regressors, excluding the constant
%   M           - Index to first element to use in the recursive residual calculation
%   SIGMA       - T by 1 vector of conditional standard deviations
%
% OUTPUTS:
%   LLF         - Minus 1 times the log likelihood
%   LIKELIHOODS - Time series of likelihoods
%   ERRORS      - Time series of model errors
%
% COMMENTS:
%
%  See also ARMAXERRORS

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 4/1/2004

e = armaxerrors(parameters,p,q,constant,y,x,m,ones(size(y)));
T = length(e);

if isempty(p)
    p = 0;
end
if isempty(q)
    q = 0;
end
if max(q)>max(p) % prepend y and x, if needed
    t = (max(q)-max(p))+1:T;
else
    t    = 1:T;
end

stde = e./sigma;
sigma2 = stde(t)'*stde(t)/length(t);

% Do not divide e by sigma since this is done in armaxerrors
likelihoods =  0.5*(2*log(sigma) + log(sigma2) + stde.^2./sigma2 + log(2*pi));

likelihoods = likelihoods(t);
errors=e(t);
LLF = sum(likelihoods);

if isnan(LLF)
    LLF=1e7;
end
