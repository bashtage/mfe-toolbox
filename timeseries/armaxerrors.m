function errors = armaxerrors(parameters,p,q,constant,y,x,m,sigma)
% PURPOSE:
%   Compute errors from an ARMAX model for use in a least squares optimizer
%
% USAGE:
%   [ERRORS] = armaxfilter_likelihood(PARAMETERS,P,Q,CONSTANT,Y,X,M,SIGMA)
%
% INPUTS:
%   PARAMETERS - A vector of GARCH process aprams of the form [constant, arch, garch]
%   P          - Vector containing lag indices of the AR component
%   Q          - Vector containing lag indices of the MA component
%   CONSTANT   - Value indicating whether the model contains a constant (1) or not (0)
%   Y          - Data augments with max(max(Q)-max(P),0) zeros
%   X          - Regressors augmented with max(max(Q)-max(P),0) zeros
%   M          - Index to first element to use in the recursive residual calculation
%   SIGMA      - Vector of conditional standard deviations with the same dimension as Y for use in
%                  GLS estimation
%
% OUTPUTS:
%   ERRORS     - Vector of errors with the same size as Y.  First M elements are 0.
%
% COMMENTS:
%
%  See also ARMAXFILTER_LIKELIHOOD

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 10/19/2009


np = length(p);
nq = length(q);
T = length(y);
errors = zeros(T,1);
k = size(x,2);

for t=m+1:T
    errors(t) = y(t);
    if constant
        errors(t) = errors(t) - parameters(1);
    end
    for i=1:np
        errors(t) = errors(t) - parameters(constant+i)*y(t-p(i));
    end
    for i=1:k
        errors(t) = errors(t) - parameters(constant+np+i)*x(t,i);
    end
    for i=1:nq
        errors(t) = errors(t) - parameters(constant+np+k+i)*errors(t-q(i));
    end
    errors(t) = errors(t);
end
errors = errors./sigma;
