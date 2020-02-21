function lambda = figarch_weights(parameters,p,q,truncLag)
% FIGARCH(Q,D,P) ARCH(oo) weight computation forfor P={0,1} and Q={0,1}
%
% USAGE:
%   [LAMBDA] = figarch_weights(PARAMETERS,P,Q)
%
% INPUTS:
%   PARAMETERS    - A (1+P+Q) by 1 column vector of parameters in the order [phi d beta] where phi
%                     is omitted if P=0 and beta is omitted if Q=0 
%   P             - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q             - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   TRUNCLAG      - Number of weights to compute in ARCH(oo) representation
%
% OUTPUTS:
%   PARAMETERS    - TRUNCLAG by 1 vector or weights for FIGARCH calculation
%
%    The conditional variance, h(t), of a FIGARCH(1,d,1) process is modeled
%    as follows:
%
%    h(t) = omega + [1-beta L - phi L (1-L)^d] epsilon^2(t) + beta * h(t-1)
%    
%    which is estimated using an ARCH(oo) representation, 
%
%    h(t) = omega + sum(lambda(i) * epsilon(t-1))
%    
%    where lambda(i) is a function of the fractional differencing parameter, phi and beta
%    
%  See also FIGARCH, FIGARCH_LIKELIHOOD, FIGARCH_PARAMETER_CHECK,
%  FIGARCH_STARTING_VALUES, FIGARCH_TRANSFORM, FIGARCH_ITRANSFORM 
%

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/13/2009


% Parse parameters
if q
    beta = parameters(2+p);
else
    beta = 0;
end

if p
    phi = parameters(1);
    d = parameters(2);
else
    phi = 0;
    d = parameters(1);
end


% Recursive weight computation
lambda = zeros(truncLag,1);
delta = zeros(truncLag,1);
lambda(1) = phi - beta + d;
delta(1) = d;
for i=2:truncLag
    delta(i) = (i-1-d)/i * delta(i-1);
    lambda(i) = beta*lambda(i-1) + (delta(i) -phi *delta(i-1));
end