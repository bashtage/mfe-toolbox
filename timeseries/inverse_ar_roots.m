function [rho, stationary]=inverse_ar_roots(phi)
% Computes the inverted roots of the characteristic equation of an AR(P)
%
% USAGE:
%  [RHO, STATIONARY] = inverse_ar_roots(PARAMETERS)
% 
% INPUTS:
%   PARAMETERS - A column vector containing the AR parameters in the order
%                  y(t) = mu + phi(1) y(t-1) + phi(2)y(t-2) + ... + phi(P) y(t-P) + e(t)
%
% OUTPUTS:
%   RHO        - A max(P) by 1 vector containing the inverted roots of the characteristic equation
%                  corresponding to the AR model input 
%   STATIONARY - Logical indicating whether the AR process is stationary
% 
% COMMENTS:
% Process is stationary if min(abs(rho)) > 1
% 
%  See also ARMAROOTS, ACF
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007

%Make sure phi is a row vector;

p=length(phi);
%Make sure phi is a row vector
if size(phi,1)>=size(phi,2) && min(size(phi))==1
    phi=phi';
else
    error('Phi should be a column vector.')
end
rho=roots([-fliplr(phi) 1]);
stationary=min(abs(rho))>1;