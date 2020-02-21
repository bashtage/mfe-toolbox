function R = z2r(z)
% Transformation of an unconstrained K(K-1)/2 vector to a correlation matrix
%
% USAGE:
%  [R] = z2r(Z)
%
% INPUTS:
%   Z - K(K-1)/2 vector of values in (-inf,inf)
%
% OUTPUTS:
%   R - A K by K correlation matrix
%
% COMMENTS:
%   The unconstrained values are mapped to the correlaiotn matrix through:
%   y = 2*exp(Z)./(1+exp(Z))-1
%   C(i,j) = y(i,j) * sqrt(1-sum(C(i,j+1:k).^2), for i=2,...,k, j=i-1:-1:1
%   R = C*C'
%
% See also R2Z, R2PHI, PHI2R

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

m = length(z);
k = ceil(sqrt(2*m));
if k*(k-1)/2~=m
    error('Incorrect number of elements in z')
end
C = zeros(k);
z = (exp(z)-1)./(1+exp(z));
count = 1;
for i=2:k
    for j=1:i-1;
        C(i,j) = z(count);
        count = count + 1;
    end
end
C(1,1) = 1;
for i=2:k
    rem = 1;
    for j=i-1:-1:1
        C(i,j) = C(i,j)*sqrt(rem);
        rem = rem - C(i,j)^2;
    end
    C(i,i) = sqrt(rem);
end

R = C*C';
r = sqrt(diag(R));
R = R ./ (r*r');