function z = r2z(R)
% Transformation of a correlation matrix to an unconstrained K(K-1)/2 vector 
%
% USAGE:
%  [R] = r2z(Z)
%
% INPUTS:
%   R - A K by K correlation matrix
%
% OUTPUTS:
%   Z - K(K-1)/2 vector of values in (-inf,inf)
%
% COMMENTS:
%   See z2r for information about the transformation from Z to R.
%
% See also Z2R, R2PHI, PHI2R

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012


k = length(R);
C = zeros(k);
C2 = chol(R)';

for i=2:k
    rem = 1;
    for j=i-1:-1:1
        C(i,j) = C2(i,j)/sqrt(rem);
        rem = rem - C2(i,j)^2;
    end
end

C =C';
z = C(~tril(true(k)));
z = log((z+1)./(1-z));