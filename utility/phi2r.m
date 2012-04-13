function R = phi2r(phi,transform)
% Transformation of a vector of angles in [0,pi] to a correlation matrix
%
% USAGE:
%  [R] = phi2r(PHI)
%
% INPUTS:
%   PHI - K(K-1)/2 vector of values in [0,pi]
%
% OUTPUTS:
%   R - A K by K correlation matrix
%
% COMMENTS:
%
% See also R2PHI, Z2R, R2Z

if nargin==2 && transform
    phi = 2*pi*exp(phi)./(1+exp(phi));
end

m = length(phi);
k = ceil(sqrt(2*m));
Csel = ~tril(ones(k));
Ssel = triu(true(k));
Ssel(1,:) = false;
C = eye(k);
S = zeros(k);
S(1,:) = 1;
c = cos(phi);
s = sin(phi);
C(Csel)=c;
S(Ssel)=s;

C = C.*cumprod(S);
R = C'*C;
r = sqrt(diag(R));
R = R./(r*r');
R = (R+R')/2;
