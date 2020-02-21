function phi = r2phi(R)
% Transformation of a correlation matrix to a set of angles in [0,pi]
%
% USAGE:
%  [PHI] = r2phi(R)
%
% INPUTS:
%   R   - A K by K correlation matrix
%
% OUTPUTS:
%   PHI - K(K-1)/2 vector of values in [0,2*pi]
%
% COMMENTS:
%   See phi2r for information about the transformation from Z to R.
%
% See also PHI2R, Z2R, R2Z


% FIXME : It is necessary to invert both cos and sin to identify where in 0, 2pi the angle is
X = chol(R);
k = length(R);
S = zeros(k);

P = zeros(k);
S(1,:) = 1;
cumS = S;
for i=1:k-1
    P(i,i+1:k) = acos(X(i,i+1:k)./cumS(i,i+1:k));
    S(i+1,i+1:k) = sin(P(i,i+1:k));
    cumS(i+1,i+1:k) = prod(S(1:i+1,i+1:k));
end
phi = P(P>0);