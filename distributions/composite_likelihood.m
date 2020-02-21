function ll = composite_likelihood(S,data,indices)
% Computes the negative of the composite normal loglikelihood (bivariate) for a K-dimensional array
%
% USAGE:
%   [LL] = composite_likelihood(S,DATA,INDICES)
%
% INPUTS:
%   S       - K by K covariance matrix
%   DATA    - Either a K by K matrix (e.g. outer-produces) or a 1 by K vector (e.g. returns)
%   INDICES - Q by 2 array of indices to use when computing the composite likleihood
%
% OUTPUTS:
%   LL      - Composite likelihood
% 
% COMMENTS:
%  This is a helper function for various multivariate GARCH models. A MEX file version of this file 
%  is available which provides a large speed-up.
%
%  See also DCC, SCALAR_VEC_VECH, BEKK, RARCH

q = size(indices,1);
[m,n] = size(data);
likConst = 3.67575413281869;
ll = 0;
if m==n
    for k=1:q
        i = indices(k,1);
        j = indices(k,2);
        s11 = S(i,i);
        s12 = S(i,j);
        s22 = S(j,j);
        det = s11*s22-s12*s12;
        x11 = data(i,i);
        x12 = data(i,j);
        x22 = data(j,j);
        ll = ll + 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/q;
    end
else
    for k=1:q
        i = indices(k,1);
        j = indices(k,2);
        s11 = S(i,i);
        s12 = S(i,j);
        s22 = S(j,j);
        det = s11*s22-s12*s12;
        x11 = data(i)*data(i);
        x12 = data(i)*data(j);
        x22 = data(j)*data(j);
        ll = ll + 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/q;
    end    
end