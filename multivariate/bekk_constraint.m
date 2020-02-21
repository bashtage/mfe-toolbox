function [c,ceq] = bekk_constraint(parameters,data,dataAsym,p,o,q,backCast,backCastAsym,type) %#ok<*INUSL>
% Non-linear constraint for estimation of BEKK(p,o,q) multivariate volatility models
%
% USAGE:
%  [C,CEQ] = bekk_constraint(PARAMETERS,DATA,DATAASYM,P,O,Q,BACKCAST,BACKCASTASYM,TYPE) 
%
% INPUTS:
%   See bekk_likelihood
%
% OUTPUTS:
%   C   - Vector of inequality constraints
%   CEQ - Empty
%
% COMMENTS:
%
%  EXAMPLES:
%
% See also BEKK, BEKK_LIKELIHOOD

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

ceq = [];
k = size(data,2);
[~,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type);

switch type
    case 1
        c = sum(A(1,1,:).^2,3) + sum(B(1,1,:).^2,3) + 0.5*sum(G(1,1,:).^2,3) - 1;
    case 2
        c = diag(sum(A.^2,3) + sum(B.^2,3) + 0.5*sum(G.^2,3) - 1);
    case 3
        m = zeros(k*k);
        for i=1:p
            m = m + kron(A(:,:,i),A(:,:,i));
        end
        for i=1:o
            m = m + 0.5*kron(G(:,:,i),G(:,:,i));
        end
        for i=1:q
            m = m + kron(B(:,:,i),B(:,:,i));
        end
        c = abs(eig(m)) - .99998;
end
